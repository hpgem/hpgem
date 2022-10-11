#include <numeric>
#include <petscsys.h>
#include "JacobiDavidsonMaxwell.h"
#include "Logger.h"

namespace hpgem {

namespace EigenSolvers {

extern "C" {
// lapack diag
void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda, double *w,
            double *work, int *lwork, int *info);
}

JacobiDavidsonMaxwellSolver::JacobiDavidsonMaxwellSolver() {}

PetscErrorCode JacobiDavidsonMaxwellSolver::clean(){
    MatDestroy(&this->AmI);
    KSPDestroy(&this->ksp);
    MatDestroy(&this->Y);
    MatDestroy(&this->H);
    BVDestroy(&this->Qt);
    BVDestroy(&this->eigenvectors);
    BVDestroy(&this->V);
    VecDestroy(&this->search_vect);
    VecDestroy(&this->residue_vect);
    return(0);
}

void JacobiDavidsonMaxwellSolver::setMaxIter(int niter) {
    this->maxIter = niter;
}
void JacobiDavidsonMaxwellSolver::setSearchSpaceMaxSize(int n) {
    this->search_space_maxsize = n;
}

void JacobiDavidsonMaxwellSolver::setSearchSpaceRestartSize(int n) {
    this->search_space_restart_size = n;
}

void JacobiDavidsonMaxwellSolver::setCorrectionNiter(int n) {
    this->correction_niter = n;
}
void JacobiDavidsonMaxwellSolver::setTolerance(PetscReal tol) {
    this->tolerance = tol;
}
void JacobiDavidsonMaxwellSolver::setTarget(PetscReal target) {
    this->ev_target = target;
}
PetscInt JacobiDavidsonMaxwellSolver::getConverged() {
    return this->eigenvectors_current_size;
}

PetscInt JacobiDavidsonMaxwellSolver::getIterationCount() { return this->iter; }

PetscErrorCode JacobiDavidsonMaxwellSolver::getEigenPair(PetscInt index,
                                                         PetscScalar &eval,
                                                         Vec &evec) {
    // return the eigen value at the position index
    PetscErrorCode ierr;
    Vec tmp_vec;
    eval = this->eigenvalues[index];
    ierr = BVGetColumn(this->eigenvectors, index, &tmp_vec);
    VecCopy(tmp_vec, evec);
    BVRestoreColumn(this->eigenvectors, index, &tmp_vec);
    VecDestroy(&tmp_vec);
    return (0);
}

void JacobiDavidsonMaxwellSolver::setMatrices(const Mat Ain, const Mat Cin) {
    // Set the A and C matrices describing the problem
    // Ax = lambda x with C x = 0
    PetscErrorCode ierr;
    PetscInt n, m;
    PetscViewer viewer;
    PetscBool ishermitian;

    this->A = Ain;
    this->C = Cin;

    // chek A hermitian
    // MatIsHermitian(Ain, 1E-6, &ishermitian);
    // MatSetOption(this->A, MAT_HERMITIAN, PETSC_TRUE);

    ierr = MatGetSize(Ain, &n, &m);
    logger(INFO, "JacobiDavidsonSolver System Size : % x %", n, m);

    ierr = MatGetSize(Cin, &n, &m);
    logger(INFO, "JacobiDavidsonSolver Contraint Size : % x %", n, m);

}

void JacobiDavidsonMaxwellSolver::setLinearSystem() {

    // compute A - I
    MatDuplicate(this->A, MAT_COPY_VALUES, &this->AmI);
    MatShift(this->AmI, -1.0);

    KSPCreate(PETSC_COMM_WORLD, &this->ksp);
    KSPSetType(this->ksp, KSPGMRES);
    KSPSetOperators(this->ksp, this->AmI, this->AmI);
    KSPSetTolerances(this->ksp, 1e-6, 1.e-12, PETSC_DEFAULT, 10);
    KSPSetFromOptions(this->ksp);
}


void JacobiDavidsonMaxwellSolver::initializeMatrices() {

    // initialize the Y and H matrix needed for the JD algorithm
    PetscInt y_nrows, y_ncols;
    PetscInt m_nrows, m_ncols;
    PetscErrorCode ierr;

    // we don't transpose this->C as it would change
    // the matrix owned by the calling code
    // Since M = I, Y = C.T
    ierr = MatHermitianTranspose(this->C, MAT_INITIAL_MATRIX, &this->Y);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);

    // Compute H = Y.T M Y = C M Y = C Y
    // M is indentity
    ierr = MatProductCreate(this->C, this->Y, NULL, &this->H);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);

    ierr = MatProductSetType(this->H, MATPRODUCT_AB);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);

    ierr = MatProductSetFromOptions(this->H);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);

    ierr = MatProductSymbolic(this->H);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);

    ierr = MatProductNumeric(this->H);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
}

void JacobiDavidsonMaxwellSolver::initializeVectors() {

    PetscInt y_nrows, y_ncols;
    PetscInt y_local_nrows, y_local_ncols;

    MatGetSize(this->Y, &y_nrows, &y_ncols);
    MatGetLocalSize(this->Y, &y_local_nrows, &y_local_ncols);

    // init the initial search vector
    VecCreate(PETSC_COMM_WORLD, &this->search_vect);
    VecSetSizes(this->search_vect, y_local_nrows, y_nrows);
    VecSetFromOptions(this->search_vect);

    // Set to 1
    VecSet(this->search_vect, 1.0);

    // normalize the search space vector
    normalizeVector(this->search_vect);

    // project / normalize
    projectCorrectionVector(this->search_vect);
    normalizeVector(this->search_vect);

    // init the residue
    VecDuplicate(this->search_vect, &this->residue_vect);
}

void JacobiDavidsonMaxwellSolver::initializeSearchSpace(int nev) {

    PetscInt ncols, nrows;
    PetscInt local_ncols, local_nrows;

    // get sizes
    MatGetSize(this->A, &nrows, &ncols);
    MatGetLocalSize(this->A, &local_nrows, &local_ncols);

    // set the current size
    this->search_space_current_size = 1;

    // allocate the search space
    this->Qt_current_size = 1;
    BVCreate(PETSC_COMM_WORLD, &this->Qt);
    BVSetSizes(this->Qt, local_nrows, nrows, this->search_space_maxsize);
    BVSetFromOptions(this->Qt);
    BVInsertVec(this->Qt, 0, this->search_vect);
    BVSetActiveColumns(this->Qt, 0, this->search_space_current_size);

    // allocate the search space
    this->V_current_size = 1;
    BVCreate(PETSC_COMM_WORLD, &this->V);
    BVSetSizes(this->V, local_nrows, nrows, this->search_space_maxsize);
    BVSetFromOptions(this->V);
    BVInsertVec(this->V, 0, this->search_vect);
    BVSetActiveColumns(this->V, 0, this->V_current_size);

    // allocate the eigenvector space
    this->eigenvectors_current_size = 0;
    BVCreate(PETSC_COMM_WORLD, &this->eigenvectors);
    BVSetSizes(this->eigenvectors, local_nrows, nrows, nev);
    BVSetFromOptions(this->eigenvectors);
    BVSetActiveColumns(this->eigenvectors, 0, this->eigenvectors_current_size);
}

PetscErrorCode JacobiDavidsonMaxwellSolver::computeRayleighQuotient(
    const Vec &x, PetscReal *out) {
    // compute x.T A x / (x.T x) = a / m
    PetscReal a, m;
    PetscScalar mcplx;

    weightedVectorDot(x, this->A, &a);
    VecDot(x, x, &mcplx);
    m = PetscRealPart(mcplx);
    *out = (abs(m) > 1E-12) ? a / m : 1.0;
    return (0);
}

PetscErrorCode JacobiDavidsonMaxwellSolver::normalizeVector(Vec &x) {

    // compute x / sqrt(x.T M X)
    PetscReal norm;
    PetscReal inv_sqrt_norm;
    PetscBool isnan, isinf;

    // create x = x.T M * x
    // with M = I
    VecNorm(x, NORM_2, &norm);

    // inverse of the norm
    inv_sqrt_norm = 1.0 / norm;

    isnan = PetscIsNanReal(inv_sqrt_norm);
    isinf = PetscIsInfReal(inv_sqrt_norm);

    if (isnan || isinf) {
        PetscPrintf(PETSC_COMM_WORLD, " Warning : Nan detected in the norm\n");
        VecSet(x, 0.0);
        exit(0);
    } else {
        // scale the vector
        VecScale(x, inv_sqrt_norm);
    }

    VecNorm(x, NORM_2, &norm);
    if (norm==0.0)
    {
        PetscPrintf(PETSC_COMM_WORLD, " Norm zero after normalization\n");
        exit(0);
    }

    return (0);
}

PetscErrorCode JacobiDavidsonMaxwellSolver::weightedVectorDot(const Vec &x,
                                                              const Mat &K,
                                                              PetscReal *val) {
    // Compute val = x.T K X
    Vec tmp;
    PetscScalar norm;

    // create tmp = K * x
    VecDuplicate(x, &tmp);
    MatMult(K, x, tmp);

    // compute the norm
    VecDot(x, tmp, &norm);
    *val = PetscRealPart(norm);

    // clean up
    VecDestroy(&tmp);

    return (0);
}

PetscErrorCode JacobiDavidsonMaxwellSolver::staticMatMultCorrPrec(Mat M, Vec x,
                                                                  Vec y) {
    PetscErrorCode error;

    JacobiDavidsonMaxwellSolver *workspace;
    error = MatShellGetContext(M, &workspace);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    workspace->correctionPreconditionerMatMult(x, y);
    return 0;
}

PetscErrorCode JacobiDavidsonMaxwellSolver::correctionPreconditionerMatMult(
    Vec x, Vec y) {

    // solve the linear system
    KSPSolve(this->ksp, x, y);

    return (0);
}

PetscErrorCode JacobiDavidsonMaxwellSolver::staticMatMultCorrOp(Mat M, Vec x,
                                                                Vec y) {
    PetscErrorCode error;

    JacobiDavidsonMaxwellSolver *workspace;
    error = MatShellGetContext(M, &workspace);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    workspace->correctionOperatorMatMult(x, y);
    return 0;
}

PetscErrorCode JacobiDavidsonMaxwellSolver::correctionOperatorMatMult(Vec x,
                                                                      Vec y) {

    // compute (I - M Q Q.T) (A - eta M) (I - Q Q.T M) y
    // here M is always the identity matrix
    Vec Ay;

    // copy x in y
    VecCopy(x, y);

    // y = y - Q Q.T y
    computeProjection(y, this->Qt);

    // tmp = A y
    VecDuplicate(x, &Ay);
    MatMult(this->A, y, Ay);

    // Compute (A - eta M) y
    // here equal A y  - eta y
    VecAXPY(Ay, -1.0 * this->eta, y);

    // y = y - Q Q.T y
    computeProjection(Ay, this->Qt);
    VecCopy(Ay, y);

    VecDestroy(&Ay);

    return (0);
}

PetscErrorCode JacobiDavidsonMaxwellSolver::solveCorrectionEquation(
    const Vec &res, Vec &sol) {
    // solve the correction equation given by
    // (I - M Q Q.T) (A - eta M) (I - Q Q.T M) sol = -res
    // with M being identity

    Mat op;
    Mat prec;

    KSP ksp;
    PC pc;
    PetscInt ncols, nrows;
    PetscInt local_ncols, local_nrows;
    PetscErrorCode ierr;
    PetscInt its;
    bool print_time_local;

    auto tstart = std::chrono::high_resolution_clock::now();

    // get the linear operator
    MatGetSize(this->A, &nrows, &ncols);
    MatGetLocalSize(this->A, &local_nrows, &local_ncols);
    MatCreateShell(PETSC_COMM_WORLD, local_nrows, local_ncols, nrows, ncols,
                   this, &op);
    MatShellSetOperation(op, MATOP_MULT, (void (*)(void))staticMatMultCorrOp);

    MatCreateShell(PETSC_COMM_WORLD, local_nrows, local_ncols, nrows, ncols,
                   this, &prec);
    MatShellSetOperation(prec, MATOP_MULT,
                         (void (*)(void))staticMatMultCorrPrec);

    // set up the linear system
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPGMRES);
    KSPSetOperators(ksp, op, op);
    ierr = KSPSetFromOptions(ksp);
    CHKERRQ(ierr);
    KSPSetTolerances(ksp, 1e-12, 1.e-12, PETSC_DEFAULT, this->correction_niter);

    // // preconditionner
    // KSPGetPC(ksp, &pc);
    // PCSetType(pc, PCGAMG);

    // solve the linear system
    KSPSolve(ksp, res, sol);
    VecScale(sol, -1.0);

    // get iteration count
    KSPGetIterationNumber(ksp, &its);

    // clean
    KSPDestroy(&ksp);
    MatDestroy(&op);

    auto tstop = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::microseconds>(tstop - tstart);

    logger(DEBUG, "  -- solveCorrectionEquation done in %d mu s [%d its]\n",
           duration.count(), its);

    MatDestroy(&op);
    MatDestroy(&prec);

    return (0);
}

PetscErrorCode JacobiDavidsonMaxwellSolver::computeResidueVector(
    const Vec &q, const PetscReal rho, Vec &res) {

    // compute the residue res = A q - rho M q
    // here M is always the identity matrix

    auto tstart = std::chrono::high_resolution_clock::now();

    // compute res  = Aq
    MatMult(this->A, q, res);

    // assemble res = A q - rho M q
    VecAXPY(res, -1.0 * rho, q);

    auto tstop = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::microseconds>(tstop - tstart);

    logger(DEBUG, "  -- computeResidueVector done in %d mu s \n",
           duration.count());

    return (0);
}

PetscErrorCode JacobiDavidsonMaxwellSolver::computeProjection(Vec &v,
                                                              const BV &Q) {
    // Compute v = v - Q Q.T v
    PetscScalar *tmp_row_vector;
    PetscInt local_size, ncol, nrow;

    auto tstart = std::chrono::high_resolution_clock::now();

    // get sizes
    BVGetSizes(Q, &local_size, &nrow, &ncol);

    // create tmp_row_vector
    tmp_row_vector = new PetscScalar[ncol];

    // compute tmp_row_vect = Q.T v
    BVDotVec(Q, v, tmp_row_vector);

    // v = v - Q Q.T v
    BVMultVec(Q, -1.0, 1.0, v, tmp_row_vector);

    auto tstop = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::microseconds>(tstop - tstart);

    logger(DEBUG, "  -- computeProjection done in %d mu s \n",
           duration.count());

    delete [] tmp_row_vector;

    return (0);
}

PetscErrorCode JacobiDavidsonMaxwellSolver::projectCorrectionVector(Vec &corr) {

    // compute corr = corr - Y H^{-1} C.T corr

    Vec rhs;
    PetscInt c_ncols, c_nrows, c_local_ncols, c_local_nrows;
    PetscInt corr_size, rhs_size;
    KSP ksp;
    PetscInt its;

    auto tstart = std::chrono::high_resolution_clock::now();

    MatGetSize(this->Y, &c_nrows, &c_ncols);
    MatGetLocalSize(this->Y, &c_local_nrows, &c_local_ncols);

    // create tmp_row vec
    VecCreate(PETSC_COMM_WORLD, &rhs);
    VecSetSizes(rhs, c_local_ncols, c_ncols);
    VecSetFromOptions(rhs);

    // compute rhs = C corr
    MatMult(this->C, corr, rhs);

    // set up the linear system
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPGMRES);
    KSPSetOperators(ksp, this->H, this->H);
    KSPSetFromOptions(ksp);
    KSPSetTolerances(ksp, 1e-12, 1.e-12, PETSC_DEFAULT, 100);

    // solve the linear system
    //  H ksp_sol = (C corr)
    // -> ksp_sol = H^{-1} C corr
    // we store the solution (ksp_sol) in rhs
    KSPSolve(ksp, rhs, rhs);
    KSPGetIterationNumber(ksp, &its);

    // compute corr - Y H^{-1} C corr
    // with rhs = H^{-1} C corr
    VecScale(rhs, -1.0);
    MatMultAdd(this->Y, rhs, corr, corr);

    VecDestroy(&rhs);
    KSPDestroy(&ksp);

    auto tstop = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::microseconds>(tstop - tstart);

    logger(DEBUG, "  -- projectCorrectionVector done in %d mu s [%d its]\n",
           duration.count(), its);

    return (0);
}


PetscReal JacobiDavidsonMaxwellSolver::computeSmallResidue(
    const Mat &A, const Vec &x, const PetscScalar lambda)
{
    PetscReal resVal;
    Vec resVec;
    Vec tmp;

    VecDuplicate(x, &resVec);
    VecDuplicate(x, &tmp);
    VecCopy(tmp, x);
    VecScale(tmp, -lambda);
    MatMultAdd(A, x, tmp, resVec);
    VecNorm(resVec, NORM_2, &resVal);

    VecDestroy(&resVec);
    VecDestroy(&tmp);

    return(resVal);
}
PetscErrorCode JacobiDavidsonMaxwellSolver::computeSmallEigenvalues(
    std::vector<PetscReal> &eigenvalues, Vec *eigenvectors) {
    // project the total matrix following
    // Ak = V.T A V
    // and computes the eigenvalues of Ak using LAPACK

    Mat Ak, AkT;
    EPS eps;
    PetscErrorCode ierr;
    PetscInt nconv;
    PetscScalar eigr;
    Vec Vr, Vi;

    PetscInt sort_idx[this->V_current_size];
    PetscScalar eval_tmp[this->V_current_size];
    PetscReal eval_tmp_real[this->V_current_size];
    Vec evec_tmp[this->V_current_size];
    PetscReal small_res;
    PetscReal norm;
    PetscInt i;

    auto tstart = std::chrono::high_resolution_clock::now();

    MatCreateSeqDense(PETSC_COMM_SELF, this->V_current_size,
                      this->V_current_size, NULL, &Ak);
    MatSetOption(Ak, MAT_HERMITIAN, PETSC_TRUE);

    BVMatProject(this->V, this->A, this->V, Ak);

    // make sur it's hermitian otherwise LAPACK crashes
    MatCreateHermitianTranspose(Ak, &AkT);
    MatAXPY(Ak, 1.0, AkT, SAME_NONZERO_PATTERN);
    MatScale(Ak, 0.5);

    // create evects
    // MatCreateVecs(Ak, NULL, &Vr);
    // MatCreateVecs(Ak, NULL, &Vi);
    VecCreateSeq(PETSC_COMM_SELF, this->V_current_size, &Vr);
    VecCreateSeq(PETSC_COMM_SELF, this->V_current_size, &Vi);

    // create eps object
    ierr = EPSCreate(PETSC_COMM_SELF, &eps);
    CHKERRQ(ierr);
    ierr = EPSSetOperators(eps, Ak, NULL);
    CHKERRQ(ierr);
    ierr = EPSSetProblemType(eps, EPS_HEP);
    CHKERRQ(ierr);
    ierr = EPSSetType(eps, EPSLAPACK);
    ierr = EPSSetFromOptions(eps);
    CHKERRQ(ierr);
    ierr = EPSSolve(eps);
    CHKERRQ(ierr);
    EPSGetConverged(eps, &nconv);

    for (i = 0; i < nconv; i++) {
        ierr = EPSGetEigenpair(eps, i, &eigr, nullptr, Vr, Vi);
        CHKERRQ(ierr);
        eval_tmp[i] = eigr;
        sort_idx[i] = i;
        // ierr = MatCreateVecs(Ak, NULL, &evec_tmp[i]);
        ierr = VecCreateSeq(PETSC_COMM_SELF, this->V_current_size, &evec_tmp[i]);
        CHKERRQ(ierr);
        ierr = VecCopy(Vr, evec_tmp[i]);
        CHKERRQ(ierr);
    }

    // extract the real part of the eigenvalues
    for (i = 0; i < this->V_current_size; i++) {
        eval_tmp_real[i] = PetscRealPart(eval_tmp[i]);
    }

    // get sort indexes eigenvalues
    PetscSortRealWithArrayInt(this->V_current_size, eval_tmp_real, sort_idx);

    // sort eigenvalues
    eigenvalues.clear();
    for (i = 0; i < nconv; i++) {

        if (eval_tmp_real[i] < 1E-6) {

            VecNorm(evec_tmp[sort_idx[i]], NORM_2, &norm);
            small_res = this->computeSmallResidue(Ak, evec_tmp[sort_idx[i]], eval_tmp_real[i] );
            logger(WARN, " Zero Eigenvalue detected: ev[%]=% (norm = %) (res = %) (size = %)", i,
                   eval_tmp_real[i], norm, small_res, this->V_current_size);
            // BVView(this->V, PETSC_VIEWER_STDOUT_SELF);
        }
        
        eigenvalues.push_back(eval_tmp_real[i]);
        ierr = VecCreateSeq(PETSC_COMM_SELF, this->V_current_size, &eigenvectors[i]);
        CHKERRQ(ierr);
        ierr = VecCopy(evec_tmp[sort_idx[i]], eigenvectors[i]);
        CHKERRQ(ierr);
        
    }

    // clean up
    EPSDestroy(&eps);
    MatDestroy(&Ak);
    MatDestroy(&AkT);
    VecDestroy(&Vr);
    VecDestroy(&Vi);

    for (i = 0; i < nconv; i++) {
        VecDestroy(&evec_tmp[i]);
    }

    auto tstop = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::microseconds>(tstop - tstart);

    if (this->print_small_evs) {
        for (int ii = 0; ii < this->V_current_size; ii++) {
            PetscPrintf(PETSC_COMM_WORLD, "eigenvalue: %6.12f\n",
                        eigenvalues[ii]);
        }
    }

    logger(DEBUG, "  -- computeSmallEigenvalues done in %d mu s \n",
           duration.count());

    return (0);
}

PetscErrorCode JacobiDavidsonMaxwellSolver::computeThreshold(Vec q, Vec r,
                                                             PetscReal *eps) {
    PetscScalar num, denom;

    // compute q.T M q
    // with M = I
    VecDot(q, q, &denom);

    // compute r.T M^{-1} r
    // with M^{-1} = I
    VecDot(r, r, &num);

    // compute the threshold
    *eps = sqrt(PetscRealPart(num) / PetscRealPart(denom));

    return (0);
}

PetscErrorCode JacobiDavidsonMaxwellSolver::solve(PetscInt nev) {

    PetscReal rho, eps, norm, norm_cq;
    std::vector<PetscReal> small_evals;
    PetscInt n,m;
    Vec small_evects[this->search_space_maxsize];
    Vec tmp_v[this->search_space_maxsize];
    Vec residue_vect_copy, correction_vect;
    Vec q;
    PetscInt k, idx_evect, ii;
    PetscBool found, stop;
    PetscScalar *vec_values;
    PetscErrorCode ierr;
    PetscMPIInt rank;

    this->nconverged = 0;

    logger(INFO, "Sovling EigenValue problem using JacobiDavidsonSolver");

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    if (this->search_space_restart_size >= search_space_maxsize) {
        logger(ERROR, "Required % eigenvalues but the maximum size of the search space is set to %",
               this->search_space_restart_size, this->search_space_maxsize);
        exit(0);
    }

    for (ii=0;ii<this->search_space_maxsize;ii++){
        small_evects[ii] = PETSC_NULL;
        tmp_v[ii] = PETSC_NULL;
    }


    // init the solver matrices and vectors
    initializeMatrices();
    initializeVectors();
    initializeSearchSpace(nev);

    // init the linear solver
    setLinearSystem();

    // rayleigh quotient
    computeRayleighQuotient(this->search_vect, &rho);

    // compute initial residue vector
    computeResidueVector(this->search_vect, rho, this->residue_vect);

    // duplicate residue vector and create the correction vector
    ierr = VecDuplicate(this->residue_vect, &residue_vect_copy);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = VecDuplicate(this->residue_vect, &correction_vect);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);

    k = 0;
    stop = PETSC_FALSE;

    for (this->iter = 0; this->iter < this->maxIter; this->iter++) {

        // determine eta
        if(rank == 0){
            eps = (PetscReal)(std::rand()) / (PetscReal)(RAND_MAX);
            this->eta = (this->iter > 0 && eps < 0.5) ? rho : this->ev_target;
        }
        MPI_Bcast(&this->eta, 1, MPIU_REAL, 0, MPI_COMM_WORLD);
        
        logger(DEBUG, "JacobiDavidsonSolver iteration : %, k = %, rho = % eta = %",
               iter, k, rho, this->eta);

        // copy res to res_new and left project
        ierr = VecCopy(this->residue_vect, residue_vect_copy);
        CHKERRABORT(PETSC_COMM_WORLD, ierr);
        computeProjection(residue_vect_copy, this->Qt);

        // solve correction equation
        solveCorrectionEquation(residue_vect_copy, correction_vect);

        // ortho correction vector
        computeProjection(correction_vect, this->Qt);

        // project correction vector
        projectCorrectionVector(correction_vect);

        // create search vector
        ierr = VecCopy(correction_vect, search_vect);
        CHKERRABORT(PETSC_COMM_WORLD, ierr);
        computeProjection(search_vect, this->V);
        normalizeVector(search_vect);

        k++;

        // Add column to V
        ierr = BVInsertVec(this->V, this->V_current_size, search_vect);
        CHKERRABORT(PETSC_COMM_WORLD, ierr);
        this->V_current_size++;
        ierr = BVSetActiveColumns(this->V, 0, this->V_current_size);
        CHKERRABORT(PETSC_COMM_WORLD, ierr);

        // compute eigenpairs on small subspace
        ierr = computeSmallEigenvalues(small_evals, small_evects);
        CHKERRABORT(PETSC_COMM_WORLD, ierr);

        // create q
        ierr = VecDuplicate(this->residue_vect, &q);
        CHKERRABORT(PETSC_COMM_WORLD, ierr);

        found = PETSC_TRUE;
        idx_evect = 0;
        while (found) {

            // rho = Lk[0]
            rho = small_evals[0];

            // q = V sk0
            ierr = VecGetArray(small_evects[0], &vec_values);
            CHKERRABORT(PETSC_COMM_WORLD, ierr);
            ierr = BVMultVec(this->V, 1.0, 0.0, q, vec_values);
            CHKERRABORT(PETSC_COMM_WORLD, ierr);
            ierr = VecRestoreArray(small_evects[0], &vec_values);
            CHKERRABORT(PETSC_COMM_WORLD, ierr);

            // r = A q - rho M q
            computeResidueVector(q, rho, this->residue_vect);

            // compute residue threshold
            computeThreshold(q, this->residue_vect, &eps);
       
            // compute covnergence criteria
            found = (eps < this->tolerance && k > 0) ? PETSC_TRUE : PETSC_FALSE;
            if (found) {
                logger(INFO,
                       "JacobiDavidsonSolver new eigenvector found with "
                       "eigenvalue %",
                       rho);
                // store the eigenvalue/eigenvector
                this->eigenvalues.push_back(rho);

                // add a new vector to the basis
                ierr = BVInsertVec(this->eigenvectors,
                                   this->eigenvectors_current_size, q);
                CHKERRABORT(PETSC_COMM_WORLD, ierr);
                this->eigenvectors_current_size++;

                if (this->eigenvectors_current_size >= nev) {
                    stop = PETSC_TRUE;
                    break;
                }

                // Compute [tmp_v1 .... tmp_vk-1] = V [s2, .... sk]
                for (ii = 0; ii < k - 1; ii++) {
                    ierr = VecGetArray(small_evects[ii + 1], &vec_values);
                    CHKERRABORT(PETSC_COMM_WORLD, ierr);
                    ierr = VecDuplicate(this->residue_vect, &tmp_v[ii]);
                    CHKERRABORT(PETSC_COMM_WORLD, ierr);
                    ierr = VecSet(tmp_v[ii], 0.0);
                    CHKERRABORT(PETSC_COMM_WORLD, ierr);
                    ierr = BVMultVec(this->V, 1.0, 0.0, tmp_v[ii], vec_values);
                    CHKERRABORT(PETSC_COMM_WORLD, ierr);
                    ierr = VecRestoreArray(small_evects[ii + 1], &vec_values);
                    CHKERRABORT(PETSC_COMM_WORLD, ierr);
                }

                // insert [tmp_v1 .... tmp_vk-1] in V
                for (ii = 0; ii < k - 1; ii++) {
                    ierr = BVInsertVec(this->V, ii, tmp_v[ii]);
                    CHKERRABORT(PETSC_COMM_WORLD, ierr);
                    ierr = VecDestroy(&tmp_v[ii]);
                    CHKERRABORT(PETSC_COMM_WORLD, ierr);
                }

                // set all other columns to 0
                for (ii = k - 1; ii < this->search_space_maxsize; ii++) {
                    ierr = BVScaleColumn(this->V, ii, 0.0);
                    CHKERRABORT(PETSC_COMM_WORLD, ierr);
                }

                // number of active cols of V
                ierr = BVSetActiveColumns(this->V, 0, k - 1);
                CHKERRABORT(PETSC_COMM_WORLD, ierr);
                this->V_current_size = k - 1;

                // remove first eigenvalue
                small_evals.erase(small_evals.begin());

                // set the eigenvectors to I(k-1)
                for (ii = 0; ii < k; ii++) {
                    VecSet(small_evects[ii], 0.0);
                    if (ii < k - 1) {
                        ierr = VecSetValue(small_evects[ii], ii, 1.0,
                                           INSERT_VALUES);
                        CHKERRABORT(PETSC_COMM_WORLD, ierr);
                        ierr = VecAssemblyBegin(small_evects[ii]);
                        CHKERRABORT(PETSC_COMM_WORLD, ierr);
                        ierr = VecAssemblyEnd(small_evects[ii]);
                        CHKERRABORT(PETSC_COMM_WORLD, ierr);
                    }
                }
                k--;
                idx_evect++;
            }
        }

        if (stop) break;

        // Qt = [Q; q]
        if (this->eigenvectors_current_size == 0) {
            ierr = BVInsertVec(this->Qt, 0, q);
            CHKERRABORT(PETSC_COMM_WORLD, ierr);
            this->Qt_current_size = 1;
            ierr = BVSetActiveColumns(this->Qt, 0, this->Qt_current_size);
            CHKERRABORT(PETSC_COMM_WORLD, ierr);
        } else {
            ierr = BVCopy(this->eigenvectors, this->Qt);
            CHKERRABORT(PETSC_COMM_WORLD, ierr);
            ierr = BVInsertVec(this->Qt, this->eigenvectors_current_size, q);
            CHKERRABORT(PETSC_COMM_WORLD, ierr);
            this->Qt_current_size = this->eigenvectors_current_size + 1;
            ierr = BVSetActiveColumns(this->Qt, 0, this->Qt_current_size);
            CHKERRABORT(PETSC_COMM_WORLD, ierr);

            for (ii = this->Qt_current_size; ii < this->search_space_maxsize;
                 ii++) {
                ierr = BVScaleColumn(this->Qt, ii, 0.0);
                CHKERRABORT(PETSC_COMM_WORLD, ierr);
            }
        }

        // restart
        if (k == this->search_space_maxsize - 1) {

            logger(INFO, "JacobiDavidsonSolver Restart [%]", this->search_space_restart_size);

            // compute tmp = V Sk[:,:size_min]
            for (ii = 0; ii < this->search_space_restart_size; ii++) {

                ierr = VecGetArray(small_evects[ii], &vec_values);
                CHKERRABORT(PETSC_COMM_WORLD, ierr);
                ierr = VecDuplicate(this->residue_vect, &tmp_v[ii]);
                CHKERRABORT(PETSC_COMM_WORLD, ierr);
                ierr = VecSet(tmp_v[ii], 0.0);
                CHKERRABORT(PETSC_COMM_WORLD, ierr);

                ierr = BVMultVec(this->V, 1.0, 0.0, tmp_v[ii], vec_values);
                CHKERRABORT(PETSC_COMM_WORLD, ierr);
                ierr = VecRestoreArray(small_evects[ii + 1], &vec_values);
                CHKERRABORT(PETSC_COMM_WORLD, ierr);
            }

            // insert [tmp_v1 .... tmp_vk-1] in V
            for (ii = 0; ii < search_space_restart_size; ii++) {
                ierr = BVInsertVec(this->V, ii, tmp_v[ii]);
                CHKERRABORT(PETSC_COMM_WORLD, ierr);
                ierr = VecDestroy(&tmp_v[ii]);
                CHKERRABORT(PETSC_COMM_WORLD, ierr);
            }

            // set all other columns to 0
            for (ii = search_space_restart_size; ii < this->search_space_maxsize;
                 ii++) {
                ierr = BVScaleColumn(this->V, ii, 0.0);
                CHKERRABORT(PETSC_COMM_WORLD, ierr);
            }

            // number of active cols of V
            ierr = BVSetActiveColumns(this->V, 0, k - 1);
            CHKERRABORT(PETSC_COMM_WORLD, ierr);
            this->V_current_size = search_space_restart_size;

            k = this->search_space_restart_size;
        }

        // remove tmp vec
        for (ii=0; ii<this->search_space_maxsize; ii++){
            VecDestroy(&small_evects[ii]);
            VecDestroy(&tmp_v[ii]);
        }
        VecDestroy(&q);

    }

    ierr = this->orderEigenvalues();
    CHKERRABORT(PETSC_COMM_WORLD, ierr);

    // clean mem    
    VecDestroy(&residue_vect_copy);
    VecDestroy(&correction_vect);
    
    return (0);
}

PetscErrorCode JacobiDavidsonMaxwellSolver::orderEigenvalues() {

    PetscInt converged;
    std::vector<PetscScalar> tmp_vec; 
    converged = this->getConverged();
    tmp_vec = this->eigenvalues;
    std::vector<std::size_t> ordering(converged);

    // set the ordering
    std::iota(ordering.begin(), ordering.end(), 0);

    auto lcomp = [this](int a, int b) {
        return PetscRealPart(this->eigenvalues[a]) <
               PetscRealPart(this->eigenvalues[b]);
    };
    std::sort(ordering.begin(), ordering.end(), lcomp);

    for(int ii=0; ii<converged;ii++){
        this->eigenvalues[ii] = tmp_vec[ordering[ii]];
    }

    return(0);

}

}  // namespace EigenSolvers

}  // namespace hpgem
#include "JacobiDavidsonMaxwell.h"
#include "Logger.h"

namespace hpgem {

namespace LinearAlgebra {

extern "C" {
// lapack diag
void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda,
           double *w, double *work, int *lwork, int *info);
}


JacobiDavidsonMaxwellSolver::JacobiDavidsonMaxwellSolver() {}

void JacobiDavidsonMaxwellSolver::set_maxIter(int niter){ this->maxIter = niter; }
void JacobiDavidsonMaxwellSolver::set_search_space_maxsize(int n){ this->search_space_maxsize = n; }
void JacobiDavidsonMaxwellSolver::set_correction_niter(int n){this->correction_niter = n;}
void JacobiDavidsonMaxwellSolver::set_tolerance(PetscReal tol){this->tolerance=tol;}

PetscInt JacobiDavidsonMaxwellSolver::getConverged(){ return this->Q_current_size;}

PetscInt JacobiDavidsonMaxwellSolver::getIterationCount(){ return this->iter; }

PetscErrorCode JacobiDavidsonMaxwellSolver::getEigenPair(PetscInt index, PetscScalar &eval, Vec &evec )
{
    // return the eigen value at the position index
    PetscErrorCode ierr;
    Vec tmp_vec;
    eval = this->eigenvalues[index];
    ierr = BVGetColumn(this->Q, index, &tmp_vec);
    VecCopy(tmp_vec, evec);
    BVRestoreColumn(this->Q, index, &tmp_vec);
    VecDestroy(&tmp_vec);
    return (0);

}

void JacobiDavidsonMaxwellSolver::setMatrices(const Mat &Ain, const Mat &Cin) 
{    
    // Set the A and C matrices describing the problem
    // Ax = lambda x with C x = 0
    PetscErrorCode ierr; 
    PetscInt n,m;
    PetscViewer viewer; 
    
    this->A = Ain;
    MatSetOption(this->A, MAT_HERMITIAN, PETSC_TRUE);
    this->C = Cin;

    ierr = MatGetSize(Ain, &n, &m);
    logger(INFO, "JacobiDavidsonSolver System Size : % x %", n,m);

}

void JacobiDavidsonMaxwellSolver::initializeMatrices() 
{

    // initialize the Y and H matrix needed for the JD algorithm
    PetscInt        y_nrows, y_ncols;
    PetscInt        m_nrows, m_ncols;
    Mat             Yh;
    PetscErrorCode  ierr;

    // transpose C
    ierr = MatHermitianTranspose(this->C, MAT_INPLACE_MATRIX, &this->C);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);

    // Since M = I, Y = C
    this->Y = C;
    
    
    // Compute H = Y.T M Y 
    // M is indentity
    // we need to explicitly construct the hermitian transpose of Y
    ierr = MatDuplicate(this->Y, MAT_COPY_VALUES, &Yh);
    ierr = MatHermitianTranspose(Yh, MAT_INPLACE_MATRIX, &Yh);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);

    ierr = MatProductCreate(Yh, this->Y, NULL, &this->H);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);

    ierr = MatProductSetType(this->H, MATPRODUCT_AB);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);

    ierr = MatProductSetFromOptions(this->H);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);

    ierr = MatProductSymbolic(this->H);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);

    ierr = MatProductNumeric(this->H); 
    CHKERRABORT(PETSC_COMM_WORLD, ierr);

    MatDestroy(&Yh);

}

void  JacobiDavidsonMaxwellSolver::initializeSearchSpace(int nev) {

    PetscInt        ncols, nrows;
    PetscInt        local_ncols, local_nrows;

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

    // allocate the search space
    this->Q_current_size = 0;
    BVCreate(PETSC_COMM_WORLD, &this->Q);
    BVSetSizes(this->Q, local_nrows, nrows, nev);
    BVSetFromOptions(this->Q);
    BVSetActiveColumns(this->Q, 0, this->Q_current_size);
}

void JacobiDavidsonMaxwellSolver::initializeVectors(){

    PetscRandom     rctx;
    PetscInt        y_nrows, y_ncols;
    PetscInt        y_local_nrows, y_local_ncols;

    MatGetSize(this->Y, &y_nrows, &y_ncols);
    MatGetLocalSize(this->Y, &y_local_nrows, &y_local_ncols);

    // init the initial search vector
    PetscRandomCreate(PETSC_COMM_WORLD, &rctx);
    PetscRandomSetFromOptions(rctx);
    VecCreate(PETSC_COMM_WORLD, &this->search_vect);
    VecSetSizes(this->search_vect, y_local_nrows, y_nrows);
    VecSetFromOptions(this->search_vect);
    // VecSetRandom(this->search_vect, rctx);

    // Set to 1 for now !!
    VecSet(this->search_vect, 1.0);

    // normalize the search space vector
    normalizeVector(this->search_vect);

    // project / normalize
    projectCorrectionVector(this->search_vect);
    normalizeVector(this->search_vect);

    // init the residue
    VecDuplicate(this->search_vect, &this->residue_vect);
}



PetscErrorCode JacobiDavidsonMaxwellSolver::computeRayleighQuotient(const Vec &x, PetscReal *out)
{
    // compute x.T A x / (x.T x)
    PetscScalar        a;
    PetscScalar         m;

    VecxTAx(x, this->A, &a);
    VecDot(x,x, &m);
    *out = (abs(PetscRealPart(m))>1E-12) ? PetscRealPart(a)/PetscRealPart(m) : 1.0;
    return(0);
}


PetscErrorCode JacobiDavidsonMaxwellSolver::normalizeVector(Vec &x)
{

    // compute x / sqrt(x.T M X)
    PetscScalar      norm; 
    PetscReal        inv_sqrt_norm;
    PetscBool        isnan, isinf;

    // create x = x.T M * x
    VecDot(x, x, &norm);

    // inverse of the norm
    inv_sqrt_norm = 1.0/sqrt(PetscRealPart(norm));

    isnan = PetscIsNanReal(inv_sqrt_norm);
    isinf = PetscIsInfReal(inv_sqrt_norm);

    if(isnan || isinf){
        PetscPrintf(PETSC_COMM_WORLD, " Warning : Nan detected in the norm\n");
        VecSet(x, 0.0);
    } else {
        // scale the vector
        VecScale(x, inv_sqrt_norm);
    }

    return(0);
}

PetscErrorCode JacobiDavidsonMaxwellSolver::VecxTAx(const Vec &x, const Mat &K, PetscScalar *val)
{
    // Compute val = x.T K X
    Vec             tmp;
    PetscInt        ncols, nrows;
    PetscInt        local_ncols, local_nrows;

    // get sizes
    MatGetSize(K, &nrows, &ncols);
    MatGetLocalSize(K, &local_nrows, &local_ncols);

    // create tmp = K * x
    VecCreate(PETSC_COMM_WORLD, &tmp);
    VecSetSizes(tmp,local_nrows,nrows);
    VecSetFromOptions(tmp);
    MatMult(K, x, tmp);

    // compute the norm
    VecDot(x, tmp, val);

    // clean up
    VecDestroy(&tmp);

    return(0);
}

PetscErrorCode JacobiDavidsonMaxwellSolver::staticMatMultCorrPrec(Mat M, Vec x, Vec y)
{
    PetscErrorCode error;

    JacobiDavidsonMaxwellSolver* workspace;
    error = MatShellGetContext(M, &workspace);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    workspace->correctionPreconditionerMatMult(x, y);
    return 0;
}

PetscErrorCode JacobiDavidsonMaxwellSolver::correctionPreconditionerMatMult(Vec x, Vec y)
{

    Mat             ApM;
    KSP             ksp;

    // create A-M
    MatDuplicate(this->A, MAT_COPY_VALUES, &ApM);
    MatShift(ApM, -1.0);

    KSPCreate(PETSC_COMM_WORLD,&ksp);
    KSPSetType(ksp,KSPGMRES);
    KSPSetOperators(ksp,ApM,ApM);
    KSPSetTolerances(ksp, 1e-6, 1.e-12, PETSC_DEFAULT, 10);
    KSPSetFromOptions(ksp);

    // solve the linear system
    KSPSolve(ksp, x, y);

    MatDestroy(&ApM);
    KSPDestroy(&ksp);

    return(0);

}



PetscErrorCode JacobiDavidsonMaxwellSolver::staticMatMultCorrOp(Mat M, Vec x, Vec y)
{
    PetscErrorCode error;

    JacobiDavidsonMaxwellSolver* workspace;
    error = MatShellGetContext(M, &workspace);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    workspace->correctionOperatorMatMult(x, y);
    return 0;
}

PetscErrorCode JacobiDavidsonMaxwellSolver::correctionOperatorMatMult(Vec x, Vec y)
{

    // compute (I - M Q Q.T) (A - eta M) (I - Q Q.T M) y
    // here M is always the identity matrix
    Vec             Ay;

    // copy x in y
    VecCopy(x, y);

    // y = y - Q Q.T y
    computeProjection(y, this->Qt);

    // tmp = A y
    VecDuplicate(x, &Ay);
    MatMult(this->A, y, Ay);

    // Compute (A - eta M) y
    // here equal A y  - eta y
    VecAXPY(Ay, -1.0*this->eta, y);

    // y = y - Q Q.T y
    computeProjection(Ay, this->Qt);
    VecCopy(Ay, y);

    VecDestroy(&Ay);

    return(0);

}


PetscErrorCode JacobiDavidsonMaxwellSolver::solveCorrectionEquation(const Vec &res, Vec &sol)
{
    // solve the correction equation given by 
    // (I - M Q Q.T) (A - eta M) (I - Q Q.T M) sol = -res
    // with M being identity

    Mat         op;
    Mat         prec;

    KSP         ksp;
    PC          pc;
    PetscInt    ncols, nrows;
    PetscInt    local_ncols, local_nrows;
    Vec         rhs;
    PetscErrorCode ierr;
    PetscInt    its;
    bool print_time_local;

    auto tstart = std::chrono::high_resolution_clock::now();

    print_time_local = this->print_time;
    this->print_time = false;

    // get the linear operator
    MatGetSize(this->A, &nrows, &ncols);
    MatGetLocalSize(this->A, &local_nrows, &local_ncols);
    MatCreateShell(PETSC_COMM_WORLD, local_nrows, local_ncols, nrows, ncols, this, &op);
    MatShellSetOperation(op, MATOP_MULT, (void (*)(void))staticMatMultCorrOp);


    MatCreateShell(PETSC_COMM_WORLD, local_nrows, local_ncols, nrows, ncols, this, &prec);
    MatShellSetOperation(prec, MATOP_MULT, (void (*)(void))staticMatMultCorrPrec);


    // set up the linear system
    KSPCreate(PETSC_COMM_WORLD,&ksp);
    KSPSetType(ksp,KSPMINRES);
    KSPSetOperators(ksp,op,op);
    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
    KSPSetTolerances(ksp, 1e-12, 1.e-12, PETSC_DEFAULT, this->correction_niter);

    // // preconditionner
    // KSPGetPC(ksp, &pc);
    // PCSetType(pc, PCGAMG);

    // prepare the rhs
    VecDuplicate(res, &rhs);
    VecCopy(res, rhs);
    VecScale(rhs, -1.0);

    // solve the linear system
    KSPSolve(ksp, rhs, sol);

    // get iteration count
    KSPGetIterationNumber(ksp, &its);

    KSPDestroy(&ksp);
    MatDestroy(&op);
    VecDestroy(&rhs);

    auto tstop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(tstop - tstart);

    this->print_time = print_time_local;

    if(this->print_time){
        PetscPrintf(PETSC_COMM_WORLD, "  -- solveCorrectionEquation done in %d mu s [%d its]\n", duration.count(),its);
    }

    return(0);

}

PetscErrorCode JacobiDavidsonMaxwellSolver::computeResidueVector(const Vec &q, const PetscReal rho, Vec &res){

    // compute the residue res = A q - rho M q
    // here M is always the identity matrix
    PetscInt        ncols, nrows;

    auto tstart = std::chrono::high_resolution_clock::now();

    // get sizes
    MatGetSize(this->A, &nrows, &ncols);

    // compute res  = Aq
    VecSet(res, 0.0);
    MatMult(this->A, q, res);

    // assemble res = A q - rho M q
    VecAXPY(res, -1.0*rho, q);

    auto tstop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(tstop - tstart);

    if(this->print_time){
        PetscPrintf(PETSC_COMM_WORLD, "  -- computeResidueVector done in %d mu s \n", duration.count());
    }

    return(0);

}


PetscErrorCode JacobiDavidsonMaxwellSolver::computeProjection(Vec &v, const BV &Q)
{
    // Compute v = v - Q Q.T v
    PetscScalar *tmp_row_vector;
    PetscInt    local_size, ncol, nrow;

    auto tstart = std::chrono::high_resolution_clock::now();

    // get sizes
    BVGetSizes(Q, &local_size, &nrow, &ncol);


    // create tmp_row_vector
    tmp_row_vector = new PetscScalar [ncol];

    // compute tmp_row_vect = Q.T v
    BVDotVec(Q, v, tmp_row_vector);

    // v = v - Q Q.T v
    BVMultVec(Q, -1.0, 1.0, v, tmp_row_vector);

    auto tstop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(tstop - tstart);

    if(this->print_time){
        PetscPrintf(PETSC_COMM_WORLD, "  -- computeProjection done in %d mu s \n", duration.count());
    }

    return (0);

}

PetscErrorCode JacobiDavidsonMaxwellSolver::projectCorrectionVector(Vec &corr)
{

     // compute corr = corr - Y H^{-1} C.T corr

    Vec         rhs, ksp_sol, tmp_vec;
    PetscInt    c_ncols, c_nrows, c_local_ncols, c_local_nrows; 
    PetscInt    corr_size, rhs_size;
    KSP         ksp;
    PetscInt    its;

    auto tstart = std::chrono::high_resolution_clock::now();

    MatGetSize(this->C, &c_nrows, &c_ncols);
    MatGetLocalSize(this->C, &c_local_nrows, &c_local_ncols);

    // create tmp_row vec
    VecCreate(PETSC_COMM_WORLD, &rhs);
    VecSetSizes(rhs,c_local_ncols,c_ncols);
    VecSetFromOptions(rhs);


    // compute rhs = C.T corr
    MatMultHermitianTranspose(this->C, corr, rhs);
    
    // set up the linear system
    KSPCreate(PETSC_COMM_WORLD,&ksp);
    KSPSetType(ksp,KSPGMRES);
    KSPSetOperators(ksp,this->H,this->H);
    KSPSetFromOptions(ksp);
    KSPSetTolerances(ksp, 1e-12, 1.e-12, PETSC_DEFAULT, 100);

    // solve the linear system
    //  H ksp_sol = (C.T corr)
    // -> ksp_sol = H^{-1} C.T corr
    VecDuplicate(rhs, &ksp_sol);
    KSPSolve(ksp, rhs, ksp_sol);
    KSPGetIterationNumber(ksp, &its);

    // compute tmp_vec = Y H^{-1} C.T corr
    VecDuplicate(corr, &tmp_vec);
    MatMult(this->Y, ksp_sol, tmp_vec);

    // compute corr - Y H^{-1} C.T corr
    VecAXPY(corr, -1.0, tmp_vec);
    
    VecDestroy(&rhs);
    VecDestroy(&tmp_vec);
    KSPDestroy(&ksp);

    auto tstop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(tstop - tstart);

    if(this->print_time){
        PetscPrintf(PETSC_COMM_WORLD, "  -- projectCorrectionVector done in %d mu s [%d its]\n", duration.count(), its);
    }

    return(0);
}



PetscErrorCode JacobiDavidsonMaxwellSolver::computeSmallEigenvalues(std::vector<PetscReal> &eigenvalues, Vec *eigenvectors)
{
    // project the total matrix following
    // Ak = V.T A V
    // and computes the eigenvalues of Ak using LAPACK

    Mat             Ak, AkT;
    EPS             eps;
    PetscErrorCode  ierr;
    PetscInt        nconv;
    PetscScalar     eigr, eigi;
    Vec             Vr, Vi;
    

    PetscInt        sort_idx[this->V_current_size];
    PetscScalar     eval_tmp[this->V_current_size];
    PetscReal       eval_tmp_real[this->V_current_size];
    Vec             evec_tmp[this->V_current_size];

    PetscInt        i;

    auto tstart = std::chrono::high_resolution_clock::now();

    MatCreateSeqDense(PETSC_COMM_SELF,
                      this->V_current_size, 
                      this->V_current_size,
                      NULL, &Ak);
    MatSetOption(Ak, MAT_HERMITIAN, PETSC_TRUE);

    BVMatProject(this->V, this->A, this->V, Ak);

    // make sur it's hermitian otherwise LAPACK crashes
    MatCreateHermitianTranspose(Ak, &AkT);
    MatAXPY(Ak, 1.0, AkT, SAME_NONZERO_PATTERN);
    MatScale(Ak, 0.5);

    // create evects
    MatCreateVecs(Ak, NULL, &Vr);
    MatCreateVecs(Ak, NULL, &Vi);

    // create eps object
    ierr = EPSCreate(PETSC_COMM_SELF,&eps);CHKERRQ(ierr);
    ierr = EPSSetOperators(eps,Ak,NULL);CHKERRQ(ierr);
    ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);
    ierr = EPSSetType(eps, EPSLAPACK);
    ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);
    ierr = EPSSolve(eps);CHKERRQ(ierr);
    EPSGetConverged(eps,&nconv);
    
    for (i=0; i<nconv; i++){
        ierr = EPSGetEigenpair(eps,  i, &eigr, &eigi, Vr, Vi);CHKERRQ(ierr);
        eval_tmp[i] = eigr;
        sort_idx[i] = i;
        ierr = MatCreateVecs(Ak, NULL, &evec_tmp[i]);CHKERRQ(ierr);
        ierr = VecCopy(Vr, evec_tmp[i]);CHKERRQ(ierr);
    }

    // extract the real part of the eigenvalues
    for(i=0; i<this->V_current_size; i++){
        eval_tmp_real[i] = PetscRealPart(eval_tmp[i]);
    }

    // get sort indexes eigenvalues
    PetscSortRealWithArrayInt(this->V_current_size, eval_tmp_real, sort_idx);

    // sort eigenvalues
    eigenvalues.clear();
    for (i=0; i<nconv; i++){

        if(eval_tmp_real[i] < 1E-6){
            PetscPrintf(PETSC_COMM_WORLD, " Warning : Zero Eigenvalue detected: %d %f\n", i, eval_tmp_real[i]);
            // return(0);
        }
        
        eigenvalues.push_back(eval_tmp_real[i]);
        ierr = MatCreateVecs(Ak, NULL, &eigenvectors[i]);CHKERRQ(ierr);
        ierr = VecCopy(evec_tmp[sort_idx[i]], eigenvectors[i]);CHKERRQ(ierr);

    }

    // clean up
    EPSDestroy(&eps);
    MatDestroy(&Ak);
    MatDestroy(&AkT);
    VecDestroy(&Vr);
    VecDestroy(&Vi);
    
    for(i=0; i<nconv; i++){
        VecDestroy(&evec_tmp[i]);
    }

    auto tstop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(tstop - tstart);

    if (this->print_small_evs){
        for (int ii=0; ii<this->V_current_size;ii++)
        {
            PetscPrintf(PETSC_COMM_WORLD, "eigenvalue: %6.12f\n", eigenvalues[ii]);
            // VecView(small_evects[ii], PETSC_VIEWER_STDOUT_WORLD);
        }
    }

    if(this->print_time){
        PetscPrintf(PETSC_COMM_WORLD, "  -- computeSmallEigenvalues done in %d mu s \n", duration.count());
    }

    return(0);
}

PetscErrorCode JacobiDavidsonMaxwellSolver::computeThreshold(Vec q, Vec r, PetscReal *eps)
{
    PetscScalar num, denom;

    // compute q.T M q 
    // with M = I
    VecDot(q,q, &denom);

    // compute r.T M^{-1} r
    // with M^{-1} = I
    VecDot(r,r, &num);

    // compute the threshold
    *eps = sqrt(PetscRealPart(num) / PetscRealPart(denom));

    return(0);
}

PetscErrorCode JacobiDavidsonMaxwellSolver::solve(PetscInt nev)
{

    PetscReal       rho, eps;
    std::vector<PetscReal> small_evals;

    Vec             small_evects[this->search_space_maxsize];
    Vec             tmp_v[this->search_space_maxsize];
    Vec             residue_vect_copy, correction_vect;
    Vec             q;
    PetscInt        k, idx_evect, ii;
    PetscBool       found, stop;
    PetscScalar       *vec_values;
    PetscErrorCode  ierr;

    this->search_space_minsize = nev+1;
    this->nconverged = 0;


    if(this->search_space_minsize >= search_space_maxsize){
        PetscPrintf(PETSC_COMM_WORLD, " Too many eigenvalues required");
        exit(0);
    }

    // init the solver matrices and vectors
    
    initializeMatrices();

    initializeVectors();
    
    initializeSearchSpace(nev);
    
    // rayleigh quotient
    computeRayleighQuotient(this->search_vect, &rho);
    
    // compute initial residue vector
    computeResidueVector(this->search_vect, rho, this->residue_vect);

    // duplicate residue vector and create the correction vector
    ierr = VecDuplicate(this->residue_vect, &residue_vect_copy); CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = VecDuplicate(this->residue_vect, &correction_vect); CHKERRABORT(PETSC_COMM_WORLD, ierr);

    k = 0;
    stop = PETSC_FALSE;


    for (this->iter=0; this->iter<this->maxIter; this->iter++)
    {

        // PetscPrintf(PETSC_COMM_WORLD, "== iteration : %d, k = %d, rho = %f\n", iter, k, rho);
        logger(INFO, "JacobiDavidsonSolver iteration : %, k = %, rho = %", iter, k, rho);

        // determine eta
        eps = (PetscReal)(std::rand())/ (PetscReal)(RAND_MAX);
        this->eta = (this->iter>0 && eps<0.5) ? rho : this->tau;

        // FOR now keep eta = tau
        // this->eta = this->tau;

        // copy res to res_new and left project
        ierr = VecCopy(this->residue_vect, residue_vect_copy); CHKERRABORT(PETSC_COMM_WORLD, ierr);
        computeProjection(residue_vect_copy, this->Qt);

        // solve correction equation
        solveCorrectionEquation(residue_vect_copy, correction_vect);

        // ortho correction vector
        computeProjection(correction_vect, this->Qt);

        // project correction vector
        projectCorrectionVector(correction_vect);


        // create search vector
        ierr = VecCopy(correction_vect, search_vect); CHKERRABORT(PETSC_COMM_WORLD, ierr);
        computeProjection(search_vect, this->V);
        normalizeVector(search_vect);


        k++;

        // Add column to V
        // VecView(search_vect, PETSC_VIEWER_STDOUT_SELF);
        ierr = BVInsertVec(this->V, this->V_current_size, search_vect); CHKERRABORT(PETSC_COMM_WORLD, ierr);
        this->V_current_size ++;
        ierr = BVSetActiveColumns(this->V, 0, this->V_current_size); CHKERRABORT(PETSC_COMM_WORLD, ierr);


        // compute eigenpairs on small subspace
        ierr = computeSmallEigenvalues(small_evals, small_evects); CHKERRABORT(PETSC_COMM_WORLD, ierr);

        // create q
        ierr = VecDuplicate(this->residue_vect, &q); CHKERRABORT(PETSC_COMM_WORLD, ierr);
       
        found = PETSC_TRUE;
        idx_evect = 0;
        while(found)
        {

            // rho = Lk[0]
            rho = small_evals[0];

            // q = V sk0
            ierr = VecGetArray(small_evects[0], &vec_values);  CHKERRABORT(PETSC_COMM_WORLD, ierr);
            ierr = BVMultVec(this->V, 1.0, 0.0, q, vec_values);  CHKERRABORT(PETSC_COMM_WORLD, ierr);
            ierr = VecRestoreArray(small_evects[0], &vec_values);  CHKERRABORT(PETSC_COMM_WORLD, ierr);

            // r = A q - rho M q
            computeResidueVector(q, rho, this->residue_vect);

            // compute residue threshold
            computeThreshold(q, this->residue_vect, &eps);


            // compute covnergence criteria
            found = (eps < this->tolerance && k>0) ? PETSC_TRUE : PETSC_FALSE;
            if(found)
            {
                logger(INFO, "JacobiDavidsonSolver new eigenvector found with eigenvalue %", rho);
                // store the eigenvalue/eigenvector
                this->eigenvalues.push_back(rho);


                // add a new vector to the basis
                ierr = BVInsertVec(Q, Q_current_size, q);  CHKERRABORT(PETSC_COMM_WORLD, ierr);
                this->Q_current_size++;

                if (this->Q_current_size >= nev) {
                    stop = PETSC_TRUE;
                    break;
                }

                // Compute [tmp_v1 .... tmp_vk-1] = V [s2, .... sk]
                for (ii=0; ii<k-1;ii++){
                    ierr = VecGetArray(small_evects[ii+1], &vec_values);  CHKERRABORT(PETSC_COMM_WORLD, ierr);
                    ierr = VecDuplicate(this->residue_vect, &tmp_v[ii]);  CHKERRABORT(PETSC_COMM_WORLD, ierr);
                    ierr = VecSet(tmp_v[ii], 0.0);  CHKERRABORT(PETSC_COMM_WORLD, ierr);
                    ierr = BVMultVec(this->V, 1.0, 0.0, tmp_v[ii], vec_values);  CHKERRABORT(PETSC_COMM_WORLD, ierr);
                    ierr = VecRestoreArray(small_evects[ii+1], &vec_values);  CHKERRABORT(PETSC_COMM_WORLD, ierr);
                }

                // insert [tmp_v1 .... tmp_vk-1] in V
                for (ii=0; ii<k-1;ii++){
                    ierr = BVInsertVec(this->V, ii, tmp_v[ii]); CHKERRABORT(PETSC_COMM_WORLD, ierr);
                    ierr = VecDestroy(&tmp_v[ii]); CHKERRABORT(PETSC_COMM_WORLD, ierr);
                }

                // set all other columns to 0
                for (ii=k-1;ii<this->search_space_maxsize;ii++){
                    ierr = BVScaleColumn(this->V, ii, 0.0); CHKERRABORT(PETSC_COMM_WORLD, ierr);
                }

                // number of active cols of V
                ierr = BVSetActiveColumns(this->V, 0, k-1); CHKERRABORT(PETSC_COMM_WORLD, ierr);
                this->V_current_size = k-1;


                // remove first eigenvalue
                small_evals.erase(small_evals.begin());

                // set the eigenvectors to I(k-1)
                for (ii=0; ii<k;ii++){
                    VecSet(small_evects[ii],0.0);
                    if(ii<k-1){
                        ierr = VecSetValue(small_evects[ii], ii , 1.0, INSERT_VALUES); CHKERRABORT(PETSC_COMM_WORLD, ierr);
                        ierr = VecAssemblyBegin(small_evects[ii]); CHKERRABORT(PETSC_COMM_WORLD, ierr);
                        ierr = VecAssemblyEnd(small_evects[ii]); CHKERRABORT(PETSC_COMM_WORLD, ierr);
                    }
                }
                k--;
                idx_evect ++;

            }

        }

        if(stop)
            break;

        // Qt = [Q; q]
        if (this->Q_current_size == 0) {
            ierr = BVInsertVec(this->Qt, 0, q); CHKERRABORT(PETSC_COMM_WORLD, ierr);
            this->Qt_current_size = 1;
            ierr = BVSetActiveColumns(this->Qt, 0, this->Qt_current_size); CHKERRABORT(PETSC_COMM_WORLD, ierr);
        }
        else {
            ierr = BVCopy(this->Q, this->Qt); CHKERRABORT(PETSC_COMM_WORLD, ierr);
            ierr = BVInsertVec(this->Qt, Q_current_size, q); CHKERRABORT(PETSC_COMM_WORLD, ierr);
            this->Qt_current_size = this->Q_current_size + 1;
            ierr = BVSetActiveColumns(this->Qt, 0, this->Qt_current_size); CHKERRABORT(PETSC_COMM_WORLD, ierr);
            for (ii=this->Qt_current_size;ii<this->search_space_maxsize;ii++){
                ierr = BVScaleColumn(this->Qt, ii, 0.0); CHKERRABORT(PETSC_COMM_WORLD, ierr);
            }
        }

        // restart
        if(k == this->search_space_maxsize-1) {

            // compute tmp = V Sk[:,:size_min]
            for (ii=0; ii<this->search_space_minsize;ii++){
                
                ierr = VecGetArray(small_evects[ii], &vec_values); CHKERRABORT(PETSC_COMM_WORLD, ierr);
                ierr = VecDuplicate(this->residue_vect, &tmp_v[ii]); CHKERRABORT(PETSC_COMM_WORLD, ierr);
                ierr = VecSet(tmp_v[ii], 0.0); CHKERRABORT(PETSC_COMM_WORLD, ierr);
                
                ierr = BVMultVec(this->V, 1.0, 0.0, tmp_v[ii], vec_values); CHKERRABORT(PETSC_COMM_WORLD, ierr);
                ierr = VecRestoreArray(small_evects[ii+1], &vec_values); CHKERRABORT(PETSC_COMM_WORLD, ierr);
            }

            // insert [tmp_v1 .... tmp_vk-1] in V
            for (ii=0; ii<search_space_minsize; ii++){
                ierr = BVInsertVec(this->V, ii, tmp_v[ii]); CHKERRABORT(PETSC_COMM_WORLD, ierr);
                ierr = VecDestroy(&tmp_v[ii]); CHKERRABORT(PETSC_COMM_WORLD, ierr);
            }


            // set all other columns to 0
            for (ii=search_space_minsize;ii<this->search_space_maxsize;ii++){
                ierr = BVScaleColumn(this->V, ii, 0.0); CHKERRABORT(PETSC_COMM_WORLD, ierr);
            }

            // number of active cols of V
            ierr = BVSetActiveColumns(this->V, 0, k-1); CHKERRABORT(PETSC_COMM_WORLD, ierr);
            this->V_current_size = search_space_minsize;

            k = this->search_space_minsize;

        }

    }

    return(0);
}

} //namespace LinearAlgebra

}  // namespace hpgem
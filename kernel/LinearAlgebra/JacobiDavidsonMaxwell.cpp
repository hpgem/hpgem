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

JacobiDavidsonMaxwellSolver::JacobiDavidsonMaxwellSolver(Mat &A, Mat &M, Mat &C)
    : A(A), M(M), C(C) {}


PetscBool JacobiDavidsonMaxwellSolver::check_vect_nan(const Vec &vec) {

    PetscBool isnan;
    PetscInt ii, low, high, local_size;
    PetscInt *idx;
    PetscScalar *vals;

    VecGetOwnershipRange(vec, &low, &high);
    local_size = high-low;
    vals = new PetscScalar [local_size];
    idx = new PetscInt [local_size];

    for(ii=0; ii<local_size; ii++){
        idx[ii] = low+ii;
    }

    VecGetValues(vec, local_size, idx, vals);

    isnan = PETSC_FALSE;
    for(ii=0;ii<local_size;ii++){
        if(PetscIsNanScalar(vals[ii]))
        {
            isnan = PETSC_TRUE;
            break;
        }
    }

    return(isnan);

}

PetscInt JacobiDavidsonMaxwellSolver::getConverged(){
    return this->Q_current_size;
}

PetscInt JacobiDavidsonMaxwellSolver::getIterationCount(){
    return this->iter;
}

PetscErrorCode JacobiDavidsonMaxwellSolver::getEigenPair(PetscInt index, 
                                                     PetscScalar &eval,
                                                     Vec &evec )
{
    PetscErrorCode ierr;
    Vec tmp_vec;
    eval = this->eigenvalues[index];
    ierr = BVGetColumn(this->Q, index, &tmp_vec);
    VecCopy(tmp_vec, evec);
    BVRestoreColumn(this->Q, index, &tmp_vec);
    VecDestroy(&tmp_vec);
    return (0);

}

void JacobiDavidsonMaxwellSolver::setMatrices(const Mat &Ain, const Mat &Min, const Mat &Cin) {
    
    PetscErrorCode ierr; 
    PetscInt n,m;
    PetscViewer viewer; 
    

    this->A = Ain;
    MatSetOption(this->A, MAT_HERMITIAN, PETSC_TRUE);
    this->C = Cin;

    // If config.useHermitian_ is true, M is identity and therefore invM as well.
    ierr = MatGetSize(Min, &n, &m);
    MatCreateConstantDiagonal(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, n,m, 1.0, &this->M);
    MatSetFromOptions(this->M);
    MatConvert(this->M, MATSEQAIJ, MAT_INITIAL_MATRIX, &this->M);

    ierr = MatGetSize(Ain, &n, &m);
    logger(INFO, "JacobiDavidsonSolver System Size : % x %", n,m);

}

void JacobiDavidsonMaxwellSolver::initializeMatrices() {

    PetscInt        y_nrows, y_ncols;
    PetscInt        m_nrows, m_ncols;
    Mat             Yh;
    PetscErrorCode  ierr;

    // transpose C
    ierr = MatHermitianTranspose(this->C, MAT_INPLACE_MATRIX, &this->C);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);

    // if we know invM we can directly use it to compute Y = M^{-1} C
    if(1) {

        ierr = MatDuplicate(C, MAT_COPY_VALUES, &this->Y);
        CHKERRABORT(PETSC_COMM_WORLD, ierr);
    }

    // if we don't know invM we have to solve M y = c 
    else {

        ierr = MatGetSize(this->C, &y_nrows, &y_ncols);
        CHKERRABORT(PETSC_COMM_WORLD, ierr);

        // create Y matrix
        MatCreate(PETSC_COMM_WORLD, &this->Y);
        MatSetSizes(this->Y,PETSC_DECIDE,PETSC_DECIDE,y_nrows,y_ncols);
        MatSetFromOptions(this->Y);
        MatSetUp(this->Y);

        // solve M * Y = C for each C vector
        solveMultipleLinearSystems();
    }


    // Compute H = Y.T M Y 
    // MatCreateHermitianTranspose(this->Y, &Yh);
    ierr = MatDuplicate(this->Y, MAT_COPY_VALUES, &Yh);
    ierr = MatHermitianTranspose(Yh, MAT_INPLACE_MATRIX, &Yh);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);


    // ierr = MatProductCreate(this->M, this->Y, NULL, &this->H);
    // ierr = MatProductCreate(this->Y, this->Y, NULL, &this->H);
    ierr = MatProductCreate(Yh, this->Y, NULL, &this->H);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);

    // ierr = MatProductSetType(this->H,MATPRODUCT_PtAP);
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

    // get sizes
    MatGetSize(this->A, &nrows, &ncols);

    // set the current size
    this->search_space_current_size = 1;

    // allocate the search space
    this->Qt_current_size = 1;
    BVCreate(PETSC_COMM_WORLD, &this->Qt);
    BVSetSizes(this->Qt, PETSC_DECIDE, nrows, this->search_space_maxsize);
    BVSetFromOptions(this->Qt);
    BVInsertVec(this->Qt, 0, this->search_vect);
    BVSetActiveColumns(this->Qt, 0, this->search_space_current_size);

    // allocate the search space
    this->V_current_size = 1;
    BVCreate(PETSC_COMM_WORLD, &this->V);
    BVSetSizes(this->V, PETSC_DECIDE, nrows, this->search_space_maxsize);
    BVSetFromOptions(this->V);
    BVInsertVec(this->V, 0, this->search_vect);
    BVSetActiveColumns(this->V, 0, this->V_current_size);

    // allocate the search space
    this->Q_current_size = 0;
    BVCreate(PETSC_COMM_WORLD, &this->Q);
    BVSetSizes(this->Q, PETSC_DECIDE, nrows, nev);
    BVSetFromOptions(this->Q);
    BVSetActiveColumns(this->Q, 0, this->Q_current_size);
}

void JacobiDavidsonMaxwellSolver::initializeVectors(){

    PetscRandom     rctx;
    PetscInt        y_nrows, y_ncols;

    MatGetSize(this->Y, &y_nrows, &y_ncols);

    // init the initial search vector
    PetscRandomCreate(PETSC_COMM_WORLD, &rctx);
    PetscRandomSetFromOptions(rctx);
    VecCreate(PETSC_COMM_WORLD, &this->search_vect);
    VecSetSizes(this->search_vect, PETSC_DECIDE, y_nrows);
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


PetscErrorCode JacobiDavidsonMaxwellSolver::solveMultipleLinearSystems(){
    KSP         ksp;
    PetscScalar *tmp_val;
    PetscInt    ncols, nrows, *idx_rows;
    Vec         rhs, sol;
    PetscInt    local_row_start, local_row_end, local_size;
    PetscInt    i;

    // get size of the matrix
    MatGetSize(this->C, &nrows, &ncols);

    // create the rhs vector
    VecCreate(PETSC_COMM_WORLD, &rhs);
    VecSetSizes(rhs,PETSC_DECIDE,nrows);
    VecSetFromOptions(rhs);

    // create the sol vector
    VecCreate(PETSC_COMM_WORLD, &sol);
    VecSetSizes(sol,PETSC_DECIDE,nrows);
    VecSetFromOptions(sol);

    // get the ownership
    VecGetOwnershipRange(sol, &local_row_start, &local_row_end);
    local_size = local_row_end - local_row_start;

    // create a range of index
    idx_rows = new PetscInt[local_size];
    for (i=0; i<local_size; i++)
        idx_rows[i] = local_row_start + i;

    // create the array
    tmp_val = new PetscScalar [local_size];

    // set up the linear system
    // MatView(M,PETSC_VIEWER_STDOUT_WORLD);
    KSPCreate(PETSC_COMM_WORLD,&ksp);
    KSPSetType(ksp,KSPCG);

    KSPSetOperators(ksp,M,M);
    KSPSetTolerances(ksp, 1e-6, 1.e-6, PETSC_DEFAULT, 100);
    KSPSetFromOptions(ksp);

    // solve each M * y = c
    for (i=0; i<ncols; i++){

        PetscPrintf(PETSC_COMM_WORLD, " Linear system [%dx%d] %d/%d\n", nrows, ncols, i, ncols);
        MatGetColumnVector(this->C, rhs, i);
        KSPSolve(ksp, rhs, sol);
        VecGetArray(sol, &tmp_val);
        MatSetValues(this->Y, local_size, idx_rows, 1 , &i, tmp_val, INSERT_VALUES);
    }
    MatAssemblyBegin(this->Y, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(this->Y, MAT_FINAL_ASSEMBLY);

    KSPDestroy(&ksp);
    VecDestroy(&rhs);
    VecDestroy(&sol);

    return(0);
}


PetscErrorCode JacobiDavidsonMaxwellSolver::computeRayleighQuotient(const Vec &x, PetscReal *out)
{
    PetscScalar       a;
    PetscScalar         m;

    VecxTAx(x, this->A, &a);
    // VecxTAx(x, this->M, &m);
    // VecNorm(x, NORM_2, &m);
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
    this->VecxTAx(x, this->M, &norm);
    // VecNorm(x, NORM_2, &norm);
    // VecDot(x, x, &norm);
    // VecView(x, PETSC_VIEWER_STDOUT_WORLD);  
    inv_sqrt_norm = 1.0/sqrt(PetscRealPart(norm));

    isnan = PetscIsNanReal(inv_sqrt_norm);
    if(PetscIsNanReal(inv_sqrt_norm)){
        logger(INFO, "Norm is Nan : %", norm);
        exit(0);
    }

    if(PetscIsInfReal(inv_sqrt_norm)){
        logger(INFO, "Norm is Nan : %", norm);
        exit(0);
    }

    isinf = PetscIsInfReal(inv_sqrt_norm);

    if(isnan || isinf){
        PetscPrintf(PETSC_COMM_WORLD, " Warning : Nan detected in the norm\n");
        VecSet(x, 0.0);
    } else {
        // scale the vector
        VecScale(x, inv_sqrt_norm);
    }


    logger.assert_debug(!check_vect_nan(x),
                        "x is nans in normalization");
    
    //  isnan = check_vect_nan(x);
    // if(isnan){
    //     PetscPrintf(PETSC_COMM_WORLD, "   ==> Error : x is nans in normalization\n");
    //     PetscPrintf(PETSC_COMM_WORLD, "inv_sqrt_norm = %f", inv_sqrt_norm);
    //     exit(0);}

    return(0);
}

PetscErrorCode JacobiDavidsonMaxwellSolver::VecxTAx(const Vec &x, const Mat &K, PetscScalar *val)
{
    // Compute val = x.T K X
    Vec             tmp;
    PetscInt        ncols, nrows;

    // get sizes
    MatGetSize(K, &nrows, &ncols);

    // create tmp = K * x
    VecCreate(PETSC_COMM_WORLD, &tmp);
    VecSetSizes(tmp,PETSC_DECIDE,nrows);
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

    MatDuplicate(this->A, MAT_COPY_VALUES, &ApM);
    // MatAXPY(ApM, 1.0, this->M, UNKNOWN_NONZERO_PATTERN);
    MatShift(ApM, -1.0);

    KSPCreate(PETSC_COMM_WORLD,&ksp);
    KSPSetType(ksp,KSPGMRES);
    KSPSetOperators(ksp,ApM,ApM);
    // KSPSetOperators(ksp,this->A, this->A);
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

    Vec             Ay, My;

    // copy x in y
    VecCopy(x, y);

    // y = y - Q Q.T M y
    computeRightProjection(y, this->Qt);


    // tmp = A y
    VecDuplicate(x, &Ay);
    MatMult(this->A, y, Ay);


    // tmp2 = M y
    // VecDuplicate(x, &My);
    // MatMult(this->M, y, My);

    // tmp = tmp - eta tmp2
    //     = A x - eta M x
    // VecAXPY(Ay, -1.0*this->eta, My);
    VecAXPY(Ay, -1.0*this->eta, y);

    // // y = y - M Q Q.T y
    computeLeftProjection(Ay);
    VecCopy(Ay, y);

    VecDestroy(&Ay);
    // VecDestroy(&My);

    return(0);

}


PetscErrorCode JacobiDavidsonMaxwellSolver::solveCorrectionEquation(const Vec &res, Vec &sol)
{

    Mat         op;
    Mat         prec;

    KSP         ksp;
    PC          pc;
    PetscInt    ncols, nrows;
    Vec         rhs;
    PetscErrorCode ierr;
    PetscInt    its;
    bool print_time_local;

    auto tstart = std::chrono::high_resolution_clock::now();

    print_time_local = this->print_time;
    this->print_time = false;

    // get the linear operator
    MatGetLocalSize(this->A, &nrows, &ncols);
    MatCreateShell(PETSC_COMM_WORLD, nrows, ncols, PETSC_DECIDE,PETSC_DECIDE, this, &op);
    MatShellSetOperation(op, MATOP_MULT, (void (*)(void))staticMatMultCorrOp);


    MatCreateShell(PETSC_COMM_WORLD, nrows, ncols, PETSC_DECIDE,PETSC_DECIDE, this, &prec);
    MatShellSetOperation(prec, MATOP_MULT, (void (*)(void))staticMatMultCorrPrec);


    // set up the linear system
    KSPCreate(PETSC_COMM_WORLD,&ksp);
    KSPSetType(ksp,KSPMINRES);
    KSPSetOperators(ksp,op,op);
    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
    KSPSetTolerances(ksp, 1e-12, 1.e-12, PETSC_DEFAULT, 10);

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
    // Vec             tmp;
    PetscInt        ncols, nrows;

    auto tstart = std::chrono::high_resolution_clock::now();

    // get sizes
    MatGetSize(this->A, &nrows, &ncols);

    // create tmp vec
    // VecCreate(PETSC_COMM_WORLD, &tmp);
    // VecSetSizes(tmp,PETSC_DECIDE,nrows);
    // VecSetFromOptions(tmp);

    // compute res  = Aq
    VecSet(res, 0.0);
    MatMult(this->A, q, res);

    // compute tmp  = M q
    // MatMult(this->M, q, tmp);
    // VecDuplicate(q, &tmp);
    // VecCopy(q, tmp);

    // assemble res = A q - rho M q
    VecAXPY(res, -1.0*rho, q);

    // clean up
    // VecDestroy(&tmp);

    auto tstop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(tstop - tstart);

    if(this->print_time){
        PetscPrintf(PETSC_COMM_WORLD, "  -- computeResidueVector done in %d mu s \n", duration.count());
    }

    return(0);

}

PetscErrorCode JacobiDavidsonMaxwellSolver::computeRightProjection(Vec &v, const BV &Q)
{
    // Compute v = v - Q Q.T M v
    Vec         tmp_column_vector;
    PetscScalar *tmp_row_vector;
    PetscInt    local_size, ncol, nrow;

    auto tstart = std::chrono::high_resolution_clock::now();

    // get sizes
    BVGetSizes(Q, &local_size, &nrow, &ncol);

    // create tmp_column_vector
    // VecCreate(PETSC_COMM_WORLD, &tmp_column_vector);
    // VecSetSizes(tmp_column_vector,PETSC_DECIDE,nrow);
    // VecSetFromOptions(tmp_column_vector);

    // create tmp_row_vector
    tmp_row_vector = new PetscScalar [ncol];

    // compute tmp_col_vect = M v
    // MatMult(this->M, v, tmp_column_vector);
    // VecDuplicate(v, &tmp_column_vector);
    // VecCopy(v, tmp_column_vector);

    // compute tmp_row_vect = Q.T M v
    // BVDotVec(Q, tmp_column_vector, tmp_row_vector);
    BVDotVec(Q, v, tmp_row_vector);

    // v = v - Q Q.T M v
    BVMultVec(Q, -1.0, 1.0, v, tmp_row_vector);

    // VecDestroy(&tmp_row_vector);
    // VecDestroy(&tmp_column_vector);

    auto tstop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(tstop - tstart);

    if(this->print_time){
        PetscPrintf(PETSC_COMM_WORLD, "  -- computeRightProjection done in %d mu s \n", duration.count());
    }

    return (0);

}


PetscErrorCode JacobiDavidsonMaxwellSolver::computeLeftProjection(Vec &v)
{
    // Compute v = v - M Q Q.T v
    // Vec         tmp_row_vector;
    PetscScalar   *tmp_row_vector;
    PetscInt    local_size, ncol, nrow;
    // BV          MQ;


    auto tstart = std::chrono::high_resolution_clock::now();

    // get sizes
    BVGetSizes(this->Qt, &local_size, &nrow, &ncol);


    // create tmp_row_vector
    tmp_row_vector = new PetscScalar [this->Qt_current_size];

    // compute tmp_row_vect = Q.T v
    BVDotVec(this->Qt, v, tmp_row_vector);

    // // compute MQ = M Q
    // BVDuplicate(this->Qt, &MQ);
    // BVSetActiveColumns(MQ, 0, this->Qt_current_size);
    // BVMatMult(this->Qt, this->M, MQ);
    // BVSetActiveColumns(MQ, 0, this->Qt_current_size);

    // v = v - M Q Q.T v
    // BVMultVec(MQ, -1.0, 1.0, v, tmp_row_vector);
    BVMultVec(this->Qt, -1.0, 1.0, v, tmp_row_vector);

    // BVDestroy(&MQ);

    auto tstop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(tstop - tstart);

    if(this->print_time){
        PetscPrintf(PETSC_COMM_WORLD, "  -- computeLeftProjection done in %d mu s \n", duration.count());
    }
    return (0);

}

PetscErrorCode JacobiDavidsonMaxwellSolver::VecxTinvAx(const Vec &x, const Mat &K, PetscScalar *val)
{
    // Compute val = x.T K^{-1} x

    Vec             ksp_sol;
    KSP             ksp;


    auto tstart = std::chrono::high_resolution_clock::now();

    // set up the linear system
    KSPCreate(PETSC_COMM_WORLD,&ksp);
    KSPSetType(ksp,KSPGMRES);
    KSPSetOperators(ksp,K,K);
    KSPSetFromOptions(ksp);
    KSPSetTolerances(ksp, 1e-6, 1.e-12, PETSC_DEFAULT, 100);

    // solve the linear system
    // compute H^{-1} C.T corr
    VecDuplicate(x, &ksp_sol);
    KSPSolve(ksp, x, ksp_sol);

    // compute the norm
    VecDot(x, ksp_sol, val);

    // clean up
    KSPDestroy(&ksp);
    VecDestroy(&ksp_sol);

    auto tstop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(tstop - tstart);

    if(this->print_time){
        PetscPrintf(PETSC_COMM_WORLD, "  -- X.T A^{-1} X  done in %d mu s \n", duration.count());
    }

    return(0);
}

PetscErrorCode JacobiDavidsonMaxwellSolver::projectCorrectionVector(Vec &corr)
{

     // compute corr - Y H^{-1} C.T corr

    Vec         rhs, ksp_sol, tmp_vec;
    PetscInt    c_ncols, c_nrows;
    KSP         ksp;
    PetscInt    its;

    auto tstart = std::chrono::high_resolution_clock::now();

    MatGetSize(this->C, &c_nrows, &c_ncols);


    // create tmp_row vec
    VecCreate(PETSC_COMM_WORLD, &rhs);
    VecSetSizes(rhs,PETSC_DECIDE,c_ncols);
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
    Mat             Ak, AkT;
    EPS             eps;
    PetscErrorCode  ierr;
    PetscInt        nconv, nc, nr;
    PetscScalar     eigr, eigi;
    Vec             Vr, Vi;
    PetscBool       flg;

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

    // MatView(Ak, PETSC_VIEWER_STDOUT_SELF);
    // exit(0);

    MatCreateHermitianTranspose(Ak, &AkT);
    MatAXPY(Ak, 1.0, AkT, SAME_NONZERO_PATTERN);
    MatScale(Ak, 0.5);

    MatGetSize(Ak, &nr, &nc);
    MatIsHermitian(Ak,1E-12, &flg);

    if(flg == PETSC_FALSE){
        // std::cout << "Small Matrix not Hermitian" << std::endl;
        PetscPrintf(PETSC_COMM_WORLD, " Error : Small matrix is not Hermitian\n");
        exit(0);
    }

    MatCreateVecs(Ak, NULL, &Vr);
    MatCreateVecs(Ak, NULL, &Vi);

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
    // VecxTAx(q, this->M, &denom);
    VecDot(q,q, &denom);

    // compute r.T M^{-1} r
    // VecxTinvAx(r, this->M, &num);
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
    Vec             q, tmp_q;
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
        computeLeftProjection(residue_vect_copy);

        // solve correction equation
        solveCorrectionEquation(residue_vect_copy, correction_vect);

        // ortho correction vector
        computeRightProjection(correction_vect, this->Qt);

        // project correction vector
        projectCorrectionVector(correction_vect);


        // create search vector
        ierr = VecCopy(correction_vect, search_vect); CHKERRABORT(PETSC_COMM_WORLD, ierr);
        computeRightProjection(search_vect, this->V);
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
        ierr = VecDuplicate(this->residue_vect, &tmp_q); CHKERRABORT(PETSC_COMM_WORLD, ierr);

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
            found = (eps < 1E-3 && k>0) ? PETSC_TRUE : PETSC_FALSE;
            if(found)
            {
                logger(INFO, "JacobiDavidsonSolver new eigenvector found with eigenvalue %", rho);
                // store the eigenvalue/eigenvector
                this->eigenvalues.push_back(rho);

                // ierr = MatCreateVecs(Ak, NULL, &eigenvectors[this->nconverged]);CHKERRQ(ierr);
                // ierr = VecDuplicate(q, &eigenvectors[this->nconverged]); CHKERRQ(ierr);
                // ierr = VecCopy(q, this->eigenvectors[this->nconverged]);CHKERRQ(ierr);
                // this->nconverged ++;

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
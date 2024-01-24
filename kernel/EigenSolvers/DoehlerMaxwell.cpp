#include <numeric>
#include <iostream>
#include <iomanip>
#include <chrono>

#include "petsc.h"
#include "slepc.h"

#include "DoehlerMaxwell.h"
#include "Logger.h"

namespace hpgem {

namespace EigenSolvers {

DoehlerMaxwellSolver::DoehlerMaxwellSolver() : eigenvectors(nullptr) {
    std::cout << "Initialize DoehlerMaxwellSolver!!!!!" << std::endl;
}

DoehlerMaxwellSolver::~DoehlerMaxwellSolver() { BVDestroy(&eigenvectors); }

void DoehlerMaxwellSolver::setMatrices(const Mat Ain, const Mat Cin) {
    // Set the A and C matrices describing the eigenvalue problem
    // Ax = lambda x with C x = 0

    // Set the A, and C, but set M to be the identity matrix (NULL is a
    // shorthand)
    this->setMatrices(Ain, Cin, NULL);
}

void DoehlerMaxwellSolver::setMatrices(const Mat Ain, const Mat Cin,
                                       const Mat Min) {
    // Set the A, C, and M matrices describing the eigenvalue problem
    // Ax = lambda M x with C x = 0.
    // If Min = NULL then M is assumed to be the identity matrix.
    this->A = Ain;
    this->C = Cin;
    this->M = Min;
}

PetscErrorCode DoehlerMaxwellSolver::solve(PetscInt nev, Mat &T_Mat_in,
                                           PetscInt n_steps_projection) {

    // Initialize the matrices for enforcing the condition this->C * x = 0
    // This is done by solving a Poisson-like problem, the matrices required for
    // solving this problem are initialized here
    this->setupProjection();

    PetscInt n_eigs = nev;  // this is done to keep the function parameter with
                            // the same names across implementations, but to
                            // keep the clearer name n_eigs internally
    PetscErrorCode err;
    // Get the rank of the process for display purposes
    PetscMPIInt rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // Show info that solver started
    if (rank == 0) {
        std::cout << "\n*******************************************************"
                     "*******"
                  << std::endl;
        std::cout << "* Doehler eingenvalue solver (PETSc) START" << std::endl;
        std::cout << "*********************************************************"
                     "*****\n"
                  << std::endl;
    }

    // Initialize the algorithm with initial guesses

    // Compute initial parameters, like system size, etc
    PetscInt A_n_rows, A_n_cols;  // number of rows and columns of matrix A
    PetscInt A_n_local_rows;
    MatGetSize(this->A, &A_n_rows, &A_n_cols);
    MatGetLocalSize(this->A, &A_n_local_rows, nullptr);

    if (rank == 0) {
        std::cout << "System size:" << std::endl;
        std::cout << "   n_rows: " << A_n_rows << std::endl;
        std::cout << "   n_cols: " << A_n_cols << std::endl;
        std::cout << "   n_eigs: " << n_eigs << std::endl;
    }

    // Initialize the eigenvector solution
    if (this->eigenvectors == nullptr) {
        BVCreate(PETSC_COMM_WORLD, &this->eigenvectors);
    }
    err = BVSetSizes(this->eigenvectors, A_n_local_rows, A_n_rows, n_eigs);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    BVSetFromOptions(this->eigenvectors);
    std::vector<PetscScalar> ritzValues(2 * n_eigs);

    // Initialize tranformation matrix T as a bv system (used to project to the
    // reduced space)
    BV searchSpace, tempBV;
    BVCreate(PETSC_COMM_WORLD, &searchSpace);

    BVSetSizes(searchSpace, A_n_local_rows, A_n_rows, 2 * n_eigs);
    BVSetFromOptions(searchSpace);

    // Set the initial (guess) values
    // (the n_eigs eigenvectors we wish to find and the
    // n_eigs vector search directions)
    PetscRandom random_context(nullptr);
    PetscRandomCreate(PETSC_COMM_WORLD, &random_context);
    PetscRandomSetInterval(random_context, 0.0, 1.0);
    PetscRandomSetFromOptions(random_context);
    BVSetRandomContext(searchSpace, random_context);
    BVSetRandom(searchSpace);

    // Apply the projector to ensure X and S satisfy the divergence free
    // constraint C S = C X = 0
    this->projectBV(searchSpace);
    // Intermediate storages
    BV largeRitzVectors;  // the search space
    BVCreate(PETSC_COMM_WORLD, &largeRitzVectors);
    BVSetSizes(largeRitzVectors, A_n_local_rows, A_n_rows, n_eigs);
    BVSetFromOptions(largeRitzVectors);
    BVDuplicate(largeRitzVectors, &tempBV);

    BV residuals;  // the residual column vectors
    BVCreate(PETSC_COMM_WORLD, &residuals);
    BVSetSizes(residuals, A_n_local_rows, A_n_rows, n_eigs);
    BVSetFromOptions(residuals);

    BV doubleResiduals;  // the BV containing the residual of the residual BV
    BVCreate(PETSC_COMM_WORLD, &doubleResiduals);
    BVSetSizes(doubleResiduals, A_n_local_rows, A_n_rows, n_eigs);
    BVSetFromOptions(doubleResiduals);

    // Temporary storage
    std::vector<PetscScalar> values(n_eigs * n_eigs);
    std::vector<PetscInt> indices(n_eigs);
    std::iota(indices.begin(), indices.end(), 0);

    // Iterate to find corrected solutions to the eigenvalue
    this->iter = 0;
    PetscReal error_max =
        1.0;  // initialize error max to determine if the loop is over or not
    do {
        this->iter++;  // update counter for number of iterations performed
        this->ritzUpdate(searchSpace, n_eigs, ritzValues);

        BVSetActiveColumns(
            searchSpace, n_eigs,
            2 * n_eigs);  // activate the columns associated to the search space
        BVCopy(searchSpace, largeRitzVectors);
        BVSetActiveColumns(searchSpace, 0,
                           2 * n_eigs);  // always return to original state

        // Compute the residual with the updated eigenvectors
        this->computeResiduals(ritzValues, this->eigenvectors, 0, n_eigs,
                               residuals, tempBV);

        // Compute convergence
        PetscReal residualNorm, evNorm;
        this->eigenvectors_current_size = 0;
        std::stringstream residual_values;
        for (PetscInt i = 0; i < n_eigs; ++i) {
            BVNormColumn(residuals, i, NORM_2, &residualNorm);
            BVNormColumn(this->eigenvectors, i, NORM_2, &evNorm);
            residual_values << " " << std::setprecision(5) << PetscRealPart(ritzValues[i])
                            << "(" << std::setprecision(2)
                            << (residualNorm / evNorm) << ")";
            if (this->eigenvectors_current_size != i) {
                // Only accept an eigenvalue as converged if all previous
                // eigenvalues also have converged. This ensures that we get the
                // wanted eigenvalues.
                continue;
            }
            // Check for convergence either in either relative or absolute sense
            if (residualNorm < evNorm * this->tolerance ||
                residualNorm <
                    evNorm * this->tolerance * std::abs(ritzValues[i])) {
                this->eigenvectors_current_size++;
            }
        }
        if (this->iter % 5 == 0 && rank == 0) {
            std::string log = residual_values.str();
            logger(INFO, "iter % converged %:%", this->iter,
                   this->eigenvectors_current_size, log);
        }
        if (this->eigenvectors_current_size == n_eigs) {
            // TODO: Improve
            error_max = 0.0;
        }
        MPI_Barrier(PETSC_COMM_WORLD);

        // Compute the new augmented solution space (the correction space) and
        // the new search space
        this->computeResiduals(ritzValues, residuals, 0, n_eigs,
                               doubleResiduals, tempBV);

        //    // T_{ij} = -v_j^H(A - lm_j B) w_i / (lm_{i+p} - lm_j)
        //    // -v_j^H(A - lm_j B)w = R(v)^H W
        Mat tmatrix;  // Make it
        MatCreateSeqDense(PETSC_COMM_SELF, n_eigs, n_eigs, NULL, &tmatrix);
        err = BVDot(doubleResiduals, largeRitzVectors, tmatrix);
        CHKERRABORT(PETSC_COMM_WORLD, err);
        // Use this matrix later on in the BVMult(R, 1, 1, W_r_bv, out);
        PetscScalar* tdata;
        MatDenseGetArray(tmatrix, &tdata);
        for (std::size_t col = 0; col < n_eigs; col++) {
            for (std::size_t row = 0; row < n_eigs; row++) {
                tdata[row + col  * n_eigs] /=
                    -(ritzValues[row + n_eigs] - std::conj(ritzValues[col]));
            }
        }
        MatDenseRestoreArray(tmatrix, &tdata);

        BVMult(residuals, 1.0, 1.0, largeRitzVectors, tmatrix);
        MatDestroy(&tmatrix);

        // Restart T_bv

        // Update X part not needed, this is already done in ritzUpdate

        // Update S part
        BVSetActiveColumns(searchSpace, n_eigs, 2*n_eigs);
        BVCopy(residuals, searchSpace);
        BVSetActiveColumns(searchSpace, 0, 2*n_eigs);
        // If nothing is done, the search directions result in vectors with very
        // small norms. This leads to very badly conditioned reduced matrices. A
        // solution to this problem is to normalize these vectors associated to
        // the new search directions. This essentially means normalizing the
        // columns of S.
        BVNormalize(searchSpace, nullptr);

        // Apply the projector to ensure X and S satisfy the divergence free
        // constraint
        //    C S = C X = 0
        // This with exact arithmetics is not necessary, but due to roundoff
        // errors, the solution is poluted, so we correct it every
        // n_steps_projection
        if (this->iter % n_steps_projection == 0) {
            projectBV(searchSpace);
        }

    } while ((this->iter <= this->maxIter) && (error_max > this->tolerance));

    // Transfer eigenvalues
    this->eigenvalues.assign(
        n_eigs, 0);  // re-initialize the eigenvalues not keep old values
    for (PetscInt eigen_v_idx = 0; eigen_v_idx < n_eigs; eigen_v_idx++) {
        this->eigenvalues[eigen_v_idx] = ritzValues[eigen_v_idx];
    }

    this->cleanupProjection();
    BVDestroy(&searchSpace);
    // Cleanup work memory
    BVDestroy(&tempBV);
    BVDestroy(&largeRitzVectors);
    BVDestroy(&residuals);
    BVDestroy(&doubleResiduals);
    PetscRandomDestroy(&random_context);

    // Show info that solver started
    if (rank == 0) {
        std::cout << "\n*******************************************************"
                     "*******"
                  << std::endl;
        std::cout << "* Doehler eingenvalue solver (PETSc) END" << std::endl;
        std::cout << "*********************************************************"
                     "*****\n"
                  << std::endl;
    }

    return 0;
}

void DoehlerMaxwellSolver::setMaxIter(int niter) { this->maxIter = niter; }

void DoehlerMaxwellSolver::setTolerance(PetscReal tol) {
    this->tolerance = tol;
}

PetscInt DoehlerMaxwellSolver::getConverged() {
    return this->eigenvectors_current_size;
}

PetscErrorCode DoehlerMaxwellSolver::getEigenPair(PetscInt index,
                                                  PetscScalar &eval,
                                                  Vec &evec) {
    logger.assert_always(index < this->getConverged(),
                         "Asking for eigenvalue % with only % converged", index,
                         this->getConverged());
    eval = this->eigenvalues[index];
    return BVCopyVec(this->eigenvectors, index, evec);
}

void DoehlerMaxwellSolver::setupProjection() {
    logger.assert_always(
        this->M == nullptr,
        "Projection unsupported with non identity mass matrix");

    // initialize the Y and H matrix needed for the projector
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

    // Set up solver
    KSPCreate(PETSC_COMM_WORLD, &this->projectionSolver_);
    KSPSetType(this->projectionSolver_, KSPGMRES);
    KSPSetOperators(this->projectionSolver_, this->H, this->H);
    KSPSetTolerances(this->projectionSolver_, 1e-12, 1.e-12, PETSC_DEFAULT,
                     100);
    KSPSetFromOptions(this->projectionSolver_);

    // Temporary vector used in the projection
    MatCreateVecs(this->H, &this->projectionTempVector_, nullptr);
}

void DoehlerMaxwellSolver::cleanupProjection() {
    VecDestroy(&this->projectionTempVector_);
    KSPDestroy(&this->projectionSolver_);
    MatDestroy(&this->H);
    MatDestroy(&this->Y);
}

void DoehlerMaxwellSolver::projectBV(BV bv) {
    PetscInt lead, active;
    BVGetActiveColumns(bv, &lead, &active);
    for (PetscInt i = lead; i < active; ++i) {
        Vec vec;
        BVGetColumn(bv, i, &vec);
        projectEigenVector(vec);
        BVRestoreColumn(bv, i, &vec);
    }
}

PetscErrorCode DoehlerMaxwellSolver::projectEigenVector(Vec &eigen_v) {

    // compute eigen_v = eigen_v - Y H^{-1} C.T eigen_v
    PetscInt its;

    auto tstart = std::chrono::high_resolution_clock::now();
    // compute rhs = C corr
    MatMult(this->C, eigen_v, this->projectionTempVector_);

    // solve the linear system
    //  H ksp_sol = (C eigen_v)
    // -> ksp_sol = H^{-1} C eigen_v
    // we store the solution (ksp_sol) in rhs
    KSPSolve(this->projectionSolver_, this->projectionTempVector_,
             this->projectionTempVector_);
    KSPGetIterationNumber(this->projectionSolver_, &its);

    // compute corr - Y H^{-1} C eigen_v
    // with rhs = H^{-1} C eigen_v
    VecScale(this->projectionTempVector_, -1.0);
    MatMultAdd(this->Y, this->projectionTempVector_, eigen_v, eigen_v);

    auto tstop = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::microseconds>(tstop - tstart);

    logger(DEBUG, "  -- projectEigenVector done in %d mu s [%d its]\n",
           duration.count(), its);

    return (0);
}

PetscErrorCode DoehlerMaxwellSolver::ritzUpdate(
    BV T_bv, PetscInt n_eigs, std::vector<PetscScalar> &ritzValues) {

    PetscErrorCode err;
    PetscInt iter_idx =
        0;  // initialize counter for number of interations performed

    // Reduced matrices obtained by projecting A_Mat and M_Mat
    // into the T_bv space (approximate eigenvectors \ocirc search space)
    Mat A_Mat_p, M_Mat_p, H_Mat_p,
        H_Mat_p1;  // H_Mat_p is a temporary hermitian matrix of either A_Mat_p
                   // or M_Mat_p
    MatCreateSeqDense(PETSC_COMM_SELF, 2 * n_eigs, 2 * n_eigs, NULL, &A_Mat_p);
    MatSetUp(A_Mat_p);

    MatCreateSeqDense(PETSC_COMM_SELF, 2 * n_eigs, 2 * n_eigs, NULL, &M_Mat_p);
    MatSetUp(M_Mat_p);

    // H_Mat_p will be created on first usage

    // Compute the reduced matrices on the space spanned by T = [X, S]
    err = BVMatProject(T_bv, this->A, T_bv, A_Mat_p);
    CHKERRABORT(PETSC_COMM_WORLD, err);

    err = BVMatProject(T_bv, this->M, T_bv, M_Mat_p);
    CHKERRABORT(PETSC_COMM_WORLD, err);

    // Make sure the resulting reduced matrices are still symmetric
    // Symmetry can be lost due to roundoff and accumulation errors

    // Force symmetry in A_Mat_p
    MatHermitianTranspose(A_Mat_p, MAT_INITIAL_MATRIX, &H_Mat_p);
    MatAXPY(A_Mat_p, 1.0, H_Mat_p, SAME_NONZERO_PATTERN);
    MatScale(A_Mat_p, 0.5);

    // Force symmetry in M_Mat_p
    MatHermitianTranspose(M_Mat_p, MAT_INITIAL_MATRIX, &H_Mat_p1);
    MatAXPY(M_Mat_p, 1.0, H_Mat_p1, SAME_NONZERO_PATTERN);
    MatScale(M_Mat_p, 0.5);

    // Setup the (small) eigensolver
    DS denseSolver;
    {
        DSCreate(PETSC_COMM_WORLD, &denseSolver);
        DSSetType(denseSolver, DSGHEP);
        DSAllocate(denseSolver, 2 * n_eigs);
        DSSetDimensions(denseSolver, 2 * n_eigs, 0, 0);
        // Comparison context for sorting the eigenvalues
        SlepcSC sc;
        DSGetSlepcSC(denseSolver, &sc);
        sc->comparison = SlepcCompareSmallestMagnitude;
        sc->comparisonctx = nullptr;
        sc->map = nullptr;
        sc->mapobj = nullptr;
        sc->rg = nullptr;
    }

    // Compute the Ritz values (L) and Ritz vectors (Q) of the reduced
    // eigenvalue problem
    {
        // Reset the eigenvalue solver
        DSSetState(denseSolver, DS_STATE_RAW);

        // Set the matrices
        Mat temp;
        DSGetMat(denseSolver, DS_MAT_A, &temp);
        MatCopy(A_Mat_p, temp, DIFFERENT_NONZERO_PATTERN);
        DSRestoreMat(denseSolver, DS_MAT_A, &temp);
        DSGetMat(denseSolver, DS_MAT_B, &temp);
        MatCopy(M_Mat_p, temp, DIFFERENT_NONZERO_PATTERN);
        DSRestoreMat(denseSolver, DS_MAT_B, &temp);

        // Solve & Sort

        std::vector<PetscScalar> evs1(2 * n_eigs);
        err = DSSolve(denseSolver, ritzValues.data(), evs1.data());
        CHKERRABORT(PETSC_COMM_WORLD, err);
        DSSort(denseSolver, ritzValues.data(), evs1.data(), nullptr, nullptr,
               nullptr);
        DSSynchronize(denseSolver, ritzValues.data(), nullptr);
        // Copy back the eigenvectors
    }

    // Now we can reconstruct the eigenvectors, i.e., compute them in the full
    // space

    // Reconstruct the eigenvectors (all at once)
    Mat ritzSmallVectors;  // get the matrix associated to the search space
                           // eigenvectors to use with mult below
    DSGetMat(denseSolver, DS_MAT_Q, &ritzSmallVectors);
    BVMultInPlace(T_bv, ritzSmallVectors, 0,
                  2 * n_eigs);  // make the multiplication T_bv_new = T_bv *
                                // Q_mat (reconstruct eigenvalues)
    DSRestoreMat(denseSolver, DS_MAT_Q, &ritzSmallVectors);

    BVSetActiveColumns(T_bv, 0, n_eigs);  // activate the columns associated
                                              // to the approximate solution
    BVCopy(T_bv, this->eigenvectors);

    BVSetActiveColumns(T_bv, 0,
                       2 * n_eigs);  // always return to original state

    DSDestroy(&denseSolver);
    MatDestroy(&A_Mat_p);
    MatDestroy(&M_Mat_p);
    MatDestroy(&H_Mat_p);
    MatDestroy(&H_Mat_p1);

    return 0;
}

void DoehlerMaxwellSolver::computeResiduals(
    const std::vector<PetscScalar> &ritzValues, BV vectors,
    PetscInt eigen_idx_start, PetscInt n_eigs, BV residuals, BV scratch) {
    PetscErrorCode err;
    PetscInt lead, active;
    BVGetActiveColumns(scratch, &lead, &active);
    BVSetActiveColumns(scratch, 0, n_eigs);
    //  Computes:
    //    R = A_Mat @ X_bv - M_Mat @ (X_bv * L)
    //
    //  Where:
    //    L_vec: are the n eigenvalues of the generalised eigenvalue problem A X
    //    = L M X X_bv: are n vectors with the same dimensions as the
    //    eigenvectors of the generalised eigenvalue problem A X = L M X A_Mat:
    //    is the A matrix of the generalised eigenvalue problem M_Mat: is the M
    //    matrix of the generalised eigenvalue problem R_bv: are the n residual
    //    vectors, associated to each eigenvalue and eigenvector pair
    //    eigen_idx_start: the start index for the eigenvalues to use
    //    n_eigs: the number of eigenvalues to use to compute the residual
    //            the eigenvalues will be in the range
    //               [eigen_idx_start, eigen_idx_start + n_eigs] (excluding the
    //               last)
    //            NOTE: X_bv and R_bv must have n_eigs columns.
    //    temp_bv: Temporary storage of the same size as X_bv

    // Compute M * X -> temp_bv
    if (this->M == nullptr) {
        BVCopy(vectors, scratch);
    } else {
        BVMatMult(vectors, this->M, scratch);
    }

    // compute in place temp_bv * L = M * X * L
    for (PetscInt i = 0; i < n_eigs; ++i) {
        BVScaleColumn(scratch, i, ritzValues[i + eigen_idx_start]);
    }

    // Compute A @ X
    // We use already R_bv, we then substract M @ (X * L) to get the final
    // residual value
    err = BVMatMult(vectors, this->A, residuals);  // make the multiplication A @ X
    CHKERRABORT(PETSC_COMM_WORLD, err);

    // Compute (A @ X) - (M @ X) * L
    BVMult(residuals, -1, 1.0, scratch, nullptr);

    // Restore active columns
    BVSetActiveColumns(scratch, lead, active);
}

}  // namespace EigenSolvers

}  // namespace hpgem

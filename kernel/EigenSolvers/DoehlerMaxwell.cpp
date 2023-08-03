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

DoehlerMaxwellSolver::DoehlerMaxwellSolver()
    : eigenvectors(nullptr) {
  std::cout << "Initialize DoehlerMaxwellSolver!!!!!" << std::endl;
}

DoehlerMaxwellSolver::~DoehlerMaxwellSolver() {
    BVDestroy(&eigenvectors);
}

void DoehlerMaxwellSolver::setMatrices(const Mat Ain, const Mat Cin) {
    // Set the A and C matrices describing the eigenvalue problem
    // Ax = lambda x with C x = 0
    
    // Set the A, and C, but set M to be the identity matrix (NULL is a shorthand)
    this->setMatrices(Ain, Cin, NULL);
}


void DoehlerMaxwellSolver::setMatrices(const Mat Ain, const Mat Cin, const Mat Min) {
  // Set the A, C, and M matrices describing the eigenvalue problem
  // Ax = lambda M x with C x = 0.
  // If Min = NULL then M is assumed to be the identity matrix.
  PetscErrorCode ierr;
  PetscInt n, m;
  PetscViewer viewer;
  PetscBool knows_is_hermitian;
  PetscBool is_hermitian;

  this->A = Ain;
  this->C = Cin;
  this->M = Min;

  // Check if A is hermitian
  
  // Directly checking if A is hermitian is expensive, therefore we check if the 
  // matrix knows it is hermitian. If it knows, inform the user. If not, inform 
  // the user, we do not know if the matrix is hermitian or not.
  ierr = MatIsHermitianKnown(A, &knows_is_hermitian, &is_hermitian);
  
  if (knows_is_hermitian) {
    if (!is_hermitian) {
        // logger(WARN, "Matrix A is not Hermitian\n");
    }
  }
  else {
    // logger(WARN, "Matrix A has no information on its Hermitian property\n");
  }
  
  ierr = MatGetSize(Ain, &n, &m);
  // logger(INFO, "DoehlerSolver System Size : % x %", n, m);
  
  if (Min != NULL) {
    ierr = MatGetSize(Min, &n, &m);
    // logger(INFO, "DoehlerSolver RHS Size : % x %", n, m);
  }
  else {
    // logger(INFO, "DoehlerSolver RHS set to identity");
  }
  
  ierr = MatGetSize(Cin, &n, &m);
  // logger(INFO, "DoehlerSolver Contraint Size : % x %", n, m);
}

PetscErrorCode DoehlerMaxwellSolver::solve(PetscInt nev, Mat &T_Mat_in, PetscInt n_steps_projection) {
  
  // Initialize the matrices for enforcing the condition this->C * x = 0
  // This is done by solving a Poisson-like problem, the matrices required for
  // solving this problem are initialized here
  this->setupProjection();

  PetscInt n_eigs = nev;  // this is done to keep the function parameter with the same 
                          // names across implementations, but to keep the clearer name 
                          // n_eigs internally
  PetscErrorCode err;
  bool indefinite_dot = true;
  bool verbose = true;
  bool normalize_S = true;
  
  // Get the rank of the process for display purposes
  PetscMPIInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  
  // Show info that solver started
  if(rank == 0)
  {
    std::cout << "\n**************************************************************" << std::endl;
    std::cout << "* Doehler eingenvalue solver (PETSc) START" << std::endl;
    std::cout << "**************************************************************\n" << std::endl;
  }
    
  // Initialize the algorithm with initial guesses
  
  // Compute initial parameters, like system size, etc 
  PetscInt A_n_rows, A_n_cols;  // number of rows and columns of matrix A
  PetscInt A_n_local_rows;
  MatGetSize(this->A, &A_n_rows, &A_n_cols);
  MatGetLocalSize(this->A, &A_n_local_rows, nullptr);
  
  if(rank == 0)
  {
    std::cout << "System size:" << std::endl;
    std::cout << "   n_rows: " << A_n_rows << std::endl;
    std::cout << "   n_cols: " << A_n_cols << std::endl;
    std::cout << "   n_eigs: " << n_eigs << std::endl;
  }
    
  // Initialize the eigenvector solution
  if(this->eigenvectors == nullptr) {
      BVCreate(PETSC_COMM_WORLD, &this->eigenvectors);
  }
  err = BVSetSizes(this->eigenvectors, A_n_local_rows, A_n_rows, n_eigs);
  CHKERRABORT(PETSC_COMM_WORLD, err);
  BVSetFromOptions(this->eigenvectors);
  std::vector<PetscScalar> ritzValues (2*n_eigs);
  
  // Initialize tranformation matrix T as a bv system (used to project to the reduced space)
  BV T_bv;
  BVCreate(PETSC_COMM_WORLD, &T_bv);
  
  BVSetSizes(T_bv, A_n_local_rows, A_n_rows, 2*n_eigs);
  BVSetFromOptions(T_bv);
  
  // Set the initial (guess) values
  // (the n_eigs eigenvectors we wish to find and the 
  // n_eigs vector search directions)
  PetscRandom random_context (nullptr);
  if(false)
    {
      // Use initial values from python for one to one comparison
      Mat T_Mat;
      MatCreateDense(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, A_n_rows, 2*n_eigs, NULL, &T_Mat);
      MatSetUp(T_Mat); 
      BVGetMat(T_bv, &T_Mat);
      MatCopy(T_Mat_in, T_Mat, DIFFERENT_NONZERO_PATTERN);
      BVRestoreMat(T_bv, &T_Mat);
      if(rank == 0) std::cout << "Generate T_bv from input" << std::endl;
    }
    else
    {
      // Use randomly generated values to initialize T_bv
      PetscRandomCreate(PETSC_COMM_WORLD, &random_context);
      PetscRandomSetInterval(random_context, 0.0, 1.0);
      PetscRandomSetFromOptions(random_context);
      BVSetRandomContext(T_bv, random_context);
      BVSetRandom(T_bv);
      if(rank == 0) std::cout << "Generate T_bv internally" << std::endl;
    }

  // Apply the projector to ensure X and S satisfy the divergence free constraint
  // C S = C X = 0
  this->projectBV(T_bv);
  
  // TODO Consider forcing all small matrices and vectors to be sequential
  
  // Iterate to find corrected solutions to the eigenvalue
  PetscInt iter_idx = 0;  // initialize counter for number of interations performed
  PetscReal error_max = 1.0;  // initialize error max to determine if the loop is over or not
  
  // Declare variables required in while loop 
  
  // Reduced matrices obtained by projecting A_Mat and M_Mat
  // into the T_bv space (approximate eigenvectors \ocirc search space) 
  Mat A_Mat_p, M_Mat_p, H_Mat_p, H_Mat_p1;  // H_Mat_p is a temporary hermitian matrix of either A_Mat_p or M_Mat_p
  MatCreateSeqDense(PETSC_COMM_SELF, 2*n_eigs, 2*n_eigs, NULL, &A_Mat_p);
  MatSetUp(A_Mat_p);
   
  MatCreateSeqDense(PETSC_COMM_SELF, 2*n_eigs, 2*n_eigs, NULL, &M_Mat_p);
  MatSetUp(M_Mat_p);

  //H_Mat_p will be created on first usage
  
//  // Setup the (small) eigensolver
//  DS denseSolver;
//  {
//      DSCreate(PETSC_COMM_WORLD, &denseSolver);
//      DSSetType(denseSolver, DSGHEP);
//      DSAllocate(denseSolver, 2 * n_eigs);
//      DSSetDimensions(denseSolver, 2 * n_eigs, 0, 0);
//      // Comparison context for sorting the eigenvalues
//      SlepcSC sc;
//      DSGetSlepcSC(denseSolver, &sc);
//      sc->comparison = SlepcCompareSmallestMagnitude;
//      sc->comparisonctx = nullptr;
//      sc->map = nullptr;
//      sc->mapobj = nullptr;
//      sc->rg = nullptr;
//  }

  
  do
  {
    iter_idx++;  // update counter for number of iterations performed
    
//    // Compute the reduced matrices on the space spanned by T = [X, S]
//    err = BVMatProject(T_bv, this->A, T_bv, A_Mat_p);
//    CHKERRABORT(PETSC_COMM_WORLD, err);
//
//    err = BVMatProject(T_bv, this->M, T_bv, M_Mat_p);
//    CHKERRABORT(PETSC_COMM_WORLD, err);
//
//    // Make sure the resulting reduced matrices are still symmetric
//    // Symmetry can be lost due to roundoff and accumulation errors
//
//    // Force symmetry in A_Mat_p
//    MatHermitianTranspose(A_Mat_p, iter_idx == 1 ? MAT_INITIAL_MATRIX : MAT_REUSE_MATRIX, &H_Mat_p);
//    MatAXPY(A_Mat_p, 1.0, H_Mat_p, SAME_NONZERO_PATTERN);
//    MatScale(A_Mat_p, 0.5);
//
//    // Force symmetry in M_Mat_p
//    MatHermitianTranspose(M_Mat_p, iter_idx == 1 ? MAT_INITIAL_MATRIX : MAT_REUSE_MATRIX, &H_Mat_p1);
//    MatAXPY(M_Mat_p, 1.0, H_Mat_p1, SAME_NONZERO_PATTERN);
//    MatScale(M_Mat_p, 0.5);
//
//    // std::cout << "\n\nThe A matrix:" << std::endl;
//    // MatView(A_Mat_p, PETSC_VIEWER_STDOUT_WORLD);
//
//    // std::cout << "\n\nThe M matrix:" << std::endl;
//    // MatView(M_Mat_p, PETSC_VIEWER_STDOUT_WORLD);
//
//    // std::cout << "\n\nThe L_Vec vector:" << std::endl;
//    // VecView(L_Vec, PETSC_VIEWER_STDOUT_WORLD);
//
//    // std::cout << "\n\nThe A matrix:" << std::endl;
//    // MatView(A_Mat_p, PETSC_VIEWER_STDOUT_WORLD);
//    //
//    // std::cout << "\n\nThe M matrix:" << std::endl;
//    // MatView(M_Mat_p, PETSC_VIEWER_STDOUT_WORLD);
//
//    // Compute the Ritz values (L) and Ritz vectors (Q) of the reduced eigenvalue problem
//    {
//        // Reset the eigenvalue solver
//        DSSetState(denseSolver, DS_STATE_RAW);
//
//        // Set the matrices
//        Mat temp;
//        DSGetMat(denseSolver, DS_MAT_A, &temp);
//        MatCopy(A_Mat_p, temp, DIFFERENT_NONZERO_PATTERN);
//        DSRestoreMat(denseSolver, DS_MAT_A, &temp);
//        DSGetMat(denseSolver, DS_MAT_B, &temp);
//        MatCopy(M_Mat_p, temp, DIFFERENT_NONZERO_PATTERN);
//        DSRestoreMat(denseSolver, DS_MAT_B, &temp);
//
//        // Solve & Sort
//
//        std::vector<PetscScalar> evs1 (2*n_eigs);
//        err = DSSolve(denseSolver, ritzValues.data(), evs1.data());
//        CHKERRABORT(PETSC_COMM_WORLD, err);
//        DSSort(denseSolver, ritzValues.data(), evs1.data(), nullptr, nullptr, nullptr);
//        DSSynchronize(denseSolver, ritzValues.data(), nullptr);
//        // Copy back the eigenvectors
//    }
//
//    // std::cout << "\n\nEigenvalues:" << std::endl;
//    // VecView(L_Vec, PETSC_VIEWER_STDOUT_WORLD);
//
//    // Now we can reconstruct the eigenvectors, i.e., compute them in the full space
//    BV T_bv_new;  // the updated reconstructed vectors, below just make a copy to start with
//    BVDuplicate(T_bv, &T_bv_new);
//    BVCopy(T_bv, T_bv_new);
//
//    // Reconstruct the eigenvectors (all at once)
//    Mat ritzSmallVectors;  // get the matrix associated to the search space eigenvectors to use with mult below
//    DSGetMat(denseSolver, DS_MAT_Q, &ritzSmallVectors);
//    BVMultInPlace(T_bv_new, ritzSmallVectors, 0, 2*n_eigs); // make the multiplication T_bv_new = T_bv * Q_mat (reconstruct eigenvalues)
//    DSRestoreMat(denseSolver, DS_MAT_Q, &ritzSmallVectors);
//
//    BVSetActiveColumns(T_bv_new, 0, n_eigs);  // activate the columns associated to the approximate solution
//    BVCopy(T_bv_new, this->eigenvectors);
//
//    BVSetActiveColumns(T_bv_new, 0, 2*n_eigs);  // always return to original state
    BV T_bv_new;  // the updated reconstructed vectors, below just make a copy to start with
    this->computeRitzValuesAndVectors(T_bv, n_eigs, ritzValues, T_bv_new);
      
    BV W_r_bv;  // the search space
    BVCreate(PETSC_COMM_WORLD, &W_r_bv);
    BVSetSizes(W_r_bv, A_n_local_rows, A_n_rows, n_eigs);
    BVSetFromOptions(W_r_bv);
    
    BVSetActiveColumns(T_bv_new, n_eigs, 2*n_eigs);  // activate the columns associated to the search space
    BVCopy(T_bv_new, W_r_bv);
    BVSetActiveColumns(T_bv_new, 0, 2*n_eigs);  // always return to original state
    
    // Compute the residual with the updated eigenvectors
    BV R_bv;  // the residual column vectors
    BVCreate(PETSC_COMM_WORLD, &R_bv);
    BVSetSizes(R_bv, A_n_local_rows, A_n_rows, n_eigs);
    BVSetFromOptions(R_bv);
    this->compute_residual_eigen_v(this->A, this->M, ritzValues, this->eigenvectors, 0, n_eigs, R_bv);

    // Compute convergence
    PetscReal residualNorm, evNorm;
    this->eigenvectors_current_size = 0;
    std::stringstream residual_values;
    for(PetscInt i = 0; i < n_eigs; ++i) {
        BVNormColumn(R_bv, i, NORM_2, &residualNorm);
        BVNormColumn(this->eigenvectors, i, NORM_2, &evNorm);
        residual_values << " " << std::setprecision(16) << ritzValues[i]
                        << "(" << std::setprecision(2) << (residualNorm/evNorm) << ")";
        if (this->eigenvectors_current_size != i) {
            // Only accept an eigenvalue as converged if all previous
            // eigenvalues also have converged. This ensures that we get the
            // wanted eigenvalues.
            continue;
        }
        // Check for convergence either in either relative or absolute sense
        if (residualNorm < evNorm * this->tolerance
            || residualNorm < evNorm * this->tolerance * std::abs(ritzValues[i])) {
            this->eigenvectors_current_size++;
        }
    }
    if (iter_idx % 5 == 0 && rank == 0) {
        std::string log = residual_values.str();
        logger(INFO, "iter % converged %:%", iter_idx,
               this->eigenvectors_current_size, log);
    }
    if (this->eigenvectors_current_size == n_eigs) {
        // TODO: Improve
        error_max = 0.0;
    }
    MPI_Barrier(PETSC_COMM_WORLD);
      
    // TODO Apply preconditioner  
    
    // Compute the new augmented solution space (the correction space) and the new search space
    BV RR_bv;  // the BV containing the residual of the residual BV
    BVCreate(PETSC_COMM_WORLD, &RR_bv);
    BVSetSizes(RR_bv, A_n_local_rows, A_n_rows, n_eigs);
    BVSetFromOptions(RR_bv);
    this->compute_residual_eigen_v(this->A, this->M, ritzValues, R_bv, 0, n_eigs, RR_bv);
    
    BVSetActiveColumns(T_bv_new, n_eigs, 2*n_eigs);
//    Mat W_r_bv_Mat;
//    BVGetMat(W_r_bv, &W_r_bv_Mat);
//
//    BV V_bv;
//    BVCreate(PETSC_COMM_WORLD, &V_bv);
//    BVSetSizes(V_bv, PETSC_DECIDE, n_eigs, n_eigs);
//    BVSetFromOptions(V_bv);

    // Part 1 of the idea: Use BVDot to create a sequentially dense matrix
//    // T_{ij} = -v_j^H(A - lm_j B) w_i / (lm_{i+p} - lm_j)
//    // -v_j^H(A - lm_j B)w = R(v)^H W
    Mat out; // Make it
    MatCreateSeqDense(PETSC_COMM_SELF, n_eigs, n_eigs, NULL, &out);
    err = BVDot(RR_bv, W_r_bv, out);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    // Use this matrix later on in the BVMult(R, 1, 1, W_r_bv, out);

    // Part 2 of the idea: To update it pull it out of petsc
    std::vector<PetscScalar> values (n_eigs*n_eigs);
    std::vector<PetscInt> indices(n_eigs);
    std::iota(indices.begin(), indices.end(), 0);
    MatGetValues(out, indices.size(), indices.data(), indices.size(), indices.data(), values.data());
    // Update
    for(std::size_t col = 0; col < n_eigs; col++) {
      for(std::size_t row = 0; row < n_eigs; row++) {
        values[row * n_eigs + col] /= -(ritzValues[row + n_eigs] - std::conj(ritzValues[col]));
      }
    }
    MatSetValues(out, indices.size(), indices.data(), indices.size(), indices.data(),values.data(), INSERT_VALUES);
    err = MatAssemblyBegin(out, MAT_FINAL_ASSEMBLY);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    err =  MatAssemblyEnd(out, MAT_FINAL_ASSEMBLY);
    CHKERRABORT(PETSC_COMM_WORLD, err);


//    BVMatMultHermitianTranspose(RR_bv, W_r_bv_Mat, V_bv);
    
//    BVRestoreMat(W_r_bv, &W_r_bv_Mat);  // restore the matrix so that we can reuse W_r_bv

//    // Divide by -(L[n_eigs:].reshape([-1, 1]) - L[:n_eigs].conjugate())
//    for(PetscInt column_idx = 0; column_idx < n_eigs; column_idx++)
//    {
//      Vec V_bv_col_Vec;
//      BVGetColumn(V_bv, column_idx, &V_bv_col_Vec);
//      // Get ownership ranges so that we only loop over the local range
//      PetscInt start_row, end_row;
//      VecGetOwnershipRange(V_bv_col_Vec, &start_row, &end_row);
//      for(PetscInt row_idx = start_row; row_idx < end_row; row_idx++)
//      {
//        PetscScalar V_value;
//        PetscScalar L_value_row_idx;  // L[row_idx + n_eigs]
//        PetscScalar L_value_column_idx;  // L[column_idx]
//
//        PetscInt row_idx_offset = row_idx + n_eigs;
//
//        err = VecGetValues(V_bv_col_Vec, 1, &row_idx, &V_value);
//        CHKERRABORT(PETSC_COMM_WORLD, err);
//        L_value_row_idx = ritzValues[row_idx_offset];
//        L_value_column_idx = ritzValues[column_idx];
//        V_value = -V_value / (L_value_row_idx - std::conj(L_value_column_idx));
//
//        err = VecSetValue(V_bv_col_Vec, row_idx, V_value, INSERT_VALUES);
////        CHKERRABORT(PETSC_COMM_WORLD, err);
////        err = VecAssemblyBegin(V_bv_col_Vec);
////        CHKERRABORT(PETSC_COMM_WORLD, err);
////        err = VecAssemblyEnd(V_bv_col_Vec);
////        CHKERRABORT(PETSC_COMM_WORLD, err);
//      }
//
//      BVRestoreColumn(V_bv, column_idx, &V_bv_col_Vec);
//    }
    
//    // Compute the new search space
//    Mat V_bv_Mat;
//    BVGetMat(V_bv, &V_bv_Mat);
//    Mat V_bv_Mat_seq;
//    PetscInt n_communicators;
//    MPI_Comm_size(PETSC_COMM_WORLD, &n_communicators);
//    MatCreateRedundantMatrix(V_bv_Mat, n_communicators, PETSC_COMM_SELF, MAT_INITIAL_MATRIX, &V_bv_Mat_seq);

//      Mat V_bv_Mat_seq;
//      MatCreateSeqDense(PETSC_COMM_SELF, n_eigs, n_eigs, NULL, &V_bv_Mat_seq);
//      MatSetRandom(V_bv_Mat_seq, random_context);
      
      BVMult(R_bv, 1.0, 1.0, W_r_bv, out);
//    BVMult(R_bv, 1.0, 1.0, W_r_bv, V_bv_Mat_seq);
//    BVRestoreMat(V_bv, &V_bv_Mat);
    
//    MatDestroy(&V_bv_Mat_seq);
//    MatDestroy(&V_bv_Mat);
      
    // Restart T_bv
      
    // Update X part
    for(PetscInt column_idx = 0; column_idx < n_eigs; column_idx++)
    {
      Vec T_bv_column_Vec;
      BVGetColumn(T_bv, column_idx, &T_bv_column_Vec);  // get the column we want to replace in T_bv and link it to vector T_bv_column_Vec
      BVCopyVec(this->eigenvectors, column_idx, T_bv_column_Vec);     // copy the column of this->eigenvectors into the T_bv_column_Vec vector, effectively 
                                                                      // copying the column of this->eigenvectors into the column of T_bv 
      BVRestoreColumn(T_bv, column_idx, &T_bv_column_Vec);
    }
    
    // Update S part
    for(PetscInt column_idx = 0; column_idx < n_eigs; column_idx++)
    {
      PetscInt column_idx_offset = column_idx + n_eigs;  // we need to offset because now we are updating the S part of T_bv
      
      Vec T_bv_column_Vec;
      BVGetColumn(T_bv, column_idx_offset, &T_bv_column_Vec);  // get the column we want to replace in T_bv and link it to vector T_bv_column_Vec
      BVCopyVec(R_bv, column_idx, T_bv_column_Vec);     // copy the column of R_bv into the T_bv_column_Vec vector, effectively 
                                                        // copying the column of R_bv into the column of T_bv 
                                                        // Note that the column index of R_bv is not offset, since R_bv 
                                                        // contains only the search space, only T_bv contains both the 
                                                        // (approximate) solution and the search space, hence the offset needed. 
      
      // If nothing is done, the search directions result in vectors with very
      // small norms. This leads to very badly conditioned reduced matrices. A
      // solution to this problem is to normalize these vectors associated to
      // the new search directions. This essentially means normalizing the columns
      // of S.
//      if(normalize_S)
//      {
//        // This is not the same as done in Matlab, in Matlab it is a pseudo-norm
//        //    VecNormalize(T_bv_column_Vec, NULL);
//        // So we make it equal to Matlab
//        Vec T_bv_column_mult;
//        VecDuplicate(T_bv_column_Vec, &T_bv_column_mult);
//        PetscScalar T_bv_column_pseudonorm;
//        VecCopy(T_bv_column_Vec, T_bv_column_mult);  // first make a copy to store the pointwise multiplication
//        VecPointwiseMult(T_bv_column_mult, T_bv_column_mult, T_bv_column_mult);  // pointwise multiplication
//        VecSum(T_bv_column_mult, &T_bv_column_pseudonorm);  // sum the components of pointwise multiplication to get a pseudonorm
//        VecScale(T_bv_column_Vec, 1.0/T_bv_column_pseudonorm);  // scale the vector with the pseudonorm
//
//        // Cleanup work memory
//        VecDestroy(&T_bv_column_mult);
//      }
//      VecNormalize(T_bv_column_Vec, nullptr);
                                                          
      BVRestoreColumn(T_bv, column_idx_offset, &T_bv_column_Vec);
    }
    BVOrthogonalize(T_bv, nullptr);
     
    
    // Apply the projector to ensure X and S satisfy the divergence free constraint
    //    C S = C X = 0
    // This with exact arithmetics is not necessary, but due to roundoff errors, 
    // the solution is poluted, so we correct it every n_steps_projection
    if(iter_idx % n_steps_projection == 0)
    {
      projectBV(T_bv);
    }
    
    // Cleanup work memory
    BVDestroy(&T_bv_new);
    BVDestroy(&W_r_bv);
    BVDestroy(&R_bv);
    BVDestroy(&RR_bv);
//    BVDestroy(&V_bv);
//    MatDestroy(&W_r_bv_Mat);
//    MatDestroy(&V_bv_Mat);
    
  } while((iter_idx <= this->maxIter) && (error_max > this->tolerance));
  
  // Transfer eigenvalues
  this->eigenvalues.assign(n_eigs, 0);  // re-initialize the eigenvalues not keep old values 
  for(PetscInt eigen_v_idx = 0; eigen_v_idx < n_eigs; eigen_v_idx++) {
    this->eigenvalues[eigen_v_idx] = ritzValues[eigen_v_idx];
  }
  iter = iter_idx;

  this->cleanupProjection();
//  DSDestroy(&denseSolver);
  BVDestroy(&T_bv);
//  MatDestroy(&A_Mat_p);
//  MatDestroy(&M_Mat_p);
//  MatDestroy(&H_Mat_p);
//  MatDestroy(&H_Mat_p1);
  PetscRandomDestroy(&random_context);
  
  // Show info that solver started
  if(rank == 0)
  {
    std::cout << "\n**************************************************************" << std::endl;
    std::cout << "* Doehler eingenvalue solver (PETSc) END" << std::endl;
    std::cout << "**************************************************************\n" << std::endl;
  }
  
  return 0;
}

void DoehlerMaxwellSolver::setMaxIter(int niter) {
    this->maxIter = niter;
}

void DoehlerMaxwellSolver::setTolerance(PetscReal tol) {
    this->tolerance = tol;
}

PetscInt DoehlerMaxwellSolver::getConverged() {
    return this->eigenvectors_current_size;
}

PetscErrorCode DoehlerMaxwellSolver::getEigenPair(PetscInt index,
                                                  PetscScalar &eval,
                                                  Vec &evec) {
    logger.assert_always(index < this->getConverged(), "Asking for eigenvalue % with only % converged",
                         index, this->getConverged());
    eval = this->eigenvalues[index];
    return BVCopyVec(this->eigenvectors, index, evec);
}

void DoehlerMaxwellSolver::setupProjection() {
    logger.assert_always(this->M == nullptr,
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
//    PC pc;
//    KSPCreate(PETSC_COMM_WORLD, &projectionSolver_);
//    KSPSetType(projectionSolver_, KSPPREONLY);
//    KSPSetOperators(projectionSolver_, this->H, this->H);
//    KSPGetPC(projectionSolver_, &pc);
//    PCSetType(pc, PCLU);
//    PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);
//    KSPSetTolerances(this->projectionSolver_, 1e-12, 1.e-12, PETSC_DEFAULT, 100);
//    KSPSetFromOptions(this->projectionSolver_);
    
    KSPCreate(PETSC_COMM_WORLD, &this->projectionSolver_);
    KSPSetType(this->projectionSolver_, KSPGMRES);
    KSPSetOperators(this->projectionSolver_, this->H, this->H);
    KSPSetTolerances(this->projectionSolver_, 1e-12, 1.e-12, PETSC_DEFAULT, 100);
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
    for(PetscInt i = lead; i < active; ++i) {
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

    // set up the linear system
    // KSPCreate(PETSC_COMM_WORLD, &ksp);
    // KSPSetType(ksp, KSPGMRES);
    // KSPSetOperators(ksp, this->H, this->H);
    // KSPSetFromOptions(ksp);
    // KSPSetTolerances(ksp, 1e-12, 1.e-12, PETSC_DEFAULT, 100);

    // solve the linear system
    //  H ksp_sol = (C eigen_v)
    // -> ksp_sol = H^{-1} C eigen_v
    // we store the solution (ksp_sol) in rhs
    KSPSolve(this->projectionSolver_, this->projectionTempVector_, this->projectionTempVector_);
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



PetscErrorCode DoehlerMaxwellSolver::computeRitzValuesAndVectors(BV &T_bv, PetscInt n_eigs, std::vector<PetscScalar> &ritzValues, BV &T_bv_new){

    PetscErrorCode err;
    PetscInt iter_idx = 0;  // initialize counter for number of interations performed
    
    // Reduced matrices obtained by projecting A_Mat and M_Mat
    // into the T_bv space (approximate eigenvectors \ocirc search space)
    Mat A_Mat_p, M_Mat_p, H_Mat_p, H_Mat_p1;  // H_Mat_p is a temporary hermitian matrix of either A_Mat_p or M_Mat_p
    MatCreateSeqDense(PETSC_COMM_SELF, 2*n_eigs, 2*n_eigs, NULL, &A_Mat_p);
    MatSetUp(A_Mat_p);

    MatCreateSeqDense(PETSC_COMM_SELF, 2*n_eigs, 2*n_eigs, NULL, &M_Mat_p);
    MatSetUp(M_Mat_p);

    //H_Mat_p will be created on first usage
    
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

    // std::cout << "\n\nThe A matrix:" << std::endl;
    // MatView(A_Mat_p, PETSC_VIEWER_STDOUT_WORLD);

    // std::cout << "\n\nThe M matrix:" << std::endl;
    // MatView(M_Mat_p, PETSC_VIEWER_STDOUT_WORLD);

    // std::cout << "\n\nThe L_Vec vector:" << std::endl;
    // VecView(L_Vec, PETSC_VIEWER_STDOUT_WORLD);

    // std::cout << "\n\nThe A matrix:" << std::endl;
    // MatView(A_Mat_p, PETSC_VIEWER_STDOUT_WORLD);
    //
    // std::cout << "\n\nThe M matrix:" << std::endl;
    // MatView(M_Mat_p, PETSC_VIEWER_STDOUT_WORLD);
    
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

    // Compute the Ritz values (L) and Ritz vectors (Q) of the reduced eigenvalue problem
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

        std::vector<PetscScalar> evs1 (2*n_eigs);
        err = DSSolve(denseSolver, ritzValues.data(), evs1.data());
        CHKERRABORT(PETSC_COMM_WORLD, err);
        DSSort(denseSolver, ritzValues.data(), evs1.data(), nullptr, nullptr, nullptr);
        DSSynchronize(denseSolver, ritzValues.data(), nullptr);
        // Copy back the eigenvectors
    }

    // std::cout << "\n\nEigenvalues:" << std::endl;
    // VecView(L_Vec, PETSC_VIEWER_STDOUT_WORLD);

    // Now we can reconstruct the eigenvectors, i.e., compute them in the full space
    BVDuplicate(T_bv, &T_bv_new);
    BVCopy(T_bv, T_bv_new);

    // Reconstruct the eigenvectors (all at once)
    Mat ritzSmallVectors;  // get the matrix associated to the search space eigenvectors to use with mult below
    DSGetMat(denseSolver, DS_MAT_Q, &ritzSmallVectors);
    BVMultInPlace(T_bv_new, ritzSmallVectors, 0, 2*n_eigs); // make the multiplication T_bv_new = T_bv * Q_mat (reconstruct eigenvalues)
    DSRestoreMat(denseSolver, DS_MAT_Q, &ritzSmallVectors);

    BVSetActiveColumns(T_bv_new, 0, n_eigs);  // activate the columns associated to the approximate solution
    BVCopy(T_bv_new, this->eigenvectors);

    BVSetActiveColumns(T_bv_new, 0, 2*n_eigs);  // always return to original state
    
    DSDestroy(&denseSolver);
    MatDestroy(&A_Mat_p);
    MatDestroy(&M_Mat_p);
    MatDestroy(&H_Mat_p);
    MatDestroy(&H_Mat_p1);
    
    return 0;
}



void DoehlerMaxwellSolver::compute_residual_eigen_v(Mat &A_Mat, Mat &M_Mat, const std::vector<PetscScalar>& ritzValues, BV &X_bv, PetscInt eigen_idx_start, PetscInt n_eigs, BV &R_bv){
    PetscErrorCode err;
  //  Computes:
  //    R = A_Mat @ X_bv - M_Mat @ (X_bv * L)
  //   
  //  Where:
  //    L_vec: are the n eigenvalues of the generalised eigenvalue problem A X = L M X
  //    X_bv: are n vectors with the same dimensions as the eigenvectors of the generalised eigenvalue problem A X = L M X
  //    A_Mat: is the A matrix of the generalised eigenvalue problem
  //    M_Mat: is the M matrix of the generalised eigenvalue problem
  //    R_bv: are the n residual vectors, associated to each eigenvalue and eigenvector pair
  //    eigen_idx_start: the start index for the eigenvalues to use 
  //    n_eigs: the number of eigenvalues to use to compute the residual
  //            the eigenvalues will be in the range 
  //               [eigen_idx_start, eigen_idx_start + n_eigs] (excluding the last)
  //            NOTE: X_bv and R_bv must have n_eigs columns.
  
  // Compute X * L
  BV XL_bv;  // to hold X * L[:2*n_eigs]
  BVDuplicate(X_bv, &XL_bv);
  BVCopy(X_bv, XL_bv);  // start by copying X and then below we multiply each column by the eigen value
  
  for(PetscInt eigen_v_idx = eigen_idx_start; eigen_v_idx < (eigen_idx_start + n_eigs); eigen_v_idx++)
  {
    // The BV contains only the eigenvectors we want from [0, n_eigs]
    // but the eigenvalue vector L_Vec contains all the eigenvalues 
    // by providing the eigen_v_idx and n_eigs parameters we can select 
    // the eigenvalues of the vector L_vec. We then just need to match 
    // the indices as done below.
    PetscInt bv_column_idx = eigen_v_idx - eigen_idx_start;  // get the column index of the eigenvector 
                                                             // of the BV associated to eigenvalue eigen_v_idx
    PetscScalar eigen_value = ritzValues[eigen_v_idx];
    BVScaleColumn(XL_bv, bv_column_idx, eigen_value);
  }
  
  // Compute M @ (X * L)
  BV MXL_bv;  // to hold M @ (X * L)
  if(M_Mat == NULL)
  {
    BVDuplicate(XL_bv, &MXL_bv);
    BVCopy(XL_bv, MXL_bv);
  }
  else
  {
    BVDuplicate(XL_bv, &MXL_bv);  
    err = BVMatMult(XL_bv, M_Mat, MXL_bv);  // make the multiplication M @ (X * L[eigen_idx_start:eigen_idx_end])
    CHKERRABORT(PETSC_COMM_WORLD, err);
  }
  
  // Compute A @ X
  // We use already R_bv, we then substract M @ (X * L) to get the final residual value
  err = BVMatMult(X_bv, A_Mat, R_bv);  // make the multiplication A @ X
  CHKERRABORT(PETSC_COMM_WORLD, err);
  
  // Compute (A @ X) - (M @ (X * L))
  BVMult(R_bv, -1, 1.0, MXL_bv, NULL);
  
  // Cleanup work memory
  BVDestroy(&XL_bv);
  BVDestroy(&MXL_bv);
}

}  // namespace EigenSolvers

}  // namespace hpgem

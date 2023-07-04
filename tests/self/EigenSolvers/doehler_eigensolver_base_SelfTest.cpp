/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2020, University of Twente
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors
 may be used to endorse or promote products derived from this software without
 specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <iostream>
#include <string>
#include <iomanip>

#include "petsc.h"
#include "slepc.h"

#include "EigenSolvers/DoehlerMaxwell.h"

// #include "Base/CommandLineOptions.h"
// 
// #include "DGMaxLogger.h"
// #include "DGMaxProgramUtils.h"
// #include "Utils/Verification/DivDGMaxEVConvergenceTest.h"
// #include "Utils/Verification/EVTestPoint.h"

static char help[] = "Solve an eigenvalue problem Ax = lMx using Doehler's algorithm.";

PetscErrorCode read_matrix(Mat &M, std::string &name, std::string &filename, std::string &read_format);
PetscErrorCode read_matrix(Mat &M, std::string &name, PetscViewer &viewer);

PetscErrorCode main(int argc, char** argv) {
    std::cout << "Testing Doehler eigensolver..." << std::endl; 
    
    /*
      PETSc and SLEPc need to be initialized, if we initialize SLEPc, PETSc is
      also initialized, so we just initalize SLEPc with the helper message.
    */
    SlepcInitialize(&argc, &argv, PETSC_NULL, help);
    
    // Computation parameters
    PetscInt n_eigenvalues = 20;  // number of eigenvalues to compute
    PetscReal tol = 1e-5;  // tolerance to use as stopping criterium for eigenvalues
    int max_iter = 3000;  // maximum number of iterations to perform before stopping
    
    /*
      Read the matrices from the data file.
    */
    
    // Declare matrices 
    Mat A, M, T, C;
    
    // The input file and data format
    std::string ifilename("/Users/apalha/work/dev/hpgem/eigensolver_doehler/petsc__eigensolver_doehler/src/python/A_M_T_C_4_matrices.dat");
    std::string data_format("binary");
    
    // The matrix names 
    std::string A_name("A");
    std::string M_name("M");
    std::string T_name("T");
    std::string C_name("C");
    
    // Load the matrices
    PetscViewer viewer;
    PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
    PetscViewerSetType(viewer, data_format.data());
    
    PetscViewerFileSetMode(viewer, FILE_MODE_READ);
    PetscViewerFileSetName(viewer, ifilename.data());

    read_matrix(A, A_name, viewer);
    read_matrix(M, M_name, viewer);
    read_matrix(T, M_name, viewer);
    read_matrix(C, C_name, viewer);

    // Clean up viewer
    PetscViewerDestroy(&viewer);
    
    // Setup eigenvalue problem solver
    hpgem::EigenSolvers::DoehlerMaxwellSolver doehler_eigensolver;
    
    doehler_eigensolver.setMatrices(A, C);
    doehler_eigensolver.setMaxIter(max_iter);
    doehler_eigensolver.setTolerance(tol);
    doehler_eigensolver.solve(n_eigenvalues, T);
    
    // Display results
    PetscMPIInt rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    
    if(rank == 0) std::cout << "\nNumber of converged eigenvalues: " << doehler_eigensolver.getConverged() << std::endl;

    for(PetscInt eigen_v_idx = 0; eigen_v_idx < n_eigenvalues; eigen_v_idx++) {
      PetscErrorCode ierr;
      Vec eigen_vector;
      MatCreateVecs(A, &eigen_vector, NULL);
      PetscScalar eigen_value;
      ierr = doehler_eigensolver.getEigenPair(eigen_v_idx, eigen_value, eigen_vector);
      if(rank == 0)
      {
        std::cout << std::setprecision(10) << std::fixed;
        std::cout << "Eigenvalue " << eigen_v_idx << ": " << eigen_value << std::endl;
      }
      VecDestroy(&eigen_vector);
    }
    
    MPI_Barrier(PETSC_COMM_WORLD);
    
    // Finalize
    return SlepcFinalize();
}



PetscErrorCode read_matrix(Mat &M, std::string &name, std::string &filename, std::string &read_format) {
  // Data viewer (for reading matrices from file)
  PetscViewer     viewer;
  PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
  PetscViewerSetType(viewer, read_format.data());
  // PetscViewerPushFormat(viewer, PETSC_VIEWER_HDF5_MAT);
  PetscViewerFileSetMode(viewer, FILE_MODE_READ);
  PetscViewerFileSetName(viewer, filename.data());

  // Read the PETSc object from viewer
  read_matrix(M, name, viewer);

  // Clean up
  PetscViewerDestroy(&viewer);

  return(0);
}


PetscErrorCode read_matrix(Mat &M, std::string &name, PetscViewer &viewer) {
  // Create the matrix to store the saved data
  MatCreate(PETSC_COMM_WORLD, &M);
  PetscObjectSetName((PetscObject)M, name.data());

  // Read the matrix from file
  MatLoad(M, viewer);

  return(0);
}

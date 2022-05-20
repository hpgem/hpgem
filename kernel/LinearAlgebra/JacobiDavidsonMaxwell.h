/*
This file forms part of hpGEM. This package has been developed over a number of
years by various people at the University of Twente and a full list of
contributors can be found at http://hpgem.org/about-the-code/team

This code is distributed using BSD 3-Clause License. A copy of which can found
below.


Copyright (c) 2018, Univesity of Twenete
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef JACOBIDAVIDSONMAXWELLSOLVER_H
#define JACOBIDAVIDSONMAXWELLSOLVER_H



#include <petsc.h>
#include <petscksp.h>

#include <slepcbv.h>
#include <slepceps.h>

#include <iostream>
#include <math.h>
#include <cstdlib>
#include <vector>

#include <chrono>  // For timing

#include <utility>
#include <valarray>


/**
 * Solver of the Maxwell equation based on the paper
 *   Multilevel preconditioned iterative eigensolvers
 *   for Maxwell eigenvalue problems
 *   P Arbenz R Geus
 *   Applied Numerical Mathematics 54 (2005) 107-121
 *
 *
 * Implementation note:
 *
 */

namespace hpgem {

namespace LinearAlgebra {

class JacobiDavidsonMaxwellSolver final {

   public:
    JacobiDavidsonMaxwellSolver();
    JacobiDavidsonMaxwellSolver(Mat &A, Mat &M, Mat &C);
    void setMatrices(const Mat &Ain, const Mat &Min, const Mat &Cin);
    PetscErrorCode solve(PetscInt nev);
    PetscInt getConverged();
    PetscErrorCode getEigenPair(PetscInt index, PetscScalar &eval, Vec &evec);
    PetscInt getIterationCount();
   private:

    void initializeMatrices();
    void initializeVectors();
    void initializeSearchSpace(int nev);

    PetscErrorCode solveMultipleLinearSystems();
    PetscErrorCode computeRayleighQuotient(const Vec &x, PetscReal *out);
    PetscErrorCode normalizeVector(Vec &x);
    PetscErrorCode VecxTAx(const Vec &x, const Mat &K, PetscScalar *val);
    PetscErrorCode getCorrectionOperator(Mat &op);
    PetscErrorCode solveCorrectionEquation(const Vec &res, Vec &sol);
    PetscErrorCode computeResidueVector(const Vec &q, const PetscReal rho, Vec &res);
    PetscErrorCode computeRightProjection(Vec &v, const BV &Q);
    PetscErrorCode computeLeftProjection(Vec &v);
    PetscErrorCode projectCorrectionVector(Vec &corr);
    PetscErrorCode addVectorToSearchSpace(const Vec &v, PetscInt idx);
    PetscErrorCode computeSmallEigenvalues(std::vector<PetscReal> &eval, Vec *evec);
    PetscErrorCode correctionOperatorMatMult(Vec x, Vec y);
    static PetscErrorCode staticMatMultCorrOp(Mat M, Vec x, Vec y);
    PetscErrorCode VecxTinvAx(const Vec &x, const Mat &K, PetscScalar *val);
    PetscErrorCode computeThreshold(Vec q, Vec r, PetscReal *eps);
    PetscErrorCode correctionPreconditionerMatMult(Vec x, Vec y);
    static PetscErrorCode staticMatMultCorrPrec(Mat M, Vec x, Vec y);
    PetscBool check_vect_nan(const Vec &vec);
    void smallev_solver(const Mat Ak, std::vector<PetscReal> &eigenvalues, Vec *eigenvectors);


    PetscReal tau = 1.2;
    PetscReal eta;
    PetscInt maxIter = 100;
    PetscInt iter = 0;
    PetscInt search_space_maxsize = 25;
    PetscInt search_space_minsize;
    PetscInt search_space_current_size = 0;
    PetscInt V_current_size = 0;
    PetscInt Q_current_size = 0;
    PetscInt Qt_current_size = 0;
    PetscInt nconverged = 0;


    Mat A, M, C;
    Mat invM;
    Mat Y, H;
    BV Qt, Q, V;
    Vec search_vect;
    Vec residue_vect;

    std::vector<PetscScalar>  eigenvalues;

    // Vec                     eigenvectors[this->search_space_maxsize];
    // Vec eigenvalues;
    // Vec *eigenvectors;

    bool print_time = false;
    bool print_small_evs = false;
};

} // namespace LinearAlgebra

}  // namespace hpgem

#endif  // JACOBIDAVIDSONMAXWELLSOLVER_H

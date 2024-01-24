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

#ifndef DOEHLERMAXWELLSOLVER_H
#define DOEHLERMAXWELLSOLVER_H

#include <vector>

#include "petsc.h"
#include "slepc.h"

namespace hpgem {

namespace EigenSolvers {

/**
 * Solver of the Maxwell equation based on the paper
 *
 * Doehler, Ein neues Gradientenverfahren zur simultanen
 * Berechnung der kleinsten oder groessten Eigenwerte des allgemeinen
 * Eigenwertproblems. Numer. Math. 40, 79--91 (1982)
 *
 *
 *  Solve
 *     A x = lambda M x
 *
 *
 * for A = A^H and with the contraint C x = 0
 *
 * Implementation note:
 *      Y = C^H (hermitian transpose of C)
 *      H = C * Y
 */
class DoehlerMaxwellSolver final {

   public:
    DoehlerMaxwellSolver();
    // PetscErrorCode clean();
    ~DoehlerMaxwellSolver();

    /**
     * Sets the matrices A and C and defaults M to the identity matrix,
     * setting up the eigenvalue problem
     *  A x = lambda x
     * constrained to
     *  C x = 0
     *
     * @param[in] Ain the A matrix to assign to the eigenvalue problem.
     * @param[in] Cin the C matrix to assign to the eigenvalue problem.
     */
    void setMatrices(const Mat Ain, const Mat Cin);

    /**
     * Sets the matrices A, C, and M, setting up the eigenvalue problem
     *  A x = lambda M x
     * constrained to
     *  C x = 0
     *
     * @param[in] Ain the A matrix to assign to the eigenvalue problem.
     * @param[in] Cin the C matrix to assign to the eigenvalue problem.
     * @param[in] Min the M matrix to assign to the eigenvalue problem,
     *            if set to NULL M is assumed the identity matrix.
     */
    void setMatrices(const Mat Ain, const Mat Cin, const Mat Min);

    /**
     * Solves the constrained eigenvalue problem
     *  A x = lambda M x
     * constrained to
     *  C x = 0
     *
     * Computes the lowest \p nev eigenvalues and associated eigenvectors.
     *
     * The computed eigenvalues are returned from getEigenPair().
     *
     * @param[in] nev the number of eigenvalues and associated eigenvectors to
     *                compute (lowest \p nenv eigenvalues).
     * @param[in] n_steps_projection the number of iterations to perform before
     *                enforcing a projection onto div(eigen_v) = 0.
     *                 <default> 10
     */
    PetscErrorCode solve(PetscInt nev, Mat &T_Mat_in,
                         PetscInt n_steps_projection = 10);

    PetscInt getConverged();
    PetscErrorCode getEigenPair(PetscInt index, PetscScalar &eval, Vec &evec);
    PetscInt getIterationCount() const { return iter; }
    void setMaxIter(int n);
    void setTolerance(PetscReal tol);

   private:
    void setupProjection();
    void cleanupProjection();
    void projectBV(BV bv);
    PetscErrorCode projectEigenVector(Vec &eigen_v);
    PetscErrorCode ritzUpdate(
        BV T_bv, PetscInt n_eigs, std::vector<PetscScalar> &L_std_vec);
    void compute_residual_eigen_v(const std::vector<PetscScalar> &ritzValues,
                                  BV &X_bv, PetscInt eigen_idx_start,
                                  PetscInt n_eigs, BV &R_bv, BV temp_bv);
    // void initializeMatrices();
    // void initializeVectors();
    // void initializeSearchSpace(int nev);
    //
    // PetscErrorCode computeRayleighQuotient(const Vec &x, PetscReal *out);
    // PetscErrorCode normalizeVector(Vec &x);
    // PetscErrorCode weightedVectorDot(const Vec &x, const Mat &K,
    //                                  PetscReal *val);
    // PetscErrorCode getCorrectionOperator(Mat &op);
    // PetscErrorCode solveCorrectionEquation(const Vec &res, Vec &sol);
    // PetscErrorCode computeResidueVector(const Vec &q, const PetscReal rho,
    //                                     Vec &res);
    //
    // PetscErrorCode gramSchmidt(Vec &v, const BV &Q);
    // PetscErrorCode modifiedGramSchmidt(Vec &v, const BV &Q);
    //

    // TODO: This is the projection function to use to enforce the constraint C
    // x = 0 PetscErrorCode projectCorrectionVector(Vec &corr);

    //
    // PetscErrorCode computeSmallEigenvalues(std::vector<PetscReal> &eval,
    //                                        Vec *evec);
    // PetscErrorCode correctionOperatorMatMult(Vec x, Vec y);
    //
    // static PetscErrorCode staticMatMultCorrOp(Mat M, Vec x, Vec y);
    // PetscErrorCode computeThreshold(Vec q, Vec r, PetscReal *eps);
    // PetscErrorCode correctionPreconditionerMatMult(Vec x, Vec y);
    // static PetscErrorCode staticMatMultCorrPrec(PC pc, Vec x, Vec y);
    //
    // PetscReal computeSmallResidue(const Mat &A, const Vec &x,
    //                               const PetscScalar lambda);
    // PetscErrorCode orderEigenvalues();
    //
    // PetscReal ev_target;
    // PetscReal eta;
    // PetscReal prec_shift;
    PetscInt maxIter;
    // PetscInt correction_niter;
    PetscInt iter = 0;
    // PetscInt search_space_maxsize;
    // PetscInt search_space_restart_size;
    // PetscInt search_space_current_size = 0;
    // PetscInt V_current_size = 0;
    //
    PetscInt eigenvectors_current_size = 0;
    // PetscInt Qt_current_size = 0;
    PetscReal tolerance;
    //

    Mat A, M, C;
    BV eigenvectors;
    std::vector<PetscScalar> eigenvalues;

    // Objects needed during the solve for the projector
    Mat Y, H;
    KSP projectionSolver_;
    Vec projectionTempVector_;
    // BV Qt;

    // BV V;
    // Vec search_vect;
    // Vec residue_vect;
    // bool print_small_evs = false;
};

}  // namespace EigenSolvers

}  // namespace hpgem

#endif  // DOEHLERMAXWELLSOLVER_H

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

#ifndef HPGEM_APP_DGMAXEIGENVALUE_H
#define HPGEM_APP_DGMAXEIGENVALUE_H

#include "ProblemTypes/AbstractEigenvalueSolver.h"

#include "DGMaxDiscretization.h"

#include "Utilities/GlobalIndexing.h"

#include <slepceps.h>

using namespace hpgem;

class DGMaxEigenvalueBase {
   public:
    enum ProjectorUse {
        /// Don't use the projector
        NONE,
        /// Only use the projector for the initial subspace
        INITIAL,
        /// Use the projector at each step
        ALL
    };

    struct SolverConfig {
        SolverConfig()
            : useHermitian_(true),
              shiftFactor_(0),
              stab_(100),
              useProjector_(NONE){};

        /// Whether to solve M^{-1}S x = omega^2 (non Hermitian) or
        /// L^{-1} S L^{-T}y = omega^2 y (Hermitian)
        bool useHermitian_;
        /// The factor 's' in shifting the basis functions
        double shiftFactor_;
        /// Stabilization parameter (will be rescaled based on facet size).
        double stab_;
        /// Use a projector to remove the kernel of the stiffness matrix
        ProjectorUse useProjector_;

        /// Whether the config uses shifts
        bool usesShifts() const {
            // Allow for arbitrary small shifts
            return shiftFactor_ != 0.0;
        }
    };
};

/**
 * Solver of the Bloch Eigenvalue problem based on DGMaxDiscretization.
 *
 * The most basic version solves the eigenvalue problem S(k) u = omega^2 Mu,
 * where u are the coefficients for the electric field, M is the mass matrix and
 * S is the stiffness matrix. The matrix S(k) depends on the k-point and is
 * Hermitian (with complex values) positive semi definite. The matrix M is an
 * mass matrix weighted by the dielectric constant, and is thus block diagonal
 * and positive definite.
 *
 * The dependence on the k-point is part of the stiffness matrix. In most basic
 * form of the problem we have the Bloch-periodic boundary condition:
 * E(x+a) = E(x) e^{ika},
 * where E is the electric field which we are solving for, x and x+a are two
 * points that are connected through the periodic boundary, i is the imaginary
 * unit and k is the kpoint for which we are solving. As result the face-matrix
 * for the faces on the periodic boundary gain a factor e^{ika} and e^{-ika} on
 * the off diagonal blocks. The implementation of DGMaxDiscretization gives the
 * correct matrix for k=0, for different k-points the matrix S needs to be
 * adjusted.
 *
 * In addition to this most basic problem there are a few extra options that can
 * be enabled in the solver:
 *  - Hermitian/Non Hermitian: To solve the problem as a standard eigenvalue
 *    problem it was originally reformulated as M^{-1}S u = omega^2 u, which is
 *    non-Hermitian. For the Hermitian variant we instead use the Cholesky
 *    decomposition of M = LL^T, and solve L^{-1} S L^{-T} y = omega^2 y, with
 *    L^T u = y.
 *  - Shifts: the idea is that the convergence between different k-points may
 *    improve if we give each basis function psi(x) a phase shift
 *    psi(x) e^{s*ik x0} where x0 is a fixed point of the basis function (we use
 *    the center of its element) and s a configurable number. Use 0 to disable
 *    it. These phase shifts result in that the stiffness matrix is changed from
 *    S(k) to D(k) S(k) D^* (k), where D(k) is a diagonal matrix with the phase
 *    shifts e^{s ik x0} on the diagonal.
 *  - Projector: The stiffness matrix S has a large kernel, corresponding to
 *    fields that are not divergence free. The projector uses that we have an
 *    analytic expression for these to remove them. For this projector an extra
 *    matrix B(k) is introduced with a non-orthogonal basis for the kernel of S,
 *    i.e. S(k) B(k) = 0. Just like S the matrix B depends on the k-point,
 *    where entries get a phase factor e^{+-ika}.
 *
 * Implementation note:
 *  - The ordering of L from the Hermitian reformulation and D from the shifts
 *    does not matter. The L matrix is block diagonal with each block
 *    corresponding to a block. The matrix D is diagonal, but as we choose the
 *    phase shifts based on the center of the element, we get that it is
 *    constant in each block. Therefore, we have L^{-1} D = D L^{-1}. This is
 *    convenient as shifting basis functions would imply that the correct order
 *    is L^{-1} D S D^* L^{-T} (shift basis functions first), while for
 *    implementation it is more practical to use the reverse order
 *    D L^{-1} S L^{-T} D^*, as this allows assembling the middle matrix and
 *    then using a diagonal rescaling.
 *
 * @tparam DIM The dimension of the problem (and thus k-vector)
 */
template <std::size_t DIM>
class DGMaxEigenvalue : public AbstractEigenvalueSolver<DIM>,
                        public DGMaxEigenvalueBase {

   public:
    DGMaxEigenvalue(Base::MeshManipulator<DIM>& mesh, std::size_t order,
                    SolverConfig config);

    void solve(AbstractEigenvalueSolverDriver<DIM>& driver) override;

   private:
    void initializeMatrices();

    Base::MeshManipulator<DIM>& mesh_;
    std::size_t order_;
    SolverConfig config_;
    DGMaxDiscretization<DIM> discretization_;

    // Implementation details, declared as inner classes to prevent name clashes
    /**
     * Workspace used to solve at a single k-point
     */
    class SolverWorkspace;
    /**
     * Extra workspace where each of the basis functions has an extra phase
     * shift exp(i akx), where x is the center of the element of the basis
     * function, k is the k-point for the current solve and a is the factor from
     * the solver config.
     */
    class ShiftWorkspace;
    /**
     * Extra workspace for a projector that removes the eigenspace corresponding
     * to the zero eigenvalues.
     */
    class ProjectorWorkspace;

    /**
     * Result class
     */
    class Result;
};

#endif  // HPGEM_APP_DGMAXEIGENVALUE_H

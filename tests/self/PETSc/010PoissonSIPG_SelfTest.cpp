/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2014, University of Twente
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

#include <cmath>

#include "Base/HpgemAPILinearSteadyState.h"
#include "petscksp.h"
#include "Utilities/GlobalMatrix.h"
#include "Utilities/GlobalVector.h"

#include "Logger.h"
#include "../TestMeshes.h"
#include "../ConvergenceTest.h"

using namespace hpgem;

/// \brief Class for solving the Poisson problem using
/// HpgemAPILinearSteadyState.
template <std::size_t DIM>
class PoissonTest : public Base::HpgemAPILinearSteadyState<DIM> {
   public:
    using typename Base::HpgemAPIBase<DIM>::PointPhysicalT;
    using typename Base::HpgemAPIBase<DIM>::PointReferenceT;
    using typename Base::HpgemAPIBase<DIM>::PointReferenceOnFaceT;

    PoissonTest(const std::string name, const std::size_t p,
                const std::size_t n)
        : Base::HpgemAPILinearSteadyState<DIM>(1, p, true, true),
          p_(p),
          totalError_(0) {
        logger(INFO, "Reading test with mesh % and p=%", name, p);
        penalty_ = 3 * n * p_ * (p_ + DIM - 1) + 1;
        this->readMesh(name);
    }

    ///\brief Compute the integrand for the stiffness matrix at the element.
    LinearAlgebra::MiddleSizeMatrix computeIntegrandStiffnessMatrixAtElement(
        Base::PhysicalElement<DIM>& element) final {
        // Obtain the number of basisfunctions that are possibly non-zero on
        // this element.
        const std::size_t numberOfBasisFunctions =
            element.getElement()->getNumberOfBasisFunctions();

        // Create the integrandVal such that it contains as many rows and
        // columns as the number of basisfunctions.
        LinearAlgebra::MiddleSizeMatrix& integrandVal =
            element.getResultMatrix();

        for (std::size_t i = 0; i < numberOfBasisFunctions; ++i) {
            for (std::size_t j = 0; j < numberOfBasisFunctions; ++j) {
                // Compute the value of gradient(phi_i).gradient(phi_j) at point
                // p and store it at the appropriate place in the matrix
                // integrandVal.
                integrandVal(i, j) = element.basisFunctionDeriv(i) *
                                     element.basisFunctionDeriv(j);
            }
        }

        return integrandVal;
    }

    /// \brief Compute the integrand for the siffness matrix at the face.
    Base::FaceMatrix computeIntegrandStiffnessMatrixAtFace(
        Base::PhysicalFace<DIM>& face) final {
        // Get the number of basis functions, first of both sides of the face
        // and then only the basis functions associated with the left and right
        // element.
        std::size_t numberOfBasisFunctions =
            face.getFace()->getNumberOfBasisFunctions();

        // Create the FaceMatrix integrandVal with the correct size.
        Base::FaceMatrix& integrandVal = face.getResultMatrix();

        // Initialize the vectors that contain gradient(phi_i), gradient(phi_j),
        // normal_i phi_i and normal_j phi_j
        LinearAlgebra::SmallVector<DIM> phiNormalI, phiNormalJ, phiDerivI,
            phiDerivJ;

        // Transform the point from the reference value to its physical value.
        // This is necessary to check at which boundary we are if we are at a
        // boundary face.
        const PointPhysicalT& pPhys = face.getPointPhysical();

        for (std::size_t i = 0; i < numberOfBasisFunctions; ++i) {
            // normal_i phi_i is computed at point p, the result is stored in
            // phiNormalI.
            phiNormalI = face.basisFunctionUnitNormal(i);
            // The gradient of basisfunction phi_i is computed at point p, the
            // result is stored in phiDerivI.
            phiDerivI = face.basisFunctionDeriv(i);

            for (std::size_t j = 0; j < numberOfBasisFunctions; ++j) {
                // normal_j phi_j is computed at point p, the result is stored
                // in phiNormalJ.
                phiNormalJ = face.basisFunctionUnitNormal(j);
                // The gradient of basisfunction phi_j is computed at point p,
                // the result is stored in phiDerivJ.
                phiDerivJ = face.basisFunctionDeriv(j);

                // Switch to the correct type of face, and compute the integrand
                // accordingly you could also compute the integrandVal by
                // directly using face->basisFunctionDeriv and
                // face->basisFunctionNormal in the following lines, but this
                // results in very long expressions Internal face:
                if (face.isInternal()) {
                    integrandVal(j, i) =
                        -(phiNormalI * phiDerivJ + phiNormalJ * phiDerivI) / 2 +
                        penalty_ * phiNormalI * phiNormalJ;
                }
                // Boundary face with Dirichlet boundary conditions:
                else if (std::abs(pPhys[0]) < 1e-9 ||
                         std::abs(pPhys[0] - 1.) < 1e-9) {
                    integrandVal(j, i) =
                        -(phiNormalI * phiDerivJ + phiNormalJ * phiDerivI) +
                        penalty_ * phiNormalI * phiNormalJ * 2;
                }
                // Boundary face with homogeneous Neumann boundary conditions:
                else {
                    integrandVal(j, i) = 0;
                }
            }
        }

        return integrandVal;
    }

    /// \brief Define the exact solution
    LinearAlgebra::MiddleSizeVector getExactSolution(
        const PointPhysicalT& p) final {
        LinearAlgebra::MiddleSizeVector exactSolution(1);

        double ret = std::sin(M_PI * p[0]);
        if (p.size() > 1) {
            ret *= std::cos(M_PI * p[1]);
        }
        if (p.size() > 2) {
            ret *= std::cos(M_PI * p[2]);
        }

        exactSolution[0] = ret;
        return exactSolution;
    }

    ///\brief Define the source term.
    LinearAlgebra::MiddleSizeVector getSourceTerm(
        const PointPhysicalT& p) final {
        LinearAlgebra::MiddleSizeVector sourceTerm(1);

        double ret = -std::sin(M_PI * p[0]) * (M_PI * M_PI);
        if (DIM > 1) {
            ret *= std::cos(M_PI * p[1]) * 2;  // 2 pi^2 prefactor
        }
        if (DIM > 2) {
            ret *= std::cos(M_PI * p[2]) * 1.5;  // 3 pi^2 prefactor
        }

        sourceTerm[0] = ret;
        return sourceTerm;
    }

    /// \brief Compute the integrals of the right-hand side associated with
    /// faces.
    LinearAlgebra::MiddleSizeVector computeIntegrandSourceTermAtFace(
        Base::PhysicalFace<DIM>& face) final {
        // Obtain the number of basisfunctions that are possibly non-zero
        const std::size_t numberOfBasisFunctions =
            face.getFace()->getNumberOfBasisFunctions();
        // Resize the integrandVal such that it contains as many rows as
        // the number of basisfunctions.
        LinearAlgebra::MiddleSizeVector& integrandVal = face.getResultVector();

        // Compute the value of the integrand
        // We have no rhs face integrals, so this is just 0.
        for (std::size_t i = 0; i < numberOfBasisFunctions; ++i) {
            integrandVal[i] = 0;
        }

        return integrandVal;
    }

    void solveSteadyStateWithPetsc(bool doComputeError) final {
#if defined(HPGEM_USE_ANY_PETSC)
        // Create and Store things before solving the problem.
        this->tasksBeforeSolving();

        // Solve the linear problem
        Utilities::GlobalIndexing indexing(this->meshes_[0]);
        // Assemble the matrix A of the system Ax = b.
        Utilities::GlobalPetscMatrix A(indexing,
                                       this->stiffnessElementMatrixID_,
                                       this->stiffnessFaceMatrixID_);
        MatScale(A, -1);
        // Declare the vectors x and b of the system Ax = b.
        Utilities::GlobalPetscVector b(indexing, this->sourceElementVectorID_,
                                       this->sourceFaceVectorID_),
            x(indexing);

        // Assemble the vector b. This is needed because Petsc assumes you don't
        // know yet whether a vector is a variable or right-hand side the moment
        // it is declared.
        b.assemble();

        // Make a solver, by default solve using a direct method. This to
        // prevent spurious errors from the solver.
        KSP ksp;
        KSPCreate(PETSC_COMM_WORLD, &ksp);

        KSPSetType(ksp, KSPPREONLY);
        PC pc;
        KSPGetPC(ksp, &pc);
        PCSetType(pc, PCLU);

        KSPSetTolerances(ksp, 1e-12, PETSC_DEFAULT, PETSC_DEFAULT,
                         PETSC_DEFAULT);
        // Tell ksp that it will solve the system Ax = b.
        KSPSetOperators(ksp, A, A);
        KSPSetFromOptions(ksp);
        KSPSolve(ksp, b, x);
        KSPType type;
        KSPGetType(ksp, &type);
        logger(INFO, "Solving using %", type);

        // Do PETSc magic, including solving.
        KSPConvergedReason converge;
        KSPGetConvergedReason(ksp, &converge);
        int iterations;
        KSPGetIterationNumber(ksp, &iterations);
        logger(INFO, "KSP solver ended because of % in % iterations.",
               KSPConvergedReasons[converge], iterations);

        x.writeTimeIntegrationVector(this->solutionVectorId_);

        if (doComputeError) {
            LinearAlgebra::MiddleSizeVector::type totalError =
                this->computeTotalError(this->solutionVectorId_, 0);
            totalError_ = totalError;
            logger(INFO, "Total error: %.", totalError);
            LinearAlgebra::MiddleSizeVector maxError =
                this->computeMaxError(this->solutionVectorId_, 0);
            logger.assert_debug(
                maxError.size() == this->configData_->numberOfUnknowns_,
                "Size of maxError (%) not equal to the number of variables (%)",
                maxError.size(), this->configData_->numberOfUnknowns_);
            for (std::size_t iV = 0; iV < this->configData_->numberOfUnknowns_;
                 iV++) {
                logger(INFO, "Maximum error %: %", this->variableNames_[iV],
                       maxError(iV));
            }
        }

        return;
#endif
    }

    LinearAlgebra::MiddleSizeVector::type getTotalError() {
        return this->totalError_;
    }

   private:
    /// polynomial order of the approximation
    int p_;

    ///\brief Penalty parameter
    ///
    /// Penalty parameter that is associated with the interior penalty
    /// discontinuous Galerkin method. This parameter is initialized in the
    /// constructor, and has to be greater than 3 * n_ * p_ * (p_ + DIM - 1) in
    /// order for the method to be stable.
    double penalty_;

    /// Weighted L2 norm of the error
    LinearAlgebra::MiddleSizeVector::type totalError_;
};

struct PoissonTestParameters {
    std::size_t p;  // Polynomial order
    std::size_t n;  // Base number of Segments
};

template <std::size_t DIM>
void runPoissonTestSeries(ConvergenceTestSet& testSet,
                          PoissonTestParameters& testParameters,
                          bool ignoreErrors) {
    runConvergenceTest(
        testSet, ignoreErrors,
        [&testParameters](std::string meshName, std::size_t level) {
            std::size_t n =
                testParameters.n * (1 << level);  // Simple way to compute 2^N
            PoissonTest<DIM> test(meshName, testParameters.p, n);
            test.solveSteadyStateWithPetsc(true);
            return std::real(test.getTotalError());
        });
}

int main(int argc, char** argv) {
    Base::parse_options(argc, argv);

    /*
     * Test solving a simple Poisson problem using SIPG on several meshes that
     * form refinements. As this gives numerical results the failure does not
     * necessarily mean a bug, but merely that the convergence of the results
     * should be checked and the expectations may need to be updated.
     *
     * The tests compute the L2 norm of the error, so expect a convergence rate
     * of 2^{p+1}. Note the meshes start at a single element, so the first 2-3
     * convergence rates may be quite off.
     */

    // For regenerating the errors
    bool ignoreErrors = false;

    PoissonTestParameters dim1P2Params = {2, 1};
    ConvergenceTestSet dim1P2Meshes = {getUnitSegmentMeshes(),
                                       {
                                           7.45240563e-03,  //------
                                           1.25626260e-02,  //  0.59
                                           1.45725197e-03,  //  8.62
                                           1.77601092e-04,  //  8.21
                                           2.20303522e-05,  //  8.06
                                           2.74709611e-06,  //  8.02
                                       }};
    runPoissonTestSeries<1>(dim1P2Meshes, dim1P2Params, ignoreErrors);

    PoissonTestParameters dim1P4Params = {4, 1};
    ConvergenceTestSet dim1P4Meshes = {getUnitSegmentMeshes(),
                                       {
                                           2.95738208e-04,  //------
                                           9.80215623e-05,  //  3.02
                                           2.99105904e-06,  // 32.77
                                           9.27117984e-08,  // 32.26
                                           2.89065244e-09,  // 32.07
                                           9.02782781e-11,  // 32.02
                                       }};
    runPoissonTestSeries<1>(dim1P4Meshes, dim1P4Params, ignoreErrors);

    PoissonTestParameters dim2P2Params = {2, 1};
    ConvergenceTestSet dim2P2Meshes = {getUnitSquareTriangleMeshes(),
                                       {
                                           2.16458877e-01,  //------
                                           3.73833650e-02,  //  5.79
                                           3.65233994e-03,  // 10.24
                                           4.66915221e-04,  //  7.82
                                           5.86899243e-05,  //  7.96
                                           7.34507877e-06,  //  7.99
                                           9.18317173e-07,  //  8.00
                                       }};
    runPoissonTestSeries<2>(dim2P2Meshes, dim2P2Params, ignoreErrors);

    // Use only p=1 to reduce computational time for this test
    PoissonTestParameters dim3P1Params{1, 1};
    ConvergenceTestSet dim3P1Meshes{getUnitCubeCubeMeshes(),
                                    {
                                        2.08032213e-01,  //------
                                        8.26794259e-02,  //  2.52
                                        2.22900949e-02,  //  3.71
                                        5.70074929e-03,  //  3.91
                                        1.43383386e-03,  //  3.98

                                    }};
    runPoissonTestSeries<3>(dim3P1Meshes, dim3P1Params, ignoreErrors);

    return 0;
}

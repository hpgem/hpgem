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

//#define hpGEM_INCLUDE_PETSC_SUPPORT//temporarily activating this definition
// makes development easier on some IDEs

#include <cmath>

#include "Base/HpgemAPILinearSteadyState.h"
#include "Output/TecplotDiscontinuousSolutionWriter.h"
#include "petscksp.h"
#include "Utilities/GlobalMatrix.h"
#include "Utilities/GlobalVector.h"
#include "../TestMeshes.h"
using namespace hpgem;
// If this test ever breaks it is not a bad thing per se.
// If the results are still readable by tecplot, and you are convinced that your
// changes improved the code, you should update the data file to reflect the
// updated result. Always confer with other developers if you do this.

/// \brief Class for solving the Poisson problem using
/// HpgemAPILinearSteadyState.
class PoissonTest : public Base::HpgemAPILinearSteadyState<2> {
   public:
    PoissonTest(const std::string name, const std::size_t p,
                const std::size_t n)
        : HpgemAPILinearSteadyState(1, p, true, true), p_(p), totalError_(0) {
        penalty_ = 3 * n * p_ * (p_ + 2 - 1) + 1;
        readMesh(name);
    }

    ///\brief Compute the integrand for the stiffness matrix at the element.
    LinearAlgebra::MiddleSizeMatrix computeIntegrandStiffnessMatrixAtElement(
        Base::PhysicalElement<2>& element) final {
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
                integrandVal(j, i) = element.basisFunctionDeriv(i) *
                                     element.basisFunctionDeriv(j);
            }
        }

        return integrandVal;
    }

    /// \brief Compute the integrand for the siffness matrix at the face.
    Base::FaceMatrix computeIntegrandStiffnessMatrixAtFace(
        Base::PhysicalFace<2>& face) final {
        // Get the number of basis functions, first of both sides of the face
        // and then only the basis functions associated with the left and right
        // element.
        std::size_t numberOfBasisFunctions =
            face.getFace()->getNumberOfBasisFunctions();

        // Create the FaceMatrix integrandVal with the correct size.
        Base::FaceMatrix& integrandVal = face.getResultMatrix();

        // Initialize the vectors that contain gradient(phi_i), gradient(phi_j),
        // normal_i phi_i and normal_j phi_j
        LinearAlgebra::SmallVector<2> phiNormalI, phiNormalJ, phiDerivI,
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

        double ret = std::sin(2 * M_PI * p[0]);
        if (p.size() > 1) {
            ret *= std::cos(2 * M_PI * p[1]) / 2.;
        }
        if (p.size() > 2) {
            ret *= std::cos(2 * M_PI * p[2]) * 2.;
        }

        exactSolution[0] = ret;
        return exactSolution;
    }

    ///\brief Define the source term.
    LinearAlgebra::MiddleSizeVector getSourceTerm(
        const PointPhysicalT& p) final {
        LinearAlgebra::MiddleSizeVector sourceTerm(1);

        double ret = -std::sin(2 * M_PI * p[0]) * (4 * M_PI * M_PI);
        if (2 > 1) {
            ret *= std::cos(2 * M_PI * p[1]);
        }
        if (2 > 2) {
            ret *= std::cos(2 * M_PI * p[2]) * 3;
        }

        sourceTerm[0] = ret;
        return sourceTerm;
    }

    /// \brief Compute the integrals of the right-hand side associated with
    /// faces.
    LinearAlgebra::MiddleSizeVector computeIntegrandSourceTermAtFace(
        Base::PhysicalFace<2>& face) final {
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
        tasksBeforeSolving();

        // Solve the linear problem
        Utilities::GlobalIndexing indexing(HpgemAPIBase::meshes_[0]);
        // Assemble the matrix A of the system Ax = b.
        Utilities::GlobalPetscMatrix A(indexing, stiffnessElementMatrixID_,
                                       stiffnessFaceMatrixID_);
        MatScale(A, -1);
        // Declare the vectors x and b of the system Ax = b.
        Utilities::GlobalPetscVector b(indexing, sourceElementVectorID_,
                                       sourceFaceVectorID_),
            x(indexing);

        // Assemble the vector b. This is needed because Petsc assumes you don't
        // know yet whether a vector is a variable or right-hand side the moment
        // it is declared.
        b.assemble();

        // Make the Krylov supspace method
        KSP ksp;
        KSPCreate(PETSC_COMM_WORLD, &ksp);
        KSPSetTolerances(ksp, 1e-12, PETSC_DEFAULT, PETSC_DEFAULT,
                         PETSC_DEFAULT);
        // Tell ksp that it will solve the system Ax = b.
        KSPSetOperators(ksp, A, A);
        KSPSetFromOptions(ksp);
        KSPSolve(ksp, b, x);
        // Do PETSc magic, including solving.
        KSPConvergedReason converge;
        KSPGetConvergedReason(ksp, &converge);
        int iterations;
        KSPGetIterationNumber(ksp, &iterations);
        logger(INFO, "KSP solver ended because of % in % iterations.",
               KSPConvergedReasons[converge], iterations);

        x.writeTimeIntegrationVector(solutionVectorId_);

        std::ofstream outFile("030TecplotOutput_SelfTest_output.dat");
        Output::TecplotDiscontinuousSolutionWriter<2> writeFunc(outFile, "test",
                                                                "01", "value");
        writeFunc.write(meshes_[0], "monomial solution", false, this);

        if (doComputeError) {
            LinearAlgebra::MiddleSizeVector::type totalError =
                computeTotalError(solutionVectorId_, 0);
            totalError_ = totalError;
            logger(INFO, "Total error: %.", totalError);
            LinearAlgebra::MiddleSizeVector maxError =
                computeMaxError(solutionVectorId_, 0);
            logger.assert_debug(
                maxError.size() == configData_->numberOfUnknowns_,
                "Size of maxError (%) not equal to the number of variables (%)",
                maxError.size(), configData_->numberOfUnknowns_);
            for (std::size_t iV = 0; iV < configData_->numberOfUnknowns_;
                 iV++) {
                logger(INFO, "Maximum error %: %", variableNames_[iV],
                       maxError(iV));
            }
        }

        return;
#endif
    }

    LinearAlgebra::MiddleSizeVector::type getTotalError() {
        return totalError_;
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

int main(int argc, char** argv) {
    Base::parse_options(argc, argv);

    PoissonTest test8(getUnitSquareTriangleMeshes()[3], 5, 8);
    test8.solveSteadyStateWithPetsc(false);
    // actual test is done by comparing output files
    return 0;
}

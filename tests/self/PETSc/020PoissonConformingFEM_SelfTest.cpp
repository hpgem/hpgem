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

#include "../ConvergenceTest.h"
#include "../TestMeshes.h"

#include "Logger.h"

using namespace hpgem;
// If this test ever breaks it is not a bad thing per se. However, once this
// breaks a thorough convergence analysis needs to be done. If the results still
// show the theoretically optimal order of convergence, and you are convinced
// that your changes improved the code, you should change the numbers in this
// test to reflect the updated result. Always confer with other developers if
// you do this.

/// \brief Class for solving the Poisson problem using HpgemAPILinearSteadyState
/// and conforming basis functions.
template <std::size_t DIM>
class PoissonTest : public Base::HpgemAPILinearSteadyState<DIM> {
   public:
    using typename Base::HpgemAPIBase<DIM>::PointPhysicalT;
    using typename Base::HpgemAPIBase<DIM>::PointReferenceT;
    using typename Base::HpgemAPIBase<DIM>::PointReferenceOnFaceT;

    PoissonTest(const std::string name, const std::size_t p)
        : Base::HpgemAPILinearSteadyState<DIM>(1, p, true, true),
          p_(p),
          totalError_(0) {
        using namespace std::string_literals;
        readMesh(name);
    }

    void readMesh(const std::string meshName) final {

        // Set the number of Element/Face Matrices/Vectors.
        std::size_t numberOfElementMatrices =
            2;  // Mass matrix and stiffness matrix
        std::size_t numberOfElementVectors = 1;  // Source term vector
        std::size_t numberOfFaceMatrices = 1;    // Stiffness matrix
        std::size_t numberOfFaceVectors = 1;  // Source term vector at boundary

        // Create mesh and set basis functions.
        this->addMesh(meshName, numberOfElementMatrices, numberOfElementVectors,
                      numberOfFaceMatrices, numberOfFaceVectors);
        this->meshes_[0]->useDefaultConformingBasisFunctions(p_);

        // Set the number of time integration vectors according to the size of
        // the Butcher tableau.
        this->setNumberOfTimeIntegrationVectorsGlobally(
            this->globalNumberOfTimeIntegrationVectors_);

        // Plot info about the mesh
        std::size_t nElements = this->meshes_[0]->getNumberOfElements();
        logger(VERBOSE, "Total number of elements: %", nElements);
    }

    ///\brief Compute the integrand for the stiffness matrix at the element.
    LinearAlgebra::MiddleSizeMatrix computeIntegrandStiffnessMatrixAtElement(
        Base::PhysicalElement<DIM>& element) final {
        // Obtain the number of basisfunctions that are possibly non-zero on
        // this element.
        const std::size_t numberOfBasisFunctions =
            element.getElement()->getNrOfBasisFunctions();

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

    // the default hpGEM solver expects to have to construct a face matrix, just
    // give it the default one
    Base::FaceMatrix computeIntegrandStiffnessMatrixAtFace(
        Base::PhysicalFace<DIM>& face) final {
        return face.getResultMatrix();
    }

    // the default hpGEM solver expects to have to construct a face vector, just
    // give it the default one
    LinearAlgebra::MiddleSizeVector computeIntegrandSourceTermAtFace(
        Base::PhysicalFace<DIM>& face) final {
        return face.getResultVector();
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
            ret *= std::cos(M_PI * p[1]) * 2;
        }
        if (DIM > 2) {
            ret *= std::cos(M_PI * p[2]) * 1.5;
        }

        sourceTerm[0] = ret;
        return sourceTerm;
    }

    // This routine alters the matrix such that it can deal with conforming
    // boundaries. It assumes correct boundary values are provided in its third
    // argument (the rest of the vector can be garbage) and that the second
    // vector will be used as the RHS of a linear system solve it clears the
    // rows corresponding to the boudary nodes to have only a 1 on the diagonal
    // and sets the RHS to the appropriate value
    void insertDirichletBoundary(Utilities::GlobalPetscMatrix& A,
                                 Utilities::GlobalPetscVector& b,
                                 Utilities::GlobalPetscVector& x) {
        std::size_t numberOfRows(0);
        std::vector<int> rows(0);
        Geometry::PointPhysical<DIM> pPhys;
        for (Base::Face* face : this->meshes_[0]->getFacesList()) {
            const PointReferenceOnFaceT& center =
                face->getReferenceGeometry()->getCenter();
            pPhys = face->referenceToPhysical(center);
            // if the face is on a dirichlet boundary
            // if(face->faceType_=(...))
            if (std::abs(pPhys[0]) < 1e-9 || std::abs(pPhys[0] - 1) < 1e-9) {
                // fetch the row numbers
                A.getMatrixBCEntries(face, numberOfRows, rows);
            }
        }
        int ierr = MatZeroRows(A, numberOfRows, &rows[0], 1.0, x, b);
        CHKERRV(ierr);
    }

    void solveSteadyStateWithPetsc(bool doComputeError) final {
#if defined(HPGEM_USE_ANY_PETSC)
        // Create and Store things before solving the problem.
        this->tasksBeforeSolving();

        // Solve the linear problem
        Utilities::GlobalIndexing indexing(this->meshes_[0]);
        // Assemble the matrix A of the system Ax = b.
        // The special value -1 is used to indicate there is no face matrix
        // (since there is no flux in the conforming case)
        Utilities::GlobalPetscMatrix A(indexing,
                                       this->stiffnessElementMatrixID_, -1);
        MatScale(A, -1);
        // Declare the vectors x and b of the system Ax = b.
        // The special value -1 is used to indicate there is no face matrix
        // (since there is no flux in the conforming case)
        Utilities::GlobalPetscVector b(indexing, this->sourceElementVectorID_,
                                       -1),
            x(indexing);

        // Assemble the vector b. This is needed because Petsc assumes you don't
        // know yet whether a vector is a variable or right-hand side the moment
        // it is declared.
        b.assemble();
        VecSet(x, 0);
        insertDirichletBoundary(A, b, x);

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
        return totalError_;
    }

   private:
    /// number of elements per cardinal direction
    int n_;

    /// polynomial order of the approximation
    int p_;

    /// Weighted L2 norm of the error
    LinearAlgebra::MiddleSizeVector::type totalError_;
};

template <std::size_t DIM>
void runConformingPoissonTest(ConvergenceTestSet& testSet, std::size_t p,
                              bool ignoreErrors) {
    runConvergenceTest(testSet, ignoreErrors,
                       [&p](std::string meshName, std::size_t level) {
                           PoissonTest<DIM> test(meshName, p);
                           test.solveSteadyStateWithPetsc(true);
                           return std::real(test.getTotalError());
                       });
}

int main(int argc, char** argv) {

    /*
     * Test solving a simple Poisson problem using conforming basis functions on
     * several meshes that form refinements. As this gives numerical results the
     * failure does not necessarily mean a bug, but merely that the convergence
     * of the results should be checked and the expectations may need to be
     * updated.
     *
     * The tests compute the L2 norm of the error, so expect a convergence rate
     * of 2^{p+1}. Note the meshes start at a single element, so the first 2-3
     * convergence rates may be quite off.
     */

    Base::parse_options(argc, argv);

    bool ignoreErrors = false;

    ConvergenceTestSet dim1P2 = {getUnitSegmentMeshes(),
                                 {
                                     2.03319622e-02,  //------
                                     1.52069546e-02,  //  1.34
                                     1.95232226e-03,  //  7.79
                                     2.45693533e-04,  //  7.95
                                     3.07637057e-05,  //  7.99
                                     3.84709137e-06,  //  8.00
                                 }};
    runConformingPoissonTest<1>(dim1P2, 2, ignoreErrors);

    ConvergenceTestSet dim1P4 = {getUnitSegmentMeshes(),
                                 {
                                     4.19454427e-04,  //------
                                     1.05581411e-04,  //  3.97
                                     3.35867926e-06,  // 31.44
                                     1.05426524e-07,  // 31.86
                                     3.29824298e-09,  // 31.96
                                     1.03098789e-10,  // 31.99
                                 }};
    runConformingPoissonTest<1>(dim1P4, 4, ignoreErrors);

    ConvergenceTestSet dim2P2 = {getUnitSquareTriangleMeshes(),
                                 {
                                     2.59545045e-01,  //------
                                     4.48506825e-02,  //  5.79
                                     4.83286994e-03,  //  9.28
                                     6.39494399e-04,  //  7.56
                                     8.12241571e-05,  //  7.87
                                     1.01948412e-05,  //  7.97
                                     1.27567489e-06,  //  7.99
                                 }};
    runConformingPoissonTest<2>(dim2P2, 2, ignoreErrors);

    ConvergenceTestSet dim3P1 = {getUnitCubeCubeMeshes(),
                                 {
                                     3.49473341e-01,  //------
                                     9.30100954e-02,  //  3.76
                                     2.29923181e-02,  //  4.05
                                     5.74619219e-03,  //  4.00
                                     1.43671073e-03,  //  4.00

                                 }};
    runConformingPoissonTest<3>(dim3P1, 1, ignoreErrors);

    // For 3D one needs 4-th order basis functions to have basis functions on
    // every part (nodes, edges, etc.). To reduce computational time we do not
    // use all available meshes.
    {
        ConvergenceTestSet dim3P4{getUnitCubeCubeMeshes(0, 3),
                                  {
                                      2.94250804e-03,  //------
                                      8.97254213e-05,  // 32.79
                                      2.89367228e-06,  // 31.01

                                  }};
        runConformingPoissonTest<3>(dim3P4, 4, ignoreErrors);
    }

    return 0;
}

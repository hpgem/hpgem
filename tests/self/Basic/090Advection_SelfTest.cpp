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
#include <functional>
#include <chrono>
#include <CMakeDefinitions.h>

#include "Base/CommandLineOptions.h"
#include "Base/ConfigurationData.h"
#include "Base/Element.h"
#include "Base/Face.h"
#include "Base/HpgemAPILinear.h"
#include "Integration/ElementIntegral.h"
#include "Integration/FaceIntegral.h"
#include "Output/TecplotDiscontinuousSolutionWriter.h"
#include "Output/TecplotSingleElementWriter.h"
#include "Logger.h"

#include "../ConvergenceTest.h"
#include "../TestMeshes.h"

// This test case is similar to that of 080AdvectionLinear_SelfTest. However
// there are two differences:
//  1. Periodic meshes are used
//  2. The local stiffness matrix is not stored, but instead the matrix vector
//  product is computed from the basis functions.
// Both choices are to diversify the number of test paths.
//
// Note that the approach of not storing matrices is common for non
// linear problems and problems with time dependent coefficients (and thus
// matrices), where such matrices would need to be updated at each timestep. To
// keep the test case simple the advection speed is constant, allowing
// comparison with the results of 080AdvectionLinear_SelfTest.

/// This class is used to test if the advection equation is solved correctly
/// when using HpgemAPISimplified. Linear advection equation: du/dt + a[0] du/dx
/// + a[1] du/dy = 0, or equivalently: du/dt = - a[0] du/dx - a[1] du/dy = 0.
using namespace hpgem;
template <std::size_t DIM>
class Advection : public Base::HpgemAPISimplified<DIM> {
   public:
    using typename Base::HpgemAPIBase<DIM>::PointPhysicalT;
    using typename Base::HpgemAPIBase<DIM>::PointReferenceT;
    using typename Base::HpgemAPIBase<DIM>::PointReferenceOnFaceT;

    /// Constructor. Assign all private variables.
    Advection(const std::string fileType, const std::size_t p)
        : Base::HpgemAPISimplified<DIM>(1, p), fileType(fileType) {
        for (std::size_t i = 0; i < DIM; ++i) {
            a[i] = 0.1 + 0.1 * i;
        }
    }

    /// \brief Compute the integrand of the right-hand side associated with
    /// elements.
    LinearAlgebra::MiddleSizeVector computeIntegrandRightHandSideAtElement(
        Base::PhysicalElement<DIM> &element,
        const LinearAlgebra::MiddleSizeVector &inputFunctionCoefficients,
        const double time) {
        std::size_t numberOfBasisFunctions =
            element.getElement()->getNumberOfBasisFunctions();
        LinearAlgebra::MiddleSizeVector &result = element.getResultVector();
        LinearAlgebra::MiddleSizeVector::type functionValue = 0;
        for (std::size_t j = 0; j < numberOfBasisFunctions; ++j) {
            functionValue +=
                inputFunctionCoefficients(j) * element.basisFunction(j);
        }
        for (std::size_t i = 0; i < numberOfBasisFunctions; ++i) {
            result(i) = functionValue * (a * element.basisFunctionDeriv(i));
        }

        return result;
    }

    /// \brief Compute the integrals of the right-hand side associated with
    /// elements.
    LinearAlgebra::MiddleSizeVector computeRightHandSideAtElement(
        Base::Element *ptrElement,
        const LinearAlgebra::MiddleSizeVector &inputFunctionCoefficients,
        const double time) final {
        // Define a function for the integrand of the right hand side at the
        // element.
        std::function<LinearAlgebra::MiddleSizeVector(
            Base::PhysicalElement<DIM> &)>
            integrandFunction = [=](Base::PhysicalElement<DIM> &element)
            -> LinearAlgebra::MiddleSizeVector {
            return this->computeIntegrandRightHandSideAtElement(
                element, inputFunctionCoefficients, time);
        };

        return this->elementIntegrator_.integrate(ptrElement,
                                                  integrandFunction);
    }

    /// \brief Compute the integrals of the right-hand side associated with
    /// faces.
    LinearAlgebra::MiddleSizeVector computeIntegrandRightHandSideAtFace(
        Base::PhysicalFace<DIM> &face, const Base::Side iSide,
        const LinearAlgebra::MiddleSizeVector &inputFunctionCoefficientsLeft,
        const LinearAlgebra::MiddleSizeVector &inputFunctionCoefficientsRight,
        const double time) {
        // Get the number of basis functions of the elements at both sides.
        std::size_t numberOfTestFunctions =
            face.getFace()->getPtrElement(iSide)->getNumberOfBasisFunctions();
        std::size_t numberOfBasisFunctionsLeft =
            face.getFace()->getPtrElementLeft()->getNumberOfBasisFunctions();
        std::size_t numberOfBasisFunctionsRight =
            face.getFace()->getPtrElementRight()->getNumberOfBasisFunctions();

        // Resize the result to the correct size and set all elements to 0.
        LinearAlgebra::MiddleSizeVector integrandVal =
            face.getResultVector(iSide);
        integrandVal *= 0;

        // Check if the outward pointing normal vector of the left element is in
        // the same direction as the advection term.
        const double A = a * face.getUnitNormalVector();

        // Compute the sign of the normal vector (1 if iSide is left, -1 if
        // iSide is right)
        int iSign = 1;
        if (iSide == Base::Side::RIGHT) {
            iSign *= -1;
        }

        // Compute the value of the jump times the advection term
        LinearAlgebra::MiddleSizeVector::type jump = 0;
        // Advection in the same direction as outward normal of left element:
        if (A > 1e-12) {
            for (std::size_t j = 0; j < numberOfBasisFunctionsLeft; ++j) {
                jump += inputFunctionCoefficientsLeft(j) *
                        face.basisFunction(Base::Side::LEFT, j);
            }
            jump *= A * iSign;
        }
        // Advection in the same direction as outward normal of right element:
        else if (A < -1e-12) {
            for (std::size_t j = 0; j < numberOfBasisFunctionsRight; ++j) {
                jump += inputFunctionCoefficientsRight(j) *
                        face.basisFunction(Base::Side::RIGHT, j);
            }
            jump *= A * iSign;
        }
        // Advection orthogonal to normal:
        else if (std::abs(A) < 1e-12) {
            for (std::size_t j = 0; j < numberOfBasisFunctionsLeft; ++j) {
                jump += inputFunctionCoefficientsLeft(j) *
                        face.basisFunction(Base::Side::LEFT, j);
            }
            for (std::size_t j = 0; j < numberOfBasisFunctionsRight; ++j) {
                jump += inputFunctionCoefficientsRight(j) *
                        face.basisFunction(Base::Side::RIGHT, j);
            }
            jump *= A * iSign / 2;
        }

        // Compute all entries of the integrand at this point:
        for (std::size_t i = 0; i < numberOfTestFunctions; ++i) {
            integrandVal(i) = -jump * face.basisFunction(iSide, i);
        }

        return integrandVal;
    }

    /// \brief Compute the integrals of the right-hand side associated with
    /// faces.
    LinearAlgebra::MiddleSizeVector computeRightHandSideAtFace(
        Base::Face *ptrFace, const Base::Side iSide,
        LinearAlgebra::MiddleSizeVector &inputFunctionCoefficientsLeft,
        LinearAlgebra::MiddleSizeVector &inputFunctionCoefficientsRight,
        const double time) final {
        // Define a function for the integrand of the right hand side at the
        // face.
        std::function<LinearAlgebra::MiddleSizeVector(
            Base::PhysicalFace<DIM> &)>
            integrandFunction = [=](Base::PhysicalFace<DIM> &face)
            -> LinearAlgebra::MiddleSizeVector {
            return this->computeIntegrandRightHandSideAtFace(
                face, iSide, inputFunctionCoefficientsLeft,
                inputFunctionCoefficientsRight, time);
        };

        return this->faceIntegrator_.integrate(ptrFace, integrandFunction);
    }

    /// Define a solution at time zero.
    double getSolutionAtTimeZero(const PointPhysicalT &point) {
        double solution;
        solution = std::sin(2 * M_PI * point[0]);
        for (std::size_t i = 1; i < DIM; i++) {
            solution *= std::sin(2 * M_PI * point[i]);
        }
        return solution;
    }

    /// Define the exact solution. In this case that is \f$
    /// u_0(\vec{x}-\vec{a}t) \f$, where \f$ u_0 \f$ is the solution at time
    /// zero.
    LinearAlgebra::MiddleSizeVector getExactSolution(
        const PointPhysicalT &point, const double &time,
        const std::size_t orderTimeDerivative) final {
        LinearAlgebra::MiddleSizeVector exactSolution(1);
        if (orderTimeDerivative == 0) {
            PointPhysicalT displacement{a * time};
            exactSolution(0) = getSolutionAtTimeZero(point - displacement);
            return exactSolution;
        }
        logger(ERROR,
               "No exact solution for order time derivative % implemented",
               orderTimeDerivative);
        exactSolution(0) = 0;
        return exactSolution;
    }

    /// Define the initial conditions. In this case it is just the exact
    /// solution at the start time.
    LinearAlgebra::MiddleSizeVector getInitialSolution(
        const PointPhysicalT &point, const double &startTime,
        const std::size_t orderTimeDerivative) final {
        return getExactSolution(point, startTime, orderTimeDerivative);
    }

    /// \brief Create a mesh, solve the problem and return the total error.
    LinearAlgebra::MiddleSizeVector::type createAndSolve(
        const double T,       // final time
        const std::size_t nT  // number of time steps
    ) {
        this->readMesh(fileType);
        this->solve(0, T, (T / nT), 1, false);
        return this->computeTotalError(this->solutionVectorId_, T);
    }

   private:
    std::string fileType;

    /// Advective vector
    LinearAlgebra::SmallVector<DIM> a;
};

struct AdvectionParameters {
    std::size_t p;
    double endTime;
    std::size_t timeSteps;
};

template <std::size_t DIM>
void runAdvectionTestSet(ConvergenceTestSet &testSet,
                         AdvectionParameters &parameters, bool ignoreFailures) {
    runConvergenceTest(testSet, ignoreFailures,
                       [&parameters](std::string meshFile, std::size_t level) {
                           Advection<DIM> test(meshFile, parameters.p);
                           return std::real(test.createAndSolve(
                               parameters.endTime, parameters.timeSteps));
                       });
}

int main(int argc, char **argv) {
    using namespace std::string_literals;
    Base::parse_options(argc, argv);

    // Note: These results should be similar to those of
    // 080AdvectionLinear_SelfTest

    // Define clocks for measuring simulation time.
    std::chrono::time_point<std::chrono::system_clock> startClock, endClock;
    std::chrono::duration<double> elapsed_seconds;
    startClock = std::chrono::system_clock::now();

    // For recomputing the error tables
    bool ignoreFailures = false;

    AdvectionParameters dim1Params = {1, 0.1, 10};
    // Minimum of 4 elements
    ConvergenceTestSet dim1Meshes = {getUnitSegmentPeriodicMeshes(2),
                                     {
                                         6.61064338e-02,  //------
                                         1.79958457e-02,  //  3.67
                                         5.18359531e-03,  //  3.47
                                         1.51019916e-03,  //  3.43
                                     }};

    runAdvectionTestSet<1>(dim1Meshes, dim1Params, ignoreFailures);

    ConvergenceTestSet dim2P1Meshes = {getUnitSquarePeriodicTriangleMeshes(),
                                       {
                                           1.03487477e-01,  //------
                                           2.28412617e-02,  //  4.53
                                           6.75474161e-03,  //  3.38
                                           2.01183288e-03,  //  3.36
                                           5.54617877e-04,  //  3.63
                                       }};
    AdvectionParameters dim2P1Params = {1, 0.1, 10};
    runAdvectionTestSet<2>(dim2P1Meshes, dim2P1Params, ignoreFailures);

    ConvergenceTestSet dim3Meshes = {getUnitCubePeriodicCubeMeshes(0, 3),
                                     {
                                         1.03487477e-01,  //------
                                         2.28412617e-02,  //  4.53
                                         6.75474161e-03,  //  3.38
                                         2.01183288e-03,  //  3.36
                                         5.54617877e-04,  //  3.63
                                     }};
    AdvectionParameters dim3Params = {1, 0.02, 2};

    endClock = std::chrono::system_clock::now();
    elapsed_seconds = endClock - startClock;
    std::cout << "Elapsed time for solving the PDE: " << elapsed_seconds.count()
              << "s\n";

    return 0;
}

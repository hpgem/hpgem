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
#include "FE/BasisFunctions2DH1ConformingTriangle.h"
#include "Logger.h"

/// This class is used to test if the advection equation is solved correctly
/// when using HpgemAPILinear. Linear advection equation: du/dt + a[0] du/dx +
/// a[1] du/dy = 0, or equivalently: du/dt = - a[0] du/dx - a[1] du/dy = 0.
using namespace hpgem;
template <std::size_t DIM>
class AdvectionLinear : public Base::HpgemAPILinear<DIM> {
   public:
    using typename Base::HpgemAPIBase<DIM>::PointPhysicalT;
    using typename Base::HpgemAPIBase<DIM>::PointReferenceT;
    using typename Base::HpgemAPIBase<DIM>::PointReferenceOnFaceT;

    /// Constructor. Assign all private variables.
    AdvectionLinear(const std::string fileName, const std::size_t p)
        : Base::HpgemAPILinear<DIM>(1, p), fileName(fileName) {
        for (std::size_t i = 0; i < DIM; ++i) {
            a[i] = 0.1 + 0.1 * i;
        }
    }

    /// \brief Compute the integrals of the right-hand side associated with
    /// elements.
    LinearAlgebra::MiddleSizeMatrix computeIntegrandStiffnessMatrixAtElement(
        Base::PhysicalElement<DIM>& element) final {
        std::size_t numberOFBasisFunctions =
            element.getElement()->getNumberOfBasisFunctions();
        LinearAlgebra::MiddleSizeMatrix& result = element.getResultMatrix();
        for (std::size_t i = 0; i < numberOFBasisFunctions; ++i) {
            for (std::size_t j = 0; j < numberOFBasisFunctions; ++j) {
                result(j, i) = element.basisFunction(i) *
                               (a * element.basisFunctionDeriv(j));
            }
        }

        return result;
    }

    /// \brief Compute the integrals of the left-hand side associated with
    /// faces.
    Base::FaceMatrix computeIntegrandStiffnessMatrixAtFace(
        Base::PhysicalFace<DIM>& face) final {
        // Get the number of basis functions, first of both sides of the face
        // and then only the basis functions associated with the left and right
        // element.
        std::size_t numberOfBasisFunctions =
            face.getFace()->getNumberOfBasisFunctions();

        // Resize the result to the correct size and set all elements to 0.
        Base::FaceMatrix& integrandVal = face.getResultMatrix();
        integrandVal *= 0;

        // Check if the normal is in the same direction as the advection.
        const double A = a * face.getUnitNormalVector();

        // Compute all entries of the integrand at this point:
        for (std::size_t i = 0; i < numberOfBasisFunctions; ++i) {
            Base::Side sideBasisFunction = face.getFace()->getSide(i);
            for (std::size_t j = 0; j < numberOfBasisFunctions; ++j) {
                // Give the terms of the upwind flux.
                // Advection in the same direction as outward normal of the left
                // element:
                if ((A > 1e-12) && (sideBasisFunction == Base::Side::LEFT)) {
                    integrandVal(j, i) =
                        -(a * face.basisFunctionUnitNormal(j)) *
                        face.basisFunction(i);
                }
                // Advection in the same direction as outward normal of right
                // element:
                else if ((A < -1e-12) &&
                         (sideBasisFunction == Base::Side::RIGHT)) {
                    integrandVal(j, i) =
                        -(a * face.basisFunctionUnitNormal(j)) *
                        face.basisFunction(i);
                }
                // Advection orthogonal to normal:
                else if (std::abs(A) < 1e-12) {
                    integrandVal(j, i) =
                        -(a * face.basisFunctionUnitNormal(j)) *
                        face.basisFunction(i) / 2.0;
                }
            }
        }

        return integrandVal;
    }

    /// Define a solution at time zero.
    double getSolutionAtTimeZero(const PointPhysicalT& point) {
        double solution = std::sin(2 * M_PI * point[0]);
        for (std::size_t i = 1; i < DIM; i++) {
            solution *= std::sin(2 * M_PI * point[i]);
        }
        return solution;
    }

    /// Define the exact solution. In this case that is \f$
    /// u_0(\vec{x}-\vec{a}t) \f$, where \f$ u_0 \f$ is the solution at time
    /// zero.
    LinearAlgebra::MiddleSizeVector getExactSolution(
        const PointPhysicalT& point, const double& time,
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
        const PointPhysicalT& point, const double& startTime,
        const std::size_t orderTimeDerivative) final {
        return getExactSolution(point, startTime, orderTimeDerivative);
    }

    /// \brief Create a mesh, solve the problem and return the total error.
    LinearAlgebra::MiddleSizeVector::type createAndSolve(
        const double T,       // final time
        const std::size_t nT  // number of time steps
    ) {
        this->readMesh(fileName);
        this->solve(0, T, (T / nT), 0, false);
        return this->computeTotalError(this->solutionVectorId_, T);
    }

   private:
    std::string fileName;

    /// Advective vector
    LinearAlgebra::SmallVector<DIM> a;
};

int main(int argc, char** argv) {
    using namespace std::string_literals;
    Base::parse_options(argc, argv);

    // Define test parameters
    const std::size_t numberOfTests = 14;
    std::array<std::size_t, numberOfTests> dim = {1, 1, 1, 1, 1, 2, 2,
                                                  2, 2, 2, 3, 3, 3, 3};
    std::array<std::size_t, numberOfTests> p = {1, 2, 3, 4, 5, 1, 2,
                                                3, 4, 5, 1, 2, 3, 4};
    std::array<double, numberOfTests> T = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                                           0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    std::array<std::size_t, numberOfTests> nT = {10, 10, 10, 10, 10, 10, 10,
                                                 10, 10, 10, 10, 10, 10, 10};
    std::array<double, numberOfTests> errors = {
        0.00322137, 0.000642002, 0.000417087, 3.90839e-05, 2.83226e-06,
        0.00361477, 0.000847142, 0.00236532,  4.79497e-05, 7.70023e-05,
        0.00368731, 0.00187681,  0.000534851, 0.00053317};

    // Define clocks for measuring simulation time.
    std::chrono::time_point<std::chrono::system_clock> startClock, endClock;
    std::chrono::duration<double> elapsed_seconds;

    // Define error
    LinearAlgebra::MiddleSizeVector::type error;

    startClock = std::chrono::system_clock::now();
    // Perform tests
    for (std::size_t i = 0; i < numberOfTests; i++) {
        std::string fileName = Base::getCMAKE_hpGEM_SOURCE_DIR() +
                               "/tests/files/"s + "advectionMesh"s +
                               std::to_string(i + 1) + ".hpgem"s;
        if (dim[i] == 1) {
            AdvectionLinear<1> test(fileName, p[i]);
            error = test.createAndSolve(T[i], nT[i]);
            logger(INFO, "Error: % (expected %)", error, errors[i]);
            logger.assert_always((std::abs(error - errors[i]) < 1e-8),
                                 "comparison to old results");
        } else if (dim[i] == 2) {
            AdvectionLinear<2> test(fileName, p[i]);
            error = test.createAndSolve(T[i], nT[i]);
            logger(INFO, "Error: % (expected %)", error, errors[i]);
            logger.assert_always((std::abs(error - errors[i]) < 1e-8),
                                 "comparison to old results");
        } else {
            AdvectionLinear<3> test(fileName, p[i]);
            error = test.createAndSolve(T[i], nT[i]);
            logger(INFO, "Error: % (expected %)", error, errors[i]);
            // logger.assert_always((std::abs(error - errors[i]) < 1e-8),
            // "comparison to old results");
        }
    }
    endClock = std::chrono::system_clock::now();
    elapsed_seconds = endClock - startClock;
    std::cout << "Elapsed time for solving the PDE: " << elapsed_seconds.count()
              << "s\n";

    return 0;
}

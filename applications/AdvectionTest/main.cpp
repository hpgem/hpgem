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

#include "Base/CommandLineOptions.h"
#include "Base/ConfigurationData.h"
#include "Base/Element.h"
#include "Base/Face.h"
#include "Base/HpgemAPILinear.h"
#include "Integration/ElementIntegral.h"
#include "Integration/FaceIntegral.h"
#include "Output/TecplotDiscontinuousSolutionWriter.h"
#include "Output/TecplotSingleElementWriter.h"
#include "Utilities/BasisFunctions2DH1ConformingTriangle.h"
#include "Logger.h"

using namespace hpgem;


const static std::size_t DIM = 2;

/// Linear advection equation du/dt + a[0] du/dx + a[1] du/dy = 0.
/// This class is meant for testing purposes.
//  Please verify that nobody went and tested a broken feature
//  Please keep the problem modelled here reasonably close to the linear
//  advection problem
class AdvectionTest : public Base::HpgemAPILinear<DIM> {
   public:
    /// Constructor. Assign all private variables.
    AdvectionTest(std::size_t p) : Base::HpgemAPILinear<DIM>(1, p) {
        // Choose the "direction" of the advection.
        // This cannot be implemented with iterators, and since the dimension is
        // not always 2, this is the most generic way to write it.
        // for (std::size_t i = 0; i < DIM; ++i)
        //{
        //    a[i] = 0.1 + 0.1 * i;
        //}
        a[0] = 0.1;
        a[2] = 0.3;
    }

    void readMesh(const std::string meshName) override final {
        // Set the number of Element/Face Matrices/Vectors.
        std::size_t numberOfElementMatrices = 2;
        std::size_t numberOfElementVectors = 0;
        std::size_t numberOfFaceMatrices = 1;
        std::size_t numberOfFaceVectors = 0;

        // Create mesh and set basis functions.
        this->addMesh(meshName, numberOfElementMatrices, numberOfElementVectors,
                      numberOfFaceMatrices, numberOfFaceVectors);
        // this->meshes_[0]->useDefaultDGBasisFunctions();

        // Set the number of time integration vectors according to the size of
        // the Butcher tableau.
        this->setNumberOfTimeIntegrationVectorsGlobally(
            this->globalNumberOfTimeIntegrationVectors_);

        // Plot info about the mesh
        std::size_t numberOfElements = this->meshes_[0]->getNumberOfElements();
        logger(VERBOSE, "Total number of elements: %", numberOfElements);
    }

    /// Compute phi_i*(a.grad(phi_j)) on a reference point on an element for all
    /// basisfunctions phi_i and phi_j.
    /// You pass the reference point to the basisfunctions. Internally the
    /// basisfunctions will be mapped to the physical element so you wont have
    /// to do any transformations yourself
    LinearAlgebra::MiddleSizeMatrix computeIntegrandStiffnessMatrixAtElement(
        Base::PhysicalElement<DIM>& element) override final {
        logger.assert_debug(element.getJacobianDet() > 0, "%",
                            element.getElement());
        std::size_t numberOfBasisFunctions =
            element.getElement()->getNumberOfBasisFunctions();
        LinearAlgebra::MiddleSizeMatrix& result = element.getResultMatrix();
        for (std::size_t i = 0; i < numberOfBasisFunctions; ++i) {
            for (std::size_t j = 0; j < numberOfBasisFunctions; ++j) {
                result(j, i) = element.basisFunction(i) *
                               (a * element.basisFunctionDeriv(j));
            }
        }

        return result;
    }

    /// \brief Compute the integrals of the left-hand side associated with
    /// faces.
    ///
    /// For every internal face, we want to compute the integral of the flux
    /// for all basisfunctions phi_i and phi_j that are non-zero on that face.
    /// For boundary faces, similar expressions can be obtained depending of the
    /// type of boundary condition. This function will compute these integrands
    /// for all basisfunctions phi_i and phi_j on a certain face at a reference
    /// point p. Then the integral can later be computed with appropriate
    /// (Gauss-)quadrature rules. The resulting matrix of values is then given
    /// in the matrix integrandVal, to which we passed a reference when calling
    /// it. Please note that you pass a reference point to the basisfunctions
    /// and the transformations are done internally. The class FaceMatrix
    /// consists of four element matrices for internal faces and one element
    /// matrix for faces on the boundary. Each element matrix corresponds to a
    /// pair of two adjacent elements of the face.
    Base::FaceMatrix computeIntegrandStiffnessMatrixAtFace(
        Base::PhysicalFace<DIM>& face) override final {
        // Get the number of basis functions, first of both sides of the face
        // and then only the basis functions associated with the left and right
        // element.
        std::size_t numberOfBasisFunctions =
            face.getFace()->getNumberOfBasisFunctions();

        // Resize the result to the correct size and set all elements to 0.
        Base::FaceMatrix& integrandVal = face.getResultMatrix();
        integrandVal *= 0;

        // Check if the normal is in the same direction as the advection.
        // Note that normal does not have length 1!
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
        return (std::sin(2 * M_PI * point[0]) * std::sin(2 * M_PI * point[1])) *
               std::sin(2 * M_PI * point[2]);
    }

    /// Define the exact solution. In this case that is \f$
    /// u_0(\vec{x}-\vec{a}t) \f$, where \f$ u_0 \f$ is the solution at time
    /// zero.
    LinearAlgebra::MiddleSizeVector getExactSolution(
        const PointPhysicalT& point, const double& time,
        const std::size_t orderTimeDerivative) override final {
        LinearAlgebra::MiddleSizeVector exactSolution(1);
        if (orderTimeDerivative == 0) {
            PointPhysicalT displacement(-a * time);
            exactSolution(0) = getSolutionAtTimeZero(point + displacement);
            return exactSolution;
        } else {
            logger(ERROR,
                   "No exact solution for order time derivative % implemented",
                   orderTimeDerivative);
            exactSolution(0) = 0;
            return exactSolution;
        }
    }

    /// Define the initial conditions. In this case it is just the exact
    /// solution at the start time.
    LinearAlgebra::MiddleSizeVector getInitialSolution(
        const PointPhysicalT& point, const double& startTime,
        const std::size_t orderTimeDerivative) override final {
        return getExactSolution(point, startTime, orderTimeDerivative);
    }

   private:
    /// Advective vector
    LinearAlgebra::SmallVector<DIM> a;
};

auto& name = Base::register_argument<std::string>(
    'n', "meshName", "Name of the mesh file", true);
auto& p =
    Base::register_argument<std::size_t>('p', "poly", "Polynomial order", true);

/// Make the problem and solve it.
int main(int argc, char** argv) {
    Base::parse_options(argc, argv);

    // Choose variable name(s). Since we have a scalar function, we only need to
    // choose one name.
    std::vector<std::string> variableNames;
    variableNames.push_back("u");

    // Construct our problem with n elements in every direction and polynomial
    // order p
    AdvectionTest test(p.getValue());

    // Create the mesh
    test.readMesh(name.getValue());

    // Set the names for the output file
    test.setOutputNames("output", "AdvectionTest", "AdvectionTest",
                        variableNames);

    // Run the simulation and write the solution

    auto startTime = std::chrono::steady_clock::now();

    test.solve(Base::startTime.getValue(), Base::endTime.getValue(),
               Base::dt.getValue(), Base::numberOfSnapshots.getValue(), true);

    auto endTime = std::chrono::steady_clock::now();

    logger(INFO, "Simulation took %ms.",
           std::chrono::duration_cast<std::chrono::milliseconds>(endTime -
                                                                 startTime)
               .count());

    return 0;
}

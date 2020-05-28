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

// this tutorial solves a 2 dimensional problem
const std::size_t DIM = 2;

/// Linear advection equation du/dt + a[0] du/dx + a[1] du/dy = 0.
/// The first self-contained (no PETSc) program to make it into the SVN
class TutorialAdvection : public Base::HpgemAPILinear<DIM> {
   public:
    /// Constructor. Assign all private variables.
    TutorialAdvection(std::size_t p) : Base::HpgemAPILinear<DIM>(1, p) {
        // Choose the "direction" of the advection.
        // This cannot be implemented with iterators, and since the dimension is
        // not always 2, this is the most generic way to write it.
        for (std::size_t i = 0; i < DIM; ++i) {
            a[i] = 0.1 + 0.1 * i;
        }
    }

    /// Compute phi_i*(a.grad(phi_j)) on an element for all
    /// basisfunctions phi_i and phi_j.
    /// hpGEM pretends the computations are done on a physical element (as
    /// opposed to a reference element), because this generally allows for
    /// easier expressions
    LinearAlgebra::MiddleSizeMatrix computeIntegrandStiffnessMatrixAtElement(
        Base::PhysicalElement<DIM>& element) final {
        // we access the actual element to find the number of basis functions
        // that are non-zero on this element
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
    /// for all basisfunctions phi_i and phi_j Then the integral can later be
    /// computed with appropriate (Gauss-)quadrature rules. The resulting matrix
    /// of values is then given in the matrix integrandVal, to which we passed a
    /// reference when calling it. hpGEM pretends the computations are done on a
    /// physical face (as opposed to a reference face), because this generally
    /// allows for easier expressions The class FaceMatrix consists of four
    /// element matrices for internal faces and one element matrix for faces on
    /// the boundary. Each element matrix corresponds to a pair of two adjacent
    /// elements of the face.
    Base::FaceMatrix computeIntegrandStiffnessMatrixAtFace(
        Base::PhysicalFace<DIM>& face) final {
        // Get the total number of basis functions of both sides of the face.
        std::size_t numberOfBasisFunctions =
            face.getFace()->getNumberOfBasisFunctions();

        // get the result with the correct size and all elements set to 0.
        Base::FaceMatrix& integrandVal = face.getResultMatrix();

        // Check if the normal is in the same direction as the advection.
        const double A = (a * face.getUnitNormalVector());

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
        return (std::sin(2 * M_PI * point[0]) * std::sin(2 * M_PI * point[1]));
    }

    /// Define the exact solution. In this case that is \f$
    /// u_0(\vec{x}-\vec{a}t) \f$, where \f$ u_0 \f$ is the solution at time
    /// zero. Note that this function can be used to compute time derivatives of
    /// the initial conditions, but for the advection equations this does not
    /// make sense
    LinearAlgebra::MiddleSizeVector getExactSolution(
        const PointPhysicalT& point, const double& time,
        const std::size_t orderTimeDerivative) final {
        logger.assert_debug(
            orderTimeDerivative == 0,
            "No exact solution for order time derivative % implemented");
        LinearAlgebra::MiddleSizeVector result(1);
        result[0] = getSolutionAtTimeZero(point - a * time);
        return result;
    }

    /// Define the initial conditions. In this case it is just the exact
    /// solution at the start time.
    LinearAlgebra::MiddleSizeVector getInitialSolution(
        const PointPhysicalT& point, const double& startTime,
        const std::size_t orderTimeDerivative) final {
        return getExactSolution(point, startTime, orderTimeDerivative);
    }

   private:
    /// Advective vector
    LinearAlgebra::SmallVector<DIM> a;
};

auto& name = Base::register_argument<std::string>(
    'n', "meshName", "name of the mesh file", true);
auto& p =
    Base::register_argument<std::size_t>('p', "poly", "Polynomial order", true);

/// Make the problem and solve it.
int main(int argc, char** argv) {
    Base::parse_options(argc, argv);

    // Choose variable name(s). Since we have a scalar function, we only need to
    // chooes one name.
    std::vector<std::string> variableNames;
    variableNames.push_back("u");

    // Construct our problem with n elements in every direction and polynomial
    // order p
    TutorialAdvection test(p.getValue());

    // read the mesh
    test.readMesh(name.getValue());

    // Set the names for the output file
    test.setOutputNames("output", "TutorialAdvection", "TutorialAdvection",
                        variableNames);

    // Run the simulation and write the solution
    test.solve(Base::startTime.getValue(), Base::endTime.getValue(),
               Base::dt.getValue(), Base::numberOfSnapshots.getValue(), true);

    return 0;
}

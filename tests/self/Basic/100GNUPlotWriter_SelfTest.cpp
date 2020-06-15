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

#include "Base/ConfigurationData.h"
#include "Base/Element.h"
#include "Base/Face.h"
#include "Base/HpgemAPILinear.h"
#include "Integration/ElementIntegral.h"
#include "Integration/FaceIntegral.h"
#include "Output/GNUPlotDiscontinuousSolutionWriter.h"
#include "Utilities/BasisFunctions2DH1ConformingTriangle.h"
#include "Logger.h"

/// This class is used to test if the advection equation is written to a GNU
/// plot file correctly when using HpgemAPISimplified. Linear advection
/// equation: du/dt + a[0] du/dx + a[1] du/dy = 0, or equivalently: du/dt = -
/// a[0] du/dx - a[1] du/dy = 0. Note that this is a copy of the
/// AdvectionSelfTest, without the error checking and with output file writing.
template <std::size_t DIM>
class Advection : public Base::HpgemAPISimplified<DIM>,
                  public Output::SingleElementWriter<DIM> {
   public:
    using typename Base::HpgemAPIBase<DIM>::PointPhysicalT;
    using typename Base::HpgemAPIBase<DIM>::PointReferenceT;
    using typename Base::HpgemAPIBase<DIM>::PointReferenceOnFaceT;

    /// Constructor. Assign all private variables.
    Advection(const std::string fileName, const std::size_t p)
        : Base::HpgemAPISimplified<DIM>(1, p), fileName(fileName) {
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
        const double time) override final {
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
        const double time) override final {
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
        const std::size_t orderTimeDerivative) override final {
        LinearAlgebra::MiddleSizeVector exactSolution(1);
        if (orderTimeDerivative == 0) {
            PointPhysicalT displacement{a * time};
            exactSolution(0) = getSolutionAtTimeZero(point - displacement);
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
        const PointPhysicalT &point, const double &startTime,
        const std::size_t orderTimeDerivative) override final {
        return getExactSolution(point, startTime, orderTimeDerivative);
    }

    /// \brief Create a mesh, solve the problem and return the total error.
    void createSolveAndWrite(const double finalTime,
                             const std::string &outputFileName) {
        using namespace std::string_literals;
        this->readMesh(Base::getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/"s +
                       fileName);
        this->solve(0, finalTime, 1e-3, 0, false);
        std::ofstream outputFile(outputFileName);
        std::string xNames;
        if (DIM == 2) {
            xNames = "01";
        } else {
            xNames = "0";
        }
        Output::GNUPlotDiscontinuousSolutionWriter<DIM>
            gnuPlotDiscontinuousSolutionWriter(
                outputFile, "Advection for testing GNUPlotWriter", xNames, "u");
        gnuPlotDiscontinuousSolutionWriter.write(this->meshes_[0], this);
    }

    void writeOutput(const Base::Element *ptrElement,
                     const PointReferenceT &pRef, std::ostream &out) override {
        const LinearAlgebra::MiddleSizeVector solution =
            ptrElement->getSolution(0, pRef);
        out << std::real(solution(0));
    }

   private:
    /// Number of elements per direction
    std::string fileName;

    /// Advective vector
    LinearAlgebra::SmallVector<DIM> a;
};

int main(int argc, char **argv) {
    Base::parse_options(argc, argv);

    Advection<1> test1R("plottingMesh1.hpgem", 1);
    test1R.createSolveAndWrite(.01, "GNUPlotWriterTestFile1DRectangular");

    Advection<2> test2T("plottingMesh2.hpgem", 1);
    test2T.createSolveAndWrite(.01, "GNUPlotWriterTestFile2DTriangular");

    Advection<2> test2R("plottingMesh3.hpgem", 1);
    test2R.createSolveAndWrite(.01, "GNUPlotWriterTestFile2DRectangular");

    return 0;
}

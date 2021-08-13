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

#include <Integration/QuadratureRules/AllGaussQuadratureRules.h>
#include "HpgemAPISimplified.h"

#include "Base/CommandLineOptions.h"
#include "Base/ConfigurationData.h"
#include "Base/Element.h"
#include "Base/Face.h"
#include "Base/MpiContainer.h"
#include "Base/TimeIntegration/AllTimeIntegrators.h"
#include "Geometry/PointReference.h"
#include "Integration/ElementIntegral.h"
#include "Integration/FaceIntegral.h"
#include "Output/TecplotDiscontinuousSolutionWriter.h"
#include "Output/TecplotSingleElementWriter.h"
#include "FE/BasisFunctions1DH1ConformingLine.h"
#include "FE/BasisFunctions2DH1ConformingSquare.h"
#include "FE/BasisFunctions2DH1ConformingTriangle.h"
#include "FE/BasisFunctions3DH1ConformingCube.h"
#include "FE/BasisFunctions3DH1ConformingTetrahedron.h"
#include "LinearAlgebra/Axpy.h"

#include "Logger.h"
namespace hpgem {
namespace Base {

///\bug Workaround for Bug 60352 in (at least) gcc 4.8.2 (should read auto&
/// numberOfSnapshots = ...)
extern CommandLineOption<std::size_t> &numberOfSnapshots;
extern CommandLineOption<double> &endTime;
extern CommandLineOption<double> &startTime;
extern CommandLineOption<double> &dt;
extern CommandLineOption<std::string> &outputName;

/// \param[in] numberOfVariables Number of variables in the PDE
/// \param[in] polynomialOrder Polynomial order of the basis functions
/// \param[in] ptrButcherTableau A butcherTableau used to solve the PDE with a
/// Runge-Kutta method. \param[in] numberOfTimeLevels Number of time levels.
/// \param[in] computeBothFaces Compute integrands for the test functions on
/// each sides of the face simultaneously (true) or seperately (false). Note:
/// requires different face integrand function
template <std::size_t DIM>
HpgemAPISimplified<DIM>::HpgemAPISimplified(
    const std::size_t numberOfVariables, const std::size_t polynomialOrder,
    const TimeIntegration::ButcherTableau *const ptrButcherTableau,
    const std::size_t numberOfTimeLevels, const bool computeBothFaces)
    : HpgemAPIBase<DIM>(
          new Base::GlobalData,
          new Base::ConfigurationData(numberOfVariables, numberOfTimeLevels)),
      ptrButcherTableau_(ptrButcherTableau),
      outputFileName_("output"),
      internalFileTitle_("output"),
      solutionTitle_("solution"),
      computeBothFaces_(computeBothFaces),
      polynomialOrder_(polynomialOrder) {
    this->globalNumberOfTimeIntegrationVectors_ =
        ptrButcherTableau->getNumberOfStages() + 1;
    solutionVectorId_ = 0;
    for (std::size_t i = 1; i < this->globalNumberOfTimeIntegrationVectors_;
         i++) {
        auxiliaryVectorIds_.push_back(i);
    }
    for (std::size_t iV = 0; iV < this->configData_->numberOfUnknowns_; iV++) {
        std::string variableName = "variable" + std::to_string(iV);
        variableNames_.push_back(variableName);
    }
}

/// \param[in] numberOfVariables Number of variables in the PDE
/// \param[in] polynomialOrder Polynomial order of the basis functions
/// \param[in] globalNumberOfTimeIntegrationVectors number of time integration
/// vectors for every element. \param[in] numberOfTimeLevels Number of time
/// levels. \param[in] computeBothFaces Compute integrands for test functions on
/// both sides of the faces simultaneously (true) or seperately (false).
template <std::size_t DIM>
HpgemAPISimplified<DIM>::HpgemAPISimplified(
    const std::size_t numberOfVariables, const std::size_t polynomialOrder,
    const std::size_t globalNumberOfTimeIntegrationVectors,
    const std::size_t numberOfTimeLevels, const bool computeBothFaces)
    : HpgemAPIBase<DIM>(
          new Base::GlobalData,
          new Base::ConfigurationData(numberOfVariables, numberOfTimeLevels)),
      ptrButcherTableau_(nullptr),
      outputFileName_("output"),
      internalFileTitle_("output"),
      solutionTitle_("solution"),
      computeBothFaces_(computeBothFaces),
      polynomialOrder_(polynomialOrder) {
    this->globalNumberOfTimeIntegrationVectors_ =
        globalNumberOfTimeIntegrationVectors;
    solutionVectorId_ = 0;
    for (std::size_t i = 1; i < this->globalNumberOfTimeIntegrationVectors_;
         i++) {
        auxiliaryVectorIds_.push_back(i);
    }
    for (std::size_t iV = 0; iV < this->configData_->numberOfUnknowns_; iV++) {
        std::string variableName = "variable" + std::to_string(iV);
        variableNames_.push_back(variableName);
    }
}

template <std::size_t DIM>
void HpgemAPISimplified<DIM>::readMesh(std::string fileName) {
    // Set the number of Element/Face Matrices/Vectors.
    std::size_t numberOfElementMatrices = 0;
    std::size_t numberOfElementVectors = 0;
    std::size_t numberOfFaceMatrices = 0;
    std::size_t numberOfFaceVectors = 0;

    // Create mesh and set basis functions.
    this->addMesh(fileName, numberOfElementMatrices, numberOfElementVectors,
                  numberOfFaceMatrices, numberOfFaceVectors);
    this->meshes_[0]->useDefaultDGBasisFunctions(this->polynomialOrder_);

    // Set the number of time integration vectors according to the size of the
    // Butcher tableau.
    this->setNumberOfTimeIntegrationVectorsGlobally(
        this->globalNumberOfTimeIntegrationVectors_);

    // Plot info about the mesh
    std::size_t numberOfElements = this->meshes_[0]->getNumberOfElements();
    logger(VERBOSE, "Total number of elements: %", numberOfElements);
}

/// \details By default this function computes the matrix of the products of all
/// basis functions corresponding to the given element.
template <std::size_t DIM>
LinearAlgebra::MiddleSizeMatrix
    HpgemAPISimplified<DIM>::computeIntegrandMassMatrix(
        Base::PhysicalElement<DIM> &element) {
    // Get a reference to the result matrix.
    LinearAlgebra::MiddleSizeMatrix &integrand = element.getResultMatrix();

    // Get the number of basis functions.
    const std::size_t numberOfBasisFunctions =
        element.getNumberOfBasisFunctions();

    // Compute the matrix of products for all basis functions.
    std::size_t iVB, jVB;  // indices for both variable and basis function.
    for (std::size_t iV = 0; iV < this->configData_->numberOfUnknowns_; ++iV) {
        for (std::size_t iB = 0; iB < numberOfBasisFunctions; ++iB) {
            for (std::size_t jB = 0; jB < numberOfBasisFunctions; ++jB) {
                iVB = element.convertToSingleIndex(iB, iV);
                jVB = element.convertToSingleIndex(jB, iV);
                integrand(iVB, jVB) =
                    element.basisFunction(iB) * element.basisFunction(jB);
            }
        }
    }
    return integrand;
}

/// \details By default this function computes the mass matrix that corresponds
/// to the integral of the inner product of the test functions on the element.
template <std::size_t DIM>
LinearAlgebra::MiddleSizeMatrix
    HpgemAPISimplified<DIM>::computeMassMatrixAtElement(
        Base::Element *ptrElement) {
    std::function<LinearAlgebra::MiddleSizeMatrix(Base::PhysicalElement<DIM> &)>
        integrandFunction = [=](Base::PhysicalElement<DIM> &element)
        -> LinearAlgebra::MiddleSizeMatrix {
        return this->computeIntegrandMassMatrix(element);
    };

    return this->elementIntegrator_.integrate(ptrElement, integrandFunction);
}

template <std::size_t DIM>
void HpgemAPISimplified<DIM>::solveMassMatrixEquationsAtElement(
    Base::Element *ptrElement,
    LinearAlgebra::MiddleSizeVector &functionCoefficients) {
    computeMassMatrixAtElement(ptrElement).solve(functionCoefficients);
}

/// \details Solve the equation \f$ Mu = r \f$ for \f$ u \f$, where \f$ r \f$ is
/// the right-hand sid and \f$ M \f$ is the mass matrix.
template <std::size_t DIM>
void HpgemAPISimplified<DIM>::solveMassMatrixEquations(
    const std::size_t timeIntegrationVectorId) {
    for (Base::Element *ptrElement : this->meshes_[0]->getElementsList()) {
        LinearAlgebra::MiddleSizeVector &functionCoefficients(
            ptrElement->getTimeIntegrationVector(timeIntegrationVectorId));

        solveMassMatrixEquationsAtElement(ptrElement, functionCoefficients);
    }

    this->synchronize(timeIntegrationVectorId);
}

template <std::size_t DIM>
LinearAlgebra::MiddleSizeVector
    HpgemAPISimplified<DIM>::computeIntegrandInitialSolution(
        Base::PhysicalElement<DIM> &element, const double startTime,
        const std::size_t orderTimeDerivative) {
    // Get a reference to the result vector.
    LinearAlgebra::MiddleSizeVector &integrand = element.getResultVector();

    // Get the physical point.
    const PointPhysicalT &pPhys = element.getPointPhysical();

    // Compute the initial solution.
    const LinearAlgebra::MiddleSizeVector initialSolution =
        getInitialSolution(pPhys, startTime, orderTimeDerivative);

    // Get the number of basis functions.
    const std::size_t numberOfBasisFunctions =
        element.getNumberOfBasisFunctions();

    // Compute the product of the initial solution and all test functions.
    for (std::size_t iV = 0; iV < this->configData_->numberOfUnknowns_; ++iV) {
        for (std::size_t iB = 0; iB < numberOfBasisFunctions; ++iB) {
            // Determine index for this combination of variable and basis
            // function.
            const std::size_t iVB = element.convertToSingleIndex(iB, iV);
            integrand(iVB) = element.basisFunction(iB) * initialSolution(iV);
        }
    }

    return integrand;
}

/// \brief By default this function computes the integral of the inner product
/// of the initial solution (for given order time derivative) and the test
/// function on the element.
template <std::size_t DIM>
LinearAlgebra::MiddleSizeVector
    HpgemAPISimplified<DIM>::integrateInitialSolutionAtElement(
        Base::Element *ptrElement, const double startTime,
        const std::size_t orderTimeDerivative) {
    // Define the integrand function for the the initial solution integral.
    const std::function<LinearAlgebra::MiddleSizeVector(
        Base::PhysicalElement<DIM> &)>
        integrandFunction = [=](Base::PhysicalElement<DIM> &element)
        -> LinearAlgebra::MiddleSizeVector {
        return this->computeIntegrandInitialSolution(element, startTime,
                                                     orderTimeDerivative);
    };

    return this->elementIntegrator_.integrate(ptrElement, integrandFunction);
}

template <std::size_t DIM>
void HpgemAPISimplified<DIM>::integrateInitialSolution(
    const std::size_t resultVectorId, const double initialTime,
    const std::size_t orderTimeDerivative) {
    for (Base::Element *ptrElement : this->meshes_[0]->getElementsList()) {
        LinearAlgebra::MiddleSizeVector &resultFunctionCoefficients =
            ptrElement->getTimeIntegrationVector(resultVectorId);
        resultFunctionCoefficients = integrateInitialSolutionAtElement(
            ptrElement, initialTime, orderTimeDerivative);
    }

    this->synchronize(resultVectorId);
}

template <std::size_t DIM>
LinearAlgebra::MiddleSizeVector::type
    HpgemAPISimplified<DIM>::computeIntegrandTotalError(
        Base::PhysicalElement<DIM> &element,
        const LinearAlgebra::MiddleSizeVector &solutionCoefficients,
        const double time) {
    // Define the integrand.
    LinearAlgebra::MiddleSizeVector::type integrand = 0;

    // Get the physical point.
    const PointPhysicalT &pPhys = element.getPointPhysical();

    // Compute the real solution.
    const LinearAlgebra::MiddleSizeVector exactSolution =
        getExactSolution(pPhys, time, 0);

    // Get the number of basis functions.
    const std::size_t numberOfBasisFunctions =
        element.getNumberOfBasisFunctions();

    // Compute the numerical solution.
    LinearAlgebra::MiddleSizeVector numericalSolution(
        this->configData_->numberOfUnknowns_);
    numericalSolution *= 0;
    for (std::size_t jB = 0; jB < numberOfBasisFunctions; ++jB) {
        for (std::size_t jV = 0; jV < this->configData_->numberOfUnknowns_;
             ++jV) {
            const std::size_t jVB = element.convertToSingleIndex(jB, jV);
            numericalSolution(jV) +=
                element.basisFunction(jB) * solutionCoefficients(jVB);
        }
    }

    // Compute the error
    const LinearAlgebra::MiddleSizeVector error =
        exactSolution - numericalSolution;

    // Compute the square of the l2 norm.
    for (std::size_t iV = 0; iV < this->configData_->numberOfUnknowns_; iV++) {
        integrand += error(iV) * error(iV);
    }

    return integrand;
}

/// \details By default the square of the standard L2 norm is integrated.
template <std::size_t DIM>
LinearAlgebra::MiddleSizeVector::type
    HpgemAPISimplified<DIM>::integrateErrorAtElement(
        Base::Element *ptrElement,
        LinearAlgebra::MiddleSizeVector &solutionCoefficients,
        const double time) {
    // Define the integrand function for the error energy.
    const std::function<LinearAlgebra::MiddleSizeVector::type(
        Base::PhysicalElement<DIM> &)>
        integrandFunction = [=](Base::PhysicalElement<DIM> &element)
        -> LinearAlgebra::MiddleSizeVector::type {
        return this->computeIntegrandTotalError(element, solutionCoefficients,
                                                time);
    };

    // Technically 2 * polynomialOrder_ would be sufficient. However, it is not
    // uncommon that this underestimates the error. Therefore we use a rule of
    // slightly higher order to ensure that we accurately compute the error.
    auto *rule = QuadratureRules::AllGaussQuadratureRules::instance().getRule(
        ptrElement->getReferenceGeometry(), 2 * this->polynomialOrder_ + 2);

    return this->elementIntegrator_.integrate(ptrElement, integrandFunction,
                                              rule);
}

/// \param[in] solutionVectorId index of the time integration vector where the
/// solution is stored. \param[in] time Time corresponding to the current
/// solution. \details The square of the total error is defined as \f[
/// e_{total}^2 := \int \|e\|^2 \,dV \f], where \f$\|e\|\f$ is some user-defined
/// norm (based on the (weighted) inner product) of the error. By default this
/// is the standard L2 norm.
template <std::size_t DIM>
LinearAlgebra::MiddleSizeVector::type
    HpgemAPISimplified<DIM>::computeTotalError(
        const std::size_t solutionVectorId, const double time) {
    LinearAlgebra::MiddleSizeVector::type localError = 0;

    for (Base::Element *ptrElement : this->meshes_[0]->getElementsList()) {
        LinearAlgebra::MiddleSizeVector &solutionCoefficients =
            ptrElement->getTimeIntegrationVector(solutionVectorId);
        localError +=
            integrateErrorAtElement(ptrElement, solutionCoefficients, time);
    }

    LinearAlgebra::MiddleSizeVector::type error;

    MPIContainer &communicator = MPIContainer::Instance();
#ifdef HPGEM_USE_MPI
    communicator.reduce(localError, MPI_SUM);
#endif
    communicator.onlyOnOneProcessor({[&]() {
        logger.assert_debug(std::abs(std::imag(localError)) < 1e-12,
                            "The computed error has an imaginary component");
        if (std::real(localError) >= 0) {
            error = std::sqrt(localError);
        } else {
            logger(WARN, "Warning: the computed total error is negative.");
            error = std::sqrt(-localError);
        }
    }});
    MPIContainer::Instance().broadcast(error);
    return error;
}

/// \details This function returns a vector of the suprema of the error of every
/// variable.
template <std::size_t DIM>
LinearAlgebra::MiddleSizeVector
    HpgemAPISimplified<DIM>::computeMaxErrorAtElement(
        Base::Element *ptrElement,
        const LinearAlgebra::MiddleSizeVector &solutionCoefficients,
        const double time) {
    // Get number of basis functions
    const std::size_t numberOfBasisFunctions =
        ptrElement->getNumberOfBasisFunctions();

    // Declare vector of maxima of the error.
    LinearAlgebra::MiddleSizeVector maxError(
        this->configData_->numberOfUnknowns_);
    maxError *= 0;

    // Get quadrature rule and number of points.
    const QuadratureRules::GaussQuadratureRule *quadratureRule =
        ptrElement->getGaussQuadratureRule();
    const std::size_t numberOfQuadPoints = quadratureRule->getNumberOfPoints();

    // For each quadrature point update the maxima of the error.
    for (std::size_t pQuad = 0; pQuad < numberOfQuadPoints; ++pQuad) {
        const Geometry::PointReference<DIM> &pRef =
            quadratureRule->getPoint(pQuad);
        const Geometry::PointPhysical<DIM> &pPhys =
            ptrElement->referenceToPhysical(pRef);

        const LinearAlgebra::MiddleSizeVector exactSolution =
            getExactSolution(pPhys, time, 0);

        LinearAlgebra::MiddleSizeVector numericalSolution(
            this->configData_->numberOfUnknowns_);
        numericalSolution *= 0;

        for (std::size_t iB = 0; iB < numberOfBasisFunctions; ++iB) {
            const double valueBasisFunction =
                ptrElement->basisFunction(iB, pRef);

            for (std::size_t iV = 0; iV < this->configData_->numberOfUnknowns_;
                 ++iV) {
                const std::size_t iVB =
                    ptrElement->convertToSingleIndex(iB, iV);

                numericalSolution(iV) +=
                    solutionCoefficients(iVB) * valueBasisFunction;
            }
        }

        for (std::size_t iV = 0; iV < this->configData_->numberOfUnknowns_;
             ++iV) {
            if (std::abs(numericalSolution(iV) - exactSolution(iV)) >
                std::real(maxError(iV))) {
                maxError(iV) =
                    std::abs(numericalSolution(iV) - exactSolution(iV));
            }
        }
    }

    return maxError;
}

#ifdef HPGEM_USE_MPI
// Note typesignature is rather liberal because it is passed to MPI_Op_create
inline void computeMPIMaximum(void *in, void *inout, int *len,
                              MPI_Datatype *type) {
    logger.assert_debug(
        *len >= 0,
        "MPI wants to create the piecewise maximum of two vectors of size %",
        len);
    for (std::size_t i = 0; i < static_cast<std::size_t>(*len); i++) {
        logger.assert_debug(
            std::abs(std::imag(
                reinterpret_cast<LinearAlgebra::MiddleSizeVector::type *>(
                    inout)[i])) < 1e-12,
            "can only do this for complex numbers");
        logger.assert_debug(
            std::abs(std::imag(
                reinterpret_cast<const LinearAlgebra::MiddleSizeVector::type *>(
                    in)[i])) < 1e-12,
            "can only do this for complex numbers");
        reinterpret_cast<LinearAlgebra::MiddleSizeVector::type *>(inout)[i] =
            std::max(
                std::real(
                    reinterpret_cast<LinearAlgebra::MiddleSizeVector::type *>(
                        inout)[i]),
                std::real(
                    reinterpret_cast<
                        const LinearAlgebra::MiddleSizeVector::type *>(in)[i]));
    }
}
#endif

/// \param[in] solutionVectorId Time level where the solution is stored.
/// \param[in] time Time corresponding to the current solution.
template <std::size_t DIM>
LinearAlgebra::MiddleSizeVector HpgemAPISimplified<DIM>::computeMaxError(
    const std::size_t solutionVectorId, const double time) {

    LinearAlgebra::MiddleSizeVector maxError(
        this->configData_->numberOfUnknowns_);
    maxError *= 0;

    for (Base::Element *ptrElement : this->meshes_[0]->getElementsList()) {
        const LinearAlgebra::MiddleSizeVector &solutionCoefficients =
            ptrElement->getTimeIntegrationVector(solutionVectorId);

        const LinearAlgebra::MiddleSizeVector maxErrorAtElement(
            computeMaxErrorAtElement(ptrElement, solutionCoefficients, time));

        for (std::size_t iV = 0; iV < this->configData_->numberOfUnknowns_;
             ++iV) {
            if (std::real(maxErrorAtElement(iV)) > std::real(maxError(iV))) {
                maxError(iV) = maxErrorAtElement(iV);
            }
        }
    }

#ifdef HPGEM_USE_MPI
    auto &comm = MPIContainer::Instance();
    MPI_Op myMax;
    MPI_Op_create(computeMPIMaximum, true, &myMax);
    comm.reduce(maxError, myMax);
    comm.broadcast(maxError);
#endif
    return maxError;
}

template <std::size_t DIM>
void HpgemAPISimplified<DIM>::computeRightHandSide(
    const std::size_t inputVectorId, const std::size_t resultVectorId,
    const double time) {
    // Apply the right hand side corresponding to integration on the elements.
    for (Base::Element *ptrElement : this->meshes_[0]->getElementsList()) {
        LinearAlgebra::MiddleSizeVector &inputFunctionCoefficients =
            ptrElement->getTimeIntegrationVector(inputVectorId);

        // Overwrite the vector at resultVectorId with the vector that is
        // returned by computeRightHandSideAtElement
        LinearAlgebra::MiddleSizeVector &resultFunctionCoefficients =
            ptrElement->getTimeIntegrationVector(resultVectorId);
        resultFunctionCoefficients = computeRightHandSideAtElement(
            ptrElement, inputFunctionCoefficients, time);
    }

    // Apply the right hand side corresponding to integration on the faces.
    for (Base::Face *ptrFace : this->meshes_[0]->getFacesList()) {
        if (ptrFace->isInternal()) {
            LinearAlgebra::MiddleSizeVector &inputFunctionCoefficientsLeft(
                ptrFace->getPtrElementLeft()->getTimeIntegrationVector(
                    inputVectorId));
            LinearAlgebra::MiddleSizeVector &inputFunctionCoefficientsRight(
                ptrFace->getPtrElementRight()->getTimeIntegrationVector(
                    inputVectorId));
            LinearAlgebra::MiddleSizeVector &resultFunctionCoefficientsLeft(
                ptrFace->getPtrElementLeft()->getTimeIntegrationVector(
                    resultVectorId));
            LinearAlgebra::MiddleSizeVector &resultFunctionCoefficientsRight(
                ptrFace->getPtrElementRight()->getTimeIntegrationVector(
                    resultVectorId));

            if (!computeBothFaces_) {
                resultFunctionCoefficientsLeft += computeRightHandSideAtFace(
                    ptrFace, Base::Side::LEFT, inputFunctionCoefficientsLeft,
                    inputFunctionCoefficientsRight, time);
                resultFunctionCoefficientsRight += computeRightHandSideAtFace(
                    ptrFace, Base::Side::RIGHT, inputFunctionCoefficientsLeft,
                    inputFunctionCoefficientsRight, time);
            } else {
                std::pair<LinearAlgebra::MiddleSizeVector,
                          LinearAlgebra::MiddleSizeVector>
                    resultFunctionCoefficients(computeBothRightHandSidesAtFace(
                        ptrFace, inputFunctionCoefficientsLeft,
                        inputFunctionCoefficientsRight, time));
                resultFunctionCoefficientsLeft +=
                    resultFunctionCoefficients.first;
                resultFunctionCoefficientsRight +=
                    resultFunctionCoefficients.second;
            }
        } else {
            LinearAlgebra::MiddleSizeVector &inputFunctionCoefficients(
                ptrFace->getPtrElementLeft()->getTimeIntegrationVector(
                    inputVectorId));
            LinearAlgebra::MiddleSizeVector &resultFunctionCoefficients(
                ptrFace->getPtrElementLeft()->getTimeIntegrationVector(
                    resultVectorId));

            resultFunctionCoefficients += computeRightHandSideAtFace(
                ptrFace, inputFunctionCoefficients, time);
        }
    }

    this->synchronize(resultVectorId);
}

template <std::size_t DIM>
LinearAlgebra::MiddleSizeVector
    HpgemAPISimplified<DIM>::getLinearCombinationOfVectors(
        const Base::Element *ptrElement,
        const std::vector<std::size_t> inputVectorIds,
        const std::vector<double> coefficientsInputVectors) {
    logger.assert_debug(
        inputVectorIds.size() == coefficientsInputVectors.size(),
        "Number of time levels and number of coefficients should be the same.");
    logger.assert_debug(inputVectorIds.size() > 0,
                        "Number of time levels should be bigger than zero.");

    LinearAlgebra::MiddleSizeVector linearCombination(
        ptrElement->getTimeIntegrationVector(inputVectorIds[0]));
    linearCombination *= coefficientsInputVectors[0];
    for (std::size_t i = 1; i < inputVectorIds.size(); i++) {
        linearCombination.axpy(
            coefficientsInputVectors[i],
            ptrElement->getTimeIntegrationVector(inputVectorIds[i]));
    }
    return linearCombination;
}

/// \details Make sure resultVectorId is different from the inputVectorIds.
/// \todo Merge this with the other computeRightHandSide to reduce code
/// duplication
template <std::size_t DIM>
void HpgemAPISimplified<DIM>::computeRightHandSide(
    const std::vector<std::size_t> inputVectorIds,
    const std::vector<double> coefficientsInputVectors,
    const std::size_t resultVectorId, const double time) {
    // Apply the right hand side corresponding to integration on the elements.
    for (Base::Element *ptrElement : this->meshes_[0]->getElementsList()) {
        LinearAlgebra::MiddleSizeVector inputFunctionCoefficients(
            getLinearCombinationOfVectors(ptrElement, inputVectorIds,
                                          coefficientsInputVectors));

        // Overwrite the vector at resultVectorId with the vector that is
        // returned by computeRightHandSideAtElement
        LinearAlgebra::MiddleSizeVector &resultFunctionCoefficients(
            ptrElement->getTimeIntegrationVector(resultVectorId));
        resultFunctionCoefficients = computeRightHandSideAtElement(
            ptrElement, inputFunctionCoefficients, time);
    }

    // Apply the right hand side corresponding to integration on the faces.
    for (Base::Face *ptrFace : this->meshes_[0]->getFacesList()) {
        if (ptrFace->isInternal()) {
            LinearAlgebra::MiddleSizeVector inputFunctionCoefficientsLeft(
                getLinearCombinationOfVectors(ptrFace->getPtrElementLeft(),
                                              inputVectorIds,
                                              coefficientsInputVectors));
            LinearAlgebra::MiddleSizeVector inputFunctionCoefficientsRight(
                getLinearCombinationOfVectors(ptrFace->getPtrElementRight(),
                                              inputVectorIds,
                                              coefficientsInputVectors));
            LinearAlgebra::MiddleSizeVector &resultFunctionCoefficientsLeft(
                ptrFace->getPtrElementLeft()->getTimeIntegrationVector(
                    resultVectorId));
            LinearAlgebra::MiddleSizeVector &resultFunctionCoefficientsRight(
                ptrFace->getPtrElementRight()->getTimeIntegrationVector(
                    resultVectorId));

            if (!computeBothFaces_) {
                resultFunctionCoefficientsLeft += computeRightHandSideAtFace(
                    ptrFace, Base::Side::LEFT, inputFunctionCoefficientsLeft,
                    inputFunctionCoefficientsRight, time);
                resultFunctionCoefficientsRight += computeRightHandSideAtFace(
                    ptrFace, Base::Side::RIGHT, inputFunctionCoefficientsLeft,
                    inputFunctionCoefficientsRight, time);
            } else {
                std::pair<LinearAlgebra::MiddleSizeVector,
                          LinearAlgebra::MiddleSizeVector>
                    resultFunctionCoefficients(computeBothRightHandSidesAtFace(
                        ptrFace, inputFunctionCoefficientsLeft,
                        inputFunctionCoefficientsRight, time));
                resultFunctionCoefficientsLeft +=
                    resultFunctionCoefficients.first;
                resultFunctionCoefficientsRight +=
                    resultFunctionCoefficients.second;
            }
        } else {
            LinearAlgebra::MiddleSizeVector inputFunctionCoefficients(
                getLinearCombinationOfVectors(ptrFace->getPtrElementLeft(),
                                              inputVectorIds,
                                              coefficientsInputVectors));
            LinearAlgebra::MiddleSizeVector &resultFunctionCoefficients(
                ptrFace->getPtrElementLeft()->getTimeIntegrationVector(
                    resultVectorId));

            resultFunctionCoefficients += computeRightHandSideAtFace(
                ptrFace, inputFunctionCoefficients, time);
        }
    }

    this->synchronize(resultVectorId);
}

template <std::size_t DIM>
void HpgemAPISimplified<DIM>::scaleVector(
    const std::size_t timeIntegrationVectorId, const double scale) {
    for (Base::Element *ptrElement : this->meshes_[0]->getElementsList()) {
        ptrElement->getTimeIntegrationVector(timeIntegrationVectorId) *= scale;
    }

    this->synchronize(timeIntegrationVectorId);
}

template <std::size_t DIM>
void HpgemAPISimplified<DIM>::scaleAndAddVector(
    const std::size_t vectorToChangeId, const std::size_t vectorToAddId,
    const double scale) {
    for (Base::Element *ptrElement : this->meshes_[0]->getElementsList()) {
        ptrElement->getTimeIntegrationVector(vectorToChangeId)
            .axpy(scale, ptrElement->getTimeIntegrationVector(vectorToAddId));
    }

    this->synchronize(vectorToChangeId);
}

template <std::size_t DIM>
std::tuple<double, double> HpgemAPISimplified<DIM>::computeErrorAndNormOfUpdate(
    double dt) {
    double error = 0.;
    double norm = 0.;
    for (Base::Element *element : this->meshes_[0]->getElementsList()) {
        const LinearAlgebra::MiddleSizeVector &v =
            element->getTimeIntegrationVector(solutionVectorId_);
        LinearAlgebra::MiddleSizeVector::type result;
        result = v.l2Norm();
        norm += std::real(result);
        LinearAlgebra::MiddleSizeVector deviation(
            element->getNumberOfBasisFunctions() *
            element->getNumberOfUnknowns());
        deviation *= 0.;
        for (std::size_t i = 0; i < ptrButcherTableau_->getNumberOfStages();
             ++i) {
            deviation +=
                element->getTimeIntegrationVector(auxiliaryVectorIds_[i]) * dt *
                ptrButcherTableau_->getErrorCoefficient(i);
        }
        LinearAlgebra::MiddleSizeVector::type result1;
        result1 = deviation.l2Norm();
        error += std::real(std::pow(result1, 2.));
    }
#ifdef HPGEM_USE_MPI
    auto &communicator = MPIContainer::Instance();
    ///\todo allreduce
    communicator.reduce(error, MPI_SUM);
    communicator.broadcast(error);
    communicator.reduce(norm, MPI_SUM);
    communicator.broadcast(norm);
#endif
    return std::make_tuple(std::sqrt(error), norm);
}

template <std::size_t DIM>
void HpgemAPISimplified<DIM>::setInitialSolution(
    const std::size_t solutionVectorId, const double initialTime,
    const std::size_t orderTimeDerivative) {
    integrateInitialSolution(solutionVectorId, initialTime,
                             orderTimeDerivative);
    solveMassMatrixEquations(solutionVectorId);
}

/// \details Computing the time derivative in this case means applying the right
/// hand side for the time integration vector with index 'inputVectorId' and
/// then solving the mass matrix equations. The result is stored at time
/// integration vector with index 'resultVectorId'.
template <std::size_t DIM>
void HpgemAPISimplified<DIM>::computeTimeDerivative(
    const std::size_t inputVectorId, const std::size_t resultVectorId,
    const double time) {
    computeRightHandSide(inputVectorId, resultVectorId, time);
    solveMassMatrixEquations(resultVectorId);
}

/// \details Computing the time derivative in this case means applying the right
/// hand side for the linear combination of time integration vectors with
/// indices 'inputVectorIds' and coefficients 'coefficientsInputVectors', and
/// then solving the mass matrix equations. The result is stored at time
/// integration vector with index 'resultVectorId'.
template <std::size_t DIM>
void HpgemAPISimplified<DIM>::computeTimeDerivative(
    const std::vector<std::size_t> inputVectorIds,
    const std::vector<double> coefficientsInputVectors,
    const std::size_t resultVectorId, const double time) {
    computeRightHandSide(inputVectorIds, coefficientsInputVectors,
                         resultVectorId, time);
    solveMassMatrixEquations(resultVectorId);
    this->synchronize(resultVectorId);
}

template <std::size_t DIM>
void HpgemAPISimplified<DIM>::computeOneTimeStep(double &time,
                                                 const double dt) {
    std::size_t numberOfStages = ptrButcherTableau_->getNumberOfStages();

    // Compute intermediate Runge-Kutta stages
    for (std::size_t iStage = 0; iStage < numberOfStages; iStage++) {
        const double stageTime = time + ptrButcherTableau_->getC(iStage) * dt;

        std::vector<std::size_t> inputVectorIds;
        std::vector<double> coefficientsInputVectors;

        inputVectorIds.push_back(solutionVectorId_);
        coefficientsInputVectors.push_back(1);
        for (std::size_t jStage = 0; jStage < iStage; jStage++) {
            inputVectorIds.push_back(auxiliaryVectorIds_[jStage]);
            coefficientsInputVectors.push_back(
                dt * ptrButcherTableau_->getA(iStage, jStage));
        }

        computeTimeDerivative(inputVectorIds, coefficientsInputVectors,
                              auxiliaryVectorIds_[iStage], stageTime);
    }

    // Update the solution
    for (std::size_t jStage = 0; jStage < numberOfStages; jStage++) {
        scaleAndAddVector(solutionVectorId_, auxiliaryVectorIds_[jStage],
                          dt * ptrButcherTableau_->getB(jStage));
    }

    // Update the time.
    time += dt;
}

template <std::size_t DIM>
void HpgemAPISimplified<DIM>::computeOneTimeStep(
    double &time, const double maximumRelativeError, const double dtMax) {
    logger.assert_debug(
        ptrButcherTableau_->hasErrorEstimate(),
        "Can only use an adaptive time step with an embedded butcher tableau");
    static double dtEstimate = dt.getValue();
    double dt = dtEstimate;
    if (dt > dtMax) dt = dtMax;
    double currentError = std::numeric_limits<double>::infinity();
    double solutionNorm = 0;
    while (maximumRelativeError * solutionNorm < currentError ||
           std::isinf(currentError)) {
        double timeCopy = time;

        // do a step, except don't update the time yet (we may need to rewind if
        // the error is too large)
        computeOneTimeStep(timeCopy, dt);
        std::tie(currentError, solutionNorm) = computeErrorAndNormOfUpdate(dt);
        if (maximumRelativeError * solutionNorm < currentError ||
            std::isinf(currentError)) {
            // should be rare enough to be more efficient than copying a back-up
            // every time step
            for (std::size_t jStage = 0;
                 jStage < ptrButcherTableau_->getNumberOfStages(); jStage++) {
                scaleAndAddVector(solutionVectorId_,
                                  auxiliaryVectorIds_[jStage],
                                  -dt * ptrButcherTableau_->getB(jStage));
            }
            dt = 0.8 * dt *
                 std::pow(maximumRelativeError / currentError * solutionNorm,
                          1. / ptrButcherTableau_->getOrder());
        }
    }
    dtEstimate = 0.9 * dt *
                 std::pow(maximumRelativeError / currentError * solutionNorm,
                          1. / ptrButcherTableau_->getOrder());
    logger(VERBOSE, "dt: %; new estimate: %", dt, dtEstimate);
    time += dt;
}

/// \param[in] outputFileName Name of the output file (minus extensions like
/// .dat). \param[in] internalFileTitle Title of the file as used by Tecplot
/// internally. \param[in] solutionTitle Title of the solution. \param[in]
/// variableNames String of variable names. The string should have the form
/// "nameVar1,nameVar2,..,nameVarN".
template <std::size_t DIM>
void HpgemAPISimplified<DIM>::setOutputNames(
    std::string outputFileName, std::string internalFileTitle,
    std::string solutionTitle, std::vector<std::string> variableNames) {
    outputFileName_ = outputFileName;
    internalFileTitle_ = internalFileTitle;
    solutionTitle_ = solutionTitle;
    logger.assert_debug(
        variableNames_.size() == this->configData_->numberOfUnknowns_,
        "Number of variable names (%) is not equal to the number of variables "
        "(%)",
        variableNames_.size(), this->configData_->numberOfUnknowns_);
    variableNames_ = variableNames;
}

template <std::size_t DIM>
void HpgemAPISimplified<DIM>::writeToTecplotFile(const Element *ptrElement,
                                                 const PointReferenceT &pRef,
                                                 std::ostream &out) {
    const std::size_t numberOfVariables = this->configData_->numberOfUnknowns_;

    const LinearAlgebra::MiddleSizeVector solution =
        ptrElement->getSolution(solutionVectorId_, pRef);

    std::size_t iV = 0;  // Index for the variable
    out << std::real(solution(iV));
    for (iV = 1; iV < numberOfVariables; iV++) {
        out << " " << std::real(solution(iV));
    }
}

template <std::size_t DIM>
void HpgemAPISimplified<DIM>::registerVTKWriteFunctions() {
    for (std::size_t iV = 0; iV < this->configData_->numberOfUnknowns_; iV++) {
        registerVTKWriteFunction(
            [=](Base::Element *element,
                const Geometry::PointReference<DIM> &pRef,
                std::size_t timeIntegrationVectorId) -> double {
                return std::real(
                    element->getSolution(timeIntegrationVectorId, pRef)[iV]);
            },
            variableNames_[iV]);
    }
}

template <std::size_t DIM>
bool HpgemAPISimplified<DIM>::checkBeforeSolving() {
    if (HpgemAPIBase<DIM>::meshes_.size() == 0) {
        logger(ERROR,
               "Error no mesh created : You need to create at least one mesh "
               "to solve a problem");
    }
    return true;
}

/// \brief Solve the PDE over the time domain [initialTime, finalTime].
/// \param[in] initialTime Initial time
/// \param[in] finalTime End time
/// \param[in] dt Size of the time step
/// \param[in] numberOutputFrames Number of times the solution is written to an
/// output file. \param[in] doComputeError Boolean to indicate if the error
/// should be computed. \todo Split in smaller functions, maybe split of
/// creating the output files?
template <std::size_t DIM>
bool HpgemAPISimplified<DIM>::solve(const double initialTime,
                                    const double finalTime, double dt,
                                    const std::size_t numberOfOutputFrames,
                                    bool doComputeError) {
    checkBeforeSolving();

    // Create output files for Paraview.
    const std::string outputFileNameVTK = outputFileName_;

    registerVTKWriteFunctions();
    Output::VTKTimeDependentWriter<DIM> VTKWriter(
        outputFileNameVTK, this->meshes_[0], this->polynomialOrder_);

    // Create output files for Tecplot.
#ifdef HPGEM_USE_MPI
    const std::string outputFileName =
        outputFileName_ + "." +
        std::to_string(Base::MPIContainer::Instance().getProcessorID());
#else

    const std::string outputFileName = outputFileName_;

#endif
    const std::string outputFileNameTecplot = outputFileName + ".dat";
    std::string dimensionsToWrite = "";
    for (std::size_t i = 0; i < DIM; i++) {
        dimensionsToWrite = dimensionsToWrite + std::to_string(i);
    }

    std::string variableString = variableNames_[0];
    for (std::size_t iV = 1; iV < variableNames_.size(); iV++) {
        variableString = variableString + "," + variableNames_[iV];
    }

    std::ofstream outputFile(outputFileNameTecplot);
    Output::TecplotDiscontinuousSolutionWriter<DIM> tecplotWriter(
        outputFile, internalFileTitle_, dimensionsToWrite, variableString);

    // Compute parameters for time integration
    const double T = finalTime - initialTime;  // Time interval
    const double outputDt = T / numberOfOutputFrames;
    double outputTime = outputDt + initialTime;

    // prefer a constant time step for non-adaptive time-stepping
    std::size_t numberOfTimeSteps = std::ceil(T / dt);
    std::size_t numberOfTimeStepsForOutput;
    if (numberOfOutputFrames > 0) {
        // Round off to above such that the number of time steps is a multiple
        // of the number of output frames.
        numberOfTimeSteps += (numberOfOutputFrames -
                              (numberOfTimeSteps % numberOfOutputFrames)) %
                             numberOfOutputFrames;

        // Recompute dt.
        dt = T / numberOfTimeSteps;

        // Compute the number of timesteps after which to create an output
        // frame.
        numberOfTimeStepsForOutput =
            (std::size_t)(numberOfTimeSteps / numberOfOutputFrames);
    }

    // if the user wants an adaptive time step, the default of 0.01 for dt is
    // silly
    if (!::hpgem::Base::dt.isUsed() && error.isUsed()) {
        dt = std::numeric_limits<double>::infinity();
    }
    double maximumRelativeError = error.getValue();

    // Set the initial time.
    double time = initialTime;

    // Create and Store things before solving the problem.
    tasksBeforeSolving();

    // Set the initial numerical solution.
    logger(INFO, "Computing and interpolating the initial solution.");
    setInitialSolution(solutionVectorId_, time, 0);
    tecplotWriter.write(this->meshes_[0], solutionTitle_, false, this, time);
    VTKWrite(VTKWriter, time, solutionVectorId_);

    // Solve the system of PDE's.
    logger(INFO, "Solving the system of PDE's.");
    logger(INFO, "Maximum time step: %.", dt);
    logger(INFO, "Maximum relative error: %.", maximumRelativeError);
    logger(INFO, "Minimum number of time steps: %.", numberOfTimeSteps);
    if (numberOfOutputFrames > 0) {
        logger(INFO, "Minimum number of time steps for output: %.",
               numberOfTimeStepsForOutput);
    }
    std::size_t actualNumberOfTimeSteps = 0;
    while (time < finalTime - 1e-14) {
        actualNumberOfTimeSteps++;
        if (error.isUsed()) {
            computeOneTimeStep(
                time, maximumRelativeError,
                std::min(std::min(dt, outputTime - time), finalTime - time));
        } else {
            computeOneTimeStep(time, dt);
        }

        if (time > outputTime - 1e-12) {
            outputTime += outputDt;
            tecplotWriter.write(this->meshes_[0], solutionTitle_, false, this,
                                time);
            VTKWrite(VTKWriter, time, solutionVectorId_);
        }
        showProgress(time, actualNumberOfTimeSteps);
    }
    logger(INFO, "Actual number of time steps: %.", actualNumberOfTimeSteps);
    if (error.isUsed()) {
        logger(
            INFO, "Maximal relative error induced by time stepping: %",
            std::pow(1. + maximumRelativeError, actualNumberOfTimeSteps) - 1.);
    }

    tasksAfterSolving();

    // Compute the energy norm of the error
    if (doComputeError) {
        const LinearAlgebra::MiddleSizeVector::type totalError =
            computeTotalError(solutionVectorId_, finalTime);
        logger(INFO, "Total error: %.", totalError);
        const LinearAlgebra::MiddleSizeVector maxError =
            computeMaxError(solutionVectorId_, finalTime);
        logger.assert_debug(
            maxError.size() == this->configData_->numberOfUnknowns_,
            "Size of maxError (%) not equal to the number of variables (%)",
            maxError.size(), this->configData_->numberOfUnknowns_);
        for (std::size_t iV = 0; iV < this->configData_->numberOfUnknowns_;
             iV++) {
            logger(INFO, "Maximum error %: %", variableNames_[iV],
                   maxError(iV));
        }
    }

    return true;
}
}  // namespace Base
}  // namespace hpgem
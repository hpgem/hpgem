/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2021, University of Twente
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

#ifndef HPGEM_L2PROJECTIONERRORQUALITYMETRIC_H
#define HPGEM_L2PROJECTIONERRORQUALITYMETRIC_H

#include "QualityMetricComputation.h"

#include "Base/CommandLineOptions.h"
#include "Base/H1ConformingTransformation.h"
#include "Base/MeshManipulator.h"
#include "Geometry/PointPhysical.h"
#include "Integration/ElementIntegral.h"
#include "Integration/QuadratureRules/AllGaussQuadratureRules.h"
#include "Output/VTKSpecificTimeWriter.h"
#include "Utilities/GlobalIndexing.h"
#include "Utilities/GlobalMatrix.h"
#include "Utilities/GlobalVector.h"

template <std::size_t dimension>
class L2ProjectionErrorQualityMetric
    : public QualityMetricComputation<dimension> {

    // Timeintegration vector to use for the solution
    const std::size_t SOLUTION_VECTOR_ID = 0;
    // Type of the value from solutions, may be complex.
    using VALUE_TYPE = hpgem::LinearAlgebra::MiddleSizeVector::type;

   public:
    // Note: We may have complex valued functions
    using TestingFunction =
        std::function<double(const hpgem::Geometry::PointPhysical<dimension>&)>;

    enum class BasisFunctionType { DISCONTINUOUS, CONFORMING };

    L2ProjectionErrorQualityMetric(BasisFunctionType bfType, std::size_t order,
                                   std::string name,
                                   std::vector<TestingFunction> functions)
        : basisFunctionType_(bfType),
          order_(order),
          name_(std::move(name)),
          testingFunctions_(functions) {}

    void computeAndPlotMetric(
        hpgem::Base::MeshManipulator<dimension>& mesh,
        hpgem::Output::VTKSpecificTimeWriter<dimension>& plotter,
        const std::string& filePrefix) final;

    void setTestingFunctions(std::vector<TestingFunction> functions) {
        testingFunctions_ = std::move(functions);
    }

   private:
    double computeElementError(
        std::size_t functionId, hpgem::Base::Element* element,
        hpgem::Integration::ElementIntegral<dimension>& integral);

    hpgem::LinearAlgebra::MiddleSizeVector computeElementVector(
        std::size_t functionId, hpgem::Base::Element* element,
        hpgem::Integration::ElementIntegral<dimension>& integral);

    VALUE_TYPE computeLocalError(
        std::size_t functionId, hpgem::Base::Element* element,
        const hpgem ::Geometry::PointReference<dimension>& p);

    std::vector<TestingFunction> testingFunctions_;
    BasisFunctionType basisFunctionType_;
    std::size_t order_;
    std::string name_;
};

template <std::size_t dimension>
void L2ProjectionErrorQualityMetric<dimension>::computeAndPlotMetric(
    hpgem::Base::MeshManipulator<dimension>& mesh,
    hpgem::Output::VTKSpecificTimeWriter<dimension>& plotter,
    const std::string& filePrefix) {

    if (testingFunctions_.empty()) {
        return;
    }
    using namespace hpgem;

    logger(INFO, "Setting up L2-projection-%", name_);

    std::map<const Base::Element*, double> elementErrorSquared;

    switch (basisFunctionType_) {
        case BasisFunctionType::CONFORMING:
            mesh.useDefaultConformingBasisFunctions(order_);
            break;
        case BasisFunctionType::DISCONTINUOUS:
            mesh.useDefaultDGBasisFunctions(order_);
            break;
        default:
            logger.assert_always(false, "Unknown basis function type");
    }

    Integration::ElementIntegral<dimension> integral;
    integral.setTransformation(
        std::shared_ptr<Base::CoordinateTransformation<dimension>>(
            new Base::H1ConformingTransformation<dimension>()));

    const std::size_t ELEM_MAT_ID = 0;
    const std::size_t ELEM_VEC_ID = 0;

    {
        // Preparation: Set up the mass matrix
        std::vector<double> phis;
        for (Base::Element* element : mesh.getElementsList()) {
            element->setNumberOfTimeIntegrationVectors(1);
            LinearAlgebra::MiddleSizeMatrix mat = integral.integrate(
                element, [&phis](Base::PhysicalElement<dimension>& pelemnt) {
                    std::size_t n = pelemnt.getNumberOfBasisFunctions();
                    LinearAlgebra::MiddleSizeMatrix mat(n, n);
                    phis.resize(n);
                    for (std::size_t i = 0; i < n; ++i) {
                        phis[i] = pelemnt.basisFunction(i);
                    }
                    for (std::size_t i = 0; i < n; ++i) {
                        for (std::size_t j = i; j < n; ++j) {
                            mat(i, j) = phis[i] * phis[j];
                            mat(j, i) = mat(i, j);
                        }
                    }
                    return mat;
                });
            element->setElementMatrix(mat, ELEM_MAT_ID);
            LinearAlgebra::MiddleSizeVector tempSource(
                element->getTotalNumberOfBasisFunctions());
            element->setElementVector(tempSource, ELEM_VEC_ID);
        }
    }

    // Preparation: Set up the global solver
    Utilities::GlobalIndexing indexing(&mesh);
    Utilities::GlobalPetscMatrix massMatrix(indexing, ELEM_MAT_ID, -1);
    Utilities::GlobalPetscVector loadVector(indexing, ELEM_VEC_ID, -1);
    Utilities::GlobalPetscVector resultVector(indexing, -1, -1);
    KSP ksp;
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPCG);
    KSPSetTolerances(ksp, 1e-12, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);

    std::stringstream detailPlotFile;
    detailPlotFile << filePrefix << name_;
    Output::VTKSpecificTimeWriter<dimension> detailPlotter(detailPlotFile.str(),
                                                           &mesh, 0, order_);

    // Compute the element error for each test function
    for (std::size_t i = 0; i < testingFunctions_.size(); ++i) {
        logger(INFO, "Creating sources for %-%", name_, i);

        for (Base::Element* element : mesh.getElementsList()) {
            element->setElementVector(
                computeElementVector(i, element, integral), ELEM_VEC_ID);
        }
        loadVector.reinit();
        logger(INFO, "Solving for %-%", name_, i);
        KSPSetOperators(ksp, massMatrix, massMatrix);
        KSPSolve(ksp, loadVector, resultVector);
        resultVector.writeTimeIntegrationVector(SOLUTION_VECTOR_ID);
        for (Base::Element* element : mesh.getElementsList()) {
            elementErrorSquared[element] +=
                computeElementError(i, element, integral);
        }

        logger(INFO, "Creating output for %-%", name_, i);

        // Plot the error of each individual function
        std::stringstream errorName;
        errorName << name_ << "-error-" << i;
        detailPlotter.write(
            [&](Base::Element* element,
                const Geometry::PointReference<dimension>& p, std::size_t) {
                return std::real(computeLocalError(i, element, p));
            },
            errorName.str());
    }

    // Plot the result
    plotter.write(
        [&elementErrorSquared](
            Base::Element* element, const Geometry::PointReference<dimension>&,
            std::size_t) { return std::sqrt(elementErrorSquared[element]); },
        name_);
}

template <std::size_t dimension>
hpgem::LinearAlgebra::MiddleSizeVector
    L2ProjectionErrorQualityMetric<dimension>::computeElementVector(
        std::size_t functionId, hpgem::Base::Element* element,
        hpgem::Integration::ElementIntegral<dimension>& integral) {
    using namespace hpgem;
    // Standard integral (f, phi_i)_element
    return integral.integrate(
        element, [&](Base::PhysicalElement<dimension>& pelement) {
            LinearAlgebra::MiddleSizeVector result(
                pelement.getTotalNumberOfBasisFunctions());
            double functionValue =
                testingFunctions_[functionId](pelement.getPointPhysical());
            for (std::size_t i = 0; i < result.size(); ++i) {
                result[i] = functionValue * pelement.basisFunction(i);
            }
            return result;
        });
}

template <std::size_t dimension>
double L2ProjectionErrorQualityMetric<dimension>::computeElementError(
    std::size_t functionId, hpgem::Base::Element* element,
    hpgem::Integration::ElementIntegral<dimension>& integral) {
    using namespace hpgem;
    double error = integral.integrate(
        element,
        [&](Base::PhysicalElement<dimension>& pelement) {
            VALUE_TYPE localError = computeLocalError(
                functionId, element, pelement.getPointReference());
            // When using complex petsc the error may include imaginary part.
            return std::real(localError) * std::real(localError) +
                   std::imag(localError) * std::imag(localError);
        },
        // Higher quadrature order than needed to improve accuracy
        QuadratureRules::AllGaussQuadratureRules::instance().getRule(
            element->getReferenceGeometry(), 2 * order_ + 4));
    return error;
}

template <std::size_t dimension>
typename L2ProjectionErrorQualityMetric<dimension>::VALUE_TYPE
    L2ProjectionErrorQualityMetric<dimension>::computeLocalError(
        std::size_t functionId, hpgem::Base::Element* element,
        const hpgem ::Geometry::PointReference<dimension>& p) {
    using namespace hpgem;

    VALUE_TYPE error =
        testingFunctions_[functionId](element->getReferenceToPhysicalMap()
                                          ->castDimension<dimension>()
                                          .transform(p));
    LinearAlgebra::MiddleSizeVector& coefficients =
        element->getTimeIntegrationVector(SOLUTION_VECTOR_ID);
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        double phi_i = element->basisFunction(i, p);

        error -= coefficients[i] * phi_i;
    }
    // In case complex valued
    return error;
}

#endif  // HPGEM_L2PROJECTIONERRORQUALITYMETRIC_H

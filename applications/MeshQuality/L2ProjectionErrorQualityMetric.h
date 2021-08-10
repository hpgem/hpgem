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

#include "GlobalSolveQualityMetric.h"

template <std::size_t dimension>
class L2ProjectionErrorQualityMetric
    : public GlobalSolveQualityMetric<dimension> {

   public:
    // Note: We may have complex valued functions
    using TestingFunction =
        std::function<double(const hpgem::Geometry::PointPhysical<dimension>&)>;

    enum class BasisFunctionType { DISCONTINUOUS, CONFORMING };

    L2ProjectionErrorQualityMetric(BasisFunctionType bfType, std::size_t order,
                                   std::string name,
                                   std::vector<TestingFunction> functions)
        : GlobalSolveQualityMetric<dimension>(order, std::move(name)),
          basisFunctionType_(bfType),
          testingFunctions_(functions) {

        integrator_.setTransformation(
            std::make_shared<
                hpgem::Base::H1ConformingTransformation<dimension>>());
    }

   private:
    bool usesFaceMatrixAndVector() const final { return false; }

    std::size_t numberOfSolves() const final {
        return testingFunctions_.size();
    }

    void computeLocalMatrices(
        hpgem::Base::MeshManipulator<dimension>& mesh) final;
    void computeLocalVectors(
        std::size_t solve, hpgem::Base::MeshManipulator<dimension>& mesh) final;

    double computeSolutionValue(
        hpgem::Base::Element* element,
        const hpgem::Geometry::PointReference<dimension>& p) const final;

    double computeFunctionValue(
        std::size_t solve,
        const hpgem::Geometry::PointPhysical<dimension>& p) const final {
        return testingFunctions_[solve](p);
    };

    hpgem::Integration::ElementIntegral<dimension>& getIntegrator() {
        return integrator_;
    }

   private:

    std::vector<TestingFunction> testingFunctions_;
    BasisFunctionType basisFunctionType_;
    hpgem::Integration::ElementIntegral<dimension> integrator_;
};

template <std::size_t dimension>
void L2ProjectionErrorQualityMetric<dimension>::computeLocalMatrices(
    hpgem::Base::MeshManipulator<dimension>& mesh) {
    using namespace hpgem;

    switch (basisFunctionType_) {
        case BasisFunctionType::CONFORMING:
            mesh.useDefaultConformingBasisFunctions(this->getOrder());
            break;
        case BasisFunctionType::DISCONTINUOUS:
            mesh.useDefaultDGBasisFunctions(this->getOrder());
            break;
        default:
            logger.assert_always(false, "Unknown basis function type");
    }

    std::vector<double> phis;
    for (Base::Element* element : mesh.getElementsList()) {
        element->setNumberOfTimeIntegrationVectors(1);
        LinearAlgebra::MiddleSizeMatrix mat = integrator_.integrate(
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
        element->setElementMatrix(mat, this->ELEMENT_MATRIX_ID);
    }
}

template <std::size_t dimension>
void L2ProjectionErrorQualityMetric<dimension>::computeLocalVectors(
    std::size_t solve, hpgem::Base::MeshManipulator<dimension>& mesh) {

    using namespace hpgem;

    for (Base::Element* element : mesh.getElementsList()) {
        // Standard integral (f, phi_i)_element
        LinearAlgebra::MiddleSizeVector vec = integrator_.integrate(
            element, [&](Base::PhysicalElement<dimension>& pelement) {
                LinearAlgebra::MiddleSizeVector result(
                    pelement.getTotalNumberOfBasisFunctions());
                double functionValue =
                    testingFunctions_[solve](pelement.getPointPhysical());
                for (std::size_t i = 0; i < result.size(); ++i) {
                    result[i] = functionValue * pelement.basisFunction(i);
                }
                return result;
            });
        element->setElementVector(vec, this->ELEMENT_VECTOR_ID);
    }
}

template <std::size_t dimension>
double L2ProjectionErrorQualityMetric<dimension>::computeSolutionValue(
    hpgem::Base::Element* element,
    const hpgem::Geometry::PointReference<dimension>& p) const {

    using namespace hpgem;
    double value;

    LinearAlgebra::MiddleSizeVector& coefficients =
        element->getTimeIntegrationVector(this->SOLUTION_VECTOR_ID);
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        double phi_i = element->basisFunction(i, p);

        value += std::real(coefficients[i]) * phi_i;
    }
    return value;
}

#endif  // HPGEM_L2PROJECTIONERRORQUALITYMETRIC_H

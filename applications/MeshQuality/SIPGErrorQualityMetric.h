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

#ifndef HPGEM_SIPGERRORQUALITYMETRIC_H
#define HPGEM_SIPGERRORQUALITYMETRIC_H

#include "QualityMetricComputation.h"

#include <Integration/FaceIntegral.h>

#include <map>
#include <memory>
#include <string>
#include <vector>

template <std::size_t dimension>
class Func {
   public:
    virtual double value(
        const hpgem::Geometry::PointPhysical<dimension>& p) const = 0;
    virtual double laplacian(
        const hpgem::Geometry::PointPhysical<dimension>& p) const = 0;
};

template <std::size_t dimension>
class SIPGErrorQualityMetric : public GlobalSolveQualityMetric<dimension> {
    // Generally 10 is sufficient, just a bit more to be safe
    const double PENALTY_CONSTANT = 10;

    using VALUE_TYPE = hpgem::LinearAlgebra::MiddleSizeVector::type;

   public:
    SIPGErrorQualityMetric(
        std::size_t order, std::string name,
        std::vector<std::shared_ptr<Func<dimension>>> functions)
        : GlobalSolveQualityMetric<dimension>(order, std::move(name)),
          functions_(functions) {

        auto transform = std::make_shared<
            hpgem::Base::H1ConformingTransformation<dimension>>();

        elementIntegrator_.setTransformation(transform);
        faceIntegrator_.setTransformation(transform);
    };

   protected:
    bool usesFaceMatrixAndVector() const final { return true; }
    std::size_t numberOfSolves() const final { return functions_.size(); }

    void computeLocalMatrices(
        hpgem::Base::MeshManipulator<dimension>& mesh) final;
    void computeLocalVectors(
        std::size_t solve,
        hpgem::Base::MeshManipulator<dimension>& mesh) final;

    virtual double computeSolutionValue(
        hpgem::Base::Element* element,
        const hpgem::Geometry::PointReference<dimension>& p) const final;
    virtual double computeFunctionValue(
        std::size_t solve,
        const hpgem::Geometry::PointPhysical<dimension>& p) const final {
        return functions_[solve]->value(p);
    }
    virtual hpgem::Integration::ElementIntegral<dimension>&
        getIntegrator() final {
        return elementIntegrator_;
    }

   private:
    // The element matrix with the term
    // integral(grad phi_i . grad phi_j)
    hpgem::LinearAlgebra::MiddleSizeMatrix computeElementMatrix(
        hpgem::Base::Element* element);

    // Element vector integral(phi_i f)
    hpgem::LinearAlgebra::MiddleSizeVector computeElementVector(
        std::size_t functionId, hpgem::Base::Element* element);

    // Face matrix
    // integral( -{{grad phi_i}} . [[phi_j]] - [[phi_i]] . {{grad phi_j}}
    //           + alpha/h [[phi_i]] [[phi_j]])
    hpgem::LinearAlgebra::MiddleSizeMatrix computeFaceMatrix(
        hpgem::Base::Face* face);

    // Face vector:
    // Internal faces: 0
    // Boundary faces integral( f(-grad phi_i . n + alpha/h phi_i))
    hpgem::LinearAlgebra::MiddleSizeVector computeFaceVector(
        std::size_t functionId, hpgem::Base::Face* face);

    std::vector<std::shared_ptr<Func<dimension>>> functions_;
    hpgem::Integration::ElementIntegral<dimension> elementIntegrator_;
    hpgem::Integration::FaceIntegral<dimension> faceIntegrator_;
};

template <std::size_t dimension>
void SIPGErrorQualityMetric<dimension>::computeLocalMatrices(
    hpgem::Base::MeshManipulator<dimension>& mesh) {
    using namespace hpgem;

    mesh.useDefaultDGBasisFunctions(this->getOrder());

    for (Base::Element* element : mesh.getElementsList()) {
        element->setElementMatrix(computeElementMatrix(element),
                                  this->ELEMENT_MATRIX_ID);
    }
    for (Base::Face* face : mesh.getFacesList()) {
        face->setFaceMatrix(computeFaceMatrix(face), this->FACE_MATRIX_ID);
    }
}

template <std::size_t dimension>
void SIPGErrorQualityMetric<dimension>::computeLocalVectors(
    std::size_t solve, hpgem::Base::MeshManipulator<dimension>& mesh) {
    using namespace hpgem;
    for (Base::Element* element : mesh.getElementsList()) {
        element->setElementVector(computeElementVector(solve, element),
                                  this->ELEMENT_VECTOR_ID);
    }
    for (Base::Face* face : mesh.getFacesList()) {
        face->setFaceVector(computeFaceVector(solve, face), this->FACE_VECTOR_ID);
    }
}

template <std::size_t dimension>
hpgem::LinearAlgebra::MiddleSizeMatrix
    SIPGErrorQualityMetric<dimension>::computeElementMatrix(
        hpgem::Base::Element* element) {
    using namespace hpgem;

    return elementIntegrator_.integrate(
        element, [](Base::PhysicalElement<dimension>& pelement) {
            std::size_t n = pelement.getNumberOfBasisFunctions();
            LinearAlgebra::MiddleSizeMatrix mat(n, n);
            for (std::size_t i = 0; i < n; ++i) {
                for (std::size_t j = 0; j < n; ++j) {
                    mat(i, j) = pelement.basisFunctionDeriv(i) *
                                pelement.basisFunctionDeriv(j);
                }
            }
            return mat;
        });
}

template <std::size_t dimension>
hpgem::LinearAlgebra::MiddleSizeVector
    SIPGErrorQualityMetric<dimension>::computeElementVector(
        std::size_t functionId, hpgem::Base::Element* element) {
    using namespace hpgem;
    return elementIntegrator_.integrate(
        element, [&](Base::PhysicalElement<dimension>& pelement) {
            std::size_t n = pelement.getNumberOfBasisFunctions();
            LinearAlgebra::MiddleSizeVector vec(n);
            double fvalue = -1.0 * functions_[functionId]->laplacian(
                                       pelement.getPointPhysical());
            for (std::size_t i = 0; i < n; ++i) {
                vec[i] = fvalue * pelement.basisFunction(i);
            }
            return vec;
        });
}

template <std::size_t dimension>
hpgem::LinearAlgebra::MiddleSizeMatrix
    SIPGErrorQualityMetric<dimension>::computeFaceMatrix(
        hpgem::Base::Face* face) {
    using namespace hpgem;

    return faceIntegrator_.integrate(
        face, [&](Base::PhysicalFace<dimension>& pface) {
            std::size_t n = face->getNumberOfBasisFunctions();
            LinearAlgebra::MiddleSizeMatrix mat(n, n);
            for (std::size_t i = 0; i < n; ++i) {
                const LinearAlgebra::SmallVector<dimension>& phiNormalI =
                    pface.basisFunctionUnitNormal(i);
                const LinearAlgebra::SmallVector<dimension>& phiDerivI =
                    pface.basisFunctionDeriv(i);

                for (std::size_t j = 0; j < n; ++j) {
                    const LinearAlgebra::SmallVector<dimension>& phiNormalJ =
                        pface.basisFunctionUnitNormal(j);
                    const LinearAlgebra::SmallVector<dimension>& phiDerivJ =
                        pface.basisFunctionDeriv(j);

                    if (face->isInternal()) {
                        mat(j, i) =
                            -(phiNormalI * phiDerivJ + phiNormalJ * phiDerivI) /
                                2 +
                            PENALTY_CONSTANT / face->getDiameter() *
                                phiNormalI * phiNormalJ;
                    } else {
                        mat(j, i) =
                            -(phiNormalI * phiDerivJ + phiNormalJ * phiDerivI) +
                            PENALTY_CONSTANT / face->getDiameter() *
                                phiNormalI * phiNormalJ;
                    }
                }
            }
            return mat;
        });
}

template <std::size_t dimension>
hpgem::LinearAlgebra::MiddleSizeVector
    SIPGErrorQualityMetric<dimension>::computeFaceVector(
        std::size_t functionId, hpgem::Base::Face* face) {
    using namespace hpgem;
    if (face->isInternal()) {
        return LinearAlgebra::MiddleSizeVector(
            face->getNumberOfBasisFunctions());
    }

    return faceIntegrator_.integrate(
        face, [&](Base::PhysicalFace<dimension>& pface) {
            std::size_t n = pface.getNumberOfBasisFunctions();
            LinearAlgebra::MiddleSizeVector vec(n);

            double fvalue =
                functions_[functionId]->value(pface.getPointPhysical());
            const LinearAlgebra::SmallVector<dimension>& normal =
                pface.getUnitNormalVector();
            for (std::size_t i = 0; i < n; ++i) {
                vec(i) += fvalue * (
                                       // Symmetry term
                                       -pface.basisFunctionDeriv(i) * normal +
                                       // Penalty term
                                       PENALTY_CONSTANT / face->getDiameter() *
                                           pface.basisFunction(i));
            }
            return vec;
        });
}

template <std::size_t dimension>
double SIPGErrorQualityMetric<dimension>::computeSolutionValue(
    hpgem::Base::Element* element,
    const hpgem::Geometry::PointReference<dimension>& p) const {
    using namespace hpgem;

    VALUE_TYPE value = 0.0;
    LinearAlgebra::MiddleSizeVector& coefficients =
        element->getTimeIntegrationVector(this->SOLUTION_VECTOR_ID);
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        double phi_i = element->basisFunction(i, p);

        value += coefficients[i] * phi_i;
    }
    // In case complex valued
    return std::real(value);
}

#endif  // HPGEM_SIPGERRORQUALITYMETRIC_H

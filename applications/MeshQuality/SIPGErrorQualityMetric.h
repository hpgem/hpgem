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
class SIPGErrorQualityMetric : public QualityMetricComputation<dimension> {
    // Timeintegration vector to use for the solution
    const std::size_t SOLUTION_VECTOR_ID = 0;
    const std::size_t ELEMENT_MATRIX = 0;
    const std::size_t ELEMENT_VECTOR = 0;
    const std::size_t FACE_MATRIX = 0;
    const std::size_t FACE_VECTOR = 0;
    // Generally 10 is sufficient, just a bit more to be safe
    const double PENALTY_CONSTANT = 12;

    using VALUE_TYPE = hpgem::LinearAlgebra::MiddleSizeVector::type;

   public:
    SIPGErrorQualityMetric(std::size_t order, std::string name,
                           std::vector<std::shared_ptr<Func<dimension>>> functions)
        : order_(order), name_(name), functions_(functions){};

    void computeAndPlotMetric(
        hpgem::Base::MeshManipulator<dimension>& mesh,
        hpgem::Output::VTKSpecificTimeWriter<dimension>& plotter) final;

   private:
    double computeElementError(
        std::size_t functionId, hpgem::Base::Element* element,
        hpgem::Integration::ElementIntegral<dimension>& integral);

    hpgem::LinearAlgebra::MiddleSizeMatrix computeElementMatrix(
        hpgem::Base::Element* element,
        hpgem::Integration::ElementIntegral<dimension>& integral);

    hpgem::LinearAlgebra::MiddleSizeVector computeElementVector(
        std::size_t functionId, hpgem::Base::Element* element,
        hpgem::Integration::ElementIntegral<dimension>& integral);

    hpgem::LinearAlgebra::MiddleSizeMatrix computeFaceMatrix(
        hpgem::Base::Face* face,
        hpgem::Integration::FaceIntegral<dimension>& integral);

    hpgem::LinearAlgebra::MiddleSizeVector computeFaceVector(
        std::size_t functionId, hpgem::Base::Face* face,
        hpgem::Integration::FaceIntegral<dimension>& integral);

    VALUE_TYPE computeLocalError(
        std::size_t functionId, hpgem::Base::Element* element,
        const hpgem ::Geometry::PointReference<dimension>& p);

    std::size_t order_;
    std::string name_;
    std::vector<std::shared_ptr<Func<dimension>>> functions_;
};

template <std::size_t dimension>
void SIPGErrorQualityMetric<dimension>::computeAndPlotMetric(
    hpgem::Base::MeshManipulator<dimension>& mesh,
    hpgem::Output::VTKSpecificTimeWriter<dimension>& plotter) {
    if (functions_.empty()) {
        return;
    }

    using namespace hpgem;

    std::map<const Base::Element*, double> elementErrorSquared;
    mesh.useDefaultDGBasisFunctions(order_);

    Integration::ElementIntegral<dimension> integral;
    integral.setTransformation(
        std::make_shared<Base::H1ConformingTransformation<dimension>>());

    // Set up matrices

    // Element matrix = grad phi_i . grad phi_j
    for (Base::Element* element : mesh.getElementsList()) {
        element->setNumberOfTimeIntegrationVectors(1);
        element->setElementMatrix(computeElementMatrix(element, integral),
                                  ELEMENT_MATRIX);
        // Temporary vector to make assembly easier
        LinearAlgebra::MiddleSizeVector tempVec(
            element->getTotalNumberOfBasisFunctions());
        element->setElementVector(tempVec, ELEMENT_VECTOR);
    }
    // Face matrix

    Integration::FaceIntegral<dimension> faceIntegral;
    faceIntegral.setTransformation(
        std::make_shared<Base::H1ConformingTransformation<dimension>>());
    // Face matrix -{{grad phi_i}}[[phi_j]] -[[phi_j]]{{grad phi_i}}
    //   + alpha/h [[phi_i]] [[phi_j]]
    for (Base::Face* face : mesh.getFacesList()) {
        face->setFaceMatrix(computeFaceMatrix(face, faceIntegral), FACE_MATRIX);
        // Temporary
        face->setFaceVector(LinearAlgebra::MiddleSizeVector(
                                face->getNumberOfBasisFunctions()),
                            FACE_VECTOR);
    }

    // Preparation: Set up the global solver
    Utilities::GlobalIndexing indexing(&mesh);
    Utilities::GlobalPetscMatrix stiffnessMatrix(indexing, ELEMENT_MATRIX,
                                                 FACE_MATRIX);
    Utilities::GlobalPetscVector loadVector(indexing, ELEMENT_VECTOR,
                                            FACE_VECTOR);
    Utilities::GlobalPetscVector resultVector(indexing, -1, -1);
    KSP ksp;
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPPREONLY);
    KSPSetOperators(ksp, stiffnessMatrix, stiffnessMatrix);

    for (std::size_t fid = 0; fid < functions_.size(); ++fid) {
        // Compute local parts of the load vector

        // Element vector = phi_i . j
        for (Base::Element* element : mesh.getElementsList()) {
            element->setElementVector(
                computeElementVector(fid, element, integral), ELEMENT_VECTOR);
        }
        // Face vector (alpha/h[[phi_j]] . n - {{grad phi_j}}) . g
        for (Base::Face* face : mesh.getFacesList()) {
            face->setFaceVector(computeFaceVector(fid, face, faceIntegral),
                                FACE_VECTOR);
        }
        // Assemble load vector & solve the system
        loadVector.reinit();

        KSPSolve(ksp, loadVector, resultVector);
        resultVector.writeTimeIntegrationVector(SOLUTION_VECTOR_ID);
        for (Base::Element* element : mesh.getElementsList()) {
            elementErrorSquared[element] +=
                computeElementError(fid, element, integral);
        }

        // Plot the error of each individual function
        std::stringstream errorName;
        errorName << name_ << "-error-" << fid;
        plotter.write(
            [&](Base::Element* element,
                const Geometry::PointReference<dimension>& p, std::size_t) {
              LinearAlgebra::MiddleSizeVector coefficients =
                  element->getTimeIntegrationVector(SOLUTION_VECTOR_ID);
              double value = 0;
              for (std::size_t i = 0; i < coefficients.size(); ++i) {
                  value += std::real(coefficients[i]) *
                           element->basisFunction(i, p);
              }
              value -= std::real(
                  functions_[fid]->value(element->referenceToPhysical(p)));
              return value;
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
hpgem::LinearAlgebra::MiddleSizeMatrix
    SIPGErrorQualityMetric<dimension>::computeElementMatrix(
        hpgem::Base::Element* element,
        hpgem::Integration::ElementIntegral<dimension>& integral) {
    using namespace hpgem;

    return integral.integrate(
        element, [](Base::PhysicalElement<dimension>& pelement) {
            std::size_t n = pelement.getNumberOfBasisFunctions();
            LinearAlgebra::MiddleSizeMatrix mat(n, n);
            for (std::size_t i = 0; i < n; ++i) {
                const LinearAlgebra::SmallVector<dimension>& gradPhii =
                    pelement.basisFunctionDeriv(i);
                for (std::size_t j = i; j < n; ++j) {
                    double value = gradPhii * pelement.basisFunctionDeriv(j);
                    mat(i, j) = value;
                    if (i != j) {
                        mat(j, i) = value;
                    }
                }
            }
            return mat;
        });
}

template <std::size_t dimension>
hpgem::LinearAlgebra::MiddleSizeVector
    SIPGErrorQualityMetric<dimension>::computeElementVector(
        std::size_t functionId, hpgem::Base::Element* element,
        hpgem::Integration::ElementIntegral<dimension>& integral) {
    using namespace hpgem;
    return integral.integrate(
        element, [&](Base::PhysicalElement<dimension>& pelement) {
            std::size_t n = pelement.getNumberOfBasisFunctions();
            LinearAlgebra::MiddleSizeVector vec(n);
            double fvalue =
                -1.0* functions_[functionId]->laplacian(pelement.getPointPhysical());
            for (std::size_t i = 0; i < n; ++i) {
                vec[i] = fvalue * pelement.basisFunction(i);
            }
            return vec;
        });
}

template <std::size_t dimension>
hpgem::LinearAlgebra::MiddleSizeMatrix
    SIPGErrorQualityMetric<dimension>::computeFaceMatrix(
        hpgem::Base::Face* face,
        hpgem::Integration::FaceIntegral<dimension>& integral) {
    using namespace hpgem;

    return integral.integrate(face, [&](Base::PhysicalFace<dimension>& pface) {
        std::size_t n = face->getNumberOfBasisFunctions();
        LinearAlgebra::MiddleSizeMatrix mat(n,n);
        for (std::size_t i = 0; i < n; ++i) {
            const LinearAlgebra::MiddleSizeVector& phin_i =
                pface.basisFunctionNormal(i);
            const LinearAlgebra::MiddleSizeVector& gradPhi_i =
                pface.basisFunctionDeriv(i);
            for (std::size_t j = 0; j < n; ++j) {
                const LinearAlgebra::MiddleSizeVector& phin_j =
                    pface.basisFunctionNormal(j);
                const LinearAlgebra::MiddleSizeVector& gradPhi_j =
                    pface.basisFunctionDeriv(j);

                // Averages over an internal face use a factor 0.5 while those
                // on the boundary face a factor 1.0.
                double factor = face->isInternal() ? 0.5 : 1.0;

                // Consistency & Symmetry terms
                mat(i, j) = -factor * (phin_i * gradPhi_j + phin_j * gradPhi_i);
                // Stability
                mat(i, j) +=
                    phin_i * phin_j * PENALTY_CONSTANT / face->getDiameter();
            }
        }
        return mat;
    });
}

template <std::size_t dimension>
hpgem::LinearAlgebra::MiddleSizeVector
    SIPGErrorQualityMetric<dimension>::computeFaceVector(
        std::size_t functionId, hpgem::Base::Face* face,
        hpgem::Integration::FaceIntegral<dimension>& integral) {
    using namespace hpgem;
    return integral.integrate(face, [&](Base::PhysicalFace<dimension>& pface) {
        std::size_t n = pface.getNumberOfBasisFunctions();
        LinearAlgebra::MiddleSizeVector vec(n);

        double fvalue = functions_[functionId]->value(pface.getPointPhysical());
        const LinearAlgebra::SmallVector<dimension>& normal =
            pface.getUnitNormalVector();
        for (std::size_t i = 0; i < n; ++i) {
            vec(i) += fvalue * (
                                   // Symmetry term
                                   0.5 * pface.basisFunctionDeriv(i) * normal +
                                   // Penalty term
                                   PENALTY_CONSTANT / face->getDiameter() *
                                       pface.basisFunction(i));
        }
        return vec;
    });
}

template <std::size_t dimension>
double SIPGErrorQualityMetric<dimension>::computeElementError(
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
typename SIPGErrorQualityMetric<dimension>::VALUE_TYPE
    SIPGErrorQualityMetric<dimension>::computeLocalError(
        std::size_t functionId, hpgem::Base::Element* element,
        const hpgem::Geometry::PointReference<dimension>& p) {
    using namespace hpgem;

    VALUE_TYPE error =
        functions_[functionId]->value(element->getReferenceToPhysicalMap()
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

#endif  // HPGEM_SIPGERRORQUALITYMETRIC_H

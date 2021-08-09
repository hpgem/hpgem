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
    const double PENALTY_CONSTANT = 10;

    using VALUE_TYPE = hpgem::LinearAlgebra::MiddleSizeVector::type;

   public:
    SIPGErrorQualityMetric(
        std::size_t order, std::string name,
        std::vector<std::shared_ptr<Func<dimension>>> functions)
        : order_(order), name_(name), functions_(functions){};

    void computeAndPlotMetric(
        hpgem::Base::MeshManipulator<dimension>& mesh,
        hpgem::Output::VTKSpecificTimeWriter<dimension>& plotter,
        const std::string& filePrefix) final;

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

    VALUE_TYPE computeSolutionValue(
        hpgem::Base::Element* element,
        const hpgem ::Geometry::PointReference<dimension>& p);

    std::size_t order_;
    std::string name_;
    std::vector<std::shared_ptr<Func<dimension>>> functions_;
};

template <std::size_t dimension>
void SIPGErrorQualityMetric<dimension>::computeAndPlotMetric(
    hpgem::Base::MeshManipulator<dimension>& mesh,
    hpgem::Output::VTKSpecificTimeWriter<dimension>& plotter,
    const std::string& filePrefix) {
    if (functions_.empty()) {
        return;
    }

    using namespace hpgem;

    logger(INFO, "Setting up %", name_);

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
        face->setFaceVector(
            LinearAlgebra::MiddleSizeVector(face->getNumberOfBasisFunctions()),
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
    KSPSetType(ksp, KSPGMRES);
    KSPSetOperators(ksp, stiffnessMatrix, stiffnessMatrix);
    KSPSetTolerances(ksp, 1e-12, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSetFromOptions(ksp);

    std::stringstream detailPlotFile;
    detailPlotFile << filePrefix << name_;
    Output::VTKSpecificTimeWriter<dimension> detailPlotter(detailPlotFile.str(),
                                                           &mesh, 0, order_);

    for (std::size_t fid = 0; fid < functions_.size(); ++fid) {
        // Compute local parts of the load vector

        logger(INFO, "Creating source for %-%", name_, fid);

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

        logger(INFO, "Solving for %-%", name_, fid);

        KSPSolve(ksp, loadVector, resultVector);
        resultVector.writeTimeIntegrationVector(SOLUTION_VECTOR_ID);

        KSPConvergedReason reason;
        KSPGetConvergedReason(ksp, &reason);
        if (reason > 0) {
            logger(INFO, "Converged with reason %", reason);
        } else {
            logger(WARN, "Unconverged with reason %", reason);
        }

        logger(INFO, "Outputting for %-%", name_, fid);

        for (Base::Element* element : mesh.getElementsList()) {
            elementErrorSquared[element] +=
                computeElementError(fid, element, integral);
        }

        // Plot the error of each individual function
        std::stringstream errorName;
        errorName << name_ << "-error-" << fid;
        detailPlotter.write(
            [&](Base::Element* element,
                const Geometry::PointReference<dimension>& p, std::size_t) {
                double error = std::real(computeSolutionValue(element, p));
                error -=
                    functions_[fid]->value(element->referenceToPhysical(p));
                return error;
            },
            errorName.str());

        // Solution
        std::stringstream solutionName;
        solutionName << name_ << "-solution-" << fid;
        detailPlotter.write(
            [&](Base::Element* element,
                const Geometry::PointReference<dimension>& p, std::size_t) {
                return std::real(computeSolutionValue(element, p));
            },
            solutionName.str());

        // Function
        std::stringstream funcName;
        funcName << name_ << "-func-" << fid;
        detailPlotter.write(
            [&](Base::Element* element,
                const Geometry::PointReference<dimension>& pref,
                std::size_t) -> double {
                Geometry::PointPhysical<dimension> phys =
                    element->referenceToPhysical(pref);
                return functions_[fid]->value(phys);
            },
            funcName.str());
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
        std::size_t functionId, hpgem::Base::Element* element,
        hpgem::Integration::ElementIntegral<dimension>& integral) {
    using namespace hpgem;
    return integral.integrate(
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
        hpgem::Base::Face* face,
        hpgem::Integration::FaceIntegral<dimension>& integral) {
    using namespace hpgem;

    return integral.integrate(face, [&](Base::PhysicalFace<dimension>& pface) {
        std::size_t n = face->getNumberOfBasisFunctions();
        LinearAlgebra::MiddleSizeMatrix mat(n, n);
        for (std::size_t i = 0; i < n; ++i) {
            // normal_i phi_i is computed at point p, the result is stored in
            // phiNormalI.
            const LinearAlgebra::SmallVector<dimension>& phiNormalI =
                pface.basisFunctionUnitNormal(i);
            // The gradient of basisfunction phi_i is computed at point p, the
            // result is stored in phiDerivI.
            const LinearAlgebra::SmallVector<dimension>& phiDerivI =
                pface.basisFunctionDeriv(i);

            for (std::size_t j = 0; j < n; ++j) {
                // normal_j phi_j is computed at point p, the result is stored
                // in phiNormalJ.
                const LinearAlgebra::SmallVector<dimension>& phiNormalJ =
                    pface.basisFunctionUnitNormal(j);
                // The gradient of basisfunction phi_j is computed at point p,
                // the result is stored in phiDerivJ.
                const LinearAlgebra::SmallVector<dimension>& phiDerivJ =
                    pface.basisFunctionDeriv(j);

                // Switch to the correct type of face, and compute the integrand
                // accordingly you could also compute the integrandVal by
                // directly using face->basisFunctionDeriv and
                // face->basisFunctionNormal in the following lines, but this
                // results in very long expressions Internal face:
                if (face->isInternal()) {
                    mat(j, i) =
                        -(phiNormalI * phiDerivJ + phiNormalJ * phiDerivI) / 2 +
                        PENALTY_CONSTANT / face->getDiameter() * phiNormalI *
                            phiNormalJ;
                } else {
                    mat(j, i) =
                        -(phiNormalI * phiDerivJ + phiNormalJ * phiDerivI) +
                        PENALTY_CONSTANT / face->getDiameter() * phiNormalI *
                            phiNormalJ;
                }
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
    if (face->isInternal()) {
        return LinearAlgebra::MiddleSizeVector(
            face->getNumberOfBasisFunctions());
    }

    return integral.integrate(face, [&](Base::PhysicalFace<dimension>& pface) {
        std::size_t n = pface.getNumberOfBasisFunctions();
        LinearAlgebra::MiddleSizeVector vec(n);

        double fvalue = functions_[functionId]->value(pface.getPointPhysical());
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
double SIPGErrorQualityMetric<dimension>::computeElementError(
    std::size_t functionId, hpgem::Base::Element* element,
    hpgem::Integration::ElementIntegral<dimension>& integral) {
    using namespace hpgem;
    double error = integral.integrate(
        element,
        [&](Base::PhysicalElement<dimension>& pelement) {
            VALUE_TYPE localError = computeSolutionValue(
                element, pelement.getPointReference());
            localError -=
                functions_[functionId]->value(pelement.getPointPhysical());
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
    SIPGErrorQualityMetric<dimension>::computeSolutionValue(
        hpgem::Base::Element* element,
        const hpgem::Geometry::PointReference<dimension>& p) {
    using namespace hpgem;

    VALUE_TYPE value = 0.0;
    LinearAlgebra::MiddleSizeVector& coefficients =
        element->getTimeIntegrationVector(SOLUTION_VECTOR_ID);
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        double phi_i = element->basisFunction(i, p);

        value += coefficients[i] * phi_i;
    }
    // In case complex valued
    return value;
}

#endif  // HPGEM_SIPGERRORQUALITYMETRIC_H

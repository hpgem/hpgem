/*
This file forms part of hpGEM. This package has been developed over a number of
years by various people at the University of Twente and a full list of
contributors can be found at http://hpgem.org/about-the-code/team

This code is distributed using BSD 3-Clause License. A copy of which can found
below.


Copyright (c) 2018, Univesity of Twenete
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "DGMaxDiscretization.h"

#include <complex>

#include "Base/HCurlConformingTransformation.h"
#include "Base/MeshManipulator.h"
#include "Integration/ElementIntegral.h"
#include "Integration/FaceIntegral.h"

#include "ElementInfos.h"

using namespace hpgem;

template <std::size_t DIM>
void DGMaxDiscretization<DIM>::initializeBasisFunctions(
    Base::MeshManipulator<DIM>& mesh, std::size_t order) {
    // We would like to configure the number of unknowns here, but this is
    // unfortunately not possible, as it is configured at the creation of
    // the mesh. The best we can do is check if it is configured correctly.
    std::size_t unknowns = mesh.getConfigData()->numberOfUnknowns_;
    if (includeProjector_) {
        logger.assert_always(unknowns == 2,
                             "DGMax+Projector expects 2 unknowns but got %",
                             unknowns);
        mesh.useNedelecDGBasisFunctions(order);
        mesh.useDefaultConformingBasisFunctions(order, 1);
    } else {
        logger.assert_always(unknowns == 1, "DGMax expects 1 unknown but got %",
                             unknowns);
        mesh.useNedelecDGBasisFunctions(order);
    }
    // TODO: This should probably also be exposed by using a constructor
    // parameter.
    // mesh.useAinsworthCoyleDGBasisFunctions();
}

template <std::size_t DIM>
void DGMaxDiscretization<DIM>::computeElementIntegrands(
    Base::MeshManipulator<DIM>& mesh, MassMatrixHandling massMatrixHandling,
    const InputFunction& sourceTerm, const InputFunction& initialCondition,
    const InputFunction& initialConditionDerivative) const {
    bool anyFuncPresent =
        sourceTerm || initialCondition || initialConditionDerivative;
    logger.assert_always(
        !anyFuncPresent || massMatrixHandling != ORTHOGONALIZE,
        "Mass matrix rescale with input functions is not implemented");
    LinearAlgebra::MiddleSizeMatrix massMatrix(1, 1), stiffnessMatrix(1, 1),
        projectorMatrix(0, 0);
    LinearAlgebra::MiddleSizeVector initialConditionVector(1),
        initialConditionDerivativeVector(1), elementVector(1);
    Integration::ElementIntegral<DIM> elIntegral;

    elIntegral.setTransformation(
        std::shared_ptr<Base::CoordinateTransformation<DIM>>(
            new Base::HCurlConformingTransformation<DIM>()));
    if (includeProjector_) {
        elIntegral.setTransformation(
            std::shared_ptr<Base::CoordinateTransformation<DIM>>(
                new Base::H1ConformingTransformation<DIM>()),
            1);
    }

    // Using the Global iterator as:
    //  - The projector needs to integration over all elements on which the
    //    conforming basis functions have support
    //  - Rescaling the stiffness face matrices (either here or later) requires
    //    the element mass matrices on both sides of the faces, including the
    //    ghost elements that share a face with a non ghost element.
    auto end = mesh.elementColEnd(Base::IteratorType::GLOBAL);
    for (typename Base::MeshManipulator<DIM>::ElementIterator it =
             mesh.elementColBegin(Base::IteratorType::GLOBAL);
         it != end; ++it) {
        Base::Element* element = *it;
        std::size_t numberOfBasisFunctions =
            element->getNumberOfBasisFunctions(0);

        // TODO: Are these resizes needed, as the content seems to be
        // overwritten by the integral.
        massMatrix.resize(numberOfBasisFunctions, numberOfBasisFunctions);
        massMatrix = elIntegral.integrate(
            element, [&](Base::PhysicalElement<DIM>& pelement) {
                LinearAlgebra::MiddleSizeMatrix res;
                elementMassMatrix(pelement, res);
                return res;
            });
        switch (massMatrixHandling) {
            case DGMaxDiscretizationBase::NORMAL:
                break;
            case DGMaxDiscretizationBase::INVERT:
                massMatrix = massMatrix.inverse();
                break;
            case DGMaxDiscretizationBase::ORTHOGONALIZE:
                massMatrix.cholesky();
                break;
            default:
                logger.assert_always(false,
                                     "Not implemented mass matrix handling");
        }
        element->setElementMatrix(massMatrix, MASS_MATRIX_ID);

        stiffnessMatrix.resize(numberOfBasisFunctions, numberOfBasisFunctions);
        stiffnessMatrix = elIntegral.integrate(
            element, [&](Base::PhysicalElement<DIM>& pelement) {
                LinearAlgebra::MiddleSizeMatrix res;
                elementStiffnessMatrix(pelement, res);
                return res;
            });
        if (massMatrixHandling == DGMaxDiscretizationBase::ORTHOGONALIZE) {
            // Compute L^{-1} S L^{-H}, where S is the stiffness matrix and
            // LL^H is the mass matrix.
            LinearAlgebra::MiddleSizeMatrix original = stiffnessMatrix;
            massMatrix.solveLowerTriangular(stiffnessMatrix,
                                            LinearAlgebra::Side::OP_LEFT,
                                            LinearAlgebra::Transpose::NOT);
            massMatrix.solveLowerTriangular(
                stiffnessMatrix, LinearAlgebra::Side::OP_RIGHT,
                LinearAlgebra::Transpose::HERMITIAN_TRANSPOSE);
            // Due to rounding errors the matrix might be slightly non
            // Hermitian, fix this by replacing S by 0.5(S + S^H).
            for (std::size_t i = 0; i < stiffnessMatrix.getNumberOfRows();
                 ++i) {
                stiffnessMatrix(i, i) = std::real(stiffnessMatrix(i, i));
                for (std::size_t j = i;
                     j < stiffnessMatrix.getNumberOfColumns(); ++j) {
                    std::complex<double> upper =
                        0.5 * (stiffnessMatrix(i, j) +
                               std::conj(stiffnessMatrix(j, i)));
                    stiffnessMatrix(i, j) = upper;
                    stiffnessMatrix(j, i) = std::conj(upper);
                }
            }
        }
        element->setElementMatrix(stiffnessMatrix, STIFFNESS_MATRIX_ID);

        if (includeProjector_) {
            std::size_t numberOfProjectorBasisFunctions =
                element->getNumberOfBasisFunctions(1);
            projectorMatrix.resize(numberOfProjectorBasisFunctions,
                                   numberOfBasisFunctions);
            projectorMatrix = elIntegral.integrate(
                element, [&](Base::PhysicalElement<DIM>& pelement) {
                    LinearAlgebra::MiddleSizeMatrix res;
                    elementProjectorMatrix(pelement, res);
                    return res;
                });

            if (massMatrixHandling == DGMaxDiscretizationBase::ORTHOGONALIZE) {
                // Compute B L^{-H}, where B is the projector matrix and L is
                // the Cholesky factor of the mass matrix.
                massMatrix.solveLowerTriangular(
                    projectorMatrix, LinearAlgebra::Side::OP_RIGHT,
                    LinearAlgebra::Transpose::HERMITIAN_TRANSPOSE);
            }

            element->setElementMatrix(projectorMatrix, PROJECTOR_MATRIX_ID);
        }

        // Note, resizes for vectors needed in case of empty functions
        initialConditionVector.resize(numberOfBasisFunctions);
        if (initialCondition) {
            initialConditionVector = elIntegral.integrate(
                element, [&](Base::PhysicalElement<DIM>& element) {
                    LinearAlgebra::MiddleSizeVector res;
                    elementInnerProduct(element, initialCondition,
                                        res);  // Initial conditions
                    return res;
                });
        }
        element->setElementVector(initialConditionVector,
                                  INITIAL_CONDITION_VECTOR_ID);

        initialConditionDerivativeVector.resize(numberOfBasisFunctions);
        if (initialConditionDerivative) {
            initialConditionDerivativeVector = elIntegral.integrate(
                element, [&](Base::PhysicalElement<DIM>& pelement) {
                    LinearAlgebra::MiddleSizeVector res;
                    elementInnerProduct(pelement, initialConditionDerivative,
                                        res);  // Initial conditions derivative
                    return res;
                });
        }
        element->setElementVector(initialConditionDerivativeVector,
                                  INITIAL_CONDITION_DERIVATIVE_VECTOR_ID);

        elementVector.resize(numberOfBasisFunctions);
        if (sourceTerm) {
            elementVector = elIntegral.integrate(
                element, [&](Base::PhysicalElement<DIM>& pelement) {
                    LinearAlgebra::MiddleSizeVector res;
                    elementInnerProduct(pelement, sourceTerm,
                                        res);  // Source  term
                    return res;
                });
        }

        element->setElementVector(elementVector, SOURCE_TERM_VECTOR_ID);
    }
}

template <std::size_t DIM>
void DGMaxDiscretization<DIM>::computeFaceIntegrals(
    Base::MeshManipulator<DIM>& mesh, MassMatrixHandling massMatrixHandling,
    const DGMaxDiscretization<DIM>::FaceInputFunction& boundaryCondition,
    double stab) const {
    logger.assert_always(
        massMatrixHandling != ORTHOGONALIZE || !boundaryCondition,
        "Rescale not implemented in combination with boundary conditions");

    LinearAlgebra::MiddleSizeMatrix stiffnessFaceMatrix(0, 0);
    // Mass matrix for the face, already Cholesky factored.
    LinearAlgebra::MiddleSizeMatrix massMatrix(0, 0);
    LinearAlgebra::MiddleSizeVector boundaryFaceVector(1);
    Integration::FaceIntegral<DIM> faIntegral;

    faIntegral.setTransformation(
        std::shared_ptr<Base::CoordinateTransformation<DIM>>(
            new Base::HCurlConformingTransformation<DIM>()));
    auto end = mesh.faceColEnd();
    for (typename Base::MeshManipulator<DIM>::FaceIterator it =
             mesh.faceColBegin();
         it != end; ++it) {
        Base::Face* face = *it;

        // Resize all the matrices and vectors;
        std::size_t numberOfBasisFunctions =
            face->getPtrElementLeft()->getNumberOfBasisFunctions(0);
        if (face->isInternal()) {
            numberOfBasisFunctions +=
                face->getPtrElementRight()->getNumberOfBasisFunctions(0);
        }

        stiffnessFaceMatrix.resize(numberOfBasisFunctions,
                                   numberOfBasisFunctions);
        boundaryFaceVector.resize(numberOfBasisFunctions);

        // Compute the actual face  integrals.
        stiffnessFaceMatrix =
            faIntegral.integrate(face, [&](Base::PhysicalFace<DIM>& pface) {
                LinearAlgebra::MiddleSizeMatrix res;
                faceMatrix(pface, res);
                LinearAlgebra::MiddleSizeMatrix temp;
                facePenaltyMatrix(pface, temp, stab);
                res += temp;
                return res;
            });
        if (massMatrixHandling == DGMaxDiscretizationBase::ORTHOGONALIZE) {
            massMatrix.resize(numberOfBasisFunctions, numberOfBasisFunctions);
            massMatrix *= 0.0;  // Clear the contents
            const LinearAlgebra::MiddleSizeMatrix& leftMassMat =
                face->getPtrElementLeft()->getElementMatrix(MASS_MATRIX_ID);
            std::size_t leftRows = leftMassMat.getNumberOfRows();
            // Copy lower triagonal part
            for (std::size_t i = 0; i < leftRows; ++i) {
                for (std::size_t j = i; j < leftRows; ++j) {
                    massMatrix(j, i) = leftMassMat(j, i);
                }
            }
            if (face->isInternal()) {
                const LinearAlgebra::MiddleSizeMatrix& rightMassMat =
                    face->getPtrElementRight()->getElementMatrix(
                        MASS_MATRIX_ID);
                std::size_t rightRows = rightMassMat.getNumberOfRows();
                // Copy lower triagonal part, now offset by leftRows.
                for (std::size_t i = 0; i < rightRows; ++i) {
                    for (std::size_t j = i; j < rightRows; ++j) {
                        massMatrix(leftRows + j, leftRows + i) =
                            rightMassMat(j, i);
                    }
                }
            }
            LinearAlgebra::MiddleSizeMatrix original = stiffnessFaceMatrix;
            // Rescaling
            massMatrix.solveLowerTriangular(stiffnessFaceMatrix,
                                            hpgem::LinearAlgebra::Side::OP_LEFT,
                                            LinearAlgebra::Transpose::NOT);
            massMatrix.solveLowerTriangular(
                stiffnessFaceMatrix, hpgem::LinearAlgebra::Side::OP_RIGHT,
                LinearAlgebra::Transpose::HERMITIAN_TRANSPOSE);
        }
        face->setFaceMatrix(stiffnessFaceMatrix, FACE_MATRIX_ID);

        if (boundaryCondition) {
            // NOTE: No check on rescaling, as the face vector is zero when
            // there is no boundaryCondition.
            boundaryFaceVector =
                faIntegral.integrate(face, [&](Base::PhysicalFace<DIM>& face) {
                    LinearAlgebra::MiddleSizeVector res;
                    faceVector(face, boundaryCondition, res, stab);
                    return res;
                });
        }

        face->setFaceVector(boundaryFaceVector, FACE_VECTOR_ID);
    }
}

template <std::size_t DIM>
void DGMaxDiscretization<DIM>::elementMassMatrix(
    Base::PhysicalElement<DIM>& el,
    LinearAlgebra::MiddleSizeMatrix& ret) const {
    const Base::Element* element = el.getElement();
    const std::size_t numberOfBasisFunctions =
        element->getNumberOfBasisFunctions(0);
    ret.resize(numberOfBasisFunctions, numberOfBasisFunctions);
    LinearAlgebra::SmallVector<DIM> phi_i, phi_j;
    double epsilon =
        static_cast<ElementInfos*>(element->getUserData())->epsilon_;
    for (std::size_t i = 0; i < numberOfBasisFunctions; ++i) {
        el.basisFunction(i, phi_i, 0);
        for (std::size_t j = 0; j < numberOfBasisFunctions; ++j) {
            el.basisFunction(j, phi_j, 0);
            ret(i, j) = phi_i * phi_j * epsilon;
            ret(j, i) = ret(i, j);
        }
    }
}

template <std::size_t DIM>
void DGMaxDiscretization<DIM>::elementStiffnessMatrix(
    Base::PhysicalElement<DIM>& el,
    LinearAlgebra::MiddleSizeMatrix& ret) const {
    const Base::Element* element = el.getElement();
    const std::size_t numberOfBasisFunctions =
        element->getNumberOfBasisFunctions(0);
    ret.resize(numberOfBasisFunctions, numberOfBasisFunctions);
    LinearAlgebra::SmallVector<DIM> phi_i, phi_j;
    for (std::size_t i = 0; i < numberOfBasisFunctions; ++i) {
        phi_i = el.basisFunctionCurl(i, 0);
        for (std::size_t j = i; j < numberOfBasisFunctions; ++j) {
            phi_j = el.basisFunctionCurl(j, 0);
            ret(i, j) = phi_i * phi_j;
            ret(j, i) = ret(i, j);
        }
    }
}

template <std::size_t DIM>
void DGMaxDiscretization<DIM>::elementProjectorMatrix(
    Base::PhysicalElement<DIM>& el,
    LinearAlgebra::MiddleSizeMatrix& ret) const {
    const Base::Element* element = el.getElement();
    double epsilon =
        static_cast<ElementInfos*>(element->getUserData())->epsilon_;
    const std::size_t dofU = element->getNumberOfBasisFunctions(0);
    const std::size_t dofP = element->getNumberOfBasisFunctions(1);
    ret.resize(dofP, dofU);
    LinearAlgebra::SmallVector<DIM> phiU;
    for (std::size_t i = 0; i < dofU; ++i) {
        el.basisFunction(i, phiU, 0);
        for (std::size_t j = 0; j < dofP; ++j) {
            ret(j, i) = epsilon * phiU * el.basisFunctionDeriv(j, 1);
        }
    }
}

template <std::size_t DIM>
void DGMaxDiscretization<DIM>::elementInnerProduct(
    Base::PhysicalElement<DIM>& el, const InputFunction& function,
    LinearAlgebra::MiddleSizeVector& ret) const {
    const Base::Element* element = el.getElement();
    const Geometry::PointReference<DIM>& p = el.getPointReference();
    const std::size_t numberOfBasisFunctions =
        element->getNumberOfBasisFunctions(0);

    ret.resize(numberOfBasisFunctions);
    const PointPhysicalT pPhys = element->referenceToPhysical(p);
    LinearAlgebra::SmallVector<DIM> val, phi;
    function(pPhys, val);
    for (std::size_t i = 0; i < numberOfBasisFunctions; ++i) {
        el.basisFunction(i, phi, 0);
        ret[i] = phi * val;
    }
}

template <std::size_t DIM>
void DGMaxDiscretization<DIM>::faceMatrix(
    Base::PhysicalFace<DIM>& fa, LinearAlgebra::MiddleSizeMatrix& ret) const {
    const Base::Face* face = fa.getFace();

    std::size_t M = face->getPtrElementLeft()->getNrOfBasisFunctions(0);
    const bool internalFace = face->isInternal();
    if (internalFace) {
        M += face->getPtrElementRight()->getNrOfBasisFunctions(0);
    }
    ret.resize(M, M);
    LinearAlgebra::SmallVector<DIM> phi_i_normal, phi_j_normal, phi_i_curl,
        phi_j_curl;
    for (std::size_t i = 0; i < M; ++i) {
        phi_i_curl = fa.basisFunctionCurl(i, 0);
        fa.basisFunctionUnitNormalCross(i, phi_i_normal, 0);

        for (std::size_t j = i; j < M; ++j) {
            phi_j_curl = fa.basisFunctionCurl(j, 0);
            fa.basisFunctionUnitNormalCross(j, phi_j_normal, 0);

            ret(i, j) = -(internalFace ? 0.5 : 1.) *
                        (phi_i_normal * phi_j_curl + phi_j_normal * phi_i_curl);
            ret(j, i) = ret(i, j);
        }
    }
}

template <std::size_t DIM>
void DGMaxDiscretization<DIM>::facePenaltyMatrix(
    Base::PhysicalFace<DIM>& fa, LinearAlgebra::MiddleSizeMatrix& ret,
    double stab) const {
    const Base::Face* face = fa.getFace();
    double diameter = face->getDiameter();

    std::size_t M = face->getPtrElementLeft()->getNrOfBasisFunctions(0);
    if (face->isInternal()) {
        M += face->getPtrElementRight()->getNrOfBasisFunctions(0);
    }

    ret.resize(M, M);
    LinearAlgebra::SmallVector<DIM> phi_i, phi_j;
    for (std::size_t i = 0; i < M; ++i) {
        fa.basisFunctionUnitNormalCross(i, phi_i, 0);
        for (std::size_t j = i; j < M; ++j) {
            fa.basisFunctionUnitNormalCross(j, phi_j, 0);

            ret(i, j) = stab / diameter * (phi_i * phi_j);
            ret(j, i) = ret(i, j);
        }
    }
}

template <std::size_t DIM>
void DGMaxDiscretization<DIM>::faceVector(
    Base::PhysicalFace<DIM>& fa, const FaceInputFunction& boundaryCondition,
    LinearAlgebra::MiddleSizeVector& ret, double stab) const {
    const Base::Face* face = fa.getFace();
    LinearAlgebra::SmallVector<DIM> normal = fa.getNormalVector();
    normal /= Base::L2Norm(normal);
    const Geometry::PointReference<DIM - 1>& p = fa.getPointReference();

    if (face->isInternal()) {
        std::size_t M = face->getPtrElementLeft()->getNrOfBasisFunctions(0) +
                        face->getPtrElementRight()->getNrOfBasisFunctions(0);
        ret.resize(M);
        for (std::size_t i = 0; i < M; ++i) ret(i) = 0;
    } else {
        double diameter = face->getDiameter();
        const Base::Element* left = face->getPtrElementLeft();
        const Geometry::PointReference<DIM>& PLeft =
            face->mapRefFaceToRefElemL(p);

        PointPhysicalT PPhys;
        PPhys = left->referenceToPhysical(PLeft);
        LinearAlgebra::SmallVector<DIM> val, phi, phi_curl, boundaryValues;

        boundaryCondition(PPhys, fa,
                          boundaryValues);  // assumes the initial conditions
                                            // and the boundary conditions match

        val = boundaryValues;
        std::size_t n = face->getPtrElementLeft()->getNrOfBasisFunctions(0);
        ret.resize(n);

        for (std::size_t i = 0; i < n; ++i) {
            fa.basisFunctionUnitNormalCross(i, phi, 0);

            phi_curl = fa.basisFunctionCurl(i, 0);

            ret(i) = -(phi_curl * val) + stab / diameter * (phi * val);
        }
    }
}

template <std::size_t DIM>
LinearAlgebra::SmallVector<DIM> DGMaxDiscretization<DIM>::computeField(
    const Base::Element* element, const Geometry::PointReference<DIM>& p,
    const LinearAlgebra::MiddleSizeVector& coefficients) const {
    logger.log(Log::WARN, "Only computing the real part of the field.");
    Base::PhysicalElement<DIM> physicalElement;
    physicalElement.setElement(element);
    std::shared_ptr<Base::CoordinateTransformation<DIM>> transform{
        new Base::HCurlConformingTransformation<DIM>()};
    physicalElement.setTransformation(transform);
    physicalElement.setPointReference(p);

    LinearAlgebra::SmallVector<DIM> result, phi;
    for (std::size_t i = 0; i < element->getNumberOfBasisFunctions(0); ++i) {
        physicalElement.basisFunction(i, phi, 0);
        result += std::real(coefficients[i]) * phi;
    }
    return result;
}

template <std::size_t DIM>
LinearAlgebra::SmallVector<DIM> DGMaxDiscretization<DIM>::computeCurlField(
    const Base::Element* element, const Geometry::PointReference<DIM>& p,
    const LinearAlgebra::MiddleSizeVector& coefficients) const {
    logger.log(Log::WARN, "Only computing the real part of the field.");
    Base::PhysicalElement<DIM> physicalElement;
    physicalElement.setElement(element);
    std::shared_ptr<Base::CoordinateTransformation<DIM>> transform{
        new Base::HCurlConformingTransformation<DIM>()};
    physicalElement.setTransformation(transform);
    physicalElement.setPointReference(p);

    LinearAlgebra::SmallVector<DIM> result;
    for (std::size_t i = 0; i < element->getNumberOfBasisFunctions(0); ++i) {
        result += std::real(coefficients[i]) *
                  physicalElement.basisFunctionCurl(i, 0);
    }
    return result;
}

// TODO: The code saves snapshots in the timeIntegrationVector, this is not
// particularly nice It might be better to pass the global vector here and
// distribute it ourselves.
template <std::size_t DIM>
std::map<typename DGMaxDiscretization<DIM>::NormType, double>
    DGMaxDiscretization<DIM>::computeError(
        Base::MeshManipulator<DIM>& mesh, std::size_t timeVector,
        DGMaxDiscretization<DIM>::InputFunction electricField,
        DGMaxDiscretization<DIM>::InputFunction electricFieldCurl,
        std::set<DGMaxDiscretization<DIM>::NormType> norms) const {
    // Note these are actually the squared norms
    double l2Norm = 0;
    double hCurlNorm = 0;
    // Setup the element integration.
    Integration::ElementIntegral<DIM> elIntegral;
    elIntegral.setTransformation(
        std::shared_ptr<Base::CoordinateTransformation<DIM>>(
            new Base::HCurlConformingTransformation<DIM>()));

    bool l2Wanted = norms.find(NormType::L2) != norms.end();
    bool hcurlWanted = norms.find(NormType::HCurl) != norms.end();
    bool dgWanted = norms.find(NormType::DG) != norms.end();

    auto end = mesh.elementColEnd();
    for (typename Base::MeshManipulator<DIM>::ElementIterator it =
             mesh.elementColBegin();
         it != end; ++it) {
        LinearAlgebra::SmallVector<2> errors =
            elIntegral.integrate((*it), [&](Base::PhysicalElement<DIM>& el) {
                return elementErrorIntegrand(el, dgWanted || hcurlWanted,
                                             timeVector, electricField,
                                             electricFieldCurl);
            });
        l2Norm += errors[0];
        hCurlNorm += errors[1];
    }
    hCurlNorm += l2Norm;

    double dgNorm = hCurlNorm;

    if (dgWanted) {
        Integration::FaceIntegral<DIM> faIntegral;
        faIntegral.setTransformation(
            std::shared_ptr<Base::CoordinateTransformation<DIM>>(
                new Base::HCurlConformingTransformation<DIM>()));
        auto end = mesh.faceColEnd();
        for (typename Base::MeshManipulator<DIM>::FaceIterator it =
                 mesh.faceColBegin();
             it != end; ++it) {
            dgNorm +=
                faIntegral.integrate(*it, [&](Base::PhysicalFace<DIM>& face) {
                    return faceErrorIntegrand(face, timeVector, electricField);
                });
        }
    }

    std::map<NormType, double> result;
    if (l2Wanted) result[L2] = sqrt(l2Norm);
    if (hcurlWanted) result[HCurl] = sqrt(hCurlNorm);
    if (dgWanted) result[DG] = sqrt(dgNorm);
    return result;
}

template <std::size_t DIM>
LinearAlgebra::SmallVector<2> DGMaxDiscretization<DIM>::elementErrorIntegrand(
    Base::PhysicalElement<DIM>& el, bool computeCurl, std::size_t timeVector,
    DGMaxDiscretization<DIM>::InputFunction exactValues,
    DGMaxDiscretization<DIM>::InputFunction curlValues) const {
    const Base::Element* element = el.getElement();
    const Geometry::PointReference<DIM>& p = el.getPointReference();
    PointPhysicalT pPhys;
    pPhys = element->referenceToPhysical(p);

    LinearAlgebra::SmallVector<DIM> phi, phiCurl, error, errorCurl;

    exactValues(pPhys, error);
    if (computeCurl) {
        curlValues(pPhys, errorCurl);
    }
    LinearAlgebra::MiddleSizeVector data;
    data = element->getTimeIntegrationVector(timeVector);
    for (std::size_t i = 0; i < element->getNrOfBasisFunctions(0); ++i) {
        el.basisFunction(i, phi, 0);
        error -= (std::real(data[i]) * phi);
        if (computeCurl) {
            phiCurl = el.basisFunctionCurl(i, 0);
            errorCurl -= (std::real(data[i]) * phiCurl);
        }
    }
    double l2Error = Base::L2Norm(error);
    l2Error *= l2Error;
    double curlError = Base::L2Norm(errorCurl);
    curlError *= curlError;
    LinearAlgebra::SmallVector<2> errors;
    errors[0] = l2Error;
    errors[1] = curlError;
    return errors;
}

template <std::size_t DIM>
double DGMaxDiscretization<DIM>::faceErrorIntegrand(
    Base::PhysicalFace<DIM>& fa, std::size_t timeVector,
    DGMaxDiscretization<DIM>::InputFunction exactSolution) const {
    // The face error part of the DG norm is given by
    // || h^0.5 [[u - u_h]]_T ||^2. Where h is the diameter of the face, u and
    // u_h are the exact and computed solutions. Further more the [[ . ]]_T is
    // the tangential jump, [[ a ]]_T = a_L x n_L + a_R x n_R for internal faces
    // and a x n for boundary faces (n is the normal).
    const Base::Face* face = fa.getFace();
    LinearAlgebra::SmallVector<DIM> normal = fa.getNormalVector();
    normal /= Base::L2Norm(normal);
    const Geometry::PointReference<2>& p = fa.getPointReference();

    Base::Element* element =
        const_cast<Base::Element*>(face->getPtrElementLeft());
    PointPhysicalT PPhys;
    const Geometry::PointReference<DIM>& pElement =
        face->mapRefFaceToRefElemL(p);

    PPhys = element->referenceToPhysical(pElement);
    LinearAlgebra::SmallVector<DIM> error, phiNormal, solutionValues;

    // Compute u_L x n_L
    exactSolution(PPhys, solutionValues);
    error = normal.crossProduct(solutionValues);
    std::size_t n = face->getPtrElementLeft()->getNrOfBasisFunctions(0);
    LinearAlgebra::MiddleSizeVector solutionCoefficients =
        element->getTimeIntegrationVector(
            timeVector);  // Issue regarding parallelisation is in this
                          // line....it goes out of bound for memory

    for (int i = 0; i < n; ++i) {
        // subtract the solution part, u_hL x n_L
        fa.basisFunctionUnitNormalCross(i, phiNormal, 0);
        // TODO: What about the complex part of the solution.
        error -= std::real(solutionCoefficients[i]) * phiNormal;
    }
    if (face->isInternal()) {
        // Note we reuse most of the vectors from the computation on the left
        // side of the face.
        LinearAlgebra::SmallVector<DIM> otherSideError;
        element = const_cast<Base::Element*>(face->getPtrElementRight());
        const Geometry::PointReference<DIM>& pElement =
            face->mapRefFaceToRefElemR(p);
        PPhys = element->referenceToPhysical(pElement);
        // Compute u_R x n_L
        exactSolution(PPhys, solutionValues);
        otherSideError = normal.crossProduct(solutionValues);
        error -= otherSideError;  // Note subtraction as n_L = - n_R.

        solutionCoefficients = element->getTimeIntegrationVector(timeVector);
        std::size_t M = face->getPtrElementLeft()->getNrOfBasisFunctions(0) +
                        face->getPtrElementRight()->getNrOfBasisFunctions(0);

        for (std::size_t i = n; i < M; ++i) {
            // Subtract u_hR x n_R, note that the normal used internally is n_R
            // and not n_L as we have with the solution on the right, so again
            // subtraction.
            fa.basisFunctionUnitNormalCross(i, phiNormal, 0);
            error -= (std::real(solutionCoefficients[i - n]) * phiNormal);
        }
    }
    double result = Base::L2Norm(error);
    result *= result;
    // To remove double contribution of flux computed on the boudary faces by
    // different processors
    if (face->getFaceType() == Geometry::FaceType::SUBDOMAIN_BOUNDARY ||
        face->getFaceType() == Geometry::FaceType::PERIODIC_SUBDOMAIN_BC) {
        result /= 2;
    }
    return result;
}

template class DGMaxDiscretization<2>;
template class DGMaxDiscretization<3>;
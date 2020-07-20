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

#include "DivDGMaxDiscretization.h"

#include "Base/MeshManipulator.h"
#include "Base/HCurlConformingTransformation.h"
#include "Integration/ElementIntegral.h"
#include "Integration/FaceIntegral.h"

#include "ElementInfos.h"

using namespace hpgem;

// Utility function
struct FaceDoFInfo {
    std::size_t leftUDoFs;
    std::size_t rightUDoFs;
    std::size_t leftPDoFs;
    std::size_t rightPDoFs;
    bool internal;

    std::size_t totalUDoFs() { return leftUDoFs + rightUDoFs; }

    std::size_t totalPDoFs() { return leftPDoFs + rightPDoFs; }

    std::size_t totalDoFs() { return totalUDoFs() + totalPDoFs(); }
};

FaceDoFInfo getFaceDoFInfo(const Base::Face& face) {
    FaceDoFInfo info;
    info.leftUDoFs = face.getPtrElementLeft()->getNumberOfBasisFunctions(0);
    info.leftPDoFs = face.getPtrElementLeft()->getNumberOfBasisFunctions(1);
    info.internal = face.isInternal();
    if (info.internal) {
        info.rightUDoFs =
            face.getPtrElementRight()->getNumberOfBasisFunctions(0);
        info.rightPDoFs =
            face.getPtrElementRight()->getNumberOfBasisFunctions(1);
    } else {
        info.rightUDoFs = 0;
        info.rightPDoFs = 0;
    }
    return info;
}

template <std::size_t DIM>
void DivDGMaxDiscretization<DIM>::initializeBasisFunctions(
    Base::MeshManipulator<DIM>& mesh, std::size_t order) {
    // We would like to configure the number of unknowns here, but this is
    // unfortunately not possible, as it is configured at the creation of
    // the mesh. The best we can do is check if it is configured correctly.
    std::size_t unknowns = mesh.getConfigData()->numberOfUnknowns_;
    logger.assert_always(unknowns == 2, "DivDGMax expects 2 unknowns but got %",
                         unknowns);
    // TODO: This needs the additional unknown id.
    mesh.useNedelecDGBasisFunctions(order);
    mesh.useDefaultDGBasisFunctions(order, 1);
}

template <std::size_t DIM>
void DivDGMaxDiscretization<DIM>::computeElementIntegrands(
    Base::MeshManipulator<DIM>& mesh, bool invertMassMatrix,
    const DivDGMaxDiscretization<DIM>::InputFunction& sourceTerm,
    const DivDGMaxDiscretization<DIM>::InputFunction& initialCondition,
    const DivDGMaxDiscretization<DIM>::InputFunction&
        initialConditionDerivative) const {
    // TODO: Add initial condition integration.
    LinearAlgebra::MiddleSizeMatrix massMatrix(2, 2), stiffnessMatrix(2, 2);

    LinearAlgebra::MiddleSizeVector vector1(2), vector2(2), sourceVector(2);

    std::size_t totalDoFs = 0, totalUDoFs = 0, totalPDoFs = 0;
    Integration::ElementIntegral<DIM> elementIntegral;

    elementIntegral.setTransformation(
        std::shared_ptr<Base::CoordinateTransformation<DIM>>(
            new Base::HCurlConformingTransformation<DIM>()),
        0);
    elementIntegral.setTransformation(
        std::shared_ptr<Base::CoordinateTransformation<DIM>>(
            new Base::H1ConformingTransformation<DIM>()),
        1);

    auto end = mesh.elementColEnd();
    for (typename Base::MeshManipulator<DIM>::ElementIterator it =
             mesh.elementColBegin();
         it != end; ++it) {
        totalUDoFs = (*it)->getNumberOfBasisFunctions(0);
        totalPDoFs = (*it)->getNumberOfBasisFunctions(1);
        totalDoFs = totalUDoFs + totalPDoFs;

        // mass matrix
        massMatrix.resize(totalDoFs, totalDoFs);
        massMatrix = elementIntegral.integrate(
            (*it), [&](Base::PhysicalElement<DIM>& element) {
                LinearAlgebra::MiddleSizeMatrix result;
                elementMassMatrix(element, result);
                return result;
            });
        if (invertMassMatrix) {
            massMatrix = massMatrix.inverse();
        }
        (*it)->setElementMatrix(massMatrix, ELEMENT_MASS_MATRIX_ID);

        stiffnessMatrix.resize(totalDoFs, totalDoFs);
        stiffnessMatrix = elementIntegral.integrate(
            (*it), [&](Base::PhysicalElement<DIM>& element) {
                LinearAlgebra::MiddleSizeMatrix result, temp;
                elementStiffnessMatrix(element, result);
                elementScalarVectorCoupling(element, temp);
                result += temp;
                return result;
            });

        (*it)->setElementMatrix(stiffnessMatrix, ELEMENT_STIFFNESS_MATRIX_ID);

        if (initialCondition) {
            logger.assert_debug(
                false, "Initial condition integration not ported yet.");
        }

        if (initialConditionDerivative) {
            logger.assert_debug(
                false,
                "Initial condition derivative integration not ported yet.");
        }

        if (sourceTerm) {
            sourceVector.resize(totalDoFs);
            sourceVector = elementIntegral.integrate(
                (*it), [&](Base::PhysicalElement<DIM>& element) {
                    LinearAlgebra::MiddleSizeVector result;
                    elementSourceVector(element, sourceTerm, result);
                    return result;
                });
            (*it)->setElementVector(sourceVector, ELEMENT_SOURCE_VECTOR_ID);
        }
    }
}

template <std::size_t DIM>
void DivDGMaxDiscretization<DIM>::computeFaceIntegrals(
    Base::MeshManipulator<DIM>& mesh,
    DivDGMaxDiscretization<DIM>::FaceInputFunction boundaryCondition,
    Stab stab) const {
    LinearAlgebra::MiddleSizeMatrix faceMatrix(2, 2);
    LinearAlgebra::MiddleSizeVector faceVector(2);
    Integration::FaceIntegral<DIM> faceIntegral;

    faceIntegral.setTransformation(
        std::shared_ptr<Base::CoordinateTransformation<DIM>>(
            new Base::HCurlConformingTransformation<DIM>()),
        0);
    faceIntegral.setTransformation(
        std::shared_ptr<Base::CoordinateTransformation<DIM>>(
            new Base::H1ConformingTransformation<DIM>()),
        1);

    auto end = mesh.faceColEnd();
    for (typename Base::MeshManipulator<DIM>::FaceIterator it =
             mesh.faceColBegin();
         it != end; ++it) {

        std::size_t totalDoFs = 0;

        totalDoFs = (*it)->getPtrElementLeft()->getNumberOfBasisFunctions(0) +
                    (*it)->getPtrElementLeft()->getNumberOfBasisFunctions(1);
        if ((*it)->isInternal()) {
            totalDoFs +=
                (*it)->getPtrElementRight()->getNumberOfBasisFunctions(0) +
                (*it)->getPtrElementRight()->getNumberOfBasisFunctions(1);
        }

        faceMatrix.resize(totalDoFs, totalDoFs);
        faceVector.resize(totalDoFs);

        faceMatrix =
            faceIntegral.integrate((*it), [&](Base::PhysicalFace<DIM>& face) {
                LinearAlgebra::MiddleSizeMatrix result, temp;
                faceStiffnessMatrix1(face, result);

                if (stab.fluxType1 == FluxType::IP) {
                    faceStiffnessMatrix2(face, temp, stab.stab1);
                    result += temp;
                    temp *= 0;  // Reset the variable;
                }
                if (stab.fluxType2 == FluxType::IP) {
                    faceStiffnessMatrix3(face, temp, stab.stab2);
                    result += temp;
                    temp *= 0;  // Reset the variable;
                }
                if (stab.fluxType3 == FluxType::IP) {
                    faceStiffnessScalarMatrix4(face, temp, stab.stab3);
                    // Note the matrix contribution is -C.
                    result -= temp;
                    temp *= 0;  // Reset the variable;
                }

                faceScalarVectorCoupling(face, temp);
                result += temp;

                // Reset no longer needed.
                return result;
            });

        if (stab.hasFlux(FluxType::BREZZI)) {
            faceMatrix += brezziFluxBilinearTerm(it, stab);
        }
        (*it)->setFaceMatrix(faceMatrix, FACE_STIFFNESS_MATRIX_ID);

        if (boundaryCondition) {
            faceVector = faceIntegral.integrate(
                (*it), [&](Base::PhysicalFace<DIM>& face) {
                    LinearAlgebra::MiddleSizeVector result;
                    faceBoundaryVector(face, boundaryCondition, result, stab);
                    return result;
                });
            if (stab.fluxType1 == FluxType::BREZZI) {
                faceVector +=
                    brezziFluxBoundaryVector(it, boundaryCondition, stab);
            }

            (*it)->setFaceVector(faceVector, FACE_BOUNDARY_VECTOR_ID);
        }
    }
}

template <std::size_t DIM>
LinearAlgebra::SmallVector<DIM> DivDGMaxDiscretization<DIM>::computeField(
    const Base::Element* element, const Geometry::PointReference<DIM>& point,
    const LinearAlgebra::MiddleSizeVector& coefficients) const {
    logger.log(Log::WARN, "Only computing the real part of the field.");
    LinearAlgebra::SmallVector<DIM> result, phiU;
    std::size_t nb0 = element->getNumberOfBasisFunctions(0),
                nb1 = element->getNumberOfBasisFunctions(1);

    // Setup the physical element with the correct transformations
    Base::PhysicalElement<DIM> physicalElement;
    std::shared_ptr<Base::CoordinateTransformation<DIM>>
        coordinateTransformationU{
            new Base::HCurlConformingTransformation<DIM>()};
    std::shared_ptr<Base::CoordinateTransformation<DIM>>
        coordinateTransformationP{new Base::H1ConformingTransformation<DIM>()};
    physicalElement.setTransformation(coordinateTransformationU, 0);
    physicalElement.setTransformation(coordinateTransformationP, 1);
    physicalElement.setElement(element);
    physicalElement.setPointReference(point);
    Geometry::PointPhysical<DIM> pointPhysical;
    pointPhysical = element->referenceToPhysical(point);

    // reconstruct the field E = u + grad q.
    for (std::size_t i = 0; i < nb0; ++i) {
        physicalElement.basisFunction(i, phiU, 0);
        result += std::real(coefficients[i]) * phiU;
    }
    for (std::size_t i = 0; i < nb1; ++i) {
        result += std::real(coefficients[nb0 + i]) *
                  physicalElement.basisFunctionDeriv(i, 1);
    }
    return result;
}

template <std::size_t DIM>
double DivDGMaxDiscretization<DIM>::computePotential(
    const Base::Element* element, const Geometry::PointReference<DIM>& point,
    const LinearAlgebra::MiddleSizeVector& coefficients) const {
    logger.log(Log::WARN, "Only computing the real part of the field.");
    LinearAlgebra::SmallVector<DIM> phiU;
    double result;
    std::size_t nb0 = element->getNumberOfBasisFunctions(0),
                nb1 = element->getNumberOfBasisFunctions(1);

    // Setup the physical element with the correct transformations
    Base::PhysicalElement<DIM> physicalElement;
    std::shared_ptr<Base::CoordinateTransformation<DIM>>
        coordinateTransformationU{
            new Base::HCurlConformingTransformation<DIM>()};
    std::shared_ptr<Base::CoordinateTransformation<DIM>>
        coordinateTransformationP{new Base::H1ConformingTransformation<DIM>()};
    physicalElement.setTransformation(coordinateTransformationU, 0);
    physicalElement.setTransformation(coordinateTransformationP, 1);
    physicalElement.setElement(element);
    physicalElement.setPointReference(point);
    Geometry::PointPhysical<DIM> pointPhysical;
    pointPhysical = element->referenceToPhysical(point);

    for (std::size_t i = 0; i < nb1; ++i) {
        result += std::real(coefficients[nb0 + i]) *
                  physicalElement.basisFunction(i, 1);
    }
    return result;
}

template <std::size_t DIM>
void DivDGMaxDiscretization<DIM>::elementMassMatrix(
    Base::PhysicalElement<DIM>& el,
    LinearAlgebra::MiddleSizeMatrix& ret) const {
    const Base::Element* element = el.getElement();
    std::size_t el_size1 = element->getNumberOfBasisFunctions(0);
    std::size_t el_size2 = element->getNumberOfBasisFunctions(1);
    ret.resize(el_size1 + el_size2, el_size1 + el_size2);
    double epsilon =
        static_cast<ElementInfos*>(element->getUserData())->epsilon_;
    LinearAlgebra::SmallVector<DIM> phi_i, phi_j;
    for (std::size_t i = 0; i < el_size1; ++i) {
        el.basisFunction(i, phi_i, 0);
        for (std::size_t j = i; j < el_size1; ++j) {
            el.basisFunction(j, phi_j, 0);
            ret(i, j) = phi_i * phi_j * epsilon;
            ret(j, i) = ret(i, j);
        }
    }
}

template <std::size_t DIM>
void DivDGMaxDiscretization<DIM>::elementStiffnessMatrix(
    Base::PhysicalElement<DIM>& el,
    LinearAlgebra::MiddleSizeMatrix& ret) const {
    const Base::Element* element = el.getElement();
    std::size_t el_size1 = element->getNumberOfBasisFunctions(0);
    std::size_t el_size2 = element->getNumberOfBasisFunctions(1);
    ret.resize(el_size1 + el_size2, el_size1 + el_size2);
    LinearAlgebra::SmallVector<DIM> phi_i, phi_j;

    for (std::size_t i = 0; i < el_size1; ++i) {
        phi_i = el.basisFunctionCurl(i, 0);
        for (std::size_t j = i; j < el_size1; ++j) {
            phi_j = el.basisFunctionCurl(j, 0);
            ret(i, j) = phi_i * phi_j;
            ret(j, i) = ret(i, j);
        }
    }
}

template <std::size_t DIM>
void DivDGMaxDiscretization<DIM>::elementScalarVectorCoupling(
    Base::PhysicalElement<DIM>& el,
    LinearAlgebra::MiddleSizeMatrix& ret) const {
    const Base::Element* element = el.getElement();
    std::size_t uDoFs = element->getNumberOfBasisFunctions(0);
    std::size_t pDoFs = element->getNumberOfBasisFunctions(1);
    ret.resize(uDoFs + pDoFs, uDoFs + pDoFs);
    double epsilon =
        static_cast<ElementInfos*>(element->getUserData())->epsilon_;
    LinearAlgebra::SmallVector<DIM> phi_i, phi_j;

    // Note, this loop only loops over the basis functions for the second
    // unknown. However, it needs the offset for the first unknown.
    for (std::size_t i = uDoFs; i < uDoFs + pDoFs; ++i) {
        phi_i = el.basisFunctionDeriv(i - uDoFs, 1);

        for (std::size_t j = 0; j < uDoFs; ++j) {
            el.basisFunction(j, phi_j, 0);
            ret(j, i) = -1.0 * (phi_i * phi_j) * epsilon;
            ret(i, j) = ret(j, i);
        }
    }
}

template <std::size_t DIM>
void DivDGMaxDiscretization<DIM>::elementSourceVector(
    Base::PhysicalElement<DIM>& el, const InputFunction& source,
    LinearAlgebra::MiddleSizeVector& ret) const {
    const Base::Element* element = el.getElement();
    const Geometry::PointReference<DIM>& p = el.getPointReference();
    std::size_t uDoFs = element->getNumberOfBasisFunctions(0);
    std::size_t pDoFs = element->getNumberOfBasisFunctions(1);
    ret.resize(uDoFs + pDoFs);

    PointPhysicalT pPhys;
    pPhys = element->referenceToPhysical(p);
    LinearAlgebra::SmallVector<DIM> sourceValue, phi;
    source(pPhys, sourceValue);
    for (std::size_t i = 0; i < (uDoFs); ++i) {
        el.basisFunction(i, phi, 0);
        ret(i) = phi * sourceValue;
    }
}

template <std::size_t DIM>
void DivDGMaxDiscretization<DIM>::faceStiffnessMatrix1(
    Base::PhysicalFace<DIM>& fa, LinearAlgebra::MiddleSizeMatrix& ret) const {
    const Base::Face* face = fa.getFace();
    std::size_t totalUDoFs =
        face->getPtrElementLeft()->getNumberOfBasisFunctions(0);
    std::size_t totalPDoFs =
        face->getPtrElementLeft()->getNumberOfBasisFunctions(1);

    const std::size_t leftUDoFs = totalUDoFs;
    const std::size_t leftPDoFs = totalPDoFs;

    if (face->isInternal()) {
        totalUDoFs = face->getPtrElementLeft()->getNumberOfBasisFunctions(0) +
                     face->getPtrElementRight()->getNumberOfBasisFunctions(0);
        totalPDoFs = face->getPtrElementLeft()->getNumberOfBasisFunctions(1) +
                     face->getPtrElementRight()->getNumberOfBasisFunctions(1);
    }
    ret.resize(totalUDoFs + totalPDoFs, totalUDoFs + totalPDoFs);

    std::vector<std::size_t> indices(totalUDoFs + totalPDoFs);
    std::vector<LinearAlgebra::SmallVector<DIM>> phiCurl(totalUDoFs);
    std::vector<LinearAlgebra::SmallVector<DIM>> phiNormal(totalUDoFs);

    // Left element
    for (std::size_t i = 0; i < leftUDoFs; ++i) {
        indices[i] = i;
        phiCurl[i] = fa.basisFunctionCurl(i, 0);
        fa.basisFunctionUnitNormalCross(i, phiNormal[i], 0);
    }
    // Right element
    for (std::size_t i = leftUDoFs; i < totalUDoFs; ++i) {
        // Note, only here we need to offset by leftPDoFs, as the indices are
        // over all unknowns and not just those for one (as in
        // fa.basisFunction*)
        indices[i] = i + leftPDoFs;
        phiCurl[i] = fa.basisFunctionCurl(i, 0);
        fa.basisFunctionUnitNormalCross(i, phiNormal[i], 0);
    }

    double factor = face->isInternal() ? -0.5 : -1;
    for (std::size_t i = 0; i < totalUDoFs; ++i) {
        const std::size_t iIndex = indices[i];
        for (std::size_t j = i; j < totalUDoFs; ++j) {
            const double entry = factor * (phiCurl[i] * phiNormal[j] +
                                           phiCurl[j] * phiNormal[i]);
            const std::size_t jIndex = indices[j];
            ret(iIndex, jIndex) = entry;
            ret(jIndex, iIndex) = entry;
        }
    }
}

template <std::size_t DIM>
void DivDGMaxDiscretization<DIM>::faceStiffnessMatrix2(
    Base::PhysicalFace<DIM>& fa, LinearAlgebra::MiddleSizeMatrix& ret,
    double stab) const {
    const Base::Face* face = fa.getFace();
    double diameter = face->getDiameter();
    std::size_t totalUDoFs =
        face->getPtrElementLeft()->getNumberOfBasisFunctions(0);
    std::size_t totalPDoFs =
        face->getPtrElementLeft()->getNumberOfBasisFunctions(1);
    std::size_t leftUDoFs = totalUDoFs;
    std::size_t leftPDoFs = totalPDoFs;

    if (face->isInternal()) {
        totalUDoFs = face->getPtrElementLeft()->getNumberOfBasisFunctions(0) +
                     face->getPtrElementRight()->getNumberOfBasisFunctions(0);
        totalPDoFs = face->getPtrElementLeft()->getNumberOfBasisFunctions(1) +
                     face->getPtrElementRight()->getNumberOfBasisFunctions(1);
    }

    ret.resize(totalUDoFs + totalPDoFs, totalUDoFs + totalPDoFs);

    std::vector<std::size_t> indices(totalUDoFs);
    std::vector<LinearAlgebra::SmallVector<DIM>> phiNormalCross(totalUDoFs);

    // Left element
    for (std::size_t i = 0; i < leftUDoFs; ++i) {
        indices[i] = i;
        fa.basisFunctionUnitNormalCross(i, phiNormalCross[i], 0);
    }
    // Right element
    for (std::size_t i = leftUDoFs; i < totalUDoFs; ++i) {
        indices[i] = i + leftPDoFs;
        fa.basisFunctionUnitNormalCross(i, phiNormalCross[i], 0);
    }

    for (std::size_t i = 0; i < totalUDoFs; ++i) {
        const std::size_t iIndex = indices[i];
        for (std::size_t j = 0; j < totalUDoFs; ++j) {
            const std::size_t jIndex = indices[j];
            // Possibly scale with mu^{-1} in the future
            double entry =
                stab / (diameter) * (phiNormalCross[i] * phiNormalCross[j]);
            ret(iIndex, jIndex) = entry;
            ret(jIndex, iIndex) = entry;
        }
    }
}

template <std::size_t DIM>
void DivDGMaxDiscretization<DIM>::faceStiffnessMatrix3(
    Base::PhysicalFace<DIM>& fa, LinearAlgebra::MiddleSizeMatrix& ret,
    double stab2) const {
    // TODO: Cleanup.
    const Base::Face* face = fa.getFace();
    double diameter = face->getDiameter();
    std::size_t totalUDoFs =
        face->getPtrElementLeft()->getNumberOfBasisFunctions(0);
    std::size_t totalPDoFs =
        face->getPtrElementLeft()->getNumberOfBasisFunctions(1);
    std::size_t leftUDoFs =
        totalUDoFs;  // to determine left/right part of totalUDoFs
    std::size_t leftPDoFs =
        totalPDoFs;  // to determine left/right part of totalUDoFs
    if (face->isInternal()) {

        totalUDoFs += face->getPtrElementRight()->getNumberOfBasisFunctions(0);
        totalPDoFs += face->getPtrElementRight()->getNumberOfBasisFunctions(1);
        ElementInfos* rightInfo = static_cast<ElementInfos*>(
            face->getPtrElementRight()->getUserData());
        //}
        ret.resize(totalUDoFs + totalPDoFs, totalUDoFs + totalPDoFs);
        ElementInfos* leftInfo = static_cast<ElementInfos*>(
            face->getPtrElementLeft()->getUserData());

        const double epsmax = std::max(leftInfo->epsilon_, rightInfo->epsilon_);

        LinearAlgebra::SmallVector<DIM> normal = fa.getUnitNormalVector();

        std::vector<std::size_t> indices(totalUDoFs);
        std::vector<double> phiNormalEpsilon(totalUDoFs);
        LinearAlgebra::SmallVector<DIM> phi;

        // Left element
        for (std::size_t i = 0; i < totalUDoFs; ++i) {
            indices[i] = i;
            fa.basisFunction(i, phi, 0);
            phiNormalEpsilon[i] = leftInfo->epsilon_ * (phi * normal);
        }
        // Right element
        for (std::size_t i = leftUDoFs; i < totalUDoFs; ++i) {
            indices[i] = i + leftPDoFs;
            fa.basisFunction(i, phi, 0);
            // - from the conversion between the left and right normal
            phiNormalEpsilon[i] = -rightInfo->epsilon_ * (phi * normal);
        }

        for (std::size_t i = 0; i < totalUDoFs; ++i) {
            const std::size_t iIndex = indices[i];
            for (std::size_t j = i; j < totalUDoFs; ++j) {
                const std::size_t jIndex = indices[j];
                // TODO: Scaling
                double entry = diameter * stab2 / epsmax * phiNormalEpsilon[i] *
                               phiNormalEpsilon[j];
                ret(iIndex, jIndex) = entry;
                ret(jIndex, iIndex) = entry;
            }
        }
    } else {
        ret.resize(totalUDoFs + totalPDoFs, totalUDoFs + totalPDoFs);
    }
}

template <std::size_t DIM>
void DivDGMaxDiscretization<DIM>::faceScalarVectorCoupling(
    Base::PhysicalFace<DIM>& fa, LinearAlgebra::MiddleSizeMatrix& ret) const {
    const Base::Face* face = fa.getFace();
    std::size_t totalUDoFs =
        face->getPtrElementLeft()->getNumberOfBasisFunctions(0);
    std::size_t totalPDoFs =
        face->getPtrElementLeft()->getNumberOfBasisFunctions(1);

    const double epsilonLeft =
        static_cast<ElementInfos*>(face->getPtrElementLeft()->getUserData())
            ->epsilon_;
    const double epsilonRight =
        face->isInternal()
            ? (static_cast<ElementInfos*>(
                   face->getPtrElementRight()->getUserData())
                   ->epsilon_)
            : 1.0;  // If no right face is present this will not be used.
    // From the averaging terms.
    double averageFactor = face->isInternal() ? 0.5 : 1;

    std::size_t leftUDofs = totalUDoFs;
    std::size_t leftPDofs = totalPDoFs;

    if (face->isInternal()) {
        totalUDoFs += face->getPtrElementRight()->getNumberOfBasisFunctions(0);
        totalPDoFs += face->getPtrElementRight()->getNumberOfBasisFunctions(1);
    }
    ret.resize(totalUDoFs + totalPDoFs, totalUDoFs + totalPDoFs);

    LinearAlgebra::SmallVector<DIM> normal = fa.getUnitNormalVector();

    std::vector<std::size_t> indicesU(totalUDoFs);
    std::vector<std::size_t> indicesP(totalPDoFs);
    // Note: while the discretization_ say {{ phiU }} * [[phiP]], which leads to
    // terms (phiP_[lr] n_[lr]) * phiU (subscripts for the left/right element).
    // We could thus multiply phiP_[lr] by the correct normal and compute the
    // innerproduct later. For better performance we use that n_l = -n_r, thus
    // the terms can also be written as (s_[lr] phiP_[lr]) (n_l * phiU), where
    // s_[lr] is the sign (1 for l, -1 for r). Note that now both terms are
    // numbers instead of vectors.
    std::vector<double> phiUNormal(totalUDoFs);
    std::vector<double> phiP(totalPDoFs);
    LinearAlgebra::SmallVector<DIM> phiU;

    // Left element U
    for (std::size_t i = 0; i < leftUDofs; ++i) {
        indicesU[i] = i;
        fa.basisFunction(i, phiU, 0);
        phiUNormal[i] = phiU * normal * epsilonLeft;
    }
    // Right element U
    // Note that for a boundary element leftUDofs == totalUDoFs and this loop is
    // skipped.
    for (std::size_t i = leftUDofs; i < totalUDoFs; ++i) {
        indicesU[i] = i + leftPDofs;
        fa.basisFunction(i, phiU, 0);
        phiUNormal[i] = phiU * normal * epsilonRight;
    }

    // Left element P
    for (std::size_t i = 0; i < leftPDofs; ++i) {
        indicesP[i] = leftUDofs + i;
        phiP[i] = fa.basisFunction(i, 1);
    }
    // Right element P, also skipped for boundary elements.
    for (std::size_t i = leftPDofs; i < totalPDoFs; ++i) {
        indicesP[i] = totalUDoFs + i;
        // - because the right normal is -1 times the left normal.
        phiP[i] = -fa.basisFunction(i, 1);
    }

    for (std::size_t i = 0; i < totalUDoFs; ++i) {
        const std::size_t iIndex = indicesU[i];
        const double& phiUi = phiUNormal[i];
        for (std::size_t j = 0; j < totalPDoFs; ++j) {
            const double entry = averageFactor * (phiUi * phiP[j]);
            const std::size_t jIndex = indicesP[j];
            ret(iIndex, jIndex) = entry;
            ret(jIndex, iIndex) = entry;
        }
    }
}

template <std::size_t DIM>
void DivDGMaxDiscretization<DIM>::faceStiffnessScalarMatrix4(
    Base::PhysicalFace<DIM>& fa, LinearAlgebra::MiddleSizeMatrix& ret,
    double stab3) const {
    // TODO: Cleanup.
    const Base::Face* face = fa.getFace();
    double diameter = face->getDiameter();
    std::size_t totalUDoFs =
        face->getPtrElementLeft()->getNumberOfBasisFunctions(0);
    std::size_t totalPDoFs =
        face->getPtrElementLeft()->getNumberOfBasisFunctions(1);

    const std::size_t leftUDoFs = totalUDoFs;
    const std::size_t leftPDoFs = totalPDoFs;

    const double epsilonLeft =
        static_cast<ElementInfos*>(face->getPtrElementLeft()->getUserData())
            ->epsilon_;

    double epsmax;
    if (face->isInternal()) {
        totalUDoFs += face->getPtrElementRight()->getNumberOfBasisFunctions(0);
        totalPDoFs += face->getPtrElementRight()->getNumberOfBasisFunctions(1);
        const double epsilonRight =
            static_cast<ElementInfos*>(
                face->getPtrElementRight()->getUserData())
                ->epsilon_;
        epsmax = std::max(epsilonLeft, epsilonRight);
    } else {
        epsmax = epsilonLeft;
    }
    ret.resize(totalUDoFs + totalPDoFs, totalUDoFs + totalPDoFs);

    std::vector<std::size_t> indices(totalPDoFs);
    std::vector<double> phiP(totalPDoFs);
    // Note we leave out the actual normal as nL * nL = +1, nR * nL = -1, etc.
    // instead we multiply the right functions by -1.

    // Left element
    for (std::size_t i = 0; i < leftPDoFs; ++i) {
        indices[i] = leftUDoFs + i;
        phiP[i] = fa.basisFunction(i, 1);
    }
    // Right element
    for (std::size_t i = leftPDoFs; i < totalPDoFs; ++i) {
        indices[i] = totalUDoFs + i;
        // - from the right normal being -1 times the left normal.
        phiP[i] = -fa.basisFunction(i, 1);
    }

    for (std::size_t i = 0; i < totalPDoFs; ++i) {
        const std::size_t iIndex = indices[i];
        const double& phiPi = phiP[i];
        for (std::size_t j = i; j < totalPDoFs; ++j) {
            // Note, positive here. The minus for the C matrix is added in
            // computeFaceIntegrals7.
            const double entry = stab3 / diameter * epsmax * phiP[j] * phiPi;
            const std::size_t jIndex = indices[j];
            ret(iIndex, jIndex) = entry;
            ret(jIndex, iIndex) = entry;
        }
    }
}

template <std::size_t DIM>
LinearAlgebra::MiddleSizeMatrix
    DivDGMaxDiscretization<DIM>::computeFaceVectorMassMatrix(
        typename Base::MeshManipulator<DIM>::FaceIterator rawFace) const {
    Integration::ElementIntegral<DIM> elementIntegral;

    elementIntegral.setTransformation(
        std::shared_ptr<Base::CoordinateTransformation<DIM>>(
            new Base::HCurlConformingTransformation<DIM>()),
        0);
    elementIntegral.setTransformation(
        std::shared_ptr<Base::CoordinateTransformation<DIM>>(
            new Base::H1ConformingTransformation<DIM>()),
        1);

    FaceDoFInfo faceInfo = getFaceDoFInfo(**rawFace);

    LinearAlgebra::MiddleSizeMatrix result = elementIntegral.integrate(
        (*rawFace)->getPtrElementLeft(),
        [&](Base::PhysicalElement<DIM>& element) {
            LinearAlgebra::MiddleSizeMatrix lmassMat(faceInfo.totalUDoFs(),
                                                     faceInfo.totalUDoFs());

            // Storage of basis functions
            LinearAlgebra::SmallVector<DIM> phii, phij;

            for (std::size_t i = 0; i < faceInfo.leftUDoFs; ++i) {
                element.basisFunction(i, phii, 0);
                for (std::size_t j = i; j < faceInfo.leftUDoFs; ++j) {
                    element.basisFunction(j, phij, 0);
                    const double value = phii * phij;
                    lmassMat(i, j) = value;
                    lmassMat(j, i) = value;
                }
            }
            return lmassMat;
        });
    if (faceInfo.internal) {
        // Build the mass matrix for the right element if it exists.
        result += elementIntegral.integrate(
            (*rawFace)->getPtrElementRight(),
            [&](Base::PhysicalElement<DIM>& element) {
                LinearAlgebra::MiddleSizeMatrix rmassMat(faceInfo.totalUDoFs(),
                                                         faceInfo.totalUDoFs());

                // Storage of basis functions
                LinearAlgebra::SmallVector<DIM> phii, phij;
                // Offset to correctly place mass matrix as second block on the
                // diagonal.
                std::size_t offset = faceInfo.leftUDoFs;
                for (std::size_t i = 0; i < faceInfo.rightUDoFs; ++i) {
                    element.basisFunction(i, phii, 0);
                    for (std::size_t j = i; j < faceInfo.rightUDoFs; ++j) {
                        element.basisFunction(j, phij, 0);
                        const double value = phii * phij;
                        rmassMat(offset + i, offset + j) = value;
                        rmassMat(offset + j, offset + i) = value;
                    }
                }
                return rmassMat;
            });
    }
    return result;
}

template <std::size_t DIM>
LinearAlgebra::MiddleSizeMatrix
    DivDGMaxDiscretization<DIM>::computeFaceScalarMassMatrix(
        typename Base::MeshManipulator<DIM>::FaceIterator rawFace) const {
    Integration::ElementIntegral<DIM> elementIntegral;

    elementIntegral.setTransformation(
        std::shared_ptr<Base::CoordinateTransformation<DIM>>(
            new Base::HCurlConformingTransformation<DIM>()),
        0);
    elementIntegral.setTransformation(
        std::shared_ptr<Base::CoordinateTransformation<DIM>>(
            new Base::H1ConformingTransformation<DIM>()),
        1);

    FaceDoFInfo faceInfo = getFaceDoFInfo(**rawFace);

    LinearAlgebra::MiddleSizeMatrix result = elementIntegral.integrate(
        (*rawFace)->getPtrElementLeft(),
        [&](Base::PhysicalElement<DIM>& element) {
            LinearAlgebra::MiddleSizeMatrix lmassMat(faceInfo.totalPDoFs(),
                                                     faceInfo.totalPDoFs());

            // Storage of basis functions
            double phii, phij;

            for (std::size_t i = 0; i < faceInfo.leftPDoFs; ++i) {
                phii = element.basisFunction(i, 1);
                for (std::size_t j = i; j < faceInfo.leftPDoFs; ++j) {
                    phij = element.basisFunction(j, 1);
                    const double value = phii * phij;
                    lmassMat(i, j) = value;
                    lmassMat(j, i) = value;
                }
            }
            return lmassMat;
        });
    if (faceInfo.internal) {
        // Build the mass matrix for the right element if it exists.
        result += elementIntegral.integrate(
            (*rawFace)->getPtrElementRight(),
            [&](Base::PhysicalElement<DIM>& element) {
                LinearAlgebra::MiddleSizeMatrix rmassMat(faceInfo.totalPDoFs(),
                                                         faceInfo.totalPDoFs());

                // Storage of basis functions
                double phii, phij;
                // Offset to correctly place mass matrix as second block on the
                // diagonal.
                std::size_t offset = faceInfo.leftPDoFs;
                for (std::size_t i = 0; i < faceInfo.rightPDoFs; ++i) {
                    phii = element.basisFunction(i, 1);
                    for (std::size_t j = i; j < faceInfo.rightPDoFs; ++j) {
                        phij = element.basisFunction(j, 1);
                        const double value = phii * phij;
                        rmassMat(offset + i, offset + j) = value;
                        rmassMat(offset + j, offset + i) = value;
                    }
                }
                return rmassMat;
            });
    }
    logger.assert_debug(result.isHermitian(1e-16), "Non hermitian mass matrix");
    return result;
}

template <std::size_t DIM>
LinearAlgebra::MiddleSizeMatrix
    DivDGMaxDiscretization<DIM>::computeScalarLiftProjector(
        typename Base::MeshManipulator<DIM>::FaceIterator rawFace) const {
    Integration::FaceIntegral<DIM> faceIntegral;
    faceIntegral.setTransformation(
        std::shared_ptr<Base::CoordinateTransformation<DIM>>(
            new Base::HCurlConformingTransformation<DIM>()),
        0);
    faceIntegral.setTransformation(
        std::shared_ptr<Base::CoordinateTransformation<DIM>>(
            new Base::H1ConformingTransformation<DIM>()),
        1);
    FaceDoFInfo faceInfo = getFaceDoFInfo(**rawFace);

    return faceIntegral.integrate(*rawFace, [&](Base::PhysicalFace<DIM>& face) {
        LinearAlgebra::MiddleSizeMatrix result;
        result.resize(faceInfo.totalUDoFs(), faceInfo.totalPDoFs());
        const LinearAlgebra::SmallVector<DIM>& normal =
            face.getUnitNormalVector();
        LinearAlgebra::SmallVector<DIM> basisV;
        double factor = faceInfo.internal ? 0.5 : 1;  // For the average
        for (std::size_t i = 0; i < faceInfo.totalPDoFs(); ++i) {
            // For the normal.
            double nfactor = i < faceInfo.leftPDoFs ? 1 : -1;
            for (std::size_t j = 0; j < faceInfo.totalUDoFs(); ++j) {
                face.basisFunction(j, basisV, 0);  // v
                // [[p]] {{v}}
                result(j, i) = face.basisFunction(i, 1) * nfactor * factor *
                               (normal * basisV);
            }
        }
        return result;
    });
}

template <std::size_t DIM>
LinearAlgebra::MiddleSizeMatrix
    DivDGMaxDiscretization<DIM>::computeVectorLiftProjector(
        typename Base::MeshManipulator<DIM>::FaceIterator rawFace) const {
    Integration::FaceIntegral<DIM> faceIntegral;
    faceIntegral.setTransformation(
        std::shared_ptr<Base::CoordinateTransformation<DIM>>(
            new Base::HCurlConformingTransformation<DIM>()),
        0);
    faceIntegral.setTransformation(
        std::shared_ptr<Base::CoordinateTransformation<DIM>>(
            new Base::H1ConformingTransformation<DIM>()),
        1);
    FaceDoFInfo faceInfo = getFaceDoFInfo(**rawFace);

    return faceIntegral.integrate(*rawFace, [&](Base::PhysicalFace<DIM>& face) {
        LinearAlgebra::MiddleSizeMatrix result;
        result.resize(faceInfo.totalUDoFs(), faceInfo.totalUDoFs());
        const LinearAlgebra::SmallVector<DIM>& normal =
            face.getUnitNormalVector();
        LinearAlgebra::SmallVector<DIM> basisUnormal, basisV;
        double factor = faceInfo.internal ? 0.5 : 1;  // For the average
        for (std::size_t i = 0; i < faceInfo.totalUDoFs(); ++i) {
            face.basisFunction(i, basisV, 0);
            for (std::size_t j = 0; j < faceInfo.totalUDoFs(); ++j) {
                face.basisFunctionUnitNormalCross(j, basisUnormal, 0);
                result(i, j) = factor * (basisUnormal * basisV);
            }
        }
        return result;
    });
}

template <std::size_t DIM>
LinearAlgebra::MiddleSizeMatrix
    DivDGMaxDiscretization<DIM>::computeVectorNormalLiftProjector(
        typename Base::MeshManipulator<DIM>::FaceIterator rawFace) const {
    Integration::FaceIntegral<DIM> faceIntegral;
    faceIntegral.setTransformation(
        std::shared_ptr<Base::CoordinateTransformation<DIM>>(
            new Base::HCurlConformingTransformation<DIM>()),
        0);
    faceIntegral.setTransformation(
        std::shared_ptr<Base::CoordinateTransformation<DIM>>(
            new Base::H1ConformingTransformation<DIM>()),
        1);
    FaceDoFInfo faceInfo = getFaceDoFInfo(**rawFace);

    double epsilonLeft = static_cast<ElementInfos*>(
                             (*rawFace)->getPtrElementLeft()->getUserData())
                             ->epsilon_;
    double epsilonRight =
        faceInfo.internal ? static_cast<ElementInfos*>(
                                (*rawFace)->getPtrElementRight()->getUserData())
                                ->epsilon_
                          : 0.0;  // Should not be used.

    return faceIntegral.integrate(*rawFace, [&](Base::PhysicalFace<DIM>& face) {
        LinearAlgebra::MiddleSizeMatrix result;
        result.resize(faceInfo.totalPDoFs(), faceInfo.totalUDoFs());
        LinearAlgebra::SmallVector<DIM> basisU, normal;
        normal = face.getUnitNormalVector();
        double factor = faceInfo.internal ? 0.5 : 1.0;  // From the average
        for (std::size_t i = 0; i < faceInfo.totalPDoFs(); ++i) {
            double basisP = face.basisFunction(i, 1);
            for (std::size_t j = 0; j < faceInfo.totalUDoFs(); ++j) {
                face.basisFunction(j, basisU, 0);
                double uvalue = basisU * normal;  // [[epsilon u]]_N
                if (j > faceInfo.leftUDoFs) {
                    // - form the change in normal
                    uvalue *= -epsilonRight;
                } else {
                    uvalue *= epsilonLeft;
                }
                result(i, j) = factor * basisP * uvalue;
            }
        }
        return result;
    });
}

/// Distribute a face matrix with only vector or scalar components to the whole
/// matrix
///
/// \param faceInfo The degree of freedom information for the face
/// \param vector Whether the matrices are for the vector (true) or scalar
/// (false) degree of freedoms. \param source The source matrix to distribute
/// \param target The target matrix to which to add the entries from the source.
void distributeFaceMatrix(FaceDoFInfo faceInfo, bool vector,
                          LinearAlgebra::MiddleSizeMatrix& source,
                          LinearAlgebra::MiddleSizeMatrix& target) {
    logger.assert_debug(source.getNumberOfRows() == source.getNumberOfColumns(),
                        "Non square source matrix.");
    logger.assert_debug(target.getNumberOfRows() == target.getNumberOfColumns(),
                        "Non square source matrix.");
    std::size_t dimension =
        vector ? faceInfo.totalUDoFs() : faceInfo.totalPDoFs();
    logger.assert_debug(dimension == source.getNumberOfColumns(),
                        "Incorrect source dimension");
    logger.assert_debug(faceInfo.totalDoFs() == target.getNumberOfColumns(),
                        "Incorrect target dimension");

    // Offsets for the left and right degrees of freedom
    std::size_t leftOffset = vector ? 0 : faceInfo.leftUDoFs;
    std::size_t rightOffset =
        vector ? faceInfo.leftPDoFs : faceInfo.totalUDoFs();
    // Degrees of freedom for the left face.
    std::size_t leftDimension =
        vector ? faceInfo.leftUDoFs : faceInfo.leftPDoFs;

    // Mapping from source to target index.
    std::vector<std::size_t> faceIndices(dimension);
    for (std::size_t i = 0; i < leftDimension; ++i) {
        faceIndices[i] = i + leftOffset;
    }
    for (std::size_t i = leftDimension; i < dimension; ++i) {
        faceIndices[i] = i + rightOffset;
    }

    for (std::size_t i = 0; i < dimension; ++i) {
        for (std::size_t j = 0; j < dimension; ++j) {
            target(faceIndices[i], faceIndices[j]) += source(i, j);
        }
    }
}

template <std::size_t DIM>
LinearAlgebra::MiddleSizeMatrix
    DivDGMaxDiscretization<DIM>::brezziFluxBilinearTerm(
        typename Base::MeshManipulator<DIM>::FaceIterator rawFace,
        Stab stab) const {
    FaceDoFInfo faceInfo = getFaceDoFInfo(**rawFace);
    LinearAlgebra::MiddleSizeMatrix result(faceInfo.totalDoFs(),
                                           faceInfo.totalDoFs());

    // Mass matrix for vector valued components
    LinearAlgebra::MiddleSizeMatrix massMat =
        computeFaceVectorMassMatrix(rawFace);
    massMat.cholesky();

    if (stab.fluxType1 == FluxType::BREZZI) {
        std::size_t leftFaces =
            (*rawFace)->getPtrElementLeft()->getNumberOfFaces();
        std::size_t rightFaces =
            faceInfo.internal
                ? (*rawFace)->getPtrElementLeft()->getNumberOfFaces()
                : 0;

        // Compute vector stabilizer, R^T L^{-T} (stab1 + numFaces * muinv)
        // L^{-1} R with R the lifting matrix for [[u_j]]_T {{u_i}} (stab1 +
        // numFace*muinv) the diagonal matrix with numFaces the number of faces
        // of the element (left on top diagonal, right bottom diagonal) L the
        // lower triangular part of the cholesky decomposition of M.
        LinearAlgebra::MiddleSizeMatrix vectorStabilizer =
            computeVectorLiftProjector(rawFace);
        // massMat.solve(vectorStabilizer); // M^{-1}R
        massMat.solveLowerTriangular(vectorStabilizer);  // L^{-1}R, with LL^{T}
                                                         // = M
        // Inplace multiplication with diagonal matrix
        double lfactor = std::sqrt(stab.stab1 + leftFaces);
        for (std::size_t i = 0; i < faceInfo.totalUDoFs(); ++i) {
            for (std::size_t j = 0; j < faceInfo.leftUDoFs; ++j) {
                vectorStabilizer(i, j) *= lfactor;
            }
        }
        if (faceInfo.internal) {
            double rfactor = std::sqrt(stab.stab1 + rightFaces);
            if (faceInfo.internal && (leftFaces != rightFaces)) {
                logger(WARN,
                       "Brezzi stabilization with unequal number of left and "
                       "right faces is untested");
            }
            for (std::size_t i = 0; i < faceInfo.totalUDoFs(); ++i) {
                for (std::size_t j = faceInfo.leftUDoFs;
                     j < faceInfo.totalUDoFs(); ++j) {
                    vectorStabilizer(i, j) *= rfactor;
                }
            }
        }
        vectorStabilizer = vectorStabilizer.transpose() * vectorStabilizer;
        // Copy stabilizer to result matrix
        distributeFaceMatrix(faceInfo, true, vectorStabilizer, result);
    }

    if (faceInfo.internal && stab.fluxType2 == FluxType::BREZZI) {
        // TODO: Make hermitian and reduce hermiticity tolerance below
        // Compute normal vector stabilizer, stab2 * S^T M^{-T} S, S = lift
        // matrix Again we are using the Cholesky decomposition off LL^T = M.
        LinearAlgebra::MiddleSizeMatrix pmassMat =
            computeFaceScalarMassMatrix(rawFace);
        pmassMat.cholesky();
        LinearAlgebra::MiddleSizeMatrix vectorStabilizer =
            computeVectorNormalLiftProjector(rawFace);
        pmassMat.solveLowerTriangular(vectorStabilizer);  // L^{-T} S
        vectorStabilizer = vectorStabilizer.transpose() * vectorStabilizer;
        vectorStabilizer *= stab.stab2;
        // Copy stabilizer to result matrix
        distributeFaceMatrix(faceInfo, true, vectorStabilizer, result);
    }

    if (stab.fluxType3 == FluxType::BREZZI) {
        // Compute scalar stabilizer, stab3 * R^T M^{-1} R, R = lift matrix
        // using stab3 * R^T L^{-T} L{-1} R, where LL^{T} = M
        LinearAlgebra::MiddleSizeMatrix scalarStabilizer =
            computeScalarLiftProjector(rawFace);
        massMat.solveLowerTriangular(scalarStabilizer);  // L^{-1}R
        scalarStabilizer = scalarStabilizer.transpose() * scalarStabilizer;
        // minus sign as it is -c_h(p,q)
        scalarStabilizer *= -stab.stab3;

        // Copy the scalar stabilizer into the result matrix
        distributeFaceMatrix(faceInfo, false, scalarStabilizer, result);
    }
    logger.assert_debug(result.isHermitian(0),
                        "Non Hermitian penalty parameter.");
    return result;
}

template <std::size_t DIM>
void DivDGMaxDiscretization<DIM>::faceBoundaryVector(
    Base::PhysicalFace<DIM>& fa,
    const DivDGMaxDiscretization<DIM>::FaceInputFunction& boundaryValue,
    LinearAlgebra::MiddleSizeVector& ret, Stab stab) const {
    const Base::Face* face = fa.getFace();
    const Geometry::PointReference<DIM - 1>& p = fa.getPointReference();

    if (face->isInternal()) {
        std::size_t totalUDoFs =
            face->getPtrElementLeft()->getNumberOfBasisFunctions(0) +
            face->getPtrElementRight()->getNumberOfBasisFunctions(0);
        std::size_t totalPDoFs =
            face->getPtrElementLeft()->getNumberOfBasisFunctions(1) +
            face->getPtrElementRight()->getNumberOfBasisFunctions(1);

        ret.resize(totalUDoFs + totalPDoFs);

        for (std::size_t i = 0; i < (totalUDoFs + totalPDoFs); ++i) {
            ret(i) = 0;
        }
    } else {
        double diameter = face->getDiameter();
        const Geometry::PointReference<DIM>& PLeft =
            face->mapRefFaceToRefElemL(p);
        const PointPhysicalT PPhys =
            face->getPtrElementLeft()->referenceToPhysical(PLeft);

        LinearAlgebra::SmallVector<DIM> val, phi_curl;
        LinearAlgebra::SmallVector<DIM> phi;
        boundaryValue(PPhys, fa, val);

        std::size_t totalUDoFs =
            face->getPtrElementLeft()->getNumberOfBasisFunctions(0);
        std::size_t totalPDoFs =
            face->getPtrElementLeft()->getNumberOfBasisFunctions(1);
        ret.resize(totalUDoFs + totalPDoFs);

        for (std::size_t i = 0; i < totalUDoFs; ++i) {
            fa.basisFunctionUnitNormalCross(i, phi, 0);
            phi_curl = fa.basisFunctionCurl(i, 0);
            double value = -(phi_curl * val);
            if (stab.fluxType1 == FluxType::IP) {
                value += stab.stab1 / (diameter) * (phi * val);
            }
            // Scale with mu^{-1} in the future
            ret(i) = value;
        }
    }
}

template <std::size_t DIM>
LinearAlgebra::MiddleSizeVector
    DivDGMaxDiscretization<DIM>::brezziFluxBoundaryVector(
        typename Base::MeshManipulator<DIM>::FaceIterator rawFace,
        const DivDGMaxDiscretization<DIM>::FaceInputFunction& boundaryValue,
        DivDGMaxDiscretization<DIM>::Stab stab) const {
    FaceDoFInfo faceInfo = getFaceDoFInfo(**rawFace);
    if (faceInfo.internal) {
        LinearAlgebra::MiddleSizeVector result;
        result.resize(faceInfo.totalDoFs());
        return result;
    }

    // Contribution R^T M^{-1} (stab + numFaces*muinv) r, with
    // r = integral_F bv . u_i dS
    // and R, M the lifting matrix and mass matrix for the vector component

    LinearAlgebra::MiddleSizeMatrix massMatrix =
        computeFaceVectorMassMatrix(rawFace);
    LinearAlgebra::MiddleSizeMatrix liftMatrix =
        computeVectorLiftProjector(rawFace);

    // Compute the r-vector
    Integration::FaceIntegral<DIM> faceIntegral;
    faceIntegral.setTransformation(
        std::shared_ptr<Base::CoordinateTransformation<DIM>>(
            new Base::HCurlConformingTransformation<DIM>()),
        0);
    faceIntegral.setTransformation(
        std::shared_ptr<Base::CoordinateTransformation<DIM>>(
            new Base::H1ConformingTransformation<DIM>()),
        1);

    LinearAlgebra::MiddleSizeVector rvector =
        faceIntegral.integrate(*rawFace, [&](Base::PhysicalFace<DIM>& face) {
            LinearAlgebra::MiddleSizeVector result;
            result.resize(faceInfo.totalUDoFs());
            LinearAlgebra::SmallVector<DIM> basisV;

            // Compute boundary value
            const Base::Face* baseFace = face.getFace();
            const Geometry::PointReference<DIM - 1>& p =
                face.getPointReference();
            const Geometry::PointReference<DIM>& PLeft =
                baseFace->mapRefFaceToRefElemL(p);
            const PointPhysicalT PPhys =
                baseFace->getPtrElementLeft()->referenceToPhysical(PLeft);
            LinearAlgebra::SmallVector<DIM> val;
            boundaryValue(PPhys, face, val);

            for (std::size_t i = 0; i < faceInfo.totalUDoFs(); ++i) {
                face.basisFunction(i, basisV, 0);
                result(i) = val * basisV;
            }
            return result;
        });
    rvector *= stab.stab1 + (*rawFace)->getPtrElementLeft()->getNumberOfFaces();
    massMatrix.transpose().solve(rvector);
    rvector = liftMatrix.transpose() * rvector;
    // Extend the  result to include degrees of freedom for the scalar part
    LinearAlgebra::MiddleSizeVector result(faceInfo.totalDoFs());
    for (std::size_t i = 0; i < rvector.size(); ++i) {
        result(i) = rvector(i);
    }
    return result;
}

template <std::size_t DIM>
double DivDGMaxDiscretization<DIM>::computeL2Error(

    Base::MeshManipulator<DIM>& mesh, std::size_t timeVector,
    const DivDGMaxDiscretization<DIM>::InputFunction& electricField) const {
    Integration::ElementIntegral<DIM> elIntegral;
    elIntegral.setTransformation(
        std::shared_ptr<Base::CoordinateTransformation<DIM>>(
            new Base::HCurlConformingTransformation<DIM>()),
        0);
    elIntegral.setTransformation(
        std::shared_ptr<Base::CoordinateTransformation<DIM>>(
            new Base::H1ConformingTransformation<DIM>()),
        1);

    double error = 0;
    auto end = mesh.elementColEnd();
    for (typename Base::MeshManipulator<DIM>::ElementIterator it =
             mesh.elementColBegin();
         it != end; ++it) {
        error +=
            elIntegral.integrate((*it), [&](Base::PhysicalElement<DIM>& el) {
                return elementErrorIntegrand(el, timeVector, electricField);
            });
    }
    return std::sqrt(error);
}

template <std::size_t DIM>
double DivDGMaxDiscretization<DIM>::elementErrorIntegrand(
    Base::PhysicalElement<DIM>& el, std::size_t timeVector,
    const DivDGMaxDiscretization<DIM>::InputFunction& exactValues) const {
    const Base::Element* element = el.getElement();
    const Geometry::PointPhysical<DIM>& pPhys = el.getPointPhysical();

    std::size_t numberOfUDoFs = element->getNumberOfBasisFunctions(0),
                numberOfPDoFs = element->getNumberOfBasisFunctions(1);
    LinearAlgebra::MiddleSizeVector data =
        element->getTimeIntegrationVector(timeVector);

    LinearAlgebra::SmallVector<DIM> error, phi;
    exactValues(pPhys, error);
    for (std::size_t i = 0; i < numberOfUDoFs; ++i) {
        el.basisFunction(i, phi, 0);
        error -= std::real(data[i]) * phi;
    }
    for (std::size_t i = 0; i < numberOfPDoFs; ++i) {
        error -=
            std::real(data[i + numberOfUDoFs]) * el.basisFunctionDeriv(i, 1);
    }

    return error.l2NormSquared();
}

template <std::size_t DIM>
char fluxName(typename DivDGMaxDiscretization<DIM>::FluxType f) {
    switch (f) {
        case DivDGMaxDiscretization<DIM>::FluxType::BREZZI:
            return 'b';
        case DivDGMaxDiscretization<DIM>::FluxType::IP:
            return 'i';
        default:
            logger.assert_always(false, "Unknown flux type.");
            return '0';
    }
}

template <std::size_t DIM>
std::ostream& printStab(
    std::ostream& os, const typename DivDGMaxDiscretization<DIM>::Stab& stab) {
    os << "Stab{" << fluxName<DIM>(stab.fluxType1) << "=" << stab.stab1 << ", "
       << fluxName<DIM>(stab.fluxType2) << "=" << stab.stab2 << ", "
       << fluxName<DIM>(stab.fluxType3) << "=" << stab.stab3 << "}";
    return os;
}

std::ostream& operator<<(std::ostream& os,
                         typename DivDGMaxDiscretization<2>::Stab& stab) {
    return printStab<2>(os, stab);
}
std::ostream& operator<<(std::ostream& os,
                         typename DivDGMaxDiscretization<3>::Stab& stab) {
    return printStab<3>(os, stab);
}

template class DivDGMaxDiscretization<2>;
template class DivDGMaxDiscretization<3>;

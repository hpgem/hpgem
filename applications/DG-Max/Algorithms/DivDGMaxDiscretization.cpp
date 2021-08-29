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

#include "Utilities/ElementLocalIndexing.h"
#include "Utilities/FaceLocalIndexing.h"

#include "ElementInfos.h"

using namespace hpgem;

// Definition of the constants to reference to.
const std::size_t DivDGMaxDiscretizationBase::ELEMENT_MASS_MATRIX_ID;
const std::size_t DivDGMaxDiscretizationBase::ELEMENT_STIFFNESS_MATRIX_ID;
const std::size_t DivDGMaxDiscretizationBase::ELEMENT_SOURCE_VECTOR_ID;
const std::size_t DivDGMaxDiscretizationBase::FACE_STIFFNESS_MATRIX_ID;
const std::size_t DivDGMaxDiscretizationBase::FACE_BOUNDARY_VECTOR_ID;

/// Helper struct to access the material information in elements adjacent to a
/// face.
struct FaceMaterialInfo {
    template <std::size_t DIM>
    FaceMaterialInfo(Base::PhysicalFace<DIM>& fa) {
        const Base::Face* face = fa.getFace();
        auto* leftInfo = dynamic_cast<ElementInfos*>(
            face->getPtrElementLeft()->getUserData());
        epsilonLeft = leftInfo->epsilon_;
        epsMax = epsilonLeft;

        if (face->isInternal()) {
            auto* rightInfo = dynamic_cast<ElementInfos*>(
                face->getPtrElementRight()->getUserData());
            epsilonRight = rightInfo->epsilon_;
            if (epsilonLeft < epsilonRight) {
                epsMax = epsilonRight;
            }
        } else {
            // Just in case
            epsilonRight = std::numeric_limits<double>::signaling_NaN();
        }
    }

    double epsilonLeft;
    double epsilonRight;
    double epsMax;
};

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
DivDGMaxDiscretization<DIM>::DivDGMaxDiscretization()
    : fieldTransform_(
          std::make_shared<Base::HCurlConformingTransformation<DIM>>()),
      potentialTransform_(
          std::make_shared<Base::H1ConformingTransformation<DIM>>()) {
    elementIntegrator_.setTransformation(fieldTransform_, 0);
    elementIntegrator_.setTransformation(potentialTransform_, 1);
    faceIntegrator_.setTransformation(fieldTransform_, 0);
    faceIntegrator_.setTransformation(potentialTransform_, 1);
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
        initialConditionDerivative) {

    Utilities::ElementLocalIndexing indexing;
    indexing.reinit({0, 1});

    // TODO: Add initial condition integration.
    LinearAlgebra::MiddleSizeMatrix massMatrix(2, 2), stiffnessMatrix(2, 2);

    LinearAlgebra::MiddleSizeVector vector1(2), vector2(2), sourceVector(2);

    auto end = mesh.elementColEnd();
    for (typename Base::MeshManipulator<DIM>::ElementIterator it =
             mesh.elementColBegin();
         it != end; ++it) {
        Base::Element* element = *it;
        indexing.reinit(element);
        computeElementMatrices(element, indexing);

        if (invertMassMatrix) {
            // Note reference to allow overwriting it
            LinearAlgebra::MiddleSizeMatrix& massMatrix =
                element->getElementMatrix(ELEMENT_MASS_MATRIX_ID);
            massMatrix = massMatrix.inverse();
        }

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
            sourceVector.resize(indexing.getNumberOfDoFs());
            sourceVector = elementIntegrator_.integrate(
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
    Stab stab) {
    LinearAlgebra::MiddleSizeMatrix faceMatrix;
    LinearAlgebra::MiddleSizeVector faceVector;

    Utilities::FaceLocalIndexing indexing;
    indexing.reinit({0, 1});

    auto end = mesh.faceColEnd();
    for (auto it = mesh.faceColBegin(); it != end; ++it) {

        const Base::Face* face = *it;
        indexing.reinit(face);

        std::size_t totalDoFs = indexing.getNumberOfDoFs();
        faceMatrix =
            faceIntegrator_.integrate(face, [&](Base::PhysicalFace<DIM>& face) {
                LinearAlgebra::MiddleSizeMatrix result, temp;
                faceStiffnessMatrix1(face, indexing, stab, result);

                //                if (stab.fluxType2 == FluxType::IP) {
                //                    faceStiffnessMatrix3(face, temp,
                //                    stab.stab2); result += temp;
                //                }

                addScalarFaceMatrixTerms(face, indexing, stab, result);

                return result;
            });

        if (stab.hasFlux(FluxType::BREZZI)) {
            faceMatrix += brezziFluxBilinearTerm(it, stab);
        }
        (*it)->setFaceMatrix(faceMatrix, FACE_STIFFNESS_MATRIX_ID);

        if (boundaryCondition) {
            faceVector = faceIntegrator_.integrate(
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
typename DivDGMaxDiscretization<DIM>::Fields
    DivDGMaxDiscretization<DIM>::computeFields(
        const Base::Element* element,
        const Geometry::PointReference<DIM>& point,
        const LinearAlgebra::MiddleSizeVector& coefficients) const {
    Fields result;
    std::size_t nPhiU = element->getNumberOfBasisFunctions(0);
    std::size_t nPhiP = element->getNumberOfBasisFunctions(1);
    // Setup transformation to real space
    Base::PhysicalElement<DIM> physicalElement;
    std::shared_ptr<Base::CoordinateTransformation<DIM>>
        coordinateTransformationU{
            new Base::HCurlConformingTransformation<DIM>()};
    std::shared_ptr<Base::CoordinateTransformation<DIM>>
        coordinateTransformationP{new Base::H1ConformingTransformation<DIM>()};
    physicalElement.setTransformation(fieldTransform_, 0);
    physicalElement.setTransformation(potentialTransform_, 1);
    physicalElement.setElement(element);
    physicalElement.setPointReference(point);
    Geometry::PointPhysical<DIM> pointPhysical;
    pointPhysical = element->referenceToPhysical(point);

    // Actual value computation
    // Compute field part
    for (std::size_t i = 0; i < nPhiU; ++i) {
        LinearAlgebra::SmallVector<DIM> phiU;
        physicalElement.basisFunction(i, phiU, 0);
        result.realEField += std::real(coefficients[i]) * phiU;
        result.imagEField += std::imag(coefficients[i]) * phiU;
    }
    // Compute potential
    for (std::size_t i = 0; i < nPhiP; ++i) {
        result.potential +=
            coefficients[nPhiU + i] * physicalElement.basisFunction(i, 1);
    }
    return result;
}

template <std::size_t DIM>
LinearAlgebra::SmallVector<DIM> DivDGMaxDiscretization<DIM>::computeField(
    const Base::Element* element, const Geometry::PointReference<DIM>& point,
    const LinearAlgebra::MiddleSizeVector& coefficients) const {

    logger.log(Log::WARN, "Only computing the real part of the field.");
    Fields fields = computeFields(element, point, coefficients);
    return fields.realEField;
}

template <std::size_t DIM>
double DivDGMaxDiscretization<DIM>::computePotential(
    const Base::Element* element, const Geometry::PointReference<DIM>& point,
    const LinearAlgebra::MiddleSizeVector& coefficients) const {

    logger.log(Log::WARN, "Only computing the real part of the potential.");
    Fields fields = computeFields(element, point, coefficients);
    return fields.potential.real();
}

template <std::size_t DIM>
void DivDGMaxDiscretization<DIM>::computeElementMatrices(
    Base::Element* element, Utilities::ElementLocalIndexing& indexing) {

    double epsilon =
        static_cast<ElementInfos*>(element->getUserData())->epsilon_;

    LinearAlgebra::MiddleSizeMatrix massMatrix = elementIntegrator_.integrate(
        element, [&indexing, &epsilon](Base::PhysicalElement<DIM>& pelem) {
            std::size_t numDoFs = indexing.getNumberOfDoFs();
            std::size_t numUDoFs = indexing.getNumberOfDoFs(0);
            LinearAlgebra::MiddleSizeMatrix ret(numDoFs, numDoFs);

            LinearAlgebra::SmallVector<DIM> phi_i, phi_j;
            for (std::size_t i = 0; i < numUDoFs; ++i) {
                pelem.basisFunction(i, phi_i, 0);
                for (std::size_t j = i; j < numUDoFs; ++j) {
                    pelem.basisFunction(j, phi_j, 0);
                    double value = epsilon * (phi_i * phi_j);
                    ret(i, j) = value;
                    if (i != j) {
                        ret(j, i) = value;
                    }
                }
            }
            return ret;
        });
    element->setElementMatrix(massMatrix, ELEMENT_MASS_MATRIX_ID);

    LinearAlgebra::MiddleSizeMatrix stiffnessMatrix =
        elementIntegrator_.integrate(
            element, [&indexing, &epsilon](Base::PhysicalElement<DIM>& pelem) {
                std::size_t numDoFs = indexing.getNumberOfDoFs();
                std::size_t numUDoFs = indexing.getNumberOfDoFs(0);
                std::size_t numPDoFs = indexing.getNumberOfDoFs(1);
                std::size_t offPDoFs = indexing.getDoFOffset(1);
                LinearAlgebra::MiddleSizeMatrix ret(numDoFs, numDoFs);

                LinearAlgebra::SmallVector<DIM> phiI;
                for (std::size_t i = 0; i < numUDoFs; ++i) {
                    const LinearAlgebra::SmallVector<DIM>& phiIC =
                        pelem.basisFunctionCurl(i, 0);
                    for (std::size_t j = i; j < numUDoFs; ++j) {
                        double value = phiIC * pelem.basisFunctionCurl(j, 0);
                        ret(i, j) = value;
                        if (i != j) {
                            ret(j, i) = value;
                        }
                    }
                    pelem.basisFunction(i, phiI, 0);

                    for (std::size_t j = 0; j < numPDoFs; ++j) {
                        double value =
                            (phiI * pelem.basisFunctionDeriv(j, 1)) * -epsilon;
                        ret(i, j + offPDoFs) = value;
                        ret(j + offPDoFs, i) = value;
                    }
                }
                return ret;
            });
    element->setElementMatrix(stiffnessMatrix, ELEMENT_STIFFNESS_MATRIX_ID);
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
    std::size_t uDoFs = element->getNumberOfBasisFunctions(0);
    std::size_t pDoFs = element->getNumberOfBasisFunctions(1);
    ret.resize(uDoFs + pDoFs);

    LinearAlgebra::SmallVector<DIM> sourceValue, phi;
    sourceValue = source(el.getPointPhysical());
    for (std::size_t i = 0; i < (uDoFs); ++i) {
        el.basisFunction(i, phi, 0);
        ret(i) = phi * sourceValue;
    }
}

template <std::size_t DIM>
void DivDGMaxDiscretization<DIM>::faceStiffnessMatrix1(
    Base::PhysicalFace<DIM>& fa, const Utilities::FaceLocalIndexing& indexing,
    const Stab& stab, LinearAlgebra::MiddleSizeMatrix& ret) const {

    const Base::Face* face = fa.getFace();
    // For IP fluxes
    const double stab1 = stab.stab1 / face->getDiameter();

    // Mapping from basis function -> face matrix entry
    std::vector<std::size_t> mapping;
    indexing.getDoFMapping(0, mapping);

    const std::size_t totalDoFs = indexing.getNumberOfDoFs();
    const std::size_t totalUDoFs = mapping.size();

    ret.resize(totalDoFs, totalDoFs);

    // Averaging factor
    double factor = face->isInternal() ? -0.5 : -1;

    LinearAlgebra::SmallVector<DIM> phiUNormali, phiUNormalj;

    for (std::size_t i = 0; i < totalUDoFs; ++i) {
        const std::size_t& iIndex = mapping[i];
        fa.basisFunctionUnitNormalCross(i, phiUNormali, 0);
        const auto& phiUCurli = fa.basisFunctionCurl(i, 0);

        for (std::size_t j = i; j < totalUDoFs; ++j) {
            const std::size_t& jIndex = mapping[j];
            fa.basisFunctionUnitNormalCross(j, phiUNormalj, 0);
            const auto& phiUCurlj = fa.basisFunctionCurl(j, 0);

            double entry =
                factor * (phiUCurli * phiUNormalj + phiUNormali * phiUCurlj);

            if (stab.fluxType1 == FluxType::IP) {
                entry += stab1 * phiUNormali * phiUNormalj;
            }

            ret(iIndex, jIndex) = entry;
            ret(jIndex, iIndex) = entry;
        }
    }

    if (face->isInternal() && stab.fluxType2 == FluxType::IP &&
        stab.stab2 != 0.0) {
        // Note checking against 0.0 exactly, as this value is set by the user
        // and often disabled by setting it to 0.

        // Vector of epsilons
        FaceMaterialInfo minfo(fa);
        // Combination of both epsilon and the sign of the normal.
        std::vector<double> signedEpsilon;
        {
            signedEpsilon.resize(totalUDoFs);
            auto leftEnd = signedEpsilon.begin() +
                           indexing.getNumberOfDoFs(0, Base::Side::LEFT);
            std::fill(signedEpsilon.begin(), leftEnd, minfo.epsilonLeft);
            std::fill(leftEnd, signedEpsilon.end(), -minfo.epsilonRight);
        }

        const double stab2 = stab.stab2 * face->getDiameter() / minfo.epsMax;
        const LinearAlgebra::SmallVector<DIM>& normal =
            fa.getUnitNormalVector();
        LinearAlgebra::SmallVector<DIM> phi;

        for (std::size_t i = 0; i < totalUDoFs; ++i) {
            const std::size_t& iIndex = mapping[i];
            fa.basisFunction(i, phi, 0);
            double phiUi = phi * normal * signedEpsilon[i];

            for (std::size_t j = i; j < totalUDoFs; ++j) {
                const std::size_t& jIndex = mapping[j];

                fa.basisFunction(j, phi, 0);

                double value =
                    stab2 * phiUi * (phi * normal) * signedEpsilon[j];
                ret(iIndex, jIndex) += value;
                if (i != j) {
                    ret(jIndex, iIndex) += value;
                }
            }
        }
    }
}

template <std::size_t DIM>
void DivDGMaxDiscretization<DIM>::addScalarFaceMatrixTerms(
    Base::PhysicalFace<DIM>& fa, const Utilities::FaceLocalIndexing& indexing,
    const Stab& stab, LinearAlgebra::MiddleSizeMatrix& ret) const {

    // Mapping from basis function -> face matrix entry
    std::vector<std::size_t> mappingU, mappingP;
    indexing.getDoFMapping(0, mappingU);
    indexing.getDoFMapping(1, mappingP);

    const std::size_t totalDoFs = indexing.getNumberOfDoFs();
    const std::size_t totalUDoFs = mappingU.size();
    const std::size_t totalPDoFs = mappingP.size();

    ret.resize(totalDoFs, totalDoFs);

    // Averaging averageFactor
    const Base::Face* face = fa.getFace();
    double averageFactor = face->isInternal() ? 0.5 : 1;

    // Vector of epsilon, indexed by the U-variable
    FaceMaterialInfo minfo(fa);
    std::vector<double> epsilons;
    {
        epsilons.resize(totalUDoFs);
        auto leftEnd =
            epsilons.begin() + indexing.getNumberOfDoFs(0, Base::Side::LEFT);
        std::fill(epsilons.begin(), leftEnd, minfo.epsilonLeft);
        std::fill(leftEnd, epsilons.end(), minfo.epsilonRight);
    }

    // Factor for the normal, indexed by P-variable;
    std::vector<int> normalSign;
    {
        normalSign.resize(totalPDoFs);
        auto leftEnd =
            normalSign.begin() + indexing.getNumberOfDoFs(1, Base::Side::LEFT);
        std::fill(normalSign.begin(), leftEnd, 1);
        std::fill(leftEnd, normalSign.end(), -1);
    }

    const auto& normal = fa.getUnitNormalVector();
    LinearAlgebra::SmallVector<DIM> phiUi;

    /// Scalar vector coupling
    /// [[p]] {{eps v}} and it's symmetric term [[q]] {{eps u}}
    for (std::size_t i = 0; i < totalUDoFs; ++i) {
        const std::size_t iIndex = mappingU[i];
        fa.basisFunction(i, phiUi, 0);
        // Minor standard optimization, the discretization requires
        // {{epsilon phiUi}} . [[phiQj]]. This translates to:
        // averageFactor * epsilon(i) phiUi . n(j) phiQj, where epsilon(i) and
        // n(j) are the value of epsilon and the normal on the same side as the
        // basis function i and j respectively. As n(right) = -n(left) this is
        // transformed into
        // [epsilon(i) phiUi . n(left)] * [phiQj(j) s(j)]
        // where s(j) = 1.0 for the left values of j, and -1 for the right ones.
        double epsUNi = (phiUi * normal) * epsilons[i];

        for (std::size_t j = 0; j < totalPDoFs; ++j) {
            double phiPj = fa.basisFunction(j, 1) * normalSign[j];

            const double entry = averageFactor * epsUNi * phiPj;
            const std::size_t jIndex = mappingP[j];
            ret(iIndex, jIndex) += entry;
            ret(jIndex, iIndex) += entry;
        }
    }
    if (stab.fluxType3 == FluxType::IP) {
        /// Stabilization of the potential term
        /// stab3/diameter * epsMax [[p]] [[q]]
        const double stab3 = stab.stab3 / face->getDiameter() * minfo.epsMax;

        for (std::size_t i = 0; i < totalPDoFs; ++i) {
            const std::size_t iIndex = mappingP[i];
            double phiPi = fa.basisFunction(i, 1) * normalSign[i];

            for (std::size_t j = i; j < totalPDoFs; ++j) {
                const std::size_t jIndex = mappingP[j];

                const double value =
                    stab3 * phiPi * fa.basisFunction(j, 1) * normalSign[j];
                // Negative contribution
                ret(iIndex, jIndex) -= value;
                if (i != j) {
                    ret(jIndex, iIndex) -= value;
                }
            }
        }
    }
}


template <std::size_t DIM>
LinearAlgebra::MiddleSizeMatrix
    DivDGMaxDiscretization<DIM>::computeFaceVectorMassMatrix(
        typename Base::MeshManipulator<DIM>::FaceIterator rawFace) {
    FaceDoFInfo faceInfo = getFaceDoFInfo(**rawFace);

    LinearAlgebra::MiddleSizeMatrix result = elementIntegrator_.integrate(
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
        result += elementIntegrator_.integrate(
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
        typename Base::MeshManipulator<DIM>::FaceIterator rawFace) {
    FaceDoFInfo faceInfo = getFaceDoFInfo(**rawFace);

    LinearAlgebra::MiddleSizeMatrix result = elementIntegrator_.integrate(
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
        result += elementIntegrator_.integrate(
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
        typename Base::MeshManipulator<DIM>::FaceIterator rawFace) {
    FaceDoFInfo faceInfo = getFaceDoFInfo(**rawFace);

    return faceIntegrator_.integrate(
        *rawFace, [&](Base::PhysicalFace<DIM>& face) {
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
        typename Base::MeshManipulator<DIM>::FaceIterator rawFace) {

    FaceDoFInfo faceInfo = getFaceDoFInfo(**rawFace);

    return faceIntegrator_.integrate(
        *rawFace, [&](Base::PhysicalFace<DIM>& face) {
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
        typename Base::MeshManipulator<DIM>::FaceIterator rawFace) {
    FaceDoFInfo faceInfo = getFaceDoFInfo(**rawFace);

    double epsilonLeft = static_cast<ElementInfos*>(
                             (*rawFace)->getPtrElementLeft()->getUserData())
                             ->epsilon_;
    double epsilonRight =
        faceInfo.internal ? static_cast<ElementInfos*>(
                                (*rawFace)->getPtrElementRight()->getUserData())
                                ->epsilon_
                          : 0.0;  // Should not be used.

    return faceIntegrator_.integrate(
        *rawFace, [&](Base::PhysicalFace<DIM>& face) {
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
        typename Base::MeshManipulator<DIM>::FaceIterator rawFace, Stab stab) {
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

        LinearAlgebra::SmallVector<DIM> val, phi_curl;
        LinearAlgebra::SmallVector<DIM> phi;
        val = boundaryValue(fa);

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
        DivDGMaxDiscretization<DIM>::Stab stab) {
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
            LinearAlgebra::SmallVector<DIM> val = boundaryValue(face);

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
    const DivDGMaxDiscretization<DIM>::InputFunction& electricField) {

    double error = 0;
    auto end = mesh.elementColEnd();
    for (typename Base::MeshManipulator<DIM>::ElementIterator it =
             mesh.elementColBegin();
         it != end; ++it) {
        error += elementIntegrator_.integrate(
            (*it), [&](Base::PhysicalElement<DIM>& el) {
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

    std::size_t numberOfUDoFs = element->getNumberOfBasisFunctions(0),
                numberOfPDoFs = element->getNumberOfBasisFunctions(1);
    LinearAlgebra::MiddleSizeVector data =
        element->getTimeIntegrationVector(timeVector);

    LinearAlgebra::SmallVector<DIM> error, phi;
    error = exactValues(el.getPointPhysical());
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

char fluxName(typename DivDGMaxDiscretizationBase::FluxType f) {
    switch (f) {
        case DivDGMaxDiscretizationBase::FluxType::BREZZI:
            return 'b';
        case DivDGMaxDiscretizationBase::FluxType::IP:
            return 'i';
        default:
            logger.assert_always(false, "Unknown flux type.");
            return '0';
    }
}

std::ostream& operator<<(std::ostream& os,
                         const DivDGMaxDiscretizationBase::Stab& stab) {
    os << "Stab{" << fluxName(stab.fluxType1) << "=" << stab.stab1 << ", "
       << fluxName(stab.fluxType2) << "=" << stab.stab2 << ", "
       << fluxName(stab.fluxType3) << "=" << stab.stab3 << "}";
    return os;
}

template class DivDGMaxDiscretization<2>;
template class DivDGMaxDiscretization<3>;

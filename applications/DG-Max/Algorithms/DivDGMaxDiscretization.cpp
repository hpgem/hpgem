/*
This file forms part of hpGEM. This package has been developed over a number of
years by various people at the University of Twente and a full list of
contributors can be found at http://hpgem.org/about-the-code/team

This code is distributed using BSD 3-Clause License. A copy of which can found
below.


Copyright (c) 2018, University of Twente
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
DivDGMaxDiscretization<DIM>::DivDGMaxDiscretization(std::size_t order,
                                                    Stab stab)
    : order_(order),
      stab_(stab),
      fieldTransform_(
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
    Base::MeshManipulator<DIM>& mesh) const {
    // We would like to configure the number of unknowns here, but this is
    // unfortunately not possible, as it is configured at the creation of
    // the mesh. The best we can do is check if it is configured correctly.
    std::size_t unknowns = mesh.getConfigData()->numberOfUnknowns_;
    logger.assert_always(unknowns == 2, "DivDGMax expects 2 unknowns but got %",
                         unknowns);
    // TODO: This needs the additional unknown id.
    mesh.useNedelecDGBasisFunctions(order_);
    mesh.useDefaultDGBasisFunctions(order_, 1);
}

template <std::size_t DIM>
void DivDGMaxDiscretization<DIM>::computeElementIntegralsImpl(
    Base::MeshManipulator<DIM>& mesh,
    const std::map<std::size_t, InputFunction>& elementVectors,
    LocalIntegrals integrals) {

    Utilities::ElementLocalIndexing indexing;
    indexing.reinit({0, 1});

    // Out of efficiency, elementColEnd may be expensive
    auto end = mesh.elementColEnd();
    for (auto it = mesh.elementColBegin(); it != end; ++it) {
        Base::Element* element = *it;
        indexing.reinit(element);
        if (integrals == LocalIntegrals::ALL) {
            computeElementMatrices(element, indexing);
        }

        for (const auto& elementVec : elementVectors) {
            LinearAlgebra::MiddleSizeVector vec;
            if (!elementVec.second) {
                vec.resize(indexing.getNumberOfDoFs());
            } else {
                vec = elementIntegrator_.integrate(
                    (*it), [&](Base::PhysicalElement<DIM>& element) {
                        LinearAlgebra::MiddleSizeVector result;
                        elementSourceVector(element, elementVec.second, result);
                        return result;
                    });
            }
            element->setElementVector(vec, elementVec.first);
        }
    }
}

template <std::size_t DIM>
void DivDGMaxDiscretization<DIM>::computeFaceIntegralsImpl(
    Base::MeshManipulator<DIM>& mesh,
    const std::map<std::size_t, FaceInputFunction>& boundaryVectors,
    DGMax::BoundaryConditionIndicator boundaryIndicator,
    LocalIntegrals integrals) {
    LinearAlgebra::MiddleSizeMatrix faceMatrix;

    Utilities::FaceLocalIndexing indexing;
    indexing.reinit({0, 1});

    // For efficiency, faceColEnd may be expensive
    auto end = mesh.faceColEnd();
    for (auto it = mesh.faceColBegin(); it != end; ++it) {

        Base::Face* face = *it;
        indexing.reinit(face);

        std::size_t totalDoFs = indexing.getNumberOfDoFs();
        using BCT = DGMax::BoundaryConditionType;
        BCT bct = BCT::INTERNAL;
        if (!face->isInternal()) {
            bct = boundaryIndicator(*face);
        }
        if (integrals == LocalIntegrals::ALL) {
            faceMatrix = faceIntegrator_.integrate(
                face, [&indexing, &bct, this](Base::PhysicalFace<DIM>& face) {
                    std::size_t numDoFs = indexing.getNumberOfDoFs();
                    LinearAlgebra::MiddleSizeMatrix result(numDoFs, numDoFs);
                    if (bct == BCT::INTERNAL || !isNaturalBoundary(bct)) {
                        faceStiffnessMatrixFieldIntegrand(face, indexing,
                                                          result);
                    }
                    addFaceMatrixPotentialIntegrand(face, indexing, bct,
                                                    result);
                    return result;
                });
            if (stab_.hasFlux(FluxType::BREZZI)) {
                faceMatrix += brezziFluxBilinearTerm(face, bct);
            }

//            if (bct != BCT::INTERNAL) {
//                faceMatrix *= 0.0;
//            }
            face->setFaceMatrix(faceMatrix, FACE_STIFFNESS_MATRIX_ID);

            if (bct == DGMax::BoundaryConditionType::SILVER_MULLER) {
                faceMatrix = faceIntegrator_.integrate(
                    face,
                    [&indexing, &bct, this](Base::PhysicalFace<DIM>& face) {
                        return computeFaceImpedanceIntegrand(face, indexing,
                                                             bct);
                    });
            } else {
                // Matrix is already of the right size, just zero it
                faceMatrix *= 0.0;
            }
            face->setFaceMatrix(faceMatrix, FACE_IMPEDANCE_MATRIX_ID);
        }

        for (const auto& boundaryVec : boundaryVectors) {
            LinearAlgebra::MiddleSizeVector vec;
            if (!boundaryVec.second) {
                vec.resize(totalDoFs);
            } else {
                vec = faceIntegrator_.integrate(
                    face, [&](Base::PhysicalFace<DIM>& face) {
                        LinearAlgebra::MiddleSizeVector result;
                        faceBoundaryVector(face, boundaryVec.second, result,
                                           bct);
                        return result;
                    });
                if (stab_.hasFlux(FluxType::BREZZI) &&
                    bct == DGMax::BoundaryConditionType::DIRICHLET) {
                    vec += brezziFluxBoundaryVector(face, boundaryVec.second);
                }
            }
            face->setFaceVector(vec, boundaryVec.first);
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

    ElementInfos& userData = ElementInfos::get(*element);
    result.material = userData.getMaterial();

    // Actual value computation
    // Compute field part
    for (std::size_t i = 0; i < nPhiU; ++i) {
        LinearAlgebra::SmallVector<DIM> phiU;
        physicalElement.basisFunction(i, phiU, 0);
        result.electricField += coefficients[i] * phiU;
        result.electricFieldCurl +=
            coefficients[i] * physicalElement.basisFunctionCurl(i, 0);
    }
    // Compute potential
    for (std::size_t i = 0; i < nPhiP; ++i) {
        result.potential +=
            coefficients[nPhiU + i] * physicalElement.basisFunction(i, 1);
    }

    // Rescale for PMLs
    const auto& pPhys = physicalElement.getPointPhysical();
    result.electricField =
        userData.getFieldRescaling(pPhys).applyDiv(result.electricField);
    result.electricFieldCurl = userData.getCurlFieldRescaling(pPhys).applyCurl(
        result.electricFieldCurl);

    return result;
}

template <std::size_t DIM>
LinearAlgebra::SmallVectorC<DIM> DivDGMaxDiscretization<DIM>::computeField(
    const Base::Element* element, const Geometry::PointReference<DIM>& point,
    const LinearAlgebra::MiddleSizeVector& coefficients) const {

    Fields fields = computeFields(element, point, coefficients);
    return fields.electricField;
}

template <std::size_t DIM>
LinearAlgebra::SmallVectorC<DIM> DivDGMaxDiscretization<DIM>::computeCurlField(
    const Base::Element* element, const PointReferenceT& point,
    const LinearAlgebra::MiddleSizeVector& coefficients) const {
    Fields fields = computeFields(element, point, coefficients);
    return fields.electricFieldCurl;
}

template <std::size_t DIM>
std::complex<double> DivDGMaxDiscretization<DIM>::computePotential(
    const Base::Element* element, const Geometry::PointReference<DIM>& point,
    const LinearAlgebra::MiddleSizeVector& coefficients) const {

    Fields fields = computeFields(element, point, coefficients);
    return fields.potential;
}

template <std::size_t DIM>
void DivDGMaxDiscretization<DIM>::writeFields(
    Output::VTKSpecificTimeWriter<DIM>& output,
    std::size_t timeIntegrationVectorId) const {
    using VecR = LinearAlgebra::SmallVector<DIM>;
    std::map<std::string, std::function<double(Fields&)>> scalars;
    std::map<std::string, std::function<VecR(Fields&)>> vectors;

    // 4 parts of the field
    vectors["E-real"] = [](Fields& fields) {
        return fields.electricField.real();
    };
    vectors["E-imag"] = [](Fields& fields) {
        return fields.electricField.imag();
    };
    scalars["p-real"] = [](Fields& fields) { return fields.potential.real(); };
    scalars["p-imag"] = [](Fields& fields) { return fields.potential.imag(); };

    // Derived quantities
    scalars["Emag"] = [](Fields& fields) {
        return fields.electricField.l2Norm();
    };
    vectors["S-kappa-real"] = [](Fields& fields) {
        // S = 1/2 Re(E x H^*)
        //   = -1/(2 omega mu) Im(E x Curl E)
        // Using i omega mu H = Curl E
        return -0.5 *
               LinearAlgebra::leftDoubledCrossProduct(
                   fields.electricField, fields.electricFieldCurl.conj())
                   .imag() /
               fields.material.getPermeability();
    };
    scalars["Energy"] = [](Fields& fields) {
        // u = 1/2(epsilon |E|^2 + mu |H|^2)
        //   = epsilon |E|^2 (via curl-curl relation)
        return fields.material.getPermittivity() *
               fields.electricField.l2NormSquared();
    };

    output.template writeMultiple<Fields>(
        [this, timeIntegrationVectorId](
            Base::Element* element, const Geometry::PointReference<DIM>& point,
            std::size_t) {
            LinearAlgebra::MiddleSizeVector coefficients =
                element->getTimeIntegrationVector(timeIntegrationVectorId);
            return computeFields(element, point, coefficients);
        },
        scalars, vectors);
    // For various use
    output.write(
        [](Base::Element* element, const Geometry::PointReference<DIM>&,
           std::size_t) {
            auto* userData = element->getUserData();
            const ElementInfos* elementInfo =
                dynamic_cast<ElementInfos*>(userData);
            if (elementInfo != nullptr) {
                return elementInfo->getPermittivity();
            } else {
                return -1.0;  // Clearly invalid value
            }
        },
        "epsilon");
}

template <std::size_t DIM>
LinearAlgebra::SmallVector<4> DivDGMaxDiscretization<DIM>::computeEnergyFluxes(
    Base::Face& face, hpgem::Base::Side side, double wavenumber,
    std::size_t timeIntegrationVectorId,
    const DGMax::FieldPattern<DIM>* background) {

    using VecC = LinearAlgebra::SmallVectorC<DIM>;

    const auto& coefficients =
        face.getTimeIntegrationVector(timeIntegrationVectorId);
    auto dofInfo = getFaceDoFInfo(face);

    double factor = face.isInternal() ? 0.5 : 1.0;

    auto flux = faceIntegrator_.integrate(
        &face, [&coefficients, &dofInfo, &factor,
                background](Base::PhysicalFace<DIM>& pface) {
            LinearAlgebra::SmallVector<4> result;
            VecC avgCurl;
            VecC avgField;
            LinearAlgebra::SmallVector<DIM> phi;
            VecC normal = pface.getUnitNormalVector();
            std::size_t leftDoFs = dofInfo.leftUDoFs + dofInfo.leftPDoFs;
            for (std::size_t i = 0; i < dofInfo.totalUDoFs(); ++i) {
                const auto& coefficient = i < dofInfo.leftUDoFs
                                              ? coefficients[i]
                                              : coefficients[i + leftDoFs];
                avgCurl += factor * coefficient * pface.basisFunctionCurl(i, 0);
                pface.basisFunction(i, phi, 0);
                avgField += factor * coefficient * normal.crossProduct(phi);
            }
            result[0] = (avgField * avgCurl).imag();
            if (background) {
                VecC backgroundE = normal.crossProduct(
                    background->field(pface.getPointPhysical()));
                VecC totalE = avgField + backgroundE;
                VecC backgroundCurl =
                    background->fieldCurl(pface.getPointPhysical());
                VecC totalCurl = avgCurl + backgroundCurl;
                result[1] = (backgroundE * backgroundCurl).imag();
                result[2] =
                    (backgroundE * avgCurl + avgField * backgroundCurl).imag();
                result[3] = (totalE * totalCurl).imag();
            }
            return result;
        });
    if (side == hpgem::Base::Side::RIGHT) {
        flux *= -1.0;
    }
    auto infos =
        dynamic_cast<ElementInfos*>(face.getPtrElement(side)->getUserData());
    logger.assert_debug(infos != nullptr, "No material information");
    return flux / (wavenumber * infos->getPermeability());
}

template <std::size_t DIM>
double DivDGMaxDiscretization<DIM>::computeFieldL2Integral(
    Base::Face& face, Base::Side side, std::size_t vector_id) {
    const LinearAlgebra::MiddleSizeVector& coefficients =
        face.getTimeIntegrationVector(vector_id);
    using VecC = LinearAlgebra::SmallVectorC<DIM>;
    std::size_t dofOffset =
        side == Base::Side::LEFT
            ? 0
            : face.getPtrElementLeft()->getTotalNumberOfBasisFunctions();
    std::size_t numDoFs =
        side == Base::Side::LEFT
            ? face.getPtrElementLeft()->getNumberOfBasisFunctions(0)
            : face.getPtrElementRight()->getNumberOfBasisFunctions(0);

    auto& elementInfo = ElementInfos::get(*face.getPtrElement(side));
    return faceIntegrator_.integrate(
        &face, [&coefficients, dofOffset, numDoFs, side,
                &elementInfo](Base::PhysicalFace<DIM>& pface) {
            VecC field;
            LinearAlgebra::SmallVector<DIM> phi;
            for (std::size_t i = 0; i < numDoFs; ++i) {
                pface.basisFunction(i + dofOffset, phi, 0);
                field += coefficients[i + dofOffset] * phi;
            }
            // Rescale fields for PMLs
            field = elementInfo
                        .getFieldRescaling(
                            pface.getPhysicalElement(side).getPointPhysical())
                        .applyDiv(field);
            return field.l2NormSquared();
        });
}

template <std::size_t DIM>
void DivDGMaxDiscretization<DIM>::computeElementMatrices(
    Base::Element* element, Utilities::ElementLocalIndexing& indexing) {

    const auto& material = ElementInfos::get(*element);

    LinearAlgebra::MiddleSizeMatrix massMatrix = elementIntegrator_.integrate(
        element, [&indexing, &material](Base::PhysicalElement<DIM>& pelem) {
            std::size_t numDoFs = indexing.getNumberOfDoFs();
            std::size_t numUDoFs = indexing.getNumberOfDoFs(0);
            LinearAlgebra::MiddleSizeMatrix ret(numDoFs, numDoFs);

            const auto materialDiv =
                material.getMaterialConstantDiv(pelem.getPointPhysical());

            LinearAlgebra::SmallVector<DIM> phi_i, phi_j;
            for (std::size_t i = 0; i < numUDoFs; ++i) {
                pelem.basisFunction(i, phi_i, 0);
                for (std::size_t j = 0; j < numUDoFs; ++j) {
                    pelem.basisFunction(j, phi_j, 0);
                    std::complex<double> value =
                        materialDiv.applyDiv(phi_j) * phi_i;

                    ret(i, j) = value;
                }
            }
            return ret;
        });
    element->setElementMatrix(massMatrix, MASS_MATRIX_ID);

    LinearAlgebra::MiddleSizeMatrix stiffnessMatrix =
        elementIntegrator_.integrate(
            element, [&indexing, &material](Base::PhysicalElement<DIM>& pelem) {
                std::size_t numDoFs = indexing.getNumberOfDoFs();
                std::size_t numUDoFs = indexing.getNumberOfDoFs(0);
                std::size_t numPDoFs = indexing.getNumberOfDoFs(1);
                std::size_t offPDoFs = indexing.getDoFOffset(1);
                LinearAlgebra::MiddleSizeMatrix ret(numDoFs, numDoFs);

                const auto& point = pelem.getPointPhysical();
                const auto materialCurl =
                    material.getMaterialConstantCurl(point);
                const auto materialDiv = material.getMaterialConstantDiv(point);

                LinearAlgebra::SmallVector<DIM> phiI;
                for (std::size_t i = 0; i < numUDoFs; ++i) {
                    LinearAlgebra::SmallVectorC<DIM> phiIC =
                        materialCurl.applyCurl(pelem.basisFunctionCurl(i, 0));
                    for (std::size_t j = 0; j < numUDoFs; ++j) {
                        std::complex<double> value =
                            phiIC * pelem.basisFunctionCurl(j, 0);
                        ret(i, j) = value;
                    }
                    pelem.basisFunction(i, phiI, 0);

                    for (std::size_t j = 0; j < numPDoFs; ++j) {
                        // The term from Maxwell is Div (epsilon E) = 0, thus
                        // for the Lagrange multiplier we have the adjoint
                        // operator of Div (epsilon .) = epsilon^H grad
                        // The (Div epsilon E, q) in the weak formulation is
                        // reworked into - (epsilon E, grad q), to ensure
                        // symmetry. We compute here -(epsilon E, grad q)
                        std::complex<double> value =
                            -materialDiv.applyDiv(phiI) *
                            pelem.basisFunctionDeriv(j, 1);
                        ret(j + offPDoFs, i) = value;
                        // -(grad p, epsilon* v) is just the conjugate value
                        ret(i, j + offPDoFs) = std::conj(value);
                    }
                }
                return ret;
            });
    element->setElementMatrix(stiffnessMatrix, STIFFNESS_MATRIX_ID);
}

template <std::size_t DIM>
void DivDGMaxDiscretization<DIM>::elementSourceVector(
    Base::PhysicalElement<DIM>& el, const InputFunction& source,
    LinearAlgebra::MiddleSizeVector& ret) const {
    const Base::Element* element = el.getElement();
    std::size_t uDoFs = element->getNumberOfBasisFunctions(0);
    std::size_t pDoFs = element->getNumberOfBasisFunctions(1);
    ret.resize(uDoFs + pDoFs);

    LinearAlgebra::SmallVectorC<DIM> sourceValue;
    LinearAlgebra::SmallVector<DIM> phi;
    sourceValue = source(*element, el.getPointPhysical());
    for (std::size_t i = 0; i < (uDoFs); ++i) {
        el.basisFunction(i, phi, 0);
        ret(i) = sourceValue * phi;
    }
}

template <std::size_t DIM>
void DivDGMaxDiscretization<DIM>::faceStiffnessMatrixFieldIntegrand(
    Base::PhysicalFace<DIM>& fa, const Utilities::FaceLocalIndexing& indexing,
    LinearAlgebra::MiddleSizeMatrix& ret) const {

    const Base::Face* face = fa.getFace();
    const auto& materialLeft = ElementInfos::get(*face->getPtrElementLeft());
    const auto* materialRight = ElementInfos::get(face->getPtrElementRight());

    const auto& p = fa.getPointPhysical();
    const auto materialCurlLeft = materialLeft.getMaterialConstantCurl(p);
    const auto materialDivLeft = materialLeft.getMaterialConstantDiv(p);
    const auto materialCurlRight =
        materialRight != nullptr ? materialRight->getMaterialConstantCurl(p)
                                 : DGMax::MaterialTensor();
    const auto materialDivRight = materialRight != nullptr
                                      ? materialRight->getMaterialConstantDiv(p)
                                      : DGMax::MaterialTensor();

    // For IP fluxes
    const double stab1 = stab_.stab1 / face->getDiameter();

    // Mapping from basis function -> face matrix entry
    std::vector<std::size_t> mapping;
    indexing.getDoFMapping(0, mapping);

    const std::size_t totalDoFs = indexing.getNumberOfDoFs();
    const std::size_t totalUDoFs = mapping.size();
    const std::size_t leftUDoFs = indexing.getNumberOfDoFs(0, Base::Side::LEFT);

    ret.resize(totalDoFs, totalDoFs);

    // Averaging factor
    double factor = face->isInternal() ? -0.5 : -1;

    LinearAlgebra::SmallVector<DIM> phiUNormali, phiUNormalj;
    LinearAlgebra::SmallVectorC<DIM> phiUCurli, phiUCurlj;

    for (std::size_t i = 0; i < totalUDoFs; ++i) {
        const std::size_t& iIndex = mapping[i];
        fa.basisFunctionUnitNormalCross(i, phiUNormali, 0);
        phiUCurli = ((i < leftUDoFs) ? materialCurlLeft : materialCurlRight)
                        .applyCurl(fa.basisFunctionCurl(i, 0));

        for (std::size_t j = 0; j < totalUDoFs; ++j) {
            const std::size_t& jIndex = mapping[j];
            fa.basisFunctionUnitNormalCross(j, phiUNormalj, 0);
            const DGMax::MaterialTensor material =
                (j < leftUDoFs) ? materialCurlLeft.adjoint()
                                : materialCurlRight.adjoint();
            phiUCurlj = material.applyCurl(fa.basisFunctionCurl(j, 0));

            // Note the ordering is important to get the conjugation of complex
            // material parameters correctly. The i=trial, j=test function and
            // the second argument of the inner product is conjugated.
            std::complex<double> entry =
                factor * (phiUCurli * phiUNormalj + phiUNormali * phiUCurlj);

            if (stab_.fluxType1 == FluxType::IP) {
                entry += stab1 * phiUNormali * phiUNormalj;
            }

            ret(jIndex, iIndex) = entry;
        }
    }

    if (face->isInternal() && stab_.fluxType2 == FluxType::IP &&
        stab_.stab2 != 0.0) {
        // Note checking against 0.0 exactly, as this value is set by the user
        // and often disabled by setting it to 0.

        // TODO: Think about the max(epsLeft, epsRight)^{-1} that was used to
        // make this independent of epsilon scaling.
        const double stab2 = stab_.stab2 * face->getDiameter();
        const LinearAlgebra::SmallVector<DIM>& normal =
            fa.getUnitNormalVector();
        LinearAlgebra::SmallVector<DIM> phi;

        for (std::size_t i = 0; i < totalUDoFs; ++i) {
            const std::size_t& iIndex = mapping[i];
            fa.basisFunction(i, phi, 0);
            const auto& materiali =
                (i < leftUDoFs ? materialDivLeft : materialDivRight);
            std::complex<double> phiUi = materiali.applyDiv(phi) * normal;
            if (i >= leftUDoFs) {
                phiUi *= -1.0;
            }

            for (std::size_t j = i; j < totalUDoFs; ++j) {
                const std::size_t& jIndex = mapping[j];

                fa.basisFunction(j, phi, 0);
                const auto materialj =
                    j < leftUDoFs ? materialDivLeft : materialDivRight;
                std::complex<double> phiUj = materialj.applyDiv(phi) * normal;
                if (j >= leftUDoFs) {
                    phiUj *= -1.0;
                }

                // Again first argument = trial function => not conjugated
                // Test functions are conjugated as an inner product is
                // computed. This is not done by * for complex<double> so we do
                // it manually.
                phiUj = std::conj(phiUi);
                std::complex<double> value = stab2 * phiUi * phiUj;
                ret(jIndex, iIndex) += value;
                if (i != j) {
                    ret(iIndex, jIndex) += std::conj(value);
                }
            }
        }
    }
}

template <std::size_t DIM>
void DivDGMaxDiscretization<DIM>::addFaceMatrixPotentialIntegrand(
    Base::PhysicalFace<DIM>& fa, const Utilities::FaceLocalIndexing& indexing,
    DGMax::BoundaryConditionType bct,
    LinearAlgebra::MiddleSizeMatrix& ret) const {

    // Mapping from basis function -> face matrix entry
    std::vector<std::size_t> mappingU, mappingP;
    indexing.getDoFMapping(0, mappingU);
    indexing.getDoFMapping(1, mappingP);

    const std::size_t totalDoFs = indexing.getNumberOfDoFs();
    const std::size_t totalUDoFs = mappingU.size();
    const std::size_t totalPDoFs = mappingP.size();
    const std::size_t leftUDoFs =
        indexing.getNumberOfDoFs(0, hpgem::Base::Side::LEFT);

    ret.resize(totalDoFs, totalDoFs);

    // Averaging averageFactor
    const Base::Face* face = fa.getFace();
    double averageFactor = face->isInternal() ? 0.5 : 1;

    // Material constant for the U-variable.
    const auto& point = fa.getPointPhysical();
    const auto materialDivLeft = ElementInfos::get(*face->getPtrElementLeft())
                                     .getMaterialConstantDiv(point);
    const auto materialDivRight =
        face->isInternal() ? ElementInfos::get(*face->getPtrElementRight())
                                 .getMaterialConstantDiv(point)
                           : DGMax::MaterialTensor();

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

        std::complex<double> epsUNi =
            (i < leftUDoFs ? materialDivLeft : materialDivRight)
                .applyDiv(phiUi) *
            normal;

        for (std::size_t j = 0; j < totalPDoFs; ++j) {
            double phiPj = fa.basisFunction(j, 1) * normalSign[j];

            // Note i = trial function -> not conjugated
            const std::complex<double> entry = averageFactor * epsUNi * phiPj;
            const std::size_t jIndex = mappingP[j];

            ret(jIndex, iIndex) += entry;
            // Hermitian counterpart => conjugate
            ret(iIndex, jIndex) += std::conj(entry);
        }
    }
    if (stab_.fluxType3 == FluxType::IP) {
        /// Stabilization of the potential term
        /// stab3/diameter * epsMax [[p]] [[q]]
        // Note: epsMax removed due to complex tensors.
        const double stab3 = stab_.stab3 / face->getDiameter();

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
    DivDGMaxDiscretization<DIM>::computeFaceImpedanceIntegrand(
        Base::PhysicalFace<DIM>& face, Utilities::FaceLocalIndexing& indexing,
        DGMax::BoundaryConditionType& bct) const {
    std::size_t nUDoFs = indexing.getNumberOfDoFs(0, Base::Side::LEFT);
    std::size_t nDoFs = indexing.getNumberOfDoFs();
    LinearAlgebra::MiddleSizeMatrix result(nDoFs, nDoFs);

    // Just to make sure
    logger.assert_debug(bct == DGMax::BoundaryConditionType::SILVER_MULLER,
                        "Function only written for SilverMuller type faces");

    std::complex<double> impedance;
    {
        auto* material = dynamic_cast<ElementInfos*>(
            face.getFace()->getPtrElementLeft()->getUserData());
        logger.assert_debug(material != nullptr, "No material information");
        impedance = std::complex<double>(0, material->getImpedance());
    }
    for (std::size_t i = 0; i < nUDoFs; ++i) {
        LinearAlgebra::SmallVector<DIM> phiUNi;
        face.basisFunctionUnitNormalCross(i, phiUNi, 0);
        for (std::size_t j = 0; j < nUDoFs; ++j) {
            LinearAlgebra::SmallVector<DIM> phiUNj;
            face.basisFunctionUnitNormalCross(j, phiUNj, 0);
            result(j, i) = -impedance * (phiUNi * phiUNj);
        }
    }
    return result;
}

template <std::size_t DIM>
LinearAlgebra::MiddleSizeMatrix
    DivDGMaxDiscretization<DIM>::computeFaceVectorMassMatrix(Base::Face* face) {
    FaceDoFInfo faceInfo = getFaceDoFInfo(*face);

    LinearAlgebra::MiddleSizeMatrix result = elementIntegrator_.integrate(
        face->getPtrElementLeft(), [&](Base::PhysicalElement<DIM>& element) {
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
            face->getPtrElementRight(),
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
    DivDGMaxDiscretization<DIM>::computeFaceScalarMassMatrix(Base::Face* face) {
    FaceDoFInfo faceInfo = getFaceDoFInfo(*face);

    LinearAlgebra::MiddleSizeMatrix result = elementIntegrator_.integrate(
        face->getPtrElementLeft(), [&](Base::PhysicalElement<DIM>& element) {
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
            face->getPtrElementRight(),
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
    DivDGMaxDiscretization<DIM>::computeScalarLiftProjector(Base::Face* face) {
    FaceDoFInfo faceInfo = getFaceDoFInfo(*face);

    return faceIntegrator_.integrate(face, [&](Base::PhysicalFace<DIM>& face) {
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
    DivDGMaxDiscretization<DIM>::computeVectorLiftProjector(Base::Face* face) {

    FaceDoFInfo faceInfo = getFaceDoFInfo(*face);

    return faceIntegrator_.integrate(face, [&](Base::PhysicalFace<DIM>& face) {
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
        Base::Face* face) {
    FaceDoFInfo faceInfo = getFaceDoFInfo(*face);

    const auto& materialLeft = ElementInfos::get(*face->getPtrElementLeft());
    const auto* materialRight = ElementInfos::get(face->getPtrElementRight());

    return faceIntegrator_.integrate(face, [&](Base::PhysicalFace<DIM>& face) {
        LinearAlgebra::MiddleSizeMatrix result;
        result.resize(faceInfo.totalPDoFs(), faceInfo.totalUDoFs());
        LinearAlgebra::SmallVector<DIM> basisU, normal;
        LinearAlgebra::SmallVectorC<DIM> matBasisU;

        const auto& point = face.getPointPhysical();
        const auto materialDivLeft = materialLeft.getMaterialConstantDiv(point);
        const auto materialDivRight =
            materialRight != nullptr
                ? materialRight->getMaterialConstantCurl(point)
                : DGMax::MaterialTensor();

        normal = face.getUnitNormalVector();
        double factor = faceInfo.internal ? 0.5 : 1.0;  // From the average
        for (std::size_t i = 0; i < faceInfo.totalPDoFs(); ++i) {
            double basisP = face.basisFunction(i, 1);
            for (std::size_t j = 0; j < faceInfo.totalUDoFs(); ++j) {
                face.basisFunction(j, basisU, 0);

                // - form the change in normal
                if (j < faceInfo.leftUDoFs) {
                    matBasisU = materialDivLeft.applyDiv(basisU);
                } else {
                    matBasisU = -materialDivRight.applyDiv(basisU);
                }

                std::complex<double> uvalue =
                    matBasisU * normal;  // [[epsilon u]]_N
                result(i, j) = factor * uvalue * basisP;
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
        Base::Face* face, DGMax::BoundaryConditionType bct) {
    FaceDoFInfo faceInfo = getFaceDoFInfo(*face);
    LinearAlgebra::MiddleSizeMatrix result(faceInfo.totalDoFs(),
                                           faceInfo.totalDoFs());

    bool isNaturalBCT = !face->isInternal() && DGMax::isNaturalBoundary(bct);

    if (isNaturalBCT && stab_.fluxType3 != FluxType::BREZZI) {
        // On a boundaries with natural boundary conditions the first two
        // stabilization terms (fluxType1,2) are not used. So prevent doing any
        // work if the last stabilization term (fluxType3) uses a different
        // flux.
        std::size_t dofs = faceInfo.totalDoFs();
        return LinearAlgebra::MiddleSizeMatrix(dofs, dofs);
    }

    // Mass matrix for vector valued components
    LinearAlgebra::MiddleSizeMatrix massMat = computeFaceVectorMassMatrix(face);
    massMat.cholesky();

    if (stab_.fluxType1 == FluxType::BREZZI && !isNaturalBCT) {
        std::size_t leftFaces = face->getPtrElementLeft()->getNumberOfFaces();
        std::size_t rightFaces =
            faceInfo.internal ? face->getPtrElementLeft()->getNumberOfFaces()
                              : 0;

        // Compute vector stabilizer, R^T L^{-T} (stab1 + numFaces * muinv)
        // L^{-1} R with R the lifting matrix for [[u_j]]_T {{u_i}} (stab1 +
        // numFace*muinv) the diagonal matrix with numFaces the number of faces
        // of the element (left on top diagonal, right bottom diagonal) L the
        // lower triangular part of the cholesky decomposition of M.
        LinearAlgebra::MiddleSizeMatrix vectorStabilizer =
            computeVectorLiftProjector(face);
        // massMat.solve(vectorStabilizer); // M^{-1}R
        massMat.solveLowerTriangular(vectorStabilizer);  // L^{-1}R, with LL^{T}
                                                         // = M
        // Inplace multiplication with diagonal matrix
        double lfactor = std::sqrt(stab_.stab1 + leftFaces);
        for (std::size_t i = 0; i < faceInfo.totalUDoFs(); ++i) {
            for (std::size_t j = 0; j < faceInfo.leftUDoFs; ++j) {
                vectorStabilizer(i, j) *= lfactor;
            }
        }
        if (faceInfo.internal) {
            double rfactor = std::sqrt(stab_.stab1 + rightFaces);
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

    if (faceInfo.internal && stab_.fluxType2 == FluxType::BREZZI &&
        !isNaturalBCT) {
        // Compute normal vector stabilizer, stab2 * S^T M^{-T} S, S = lift
        // matrix Again we are using the Cholesky decomposition off LL^T = M.
        LinearAlgebra::MiddleSizeMatrix pmassMat =
            computeFaceScalarMassMatrix(face);
        pmassMat.cholesky();
        LinearAlgebra::MiddleSizeMatrix vectorStabilizer =
            computeVectorNormalLiftProjector(face);
        pmassMat.solveLowerTriangular(vectorStabilizer);  // L^{-T} S
        vectorStabilizer = vectorStabilizer.transpose() * vectorStabilizer;
        vectorStabilizer *= stab_.stab2;
        // Copy stabilizer to result matrix
        distributeFaceMatrix(faceInfo, true, vectorStabilizer, result);
    }

    if (stab_.fluxType3 == FluxType::BREZZI) {
        // Compute scalar stabilizer, stab3 * R^T M^{-1} R, R = lift matrix
        // using stab3 * R^T L^{-T} L{-1} R, where LL^{T} = M
        LinearAlgebra::MiddleSizeMatrix scalarStabilizer =
            computeScalarLiftProjector(face);
        massMat.solveLowerTriangular(scalarStabilizer);  // L^{-1}R
        scalarStabilizer = scalarStabilizer.transpose() * scalarStabilizer;
        // minus sign as it is -c_h(p,q)
        scalarStabilizer *= -stab_.stab3;

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
    LinearAlgebra::MiddleSizeVector& ret,
    DGMax::BoundaryConditionType bct) const {
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
    } else if (bct == DGMax::BoundaryConditionType::DIRICHLET) {
        double diameter = face->getDiameter();

        LinearAlgebra::SmallVectorC<DIM> val;
        LinearAlgebra::SmallVector<DIM> phi;
        LinearAlgebra::SmallVectorC<DIM> phi_curlMat;
        val = boundaryValue(fa);

        std::size_t totalUDoFs =
            face->getPtrElementLeft()->getNumberOfBasisFunctions(0);
        std::size_t totalPDoFs =
            face->getPtrElementLeft()->getNumberOfBasisFunctions(1);
        ret.resize(totalUDoFs + totalPDoFs);

        const auto& material = ElementInfos::get(*face->getPtrElementLeft());
        const auto materialCurl =
            material.getMaterialConstantCurl(fa.getPointPhysical());

        const DGMax::MaterialTensor& material1 = materialCurl.adjoint();
        for (std::size_t i = 0; i < totalUDoFs; ++i) {
            fa.basisFunctionUnitNormalCross(i, phi, 0);
            phi_curlMat = material1.applyCurl(fa.basisFunctionCurl(i, 0));
            std::complex<double> value = -(val * phi_curlMat);
            if (stab_.fluxType1 == FluxType::IP) {
                value += stab_.stab1 / (diameter) * (val * phi);
            }
            ret(i) = value;
        }
    } else if (bct == DGMax::BoundaryConditionType::NEUMANN ||
               bct == DGMax::BoundaryConditionType::SILVER_MULLER) {
        // Valid for both IP and Brezzi fluxes
        std::size_t totalUDoFs =
            face->getPtrElementLeft()->getNumberOfBasisFunctions(0);
        std::size_t totalPDoFs =
            face->getPtrElementLeft()->getNumberOfBasisFunctions(1);
        ret.resize(totalUDoFs + totalPDoFs);

        auto val = boundaryValue(fa);
        LinearAlgebra::SmallVector<DIM> phiUN;
        for (std::size_t i = 0; i < totalUDoFs; ++i) {
            fa.basisFunctionUnitNormalCross(i, phiUN, 0);
            ret(i) = val * phiUN;
        }
    } else {
        logger(ERROR,
               "Boundary condition face vector not implemented in DivDGMax");
    }
}

template <std::size_t DIM>
LinearAlgebra::MiddleSizeVector
    DivDGMaxDiscretization<DIM>::brezziFluxBoundaryVector(
        Base::Face* face,
        const DivDGMaxDiscretization<DIM>::FaceInputFunction& boundaryValue) {
    FaceDoFInfo faceInfo = getFaceDoFInfo(*face);
    if (faceInfo.internal) {
        LinearAlgebra::MiddleSizeVector result;
        result.resize(faceInfo.totalDoFs());
        return result;
    }

    // Contribution R^T M^{-1} (stab + numFaces*muinv) r, with
    // r = integral_F bv . u_i dS
    // and R, M the lifting matrix and mass matrix for the vector component

    LinearAlgebra::MiddleSizeMatrix massMatrix =
        computeFaceVectorMassMatrix(face);
    LinearAlgebra::MiddleSizeMatrix liftMatrix =
        computeVectorLiftProjector(face);

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
        faceIntegral.integrate(face, [&](Base::PhysicalFace<DIM>& face) {
            LinearAlgebra::MiddleSizeVector result;
            result.resize(faceInfo.totalUDoFs());
            LinearAlgebra::SmallVector<DIM> basisV;

            // Compute boundary value
            LinearAlgebra::SmallVectorC<DIM> val = boundaryValue(face);

            for (std::size_t i = 0; i < faceInfo.totalUDoFs(); ++i) {
                face.basisFunction(i, basisV, 0);
                result(i) = val * basisV;
            }
            return result;
        });
    rvector *= stab_.stab1 + face->getPtrElementLeft()->getNumberOfFaces();
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
    Base::MeshManipulator<DIM>& mesh, std::size_t timeIntegrationVectorId,
    InputFunction electricField) {

    double error = 0;
    auto end = mesh.elementColEnd();
    for (typename Base::MeshManipulator<DIM>::ElementIterator it =
             mesh.elementColBegin();
         it != end; ++it) {
        error += elementIntegrator_.integrate(
            (*it), [&](Base::PhysicalElement<DIM>& el) {
                return elementErrorIntegrand(el, timeIntegrationVectorId,
                                             electricField);
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

    auto fieldRescale =
        ElementInfos::get(*element).getFieldRescaling(el.getPointPhysical());

    LinearAlgebra::SmallVectorC<DIM> error;
    std::complex<double> potentialError = 0.0;  // Should be zero
    LinearAlgebra::SmallVector<DIM> phi;
    error = exactValues(*element, el.getPointPhysical());
    for (std::size_t i = 0; i < numberOfUDoFs; ++i) {
        el.basisFunction(i, phi, 0);
        error -= fieldRescale.applyDiv(data[i] * phi);
    }
    for (std::size_t i = 0; i < numberOfPDoFs; ++i) {
        potentialError += data[i + numberOfUDoFs] * el.basisFunction(i, 1);
    }

    return error.l2NormSquared() + std::norm(potentialError);
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

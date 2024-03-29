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
#include "Base/H1ConformingTransformation.h"
#include "Base/MeshManipulator.h"

#include "ElementInfos.h"

using namespace hpgem;

template <std::size_t DIM>
DGMaxDiscretization<DIM>::DGMaxDiscretization(std::size_t order, double stab,
                                              bool includeProjector)
    : order_(order),
      stab_(stab),
      includeProjector_(includeProjector),
      matrixHandling_(NORMAL) {

    transforms_.emplace_back(new Base::HCurlConformingTransformation<DIM>());
    if (includeProjector_) {
        transforms_.emplace_back(new Base::H1ConformingTransformation<DIM>());
    }
    for (std::size_t i = 0; i < transforms_.size(); ++i) {
        elementIntegrator_.setTransformation(transforms_[i], i);
        faceIntegrator_.setTransformation(transforms_[i], i);
    }
}

template <std::size_t DIM>
void DGMaxDiscretization<DIM>::initializeBasisFunctions(
    Base::MeshManipulator<DIM>& mesh) const {
    // We would like to configure the number of unknowns here, but this is
    // unfortunately not possible, as it is configured at the creation of
    // the mesh. The best we can do is check if it is configured correctly.
    std::size_t unknowns = mesh.getConfigData()->numberOfUnknowns_;
    if (includeProjector_) {
        logger.assert_always(unknowns == 2,
                             "DGMax+Projector expects 2 unknowns but got %",
                             unknowns);
        mesh.useNedelecDGBasisFunctions(order_);
        mesh.useDefaultConformingBasisFunctions(order_, 1);
    } else {
        logger.assert_always(unknowns == 1, "DGMax expects 1 unknown but got %",
                             unknowns);
        mesh.useNedelecDGBasisFunctions(order_);
    }
    // TODO: This should probably also be exposed by using a constructor
    // parameter.
    // mesh.useAinsworthCoyleDGBasisFunctions();
}

template <std::size_t DIM>
void DGMaxDiscretization<DIM>::computeElementIntegralsImpl(
    Base::MeshManipulator<DIM>& mesh,
    const std::map<std::size_t, InputFunction>& elementVectors,
    LocalIntegrals integrals) {
    logger.assert_always(
        !(matrixHandling_ == ORTHOGONALIZE && !elementVectors.empty()),
        "Mass matrix rescale with input functions is not implemented");
    LinearAlgebra::MiddleSizeVector tempElementVector;

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

        if (integrals == LocalIntegrals::ALL) {
            computeElementMatrices(element);
            postProcessElementMatrices(element);
        }

        for (auto const& elementVectorDef : elementVectors) {
            std::size_t numberOfBasisFunctions =
                element->getNumberOfBasisFunctions(0);
            tempElementVector.resize(numberOfBasisFunctions);
            if (elementVectorDef.second) {
                tempElementVector = elementIntegrator_.integrate(
                    element, [&](Base::PhysicalElement<DIM>& element) {
                        LinearAlgebra::MiddleSizeVector res;
                        elementInnerProduct(element, elementVectorDef.second,
                                            res);  // Initial conditions
                        return res;
                    });
            }
            element->setElementVector(tempElementVector,
                                      elementVectorDef.first);
        }
    }
}

template <std::size_t DIM>
void DGMaxDiscretization<DIM>::computeFaceIntegralsImpl(
    Base::MeshManipulator<DIM>& mesh,
    const std::map<std::size_t, FaceInputFunction>& boundaryVectors,
    DGMax::BoundaryConditionIndicator boundaryIndicator,
    LocalIntegrals integrals) {
    MassMatrixHandling massMatrixHandling = matrixHandling_;
    logger.assert_always(
        !(massMatrixHandling == ORTHOGONALIZE && !boundaryVectors.empty()),
        "Rescale not implemented in combination with boundary conditions");

    LinearAlgebra::MiddleSizeMatrix stiffnessFaceMatrix(0, 0);
    // Mass matrix for the face, already Cholesky factored.
    LinearAlgebra::MiddleSizeMatrix massMatrix(0, 0);
    LinearAlgebra::MiddleSizeVector tempFaceVector;
    Integration::FaceIntegral<DIM> faIntegral;

    faIntegral.setTransformation(
        std::shared_ptr<Base::CoordinateTransformation<DIM>>(
            new Base::HCurlConformingTransformation<DIM>()));
    auto end = mesh.faceColEnd();
    for (auto it = mesh.faceColBegin(); it != end; ++it) {
        Base::Face* face = *it;

        if (integrals == LocalIntegrals::ALL) {
            computeFaceMatrix(face, boundaryIndicator);
            postProcessFaceMatrices(face);
        }

        for (auto const& faceVectorDef : boundaryVectors) {
            std::size_t numberOfBasisFunctions =
                face->getPtrElementLeft()->getNumberOfBasisFunctions(0);
            if (face->isInternal()) {
                numberOfBasisFunctions +=
                    face->getPtrElementRight()->getNumberOfBasisFunctions(0);
            }

            tempFaceVector.resize(numberOfBasisFunctions);
            using BCT = DGMax::BoundaryConditionType;
            BCT bct =
                face->isInternal() ? BCT::INTERNAL : boundaryIndicator(*face);
            if (faceVectorDef.second) {
                tempFaceVector = faIntegral.integrate(
                    face, [&](Base::PhysicalFace<DIM>& face) {
                        LinearAlgebra::MiddleSizeVector res;
                        faceVector(face, faceVectorDef.second, res, bct);
                        return res;
                    });
            }
            face->setFaceVector(tempFaceVector, faceVectorDef.first);
        }
    }
}

template <std::size_t DIM>
void DGMaxDiscretization<DIM>::computeElementMatrices(Base::Element* element) {
    // Shared information
    const auto& material = ElementInfos::get(*element);

    const std::size_t dofU = element->getNumberOfBasisFunctions(0);

    // Mass matrix
    LinearAlgebra::MiddleSizeMatrix massMatrix = elementIntegrator_.integrate(
        element, [&material, dofU = dofU](Base::PhysicalElement<DIM>& pel) {
            const auto& point = pel.getPointPhysical();
            const auto materialDiv = material.getMaterialConstantDiv(point);

            LinearAlgebra::MiddleSizeMatrix ret;
            ret.resize(dofU, dofU);
            LinearAlgebra::SmallVector<DIM> phi_i, phi_j;
            for (std::size_t i = 0; i < dofU; ++i) {
                pel.basisFunction(i, phi_i, 0);
                for (std::size_t j = 0; j < dofU; ++j) {
                    pel.basisFunction(j, phi_j, 0);
                    std::complex<double> val =
                        materialDiv.applyDiv(phi_i) * phi_j;
                    ret(j, i) = val;
                }
            }
            return ret;
        });
    element->setElementMatrix(massMatrix, MASS_MATRIX_ID);
    // Stiffness matrix
    LinearAlgebra::MiddleSizeMatrix stiffnessMatrix =
        elementIntegrator_.integrate(
            element, [&material, dofU = dofU](Base::PhysicalElement<DIM>& pel) {
                LinearAlgebra::MiddleSizeMatrix ret;
                ret.resize(dofU, dofU);

                const auto& point = pel.getPointPhysical();
                const auto materialCurl =
                    material.getMaterialConstantCurl(point);

                LinearAlgebra::SmallVector<DIM> phi_i, phi_j;
                for (std::size_t i = 0; i < dofU; ++i) {
                    phi_i = pel.basisFunctionCurl(i, 0);
                    for (std::size_t j = 0; j < dofU; ++j) {
                        phi_j = pel.basisFunctionCurl(j, 0);
                        std::complex<double> val =
                            materialCurl.applyCurl(phi_i) * phi_j;
                        ret(j, i) = val;
                    }
                }
                return ret;
            });
    element->setElementMatrix(stiffnessMatrix, STIFFNESS_MATRIX_ID);
    if (includeProjector_) {
        LinearAlgebra::MiddleSizeMatrix projectorMatrix =
            elementIntegrator_.integrate(
                element, [&](Base::PhysicalElement<DIM>& pel) {
                    LinearAlgebra::MiddleSizeMatrix ret;
                    const std::size_t dofP =
                        element->getNumberOfBasisFunctions(1);
                    ret.resize(dofP, dofU);

                    const auto& point = pel.getPointPhysical();
                    const auto materialDiv =
                        material.getMaterialConstantDiv(point);

                    LinearAlgebra::SmallVector<DIM> phiU;
                    for (std::size_t i = 0; i < dofU; ++i) {
                        pel.basisFunction(i, phiU, 0);
                        for (std::size_t j = 0; j < dofP; ++j) {
                            ret(j, i) = materialDiv.applyDiv(phiU) *
                                        pel.basisFunctionDeriv(j, 1);
                        }
                    }
                    return ret;
                });
        element->setElementMatrix(projectorMatrix, PROJECTOR_MATRIX_ID);
    }
}

template <std::size_t DIM>
void DGMaxDiscretization<DIM>::postProcessElementMatrices(
    Base::Element* element) const {
    if (matrixHandling_ == NORMAL) {
        return;
    } else if (matrixHandling_ == INVERT) {
        // Note reference to allow overwriting
        LinearAlgebra::MiddleSizeMatrix& massMat =
            element->getElementMatrix(MASS_MATRIX_ID);
        // Note: Overwrites it.
        massMat = massMat.inverse();
    } else if (matrixHandling_ == ORTHOGONALIZE) {
        // Note reference to allow overwriting
        LinearAlgebra::MiddleSizeMatrix& massMatrix =
            element->getElementMatrix(MASS_MATRIX_ID);
        // Inplace cholesky
        massMatrix.cholesky();

        // NOTE: All solves happen in place, updating the matrix in place!
        LinearAlgebra::MiddleSizeMatrix& stiffnessMatrix =
            element->getElementMatrix(STIFFNESS_MATRIX_ID);

        // Compute L^{-1} S L^{-H}, where S is the stiffness matrix and
        // LL^H is the mass matrix.
        massMatrix.solveLowerTriangular(stiffnessMatrix,
                                        LinearAlgebra::Side::OP_LEFT,
                                        LinearAlgebra::Transpose::NOT);
        massMatrix.solveLowerTriangular(
            stiffnessMatrix, LinearAlgebra::Side::OP_RIGHT,
            LinearAlgebra::Transpose::HERMITIAN_TRANSPOSE);
        // Due to rounding errors the matrix might be slightly non
        // Hermitian, fix this by replacing S by 0.5(S + S^H).
        for (std::size_t i = 0; i < stiffnessMatrix.getNumberOfRows(); ++i) {
            stiffnessMatrix(i, i) = std::real(stiffnessMatrix(i, i));
            for (std::size_t j = i; j < stiffnessMatrix.getNumberOfColumns();
                 ++j) {
                std::complex<double> upper =
                    0.5 *
                    (stiffnessMatrix(i, j) + std::conj(stiffnessMatrix(j, i)));
                stiffnessMatrix(i, j) = upper;
                stiffnessMatrix(j, i) = std::conj(upper);
            }
        }

        if (includeProjector_) {
            LinearAlgebra::MiddleSizeMatrix& projectorMatrix =
                element->getElementMatrix(PROJECTOR_MATRIX_ID);
            // Compute B L^{-H}, where B is the projector matrix and L is
            // the Cholesky factor of the mass matrix.
            massMatrix.solveLowerTriangular(
                projectorMatrix, LinearAlgebra::Side::OP_RIGHT,
                LinearAlgebra::Transpose::HERMITIAN_TRANSPOSE);
        }
    }
}

template <std::size_t DIM>
void DGMaxDiscretization<DIM>::elementInnerProduct(
    Base::PhysicalElement<DIM>& el, const InputFunction& function,
    LinearAlgebra::MiddleSizeVector& ret) const {
    const std::size_t numberOfBasisFunctions = el.getNumberOfBasisFunctions(0);

    ret.resize(numberOfBasisFunctions);
    LinearAlgebra::SmallVectorC<DIM> val;
    LinearAlgebra::SmallVector<DIM> phi;
    val = function(*el.getElement(), el.getPointPhysical());
    for (std::size_t i = 0; i < numberOfBasisFunctions; ++i) {
        el.basisFunction(i, phi, 0);
        ret[i] = val * phi;
    }
}

template <std::size_t DIM>
void DGMaxDiscretization<DIM>::computeFaceMatrix(
    Base::Face* face, DGMax::BoundaryConditionIndicator boundaryIndicator) {
    std::size_t numDoFs = face->getNumberOfBasisFunctions(0);
    std::size_t leftDoFs =
        face->getPtrElementLeft()->getNumberOfBasisFunctions(0);
    const bool internalFace = face->isInternal();

    using BCT = DGMax::BoundaryConditionType;
    BCT bct = internalFace ? BCT::INTERNAL : boundaryIndicator(*face);

    Base::FaceMatrix stiffnessMatrix(0, 0);
    if (!internalFace && DGMax::isNaturalBoundary(bct)) {
        // The Neumann boundary condition is a natural boundary condition and
        // such a face does not have a contribution to the stiffness matrix.
        stiffnessMatrix.resize(leftDoFs, numDoFs - leftDoFs);
    } else if (bct == BCT::INTERNAL || bct == BCT::DIRICHLET) {
        // Internal and Dirichlet faces have SIPG like face integrals:
        //  -[[u]]_T {{curl v}} - {{curl u}} [[v]]_T
        //  + stab/h [[u]]_T [[v]]_T

        auto& materialLeft = ElementInfos::get(*face->getPtrElementLeft());
        auto* materialRight = ElementInfos::get(face->getPtrElementRight());

        // Factor for averaging. Negative sign is from the weak formulation
        double factor = -(internalFace ? 0.5 : 1.);
        // Standard rescaling of the stability parameter so that it does not
        // need to depend on the mesh size.
        double localStab = stab_ / face->getDiameter();
        stiffnessMatrix =
            faceIntegrator_.integrate(face, [&](Base::PhysicalFace<DIM>& pfa) {
                Base::FaceMatrix ret(leftDoFs, numDoFs - leftDoFs);

                auto p = pfa.getPointPhysical();
                const auto materialLeftCurl =
                    materialLeft.getMaterialConstantCurl(p);
                const auto materialRightCurl =
                    materialRight != nullptr
                        ? materialRight->getMaterialConstantCurl(p)
                        : DGMax::MaterialTensor();

                LinearAlgebra::SmallVector<DIM> phiNi, phiNj;
                LinearAlgebra::SmallVectorC<DIM> phiCi, phiCj;

                for (std::size_t i = 0; i < numDoFs; ++i) {
                    phiCi =
                        ((i < leftDoFs) ? materialLeftCurl : materialRightCurl)
                            .applyCurl(pfa.basisFunctionCurl(i, 0));
                    pfa.basisFunctionUnitNormalCross(i, phiNi, 0);

                    for (std::size_t j = 0; j < numDoFs; ++j) {
                        phiCj = ((j < leftDoFs) ? materialLeftCurl.adjoint()
                                                : materialRightCurl.adjoint())
                                    .applyCurl(pfa.basisFunctionCurl(j, 0));
                        pfa.basisFunctionUnitNormalCross(j, phiNj, 0);
                        std::complex<double> value =
                            factor * (phiCi * phiNj + phiNi * phiCj) +
                            localStab * phiNi * phiNj;
                        ret(j, i) = value;
                    }
                }

                return ret;
            });
    }
    face->setFaceMatrix(stiffnessMatrix, FACE_STIFFNESS_MATRIX_ID);

    // Impedance contribution to the stiffness matrix
    if (bct == DGMax::BoundaryConditionType::SILVER_MULLER) {
        auto* material = dynamic_cast<ElementInfos*>(
            face->getPtrElementLeft()->getUserData());
        logger.assert_debug(material != nullptr, "No material");

        std::complex<double> impedance =
            std::complex<double>(0, material->getImpedance());
        stiffnessMatrix = faceIntegrator_.integrate(
            face, [numDoFs, impedance](Base::PhysicalFace<DIM>& pfa) {
                Base::FaceMatrix result(numDoFs, 0);
                for (std::size_t i = 0; i < numDoFs; ++i) {
                    LinearAlgebra::SmallVector<DIM> phiUNi;
                    pfa.basisFunctionUnitNormalCross(i, phiUNi, 0);
                    for (std::size_t j = 0; j < numDoFs; ++j) {
                        LinearAlgebra::SmallVector<DIM> phiUNj;
                        pfa.basisFunctionUnitNormalCross(j, phiUNj, 0);
                        result(j, i) = -impedance * (phiUNi * phiUNj);
                    }
                }
                return result;
            });
    } else {
        stiffnessMatrix *= 0.0;
    }
    face->setFaceMatrix(stiffnessMatrix, FACE_IMPEDANCE_MATRIX_ID);
}

template <std::size_t DIM>
void DGMaxDiscretization<DIM>::postProcessFaceMatrices(Base::Face* face) const {
    if (matrixHandling_ == ORTHOGONALIZE) {
        // A reference to the matrix that will be updated in place.
        Base::FaceMatrix& faceMatrix =
            face->getFaceMatrix(FACE_STIFFNESS_MATRIX_ID);

        bool isInternal = face->isInternal();
        std::size_t sideCount = isInternal ? 2 : 1;
        // The mass matrix/matrices that form a block diagonal matrix
        std::array<LinearAlgebra::MiddleSizeMatrix*, 2> massMatrices = {
            nullptr, nullptr};
        massMatrices[0] =
            &face->getPtrElementLeft()->getElementMatrix(MASS_MATRIX_ID);
        std::array<hpgem::Base::Side, 2> sides = {Base::Side::LEFT,
                                                  Base::Side::RIGHT};
        if (isInternal) {
            massMatrices[1] =
                &face->getPtrElementRight()->getElementMatrix(MASS_MATRIX_ID);
        }
        for (std::size_t i = 0; i < sideCount; ++i) {
            for (std::size_t j = 0; j < sideCount; ++j) {
                LinearAlgebra::MiddleSizeMatrix& subMatrix =
                    faceMatrix.getElementMatrix(sides[i], sides[j]);
                // Rescale the submatrix, note that different mass matrices may
                // be used for the left and right side.
                massMatrices[i]->solveLowerTriangular(
                    subMatrix, hpgem::LinearAlgebra::Side::OP_LEFT,
                    hpgem::LinearAlgebra::Transpose::NOT);
                massMatrices[j]->solveLowerTriangular(
                    subMatrix, LinearAlgebra::Side::OP_RIGHT,
                    LinearAlgebra::Transpose::HERMITIAN_TRANSPOSE);
            }
        }
    }
}

template <std::size_t DIM>
void DGMaxDiscretization<DIM>::faceVector(
    Base::PhysicalFace<DIM>& fa, const FaceInputFunction& boundaryCondition,
    LinearAlgebra::MiddleSizeVector& ret,
    DGMax::BoundaryConditionType bct) const {

    std::size_t numDoFs = fa.getNumberOfBasisFunctions();
    ret.resize(numDoFs);

    const Base::Face* face = fa.getFace();

    if (face->isInternal()) {
        // Set the vector to zero for the internal contribution
        ret.set(0.0);
    } else if (bct == DGMax::BoundaryConditionType::DIRICHLET) {
        double diameter = face->getDiameter();
        LinearAlgebra::SmallVectorC<DIM> val, phi_curl;
        LinearAlgebra::SmallVector<DIM> phi;

        const auto point = fa.getPointPhysical();
        const auto& materialCurlAdj =
            ElementInfos::get(*face->getPtrElementLeft())
                .getMaterialConstantCurl(point)
                .adjoint();

        val = boundaryCondition(fa);

        for (std::size_t i = 0; i < numDoFs; ++i) {
            fa.basisFunctionUnitNormalCross(i, phi, 0);
            phi_curl = materialCurlAdj.applyCurl(fa.basisFunctionCurl(i, 0));
            ret(i) = -(val * phi_curl) + stab_ / diameter * (val * phi);
        }
    } else if (bct == DGMax::BoundaryConditionType::NEUMANN ||
               bct == DGMax::BoundaryConditionType::SILVER_MULLER) {
        // Compute g_N (n x phi_i)
        LinearAlgebra::SmallVector<DIM> phiN;
        LinearAlgebra::SmallVectorC<DIM> val = boundaryCondition(fa);
        for (std::size_t i = 0; i < numDoFs; ++i) {
            fa.basisFunctionUnitNormalCross(i, phiN, 0);
            ret(i) = val * phiN;
        }
    } else {
        logger(ERROR,
               "No boundary source term implemented for this boundary "
               "condition type.");
    }
}

template <std::size_t DIM>
typename DGMaxDiscretization<DIM>::Fields
    DGMaxDiscretization<DIM>::computeFields(
        const Base::Element* element, const Geometry::PointReference<DIM>& p,
        const LinearAlgebra::MiddleSizeVector& coefficients) const {
    Base::PhysicalElement<DIM> physicalElement;
    physicalElement.setElement(element);
    std::shared_ptr<Base::CoordinateTransformation<DIM>> transform{
        new Base::HCurlConformingTransformation<DIM>()};
    physicalElement.setTransformation(transform);
    physicalElement.setPointReference(p);

    Fields result;

    for (std::size_t i = 0; i < element->getNumberOfBasisFunctions(0); ++i) {
        LinearAlgebra::SmallVector<DIM> phi;
        physicalElement.basisFunction(i, phi, 0);
        result.electricField += coefficients[i] * phi;
        result.electricFieldCurl +=
            coefficients[i] * physicalElement.basisFunctionCurl(i);
    }
    const ElementInfos& userData = ElementInfos::get(*element);

    result.material = userData.getMaterial();

    // Rescale for PMLs
    const auto& pPhys = physicalElement.getPointPhysical();
    result.electricField =
        userData.getFieldRescaling(pPhys).applyDiv(result.electricField);
    result.electricFieldCurl = userData.getCurlFieldRescaling(pPhys).applyCurl(
        result.electricFieldCurl);
    return result;
}

template <std::size_t DIM>
void DGMaxDiscretization<DIM>::writeFields(
    Output::VTKSpecificTimeWriter<DIM>& writer,
    std::size_t timeIntegrationVectorId) const {
    using VecR = LinearAlgebra::SmallVector<DIM>;
    std::map<std::string, std::function<double(Fields&)>> scalars;
    std::map<std::string, std::function<VecR(Fields&)>> vectors;

    vectors["Ereal"] = [](Fields& fields) {
        return fields.electricField.real();
    };
    vectors["Eimag"] = [](Fields& fields) {
        return fields.electricField.imag();
    };

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

    writer.template writeMultiple<Fields>(
        [this](Base::Element* element,
               const Geometry::PointReference<DIM>& point, std::size_t) {
            LinearAlgebra::MiddleSizeVector coefficients =
                element->getTimeIntegrationVector(0);
            // When using the Hermitian system we applied a rescaling of the
            // solution coefficients to use y = L^H x (LL^H = M is the Cholesky
            // decomposition of the mass matrix and x the actual coefficients).
            // Undo this transformation to correctly compute the fields.
            if (matrixHandling_ == ORTHOGONALIZE) {
                // In place solve
                element->getElementMatrix(MASS_MATRIX_ID)
                    .solveLowerTriangular(
                        coefficients, hpgem::LinearAlgebra::Side::OP_LEFT,
                        hpgem::LinearAlgebra::Transpose::HERMITIAN_TRANSPOSE);
            }
            return computeFields(element, point, coefficients);
        },
        scalars, vectors);

    // Output the epsilon separately
    writer.write(
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
LinearAlgebra::SmallVector<4> DGMaxDiscretization<DIM>::computeEnergyFluxes(
    Base::Face& face, hpgem::Base::Side side, double wavenumber,
    std::size_t timeIntegrationVectorId,
    const DGMax::FieldPattern<DIM>* background) {

    using VecC = LinearAlgebra::SmallVectorC<DIM>;

    auto coefficients = face.getTimeIntegrationVector(timeIntegrationVectorId);
    auto leftDoFs = face.getPtrElementLeft()->getNumberOfBasisFunctions();
    if (matrixHandling_ == ORTHOGONALIZE) {
        // We need to undo the orthogonalization
        // for boundary faces this is easy (see computing of fields). For
        // internal faces we need to solve the left and right coefficients with
        // different mass matrices.
        LinearAlgebra::MiddleSizeVector temp;
        auto& solveCoefs = face.isInternal() ? temp : coefficients;
        if (face.isInternal()) {
            // Copy for the left side
            temp.resize(leftDoFs);
            for (int i = 0; i < leftDoFs; ++i) {
                temp[i] = coefficients[i];
            }
        }

        // In place solve
        face.getPtrElementLeft()
            ->getElementMatrix(MASS_MATRIX_ID)
            .solveLowerTriangular(
                coefficients, hpgem::LinearAlgebra::Side::OP_LEFT,
                hpgem::LinearAlgebra::Transpose::HERMITIAN_TRANSPOSE);

        if (face.isInternal()) {
            // Copy updated values back to the coefficient vector
            for (std::size_t i = 0; i < leftDoFs; ++i) {
                coefficients[i] = temp[i];
            }

            // Same process but now for the coefficients on the right side.
            std::size_t rightDoFs =
                face.getPtrElementRight()->getNumberOfBasisFunctions();
            temp.resize(rightDoFs);
            for (std::size_t i = 0; i < rightDoFs; ++i) {
                temp[i] = coefficients[i + rightDoFs];
            }
            face.getPtrElementRight()
                ->getElementMatrix(MASS_MATRIX_ID)
                .solveLowerTriangular(
                    temp, hpgem::LinearAlgebra::Side::OP_LEFT,
                    hpgem::LinearAlgebra::Transpose::HERMITIAN_TRANSPOSE);
            for (std::size_t i = 0; i < rightDoFs; ++i) {
                coefficients[i + rightDoFs] = temp[i];
            }
        }
    }

    double factor = face.isInternal() ? 0.5 : 1.0;

    auto flux = faceIntegrator_.integrate(
        &face,
        [&coefficients, &factor, background](Base::PhysicalFace<DIM>& pface) {
            LinearAlgebra::SmallVector<4> result;
            // Average curl of the field
            VecC avgCurl;
            // n cross E for the two sides
            VecC avgE;
            LinearAlgebra::SmallVector<DIM> phi;
            VecC normal = pface.getUnitNormalVector();
            for (int i = 0; i < pface.getNumberOfBasisFunctions(); ++i) {
                avgCurl +=
                    factor * coefficients[i] * pface.basisFunctionCurl(i);
                pface.basisFunction(i, phi);
                avgE += factor * coefficients[i] * normal.crossProduct(phi);
            }
            // Compute n . Re(S) = n . Re(E x flux[H]^*)
            // with H = i omega muinv Curl E
            //   = omega n . Im(E x muinv flux[Curl E]^*)
            // using complex inner product
            //   = omega Im((n x E) . muinv flux[Curl E])
            // return (avgE * avgCurl).imag();
            result[0] = (avgE * avgCurl).imag();
            if (background) {
                VecC backgroundE = normal.crossProduct(
                    background->field(pface.getPointPhysical()));
                VecC totalE = avgE + backgroundE;
                VecC backgroundCurl =
                    background->fieldCurl(pface.getPointPhysical());
                VecC totalCurl = avgCurl + backgroundCurl;
                result[1] = (backgroundE * backgroundCurl).imag();
                result[2] =
                    (backgroundE * avgCurl + avgE * backgroundCurl).imag();
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
double DGMaxDiscretization<DIM>::computeFieldL2Integral(Base::Face& face,
                                                        Base::Side side,
                                                        std::size_t vector_id) {
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

// TODO: The code saves snapshots in the timeIntegrationVector, this is not
// particularly nice It might be better to pass the global vector here and
// distribute it ourselves.
template <std::size_t DIM>
std::map<typename DGMaxDiscretizationBase::NormType, double>
    DGMaxDiscretization<DIM>::computeError(
        Base::MeshManipulator<DIM>& mesh, std::size_t timeVector,
        DGMaxDiscretization<DIM>::InputFunction electricField,
        DGMaxDiscretization<DIM>::InputFunction electricFieldCurl,
        std::set<DGMaxDiscretizationBase::NormType> norms) {
    // Note these are actually the squared norms
    double l2Norm = 0;
    double hCurlNorm = 0;
    // Setup the element integration.
    bool l2Wanted = norms.find(NormType::L2) != norms.end();
    bool hcurlWanted = norms.find(NormType::HCurl) != norms.end();
    bool dgWanted = norms.find(NormType::DG) != norms.end();

    auto end = mesh.elementColEnd();
    for (typename Base::MeshManipulator<DIM>::ElementIterator it =
             mesh.elementColBegin();
         it != end; ++it) {
        LinearAlgebra::SmallVector<2> errors = elementIntegrator_.integrate(
            (*it), [&](Base::PhysicalElement<DIM>& el) {
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
        auto end = mesh.faceColEnd();
        for (typename Base::MeshManipulator<DIM>::FaceIterator it =
                 mesh.faceColBegin();
             it != end; ++it) {
            dgNorm += faceIntegrator_.integrate(
                *it, [&](Base::PhysicalFace<DIM>& face) {
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

    LinearAlgebra::SmallVector<DIM> phi, phiCurl;
    LinearAlgebra::SmallVectorC<DIM> error, errorCurl;
    auto fieldRescale =
        ElementInfos::get(*element).getFieldRescaling(el.getPointPhysical());
    auto curlRescale = ElementInfos::get(*element).getCurlFieldRescaling(
        el.getPointPhysical());

    error = exactValues(*element, el.getPointPhysical());
    if (computeCurl) {
        errorCurl = curlValues(*element, el.getPointPhysical());
    }
    LinearAlgebra::MiddleSizeVector data;
    data = element->getTimeIntegrationVector(timeVector);
    for (std::size_t i = 0; i < element->getNrOfBasisFunctions(0); ++i) {
        el.basisFunction(i, phi, 0);
        error -= fieldRescale.applyDiv(data[i] * phi);
        if (computeCurl) {
            phiCurl = el.basisFunctionCurl(i, 0);
            errorCurl -= curlRescale.applyCurl(data[i] * phiCurl);
        }
    }
    double l2Error = error.l2NormSquared();
    double curlError = errorCurl.l2NormSquared();
    LinearAlgebra::SmallVector<2> errors;
    errors[0] = l2Error;
    errors[1] = curlError;
    return errors;
}

template <std::size_t DIM>
double DGMaxDiscretization<DIM>::faceErrorIntegrand(
    Base::PhysicalFace<DIM>& fa, std::size_t timeVector,
    DGMaxDiscretization<DIM>::InputFunction exactSolution) const {
    logger(WARN, "DG-Error does not work with PMLs");

    // The face error part of the DG norm is given by
    // || h^0.5 [[u - u_h]]_T ||^2. Where h is the diameter of the face, u and
    // u_h are the exact and computed solutions. Further more the [[ . ]]_T is
    // the tangential jump, [[ a ]]_T = a_L x n_L + a_R x n_R for internal faces
    // and a x n for boundary faces (n is the normal).
    const Base::Face* face = fa.getFace();
    LinearAlgebra::SmallVector<DIM> normal = fa.getNormalVector();
    normal /= normal.l2Norm();
    const Geometry::PointReference<2>& p = fa.getPointReference();

    Base::Element* element =
        const_cast<Base::Element*>(face->getPtrElementLeft());
    PointPhysicalT PPhys;
    const Geometry::PointReference<DIM>& pElement =
        face->mapRefFaceToRefElemL(p);

    PPhys = element->referenceToPhysical(pElement);
    LinearAlgebra::SmallVector<DIM> error, phiNormal, solutionValues;

    // Compute u_L x n_L
    solutionValues = exactSolution(*element, PPhys);
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
        solutionValues = exactSolution(*element, PPhys);
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
    double result = error.l2Norm();
    result *= result;
    // To remove double contribution of flux computed on the boundary faces by
    // different processors
    if (face->getFaceType() == Geometry::FaceType::SUBDOMAIN_BOUNDARY ||
        face->getFaceType() == Geometry::FaceType::PERIODIC_SUBDOMAIN_BC) {
        result /= 2;
    }
    return result;
}

template class DGMaxDiscretization<2>;
template class DGMaxDiscretization<3>;
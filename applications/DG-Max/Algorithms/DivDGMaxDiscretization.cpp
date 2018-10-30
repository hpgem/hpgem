/*
This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
http://hpgem.org/about-the-code/team

This code is distributed using BSD 3-Clause License. A copy of which can found below.


Copyright (c) 2018, Univesity of Twenete
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "DivDGMaxDiscretization.h"

#include "Base/HCurlConformingTransformation.h"
#include "Integration/ElementIntegral.h"
#include "Integration/FaceIntegral.h"


void DivDGMaxDiscretization::initializeBasisFunctions(Base::MeshManipulator<DIM> &mesh)
{
    //TODO: This needs the additional unknown id.
    mesh.useNedelecDGBasisFunctions();
    mesh.useDefaultDGBasisFunctions(1);
}

void DivDGMaxDiscretization::computeElementIntegrands(
        Base::MeshManipulator<DIM> &mesh, bool invertMassMatrix,
        const DivDGMaxDiscretization::InputFunction &sourceTerm,
        const DivDGMaxDiscretization::InputFunction &initialCondition,
        const DivDGMaxDiscretization::InputFunction &initialConditionDerivative) const
{
    //TODO: Add initial condition integration.
    LinearAlgebra::MiddleSizeMatrix massMatrix(2, 2),  stiffnessMatrix(2, 2);

    LinearAlgebra::MiddleSizeVector vector1(2), vector2(2), sourceVector(2);

    std::size_t totalDoFs = 0, totalUDoFs = 0, totalPDoFs = 0;
    Integration::ElementIntegral<DIM> elementIntegral(false);

    elementIntegral.setTransformation(std::shared_ptr<Base::CoordinateTransformation<DIM>> (new Base::HCurlConformingTransformation<DIM>()), 0);
    elementIntegral.setTransformation(std::shared_ptr<Base::CoordinateTransformation<DIM>> (new Base::H1ConformingTransformation<DIM>()), 1);

    for (hpGemUIExtentions::ElementIterator it = mesh.elementColBegin(); it != mesh.elementColEnd(); ++it)
    {
        totalUDoFs = (*it)->getNumberOfBasisFunctions(0);
        totalPDoFs = (*it)->getNumberOfBasisFunctions(1);
        totalDoFs = totalUDoFs + totalPDoFs;

        //mass matrix
        massMatrix.resize(totalDoFs, totalDoFs);
        massMatrix = elementIntegral.integrate<LinearAlgebra::MiddleSizeMatrix>((*it),
                [&](Base::PhysicalElement<DIM>& element) {
            LinearAlgebra::MiddleSizeMatrix result;
            elementMassMatrix(element, result);
            return result;
        });
        if(invertMassMatrix)
        {
            massMatrix = massMatrix.inverse();
        }
        (*it)->setElementMatrix(massMatrix, ELEMENT_MASS_MATRIX_ID);

        stiffnessMatrix.resize(totalDoFs, totalDoFs);
        stiffnessMatrix = elementIntegral.integrate<LinearAlgebra::MiddleSizeMatrix>((*it),
                [&](Base::PhysicalElement<DIM>& element) {
            LinearAlgebra::MiddleSizeMatrix result, temp;
            elementStiffnessMatrix(element, result);
            elementScalarVectorCoupling(element, temp);
            result += temp;
            return result;
        });

        (*it)->setElementMatrix(stiffnessMatrix, ELEMENT_STIFFNESS_MATRIX_ID);

        if (initialCondition)
        {
            logger.assert(false, "Initial condition integration not ported yet.");
        }

        if (initialConditionDerivative)
        {
            logger.assert(false, "Initial condition derivative integration not ported yet.");
        }

        if (sourceTerm)
        {
            sourceVector.resize(totalDoFs);
            sourceVector = elementIntegral.integrate<LinearAlgebra::MiddleSizeVector>((*it),
                    [&](Base::PhysicalElement<DIM>& element) {
                LinearAlgebra::MiddleSizeVector result;
                elementSourceVector(element, sourceTerm, result);
                return result;
            });
            (*it)->setElementVector(sourceVector, ELEMENT_SOURCE_VECTOR_ID);
        }

    }
}

void DivDGMaxDiscretization::computeFaceIntegrals(
        Base::MeshManipulator<DIM> &mesh,
        DivDGMaxDiscretization::FaceInputFunction boundaryCondition,
        Stab stab) const
{
    LinearAlgebra::MiddleSizeMatrix faceMatrix(2, 2);
    LinearAlgebra::MiddleSizeVector faceVector(2);
    Integration::FaceIntegral<DIM> faceIntegral(false);

    faceIntegral.setTransformation(std::shared_ptr<Base::CoordinateTransformation<DIM>> (new Base::HCurlConformingTransformation<DIM>()), 0);
    faceIntegral.setTransformation(std::shared_ptr<Base::CoordinateTransformation<DIM>> (new Base::H1ConformingTransformation<DIM>()), 1);

    for (hpGemUIExtentions::FaceIterator it = mesh.faceColBegin(); it != mesh.faceColEnd(); ++it)
    {

        std::size_t totalDoFs = 0;

        totalDoFs = (*it)->getPtrElementLeft()->getNumberOfBasisFunctions(0) + (*it)->getPtrElementLeft()->getNumberOfBasisFunctions(1);
        if ((*it)->isInternal())
        {
            totalDoFs += (*it)->getPtrElementRight()->getNumberOfBasisFunctions(0)
                      + (*it)->getPtrElementRight()->getNumberOfBasisFunctions(1);
        }

        faceMatrix.resize(totalDoFs, totalDoFs);
        faceVector.resize(totalDoFs);

        faceMatrix = faceIntegral.integrate<LinearAlgebra::MiddleSizeMatrix>((*it),
                [&] (Base::PhysicalFace<DIM>& face) {
            LinearAlgebra::MiddleSizeMatrix result, temp;
            faceStiffnessMatrix1(face, result);

            faceStiffnessMatrix2(face, temp, stab.stab1);
            result += temp;

            temp *= 0; // Reset the variable;
            faceStiffnessMatrix3(face, temp, stab.stab2);
            result += temp;

            temp *= 0; // Reset the variable;
            faceScalarVectorCoupling(face, temp);
            result += temp;

            temp *= 0; // Reset the variable;
            faceStiffnessScalarMatrix4(face, temp, stab.stab3);
            result += temp;
            // Reset no longer needed.
            return result;
        });
        (*it)->setFaceMatrix(faceMatrix, FACE_STIFFNESS_MATRIX_ID);

        if (boundaryCondition)
        {
            faceVector = faceIntegral.integrate<LinearAlgebra::MiddleSizeVector>((*it),
                    [&](Base::PhysicalFace<DIM> &face) {
                LinearAlgebra::MiddleSizeVector result;
                faceBoundaryVector(face, boundaryCondition, result, stab.stab1);
                return result;
            });
            (*it)->setFaceVector(faceVector, FACE_BOUNDARY_VECTOR_ID);
        }
    }

}

LinearAlgebra::SmallVector<DIM> DivDGMaxDiscretization::computeField (
        const Base::Element *element, const Geometry::PointReference<DIM> &point,
        const LinearAlgebra::MiddleSizeVector &coefficients) const
{
    logger.log(Log::WARN, "Only computing the real part of the field.");
    LinearAlgebra::SmallVector<DIM> result, phiU;
    std::size_t nb0 = element->getNumberOfBasisFunctions(0), nb1 = element->getNumberOfBasisFunctions(1);

    // Setup the physical element with the correct transformations
    Base::PhysicalElement<DIM> physicalElement;
    physicalElement.setElement(element);
    std::shared_ptr<Base::CoordinateTransformation<DIM>> coordinateTransformationU {new Base::HCurlConformingTransformation<DIM>()};
    std::shared_ptr<Base::CoordinateTransformation<DIM>> coordinateTransformationP {new Base::H1ConformingTransformation<DIM>()};
    physicalElement.setTransformation(coordinateTransformationU, 0);
    physicalElement.setTransformation(coordinateTransformationP, 1);
    physicalElement.setPointReference(point);
    Geometry::PointPhysical<DIM> pointPhysical;
    pointPhysical = element->referenceToPhysical(point);

    // reconstruct the field E = u + grad q.
    for(std::size_t i = 0; i < nb0; ++i)
    {
        physicalElement.basisFunction(i, phiU, 0);
        result += std::real(coefficients[i]) * phiU;
    }
    for (std::size_t i = 0; i < nb1; ++i)
    {
        result += std::real(coefficients[nb0+i]) * physicalElement.basisFunctionDeriv(i, 1);
    }
    return result;
}

void DivDGMaxDiscretization::elementMassMatrix(Base::PhysicalElement<DIM> &el,
                                               LinearAlgebra::MiddleSizeMatrix &ret) const
{
    const Base::Element* element = el.getElement();
    std::size_t el_size1 = element->getNumberOfBasisFunctions(0);
    std::size_t el_size2 = element->getNumberOfBasisFunctions(1);
    ret.resize(el_size1 + el_size2, el_size1 + el_size2);
    double epsilon = static_cast<ElementInfos*>(element->getUserData())->epsilon_;
    LinearAlgebra::SmallVector<DIM> phi_i, phi_j;
    for (std::size_t i = 0; i < el_size1; ++i)
    {
        el.basisFunction(i, phi_i, 0);
        for (std::size_t j = i; j < el_size1; ++j)
        {
            el.basisFunction(j, phi_j, 0);
            ret(i, j) = phi_i * phi_j * epsilon;
            ret(j, i) = ret(i, j);
        }
    }
}

void DivDGMaxDiscretization::elementStiffnessMatrix(Base::PhysicalElement<DIM> &el,
                                                    LinearAlgebra::MiddleSizeMatrix &ret) const
{
    const Base::Element* element = el.getElement();
    std::size_t el_size1 = element->getNumberOfBasisFunctions(0);
    std::size_t el_size2 = element->getNumberOfBasisFunctions(1);
    ret.resize(el_size1 + el_size2, el_size1 + el_size2);
    LinearAlgebra::SmallVector<DIM> phi_i, phi_j;

    for (std::size_t i = 0; i < el_size1; ++i)
    {
        phi_i = el.basisFunctionCurl(i, 0);
        for (std::size_t j = i; j < el_size1; ++j)
        {
            phi_j = el.basisFunctionCurl(j, 0);
            ret(i, j) = phi_i * phi_j;
            ret(j, i) = ret(i, j);
        }
    }
}

void DivDGMaxDiscretization::elementScalarVectorCoupling(Base::PhysicalElement<DIM> &el,
                                                         LinearAlgebra::MiddleSizeMatrix &ret) const
{
    const Base::Element* element = el.getElement();
    std::size_t uDoFs = element->getNumberOfBasisFunctions(0);
    std::size_t pDoFs = element->getNumberOfBasisFunctions(1);
    ret.resize(uDoFs + pDoFs, uDoFs + pDoFs);
    double epsilon = static_cast<ElementInfos*>(element->getUserData())->epsilon_;
    LinearAlgebra::SmallVector<DIM> phi_i, phi_j;

    // Note, this loop only loops over the basis functions for the second
    // unknown. However, it needs the offset for the first unknown.
    for (std::size_t i = uDoFs; i < uDoFs + pDoFs; ++i)
    {
        phi_i = el.basisFunctionDeriv(i - uDoFs, 1);

        for (std::size_t j = 0; j < uDoFs; ++j)
        {
            el.basisFunction(j, phi_j, 0);
            ret(j, i) = -1.0 * (phi_i * phi_j) * epsilon;
            ret(i, j) = ret(j, i);

        }
    }
}

void DivDGMaxDiscretization::elementSourceVector(
        Base::PhysicalElement<DIM> &el,
        const InputFunction &source,
        LinearAlgebra::MiddleSizeVector &ret) const
{
    const Base::Element* element = el.getElement();
    const Geometry::PointReference<DIM>& p = el.getPointReference();
    std::size_t uDoFs = element->getNumberOfBasisFunctions(0);
    std::size_t pDoFs = element->getNumberOfBasisFunctions(1);
    ret.resize(uDoFs + pDoFs);

    PointPhysicalT pPhys;
    pPhys = element->referenceToPhysical(p);
    LinearAlgebra::SmallVector<DIM> sourceValue, phi;
    source(pPhys, sourceValue);
    for (std::size_t i = 0; i < (uDoFs); ++i)
    {
        el.basisFunction(i, phi, 0);
        ret(i) = phi * sourceValue;
    }
}

void DivDGMaxDiscretization::faceStiffnessMatrix1(
        Base::PhysicalFace<DIM> &fa,
        LinearAlgebra::MiddleSizeMatrix &ret) const
{
    const Base::Face* face = fa.getFace();
    std::size_t totalUDoFs = face->getPtrElementLeft()->getNumberOfBasisFunctions(0);
    std::size_t totalPDoFs = face->getPtrElementLeft()->getNumberOfBasisFunctions(1);

    const std::size_t leftUDoFs = totalUDoFs;
    const std::size_t leftPDoFs = totalPDoFs;


    if(face->isInternal())
    {
        totalUDoFs = face->getPtrElementLeft()->getNumberOfBasisFunctions(0) + face->getPtrElementRight()->getNumberOfBasisFunctions(0);
        totalPDoFs = face->getPtrElementLeft()->getNumberOfBasisFunctions(1) + face->getPtrElementRight()->getNumberOfBasisFunctions(1);
    }
    ret.resize(totalUDoFs + totalPDoFs, totalUDoFs + totalPDoFs);

    std::vector<std::size_t> indices (totalUDoFs + totalPDoFs);
    std::vector<LinearAlgebra::SmallVector<DIM>> phiCurl (totalUDoFs);
    std::vector<LinearAlgebra::SmallVector<DIM>> phiNormal (totalUDoFs);

    // Left element
    for(std::size_t i = 0; i < leftUDoFs; ++i)
    {
        indices[i] = i;
        phiCurl[i] = fa.basisFunctionCurl(i, 0);
        fa.basisFunctionUnitNormalCross(i, phiNormal[i], 0);
    }
    // Right element
    for (std::size_t i = leftUDoFs; i < totalUDoFs; ++i)
    {
        // Note, only here we need to offset by leftPDoFs, as the indices are over all
        // unknowns and not just those for one (as in fa.basisFunction*)
        indices[i] = i + leftPDoFs;
        phiCurl[i] = fa.basisFunctionCurl(i, 0);
        fa.basisFunctionUnitNormalCross(i, phiNormal[i], 0);
    }


    double factor = face->isInternal() ? -0.5 : -1;
    for (std::size_t i = 0; i < totalUDoFs; ++i)
    {
        const std::size_t iIndex = indices[i];
        for (std::size_t j = i; j < totalUDoFs; ++j)
        {
            const double entry = factor * (phiCurl[i] * phiNormal[j] + phiCurl[j] * phiNormal[i]);
            const std::size_t jIndex = indices[j];
            ret(iIndex, jIndex) = entry;
            ret(jIndex, iIndex) = entry;
        }
    }
}

void DivDGMaxDiscretization::faceStiffnessMatrix2(
        Base::PhysicalFace<DIM> &fa, LinearAlgebra::MiddleSizeMatrix &ret,
        double stab) const
{
    const Base::Face* face = fa.getFace();
    std::size_t totalUDoFs = face->getPtrElementLeft()->getNumberOfBasisFunctions(0);
    std::size_t totalPDoFs = face->getPtrElementLeft()->getNumberOfBasisFunctions(1);
    std::size_t leftUDoFs = totalUDoFs;
    std::size_t leftPDoFs = totalPDoFs;
    if(face->isInternal())
    {
        totalUDoFs = face->getPtrElementLeft()->getNumberOfBasisFunctions(0) + face->getPtrElementRight()->getNumberOfBasisFunctions(0);
        totalPDoFs = face->getPtrElementLeft()->getNumberOfBasisFunctions(1) + face->getPtrElementRight()->getNumberOfBasisFunctions(1);
    }

    ret.resize(totalUDoFs + totalPDoFs, totalUDoFs + totalPDoFs);

    std::vector<std::size_t> indices (totalUDoFs);
    std::vector<LinearAlgebra::SmallVector<DIM>> phiNormalCross (totalUDoFs);

    // Left element
    for (std::size_t i = 0; i < leftUDoFs; ++i)
    {
        indices[i] = i;
        fa.basisFunctionUnitNormalCross(i, phiNormalCross[i], 0);
    }
    // Right element
    for (std::size_t i = leftUDoFs; i < totalUDoFs; ++i)
    {
        indices[i] = i + leftPDoFs;
        fa.basisFunctionUnitNormalCross(i, phiNormalCross[i], 0);
    }

    for (std::size_t i = 0; i < totalUDoFs; ++i)
    {
        const std::size_t iIndex = indices[i];
        for (std::size_t j = 0; j < totalUDoFs; ++j)
        {
            const std::size_t jIndex = indices[j];
            double entry = stab * (phiNormalCross[i] * phiNormalCross[j]);
            ret(iIndex, jIndex) = entry;
            ret(jIndex, iIndex) = entry;
        }
    }
}

void DivDGMaxDiscretization::faceStiffnessMatrix3(Base::PhysicalFace<DIM> &fa, LinearAlgebra::MiddleSizeMatrix &ret,
                                                  double stab2) const
{
    // TODO: Cleanup.
    const Base::Face* face = fa.getFace();
    std::size_t totalUDoFs = face->getPtrElementLeft()->getNumberOfBasisFunctions(0);
    std::size_t totalPDoFs = face->getPtrElementLeft()->getNumberOfBasisFunctions(1);
    std::size_t leftUDoFs = totalUDoFs; // to determine left/right part of totalUDoFs
    std::size_t leftPDoFs = totalPDoFs; // to determine left/right part of totalUDoFs
    if (face->isInternal())
    {

        totalUDoFs += face->getPtrElementRight()->getNumberOfBasisFunctions(0);
        totalPDoFs += face->getPtrElementRight()->getNumberOfBasisFunctions(1);
        ElementInfos* rightInfo = static_cast<ElementInfos*>(face->getPtrElementRight()->getUserData());
        //}
        ret.resize(totalUDoFs + totalPDoFs, totalUDoFs + totalPDoFs);
        ElementInfos* leftInfo = static_cast<ElementInfos*>(face->getPtrElementLeft()->getUserData());


        LinearAlgebra::SmallVector<DIM> normal = fa.getUnitNormalVector();

        std::vector<std::size_t > indices (totalUDoFs);
        std::vector<double> phiNormalEpsilon (totalUDoFs);
        LinearAlgebra::SmallVector<DIM> phi;

        // Left element
        for (std::size_t i = 0; i < totalUDoFs; ++i)
        {
            indices[i] = i;
            fa.basisFunction(i, phi, 0);
            phiNormalEpsilon[i] = leftInfo->epsilon_ * (phi * normal);
        }
        // Right element
        for (std::size_t i = leftUDoFs; i < totalUDoFs; ++i)
        {
            indices[i] = i + leftPDoFs;
            fa.basisFunction(i, phi, 0);
            // - from the conversion between the left and right normal
            phiNormalEpsilon[i] = -rightInfo->epsilon_ * (phi * normal);
        }

        for (std::size_t i = 0; i < totalUDoFs; ++i)
        {
            const std::size_t iIndex = indices[i];
            for (std::size_t j = i; j < totalUDoFs; ++j)
            {
                const std::size_t jIndex = indices[j];
                double entry = stab2 * phiNormalEpsilon[i] * phiNormalEpsilon[j];
                ret(iIndex, jIndex) = entry;
                ret(jIndex, iIndex) = entry;
            }
        }
    }
    else
    {
        ret.resize(totalUDoFs + totalPDoFs, totalUDoFs + totalPDoFs);
    }
}

void DivDGMaxDiscretization::faceScalarVectorCoupling(
        Base::PhysicalFace<DIM> &fa,
        LinearAlgebra::MiddleSizeMatrix &ret) const
{
    const Base::Face* face = fa.getFace();
    std::size_t totalUDoFs = face->getPtrElementLeft()->getNumberOfBasisFunctions(0);
    std::size_t totalPDoFs = face->getPtrElementLeft()->getNumberOfBasisFunctions(1);

    const double epsilonLeft = static_cast<ElementInfos*>(face->getPtrElementLeft()->getUserData())->epsilon_;
    const double epsilonRight = face->isInternal()
            ? (static_cast<ElementInfos*>(face->getPtrElementRight()->getUserData())->epsilon_)
            : 0; // If no right face is present this will not be used.
    // From the averaging terms.
    double averageFactor = face->isInternal() ? 0.5 : 1;

    std::size_t leftUDofs = totalUDoFs;
    std::size_t leftPDofs = totalPDoFs;

    if (face->isInternal())
    {
        totalUDoFs += face->getPtrElementRight()->getNumberOfBasisFunctions(0);
        totalPDoFs += face->getPtrElementRight()->getNumberOfBasisFunctions(1);
    }
    ret.resize(totalUDoFs + totalPDoFs, totalUDoFs + totalPDoFs);

    LinearAlgebra::SmallVector<DIM> normal = fa.getUnitNormalVector();

    std::vector<std::size_t> indicesU (totalUDoFs);
    std::vector<std::size_t> indicesP (totalPDoFs);
    // Note: while the discretization_ say {{ phiU }} * [[phiP]], which leads to
    // terms (phiP_[lr] n_[lr]) * phiU (subscripts for the left/right element). We
    // could thus multiply phiP_[lr] by the correct normal and compute the
    // innerproduct later. For better performance we use that n_l = -n_r, thus
    // the terms can also be written as (s_[lr] phiP_[lr]) (n_l * phiU), where
    // s_[lr] is the sign (1 for l, -1 for r). Note that now both terms are
    // numbers instead of vectors.
    std::vector<double> phiUNormal (totalUDoFs);
    std::vector<double> phiP (totalPDoFs);
    LinearAlgebra::SmallVector<DIM> phiU;

    // Left element U
    for (std::size_t i = 0; i < leftUDofs; ++i)
    {
        indicesU[i] = i;
        fa.basisFunction(i, phiU, 0);
        phiUNormal[i] = phiU * normal * epsilonLeft;
    }
    // Right element U
    // Note that for a boundary element leftUDofs == totalUDoFs and this loop is skipped.
    for (std::size_t i = leftUDofs; i < totalUDoFs; ++i)
    {
        indicesU[i] = i + leftPDofs;
        fa.basisFunction(i, phiU, 0);
        phiUNormal[i] = phiU * normal * epsilonRight;
    }

    // Left element P
    for (std::size_t i = 0; i < leftPDofs; ++i)
    {
        indicesP[i] = leftUDofs + i;
        phiP[i] = fa.basisFunction(i, 1);
    }
    // Right element P, also skipped for boundary elements.
    for (std::size_t i = leftPDofs; i < totalPDoFs; ++i)
    {
        indicesP[i] = totalUDoFs + i;
        // - because the right normal is -1 times the left normal.
        phiP[i] = -fa.basisFunction(i, 1);
    }

    for (std::size_t i = 0; i < totalUDoFs; ++i)
    {
        const std::size_t iIndex = indicesU[i];
        const double& phiUi = phiUNormal[i];
        for (std::size_t j = 0; j < totalPDoFs; ++j)
        {
            const double entry = averageFactor * (phiUi * phiP[j]);
            const std::size_t jIndex = indicesP[j];
            ret(iIndex, jIndex) = entry;
            ret(jIndex, iIndex) = entry;
        }
    }
}

void DivDGMaxDiscretization::faceStiffnessScalarMatrix4(
        Base::PhysicalFace<DIM> &fa,
        LinearAlgebra::MiddleSizeMatrix &ret, double stab3) const
{
    //TODO: Cleanup.
    const Base::Face* face = fa.getFace();
    std::size_t totalUDoFs = face->getPtrElementLeft()->getNumberOfBasisFunctions(0);
    std::size_t totalPDoFs = face->getPtrElementLeft()->getNumberOfBasisFunctions(1);

    const std::size_t leftUDoFs = totalUDoFs;
    const std::size_t leftPDoFs = totalPDoFs;

    if (face->isInternal())
    {
        totalUDoFs += face->getPtrElementRight()->getNumberOfBasisFunctions(0);
        totalPDoFs += face->getPtrElementRight()->getNumberOfBasisFunctions(1);
    }
    ret.resize(totalUDoFs + totalPDoFs, totalUDoFs + totalPDoFs);

    std::vector<std::size_t> indices (totalPDoFs);
    std::vector<double> phiP (totalPDoFs);
    // Note we leave out the actual normal as nL * nL = +1, nR * nL = -1, etc.
    // instead we multiply the right functions by -1.

    // Left element
    for (std::size_t i = 0; i < leftPDoFs; ++i)
    {
        indices[i] = leftUDoFs + i;
        phiP[i] = fa.basisFunction(i, 1);
    }
    // Right element
    for (std::size_t i = leftPDoFs; i < totalPDoFs; ++i)
    {
        indices[i] = totalUDoFs + i;
        // - from the right normal being -1 times the left normal.
        phiP[i] = - fa.basisFunction(i, 1);
    }

    for (std::size_t i = 0; i < totalPDoFs; ++i)
    {
        const std::size_t iIndex = indices[i];
        const double& phiPi = phiP[i];
        for (std::size_t j = i; j < totalPDoFs; ++j)
        {
            const double entry = stab3 * phiP[j] * phiPi;
            const std::size_t jIndex = indices[j];
            ret(iIndex, jIndex) = entry;
            ret(jIndex, iIndex) = entry;
        }
    }
}

void DivDGMaxDiscretization::faceBoundaryVector(
        Base::PhysicalFace<DIM>& fa,
        const DivDGMaxDiscretization::FaceInputFunction &boundaryValue,
        LinearAlgebra::MiddleSizeVector &ret, double stab1) const
{
    const Base::Face* face = fa.getFace();
    const Geometry::PointReference<DIM - 1>& p = fa.getPointReference();


    if(face->isInternal())
    {
        std::size_t totalUDoFs = face->getPtrElementLeft()->getNumberOfBasisFunctions(0) + face->getPtrElementRight()->getNumberOfBasisFunctions(0);
        std::size_t totalPDoFs = face->getPtrElementLeft()->getNumberOfBasisFunctions(1) + face->getPtrElementRight()->getNumberOfBasisFunctions(1);

        ret.resize(totalUDoFs + totalPDoFs);

        for(std::size_t i = 0 ; i < (totalUDoFs + totalPDoFs); ++i)
        {
            ret(i) = 0;
        }
    }
    else
    {
        const PointElementReferenceT& PLeft = face->mapRefFaceToRefElemL(p);
        const PointPhysicalT PPhys = face->getPtrElementLeft()->referenceToPhysical(PLeft);

        LinearAlgebra::SmallVector<DIM> val, phi_curl;
        LinearAlgebra::SmallVector<DIM> phi;
        boundaryValue(PPhys, fa, val);

        std::size_t totalUDoFs = face->getPtrElementLeft()->getNumberOfBasisFunctions(0);
        std::size_t totalPDoFs = face->getPtrElementLeft()->getNumberOfBasisFunctions(1);
        ret.resize(totalUDoFs + totalPDoFs);

        for (std::size_t i = 0; i < totalUDoFs; ++i)
        {
            fa.basisFunctionUnitNormalCross(i, phi, 0);
            phi_curl = fa.basisFunctionCurl(i, 0);
            ret(i) = stab1 * (phi * val) - (phi_curl * val);
        }

    }
}

double DivDGMaxDiscretization::computeL2Error(
        Base::MeshManipulator<DIM> &mesh, std::size_t timeVector,
        const DivDGMaxDiscretization::InputFunction &electricField) const
{
    Integration::ElementIntegral<DIM> elIntegral(false);
    elIntegral.setTransformation(std::shared_ptr<Base::CoordinateTransformation<DIM>>(
            new Base::HCurlConformingTransformation<DIM>()), 0);
    elIntegral.setTransformation(std::shared_ptr<Base::CoordinateTransformation<DIM>>(
            new Base::H1ConformingTransformation<DIM>()), 1);

    double error = 0;
    for (Base::MeshManipulator<DIM>::ElementIterator it = mesh.elementColBegin(); it != mesh.elementColEnd(); ++it)
    {
        error += elIntegral.integrate<double> ((*it),
                [&](Base::PhysicalElement<DIM> &el) {
            return elementErrorIntegrand(el, timeVector, electricField);
        });
    }
    return std::sqrt(error);
}

double DivDGMaxDiscretization::elementErrorIntegrand(
        Base::PhysicalElement<DIM> &el, std::size_t timeVector,
        const DivDGMaxDiscretization::InputFunction &exactValues) const
{
    const Base::Element* element = el.getElement();
    const Geometry::PointPhysical<DIM>& pPhys = el.getPointPhysical();

    std::size_t numberOfUDoFs = element->getNumberOfBasisFunctions(0), numberOfPDoFs = element->getNumberOfBasisFunctions(1);
    LinearAlgebra::MiddleSizeVector data = element->getTimeIntegrationVector(timeVector);

    LinearAlgebra::SmallVector<DIM> error, phi;
    exactValues(pPhys, error);
    for (std::size_t i = 0; i < numberOfUDoFs; ++i)
    {
        el.basisFunction(i, phi, 0);
        error -= std::real(data[i]) * phi;
    }
    for (std::size_t i = 0; i < numberOfPDoFs; ++i)
    {
        error -= std::real(data[i + numberOfUDoFs]) * el.basisFunctionDeriv(i, 1);
    }

    return error.l2NormSquared();
}
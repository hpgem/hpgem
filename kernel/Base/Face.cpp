/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, University of Twente
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#include "Face.h"
#include "Element.h"

#include "Logger.h"
#include <iostream>
#include "Integration/QuadratureRules/GaussQuadratureRule.h"
#include "Geometry/ReferenceGeometry.h"
#include "Geometry/PointReference.h"
#include "LinearAlgebra/NumericalVector.h"
#include "L2Norm.h"
#include "FaceCacheData.h"
#include "ElementCacheData.h"
#include "Node.h"
#include "Geometry/PhysicalGeometry.h"
#include "Geometry/PointPhysical.h"

namespace Base
{
    
    class Face;
    Face::Face(Element* ptrElemL, const LocalFaceNrTypeT& localFaceNumL, Element* ptrElemR, const LocalFaceNrTypeT& localFaceNumR, std::size_t faceID, std::size_t numberOfFaceMatrixes, std::size_t numberOfFaceVectors)
            : FaceGeometryT((ElementGeometryT*) ptrElemL, localFaceNumL, (ElementGeometryT*) ptrElemR, localFaceNumR),
            FaceData(ptrElemL->getNrOfBasisFunctions() * ptrElemL->getNrOfUnknows() + ptrElemR->getNrOfBasisFunctions() * ptrElemR->getNrOfUnknows(), numberOfFaceMatrixes, numberOfFaceVectors), 
            elementLeft_(ptrElemL), elementRight_(ptrElemR), nrOfConformingDOFOnTheFace_(0), faceID_(faceID)
    {
        logger.assert(ptrElemL != nullptr, "Invalid element passed");
        logger.assert(ptrElemR != nullptr, "Error: passing a boundary face to the constructor for internal faces!");
        createQuadratureRules();
        ptrElemL->setFace(localFaceNumL, this);
        ptrElemR->setFace(localFaceNumR, this);
        
        std::vector<std::size_t> leftVertices, rightVertices;
        std::vector<std::size_t> localLeftVertices = ptrElemL->getPhysicalGeometry()->getLocalFaceNodeIndices(localFaceNumL);
        std::vector<std::size_t> localRightVertices = ptrElemR->getPhysicalGeometry()->getLocalFaceNodeIndices(localFaceNumR);
        for (std::size_t i = 0; i < getReferenceGeometry()->getNumberOfNodes(); ++i)
        {
            leftVertices.push_back(ptrElemL->getNode(localLeftVertices[i])->getID());
            rightVertices.push_back(ptrElemR->getNode(localRightVertices[i])->getID());
        }
        initialiseFaceToFaceMapIndex(leftVertices, rightVertices);
    }
    
    Face::Face(Element* ptrElemL, const LocalFaceNrTypeT& localFaceNumL, const Geometry::FaceType& faceType, std::size_t faceID, std::size_t numberOfFaceMatrixes, std::size_t numberOfFaceVectors)
            : FaceGeometryT((ElementGeometryT*) ptrElemL, localFaceNumL, faceType), FaceData(ptrElemL->getNrOfBasisFunctions() * ptrElemL->getNrOfUnknows(), numberOfFaceMatrixes, numberOfFaceVectors), elementLeft_(ptrElemL), elementRight_(nullptr), nrOfConformingDOFOnTheFace_(0), faceID_(faceID)
    {
        logger.assert(ptrElemL != nullptr, "Invalid element passed");
        createQuadratureRules();
        ptrElemL->setFace(localFaceNumL, this);
    }
    
    Face::Face(const Face& other, Element* elementL, const std::size_t localFaceL, Element* elementR, const std::size_t localFaceR)
        : FaceGeometry(other, elementL, localFaceL, elementR, localFaceR), 
        FaceData(other),
        elementLeft_(elementL), elementRight_(elementR),
        quadratureRule_(other.quadratureRule_), faceID_(other.faceID_),
        nrOfConformingDOFOnTheFace_(other.nrOfConformingDOFOnTheFace_)
    {        
        logger.assert(elementL != nullptr, "Invalid element passed");
        logger(DEBUG, "Coupling (left) face % to element %", faceID_, elementL->getID());
        elementL->setFace(localFaceL, this);
        logger.assert(elementL->getNrOfFaces() > 0, "Element does not contain any face!");
        if (elementR != nullptr)
        {
            elementR->setFace(localFaceR, this);
            logger(DEBUG, "Coupling (right) face % to element %", faceID_, elementR->getID());
        }
    }
    
    void Face::createQuadratureRules()
    {
        //order of quadrature rules:
        std::size_t rightOrder = (elementRight_ == nullptr ? 0 : elementRight_->getGaussQuadratureRule()->order());
        std::size_t leftOrder = elementLeft_->getGaussQuadratureRule()->order();
        if (leftOrder >= rightOrder)
        {
            quadratureRule_ = elementLeft_->getReferenceGeometry()->getCodim1ReferenceGeometry(FaceGeometryT::localFaceNumberLeft_)->getGaussQuadratureRule(leftOrder);
        }
        else
        {
            logger(DEBUG, "again..... Face<DIM>::createQuadratureRules(): % %.", leftOrder, rightOrder);
            quadratureRule_ = elementRight_->getReferenceGeometry()->getCodim1ReferenceGeometry(FaceGeometryT::localFaceNumberRight_)->getGaussQuadratureRule(rightOrder);
        }
        
    }
    
    std::size_t Face::getNrOfBasisFunctions() const
    {
        if (isInternal())
        {
            return getPtrElementLeft()->getNrOfBasisFunctions() + getPtrElementRight()->getNrOfBasisFunctions();
        }
        else
        {
            return getPtrElementLeft()->getNrOfBasisFunctions();
        }
    }
    
    double Face::basisFunction(std::size_t i, const Geometry::PointReference& p) const
    {
        logger.assert(i<getNrOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNrOfBasisFunctions());
        std::size_t n(getPtrElementLeft()->getNrOfBasisFunctions());
        if (i < n)
        {
            return getPtrElementLeft()->basisFunction(i, mapRefFaceToRefElemL(p));
        }
        else
        {
            return getPtrElementRight()->basisFunction(i - n, mapRefFaceToRefElemR(p));
        }
    }
    
    void Face::basisFunction(std::size_t i, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret) const
    {
        logger.assert(i<getNrOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNrOfBasisFunctions());
        Geometry::PointReference pElement(p.size() + 1);
        std::size_t n(getPtrElementLeft()->getNrOfBasisFunctions());
        if (i < n)
        {
            pElement = mapRefFaceToRefElemL(p);
            getPtrElementLeft()->basisFunction(i, pElement, ret);
        }
        else
        {
            pElement = mapRefFaceToRefElemR(p);
            getPtrElementRight()->basisFunction(i - n, pElement, ret);
        }
    }
    
    /// \param[in] iSide The index corresponding to the side of the face.
    /// \param[in] iBasisFunction The index corresponding to the basis function.
    /// \param[in] p The reference point on the reference element.
    double Face::basisFunction(Side iSide, std::size_t iBasisFunction, const Geometry::PointReference& p) const
    {
        if (iSide == Side::LEFT)
        {
            logger.assert(iBasisFunction < getPtrElementLeft()->getNrOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", iBasisFunction, getPtrElementLeft()->getNrOfBasisFunctions());
            return getPtrElementLeft()->basisFunction(iBasisFunction, mapRefFaceToRefElemL(p));
        }
        else
        {
            logger.assert(isInternal(), "boundary faces only have a \"left\" element");
            logger.assert(iBasisFunction < getPtrElementRight()->getNrOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", iBasisFunction, getPtrElementRight()->getNrOfBasisFunctions());
            return getPtrElementRight()->basisFunction(iBasisFunction, mapRefFaceToRefElemR(p));
        }
    }
    
    LinearAlgebra::NumericalVector Face::basisFunctionNormal(std::size_t i, const LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p) const
    {
        logger.assert(i<getNrOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNrOfBasisFunctions());
        LinearAlgebra::NumericalVector ret;
        Geometry::PointReference pElement(p.size() + 1);
        std::size_t n = getPtrElementLeft()->getNrOfBasisFunctions();
        if (i < n)
        {
            pElement = mapRefFaceToRefElemL(p);
            ret = normal;
            ret *= getPtrElementLeft()->basisFunction(i, pElement) / Base::L2Norm(normal);
        }
        else
        {
            pElement = mapRefFaceToRefElemR(p);
            ret = normal;
            ret *= -getPtrElementRight()->basisFunction(i - n, pElement) / Base::L2Norm(normal);
        }
        return ret;
    }
    
    /// \param[in] iSide The index corresponding to the side of the face.
    /// \param[in] iBasisFunction The index corresponding to the basis function.
    /// \param[in] normal The normal vector (pointing outwards with respect to the element on the left side).
    /// \param[in] p The reference point on the reference element.
    LinearAlgebra::NumericalVector Face::basisFunctionNormal(Side iSide, std::size_t iBasisFunction, const LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p) const
    {
        Geometry::PointReference pElement(p.size() + 1);
        if (iSide == Side::LEFT)
        {
            logger.assert(iBasisFunction < getPtrElementLeft()->getNrOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", iBasisFunction, getPtrElementLeft()->getNrOfBasisFunctions());
            pElement = mapRefFaceToRefElemL(p);
            return getPtrElementLeft()->basisFunction(iBasisFunction, pElement) * normal / Base::L2Norm(normal);
        }
        else
        {
            logger.assert(isInternal(), "boundary faces only have a \"left\" element");
            logger.assert(iBasisFunction < getPtrElementRight()->getNrOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", iBasisFunction, getPtrElementRight()->getNrOfBasisFunctions());
            pElement = mapRefFaceToRefElemR(p);
            return -getPtrElementRight()->basisFunction(iBasisFunction, pElement) * normal / Base::L2Norm(normal);
        }
    }
    
    double Face::basisFunctionDeriv(std::size_t i, std::size_t jDir, const Geometry::PointReference& p) const
    {
        logger.assert(i<getNrOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNrOfBasisFunctions());
        Geometry::PointReference pElement(p.size() + 1);
        std::size_t n = getPtrElementLeft()->getNrOfBasisFunctions();
        if (i < n)
        {
            pElement = mapRefFaceToRefElemL(p);
            return getPtrElementLeft()->basisFunctionDeriv(i, jDir, pElement);
        }
        else
        {
            pElement = mapRefFaceToRefElemR(p);
            return getPtrElementRight()->basisFunctionDeriv(i - n, jDir, pElement);
        }
    }
    
    LinearAlgebra::NumericalVector Face::basisFunctionDeriv(std::size_t i, const Geometry::PointReference& p) const
    {
        logger.assert(i<getNrOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNrOfBasisFunctions());
        Geometry::PointReference pElement(p.size() + 1);
        std::size_t n = getPtrElementLeft()->getNrOfBasisFunctions();
        if (i < n)
        {
            pElement = mapRefFaceToRefElemL(p);
            return getPtrElementLeft()->basisFunctionDeriv(i, pElement);
        }
        else
        {
            pElement = mapRefFaceToRefElemR(p);
            return getPtrElementRight()->basisFunctionDeriv(i - n, pElement);
        }
    }
    
    /// \param[in] iSide The index corresponding to the side of the face.
    /// \param[in] iBasisFunction The index corresponding to the basis function.
    /// \param[in] p The reference point on the reference element.
    LinearAlgebra::NumericalVector Face::basisFunctionDeriv(Side iSide, std::size_t iBasisFunction, const Geometry::PointReference& p) const
    {
        Geometry::PointReference pElement(p.size() + 1);
        if (iSide == Side::LEFT)
        {
            logger.assert(iBasisFunction < getPtrElementLeft()->getNrOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", iBasisFunction, getPtrElementLeft()->getNrOfBasisFunctions());
            pElement = mapRefFaceToRefElemL(p);
            return getPtrElementLeft()->basisFunctionDeriv(iBasisFunction, pElement);
        }
        else
        {
            logger.assert(isInternal(), "boundary faces only have a \"left\" element");
            logger.assert(iBasisFunction < getPtrElementRight()->getNrOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", iBasisFunction, getPtrElementRight()->getNrOfBasisFunctions());
            pElement = mapRefFaceToRefElemR(p);
            return getPtrElementRight()->basisFunctionDeriv(iBasisFunction, pElement);
        }
    }
    
    LinearAlgebra::NumericalVector Face::basisFunctionCurl(std::size_t i, const Geometry::PointReference& p) const
    {
        logger.assert(i<getNrOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNrOfBasisFunctions());
        Geometry::PointReference pElement(p.size() + 1);
        std::size_t numBasisFuncsLeft = getPtrElementLeft()->getNrOfBasisFunctions();
        if (i < numBasisFuncsLeft)
        {
            pElement = mapRefFaceToRefElemL(p);
            return getPtrElementLeft()->basisFunctionCurl(i, pElement);
        }
        else
        {
            pElement = mapRefFaceToRefElemR(p);
            return getPtrElementRight()->basisFunctionCurl(i - numBasisFuncsLeft, pElement);
        }
    }
    
    ///Get the time level data from both elements and concatenate them. 
    ///Note that we assume that the data is stored as column "vectors".
    LinearAlgebra::NumericalVector Face::getTimeLevelData(std::size_t timeLevel, std::size_t unknown) const
    {
        LinearAlgebra::NumericalVector resLeft = getPtrElementLeft()->getTimeLevelData(timeLevel, unknown);
        if (isInternal())
        {
            std::size_t numBasisFuncs = getNrOfBasisFunctions();
            std::size_t numBasisFuncsLeft = getPtrElementLeft()->getNrOfBasisFunctions();
            resLeft.resize(numBasisFuncs);
            LinearAlgebra::NumericalVector resRight = getPtrElementRight()->getTimeLevelData(timeLevel, unknown);
            for (std::size_t i = numBasisFuncsLeft; i < numBasisFuncs; ++i)
            {
                resLeft[i] = resRight[i - numBasisFuncsLeft];
            }
        }
        return resLeft;
    }
    
    LinearAlgebra::NumericalVector Face::getCurrentData()
    {
        LinearAlgebra::NumericalVector dataLeft = getPtrElementLeft()->getCurrentData();
        if (isInternal())
        {
            std::size_t numBasisFuncs = getNrOfBasisFunctions();
            std::size_t numBasisFuncsLeft = getPtrElementLeft()->getNrOfBasisFunctions();
            dataLeft.resize(numBasisFuncs);
            LinearAlgebra::NumericalVector dataRight = getPtrElementRight()->getCurrentData();
            for (std::size_t i = numBasisFuncsLeft; i < numBasisFuncs; ++i)
            {
                dataLeft[i] = dataRight[i - numBasisFuncsLeft];
            }
        }
        return dataLeft;
    }
    
    /// \param[in] side The side of the face.
    /// \param[in] varId The index corresponding to the variable.
    /// \param[in] scalarBasisFunctionId The index corresponding to the basisfunction.
    const std::size_t Face::convertToSingleIndex(Side side, std::size_t scalarBasisFunctionId, std::size_t varId) const
    {
        logger.assert(varId < getPtrElementLeft()->getNrOfUnknows(), "Asked for unknown %, but there are only % unknowns", varId, getPtrElementLeft()->getNrOfUnknows());
        if (side == Side::LEFT)
        {
            logger.assert(scalarBasisFunctionId < getPtrElementLeft()->getNrOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", scalarBasisFunctionId, getPtrElementLeft()->getNrOfBasisFunctions());
            return varId * getPtrElementLeft()->getNrOfBasisFunctions() + scalarBasisFunctionId;
        }
        else
        {
            logger.assert(isInternal(), "boundary faces only have a \"left\" element");
            logger.assert(scalarBasisFunctionId < getPtrElementLeft()->getNrOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", scalarBasisFunctionId, getPtrElementLeft()->getNrOfBasisFunctions());
            std::size_t nDOFLeft = getPtrElementLeft()->getNrOfUnknows() * getPtrElementLeft()->getNrOfBasisFunctions();
            return nDOFLeft + varId * getPtrElementRight()->getNrOfBasisFunctions() + scalarBasisFunctionId;
        }
    }
    
    const Side Face::getSide(std::size_t faceBasisFunctionId) const
    {
        std::size_t nDOFLeft = getPtrElementLeft()->getNrOfUnknows() * getPtrElementLeft()->getNrOfBasisFunctions();
        if(faceBasisFunctionId < nDOFLeft)
        {
            return Side::LEFT;
        }
        else
        {
            logger.assert(faceBasisFunctionId < nDOFLeft + (isInternal() ? getPtrElementRight()->getNrOfUnknows() * getPtrElementRight()->getNrOfBasisFunctions() : 0), "The index for the face basis (vector)function (%) is larger than the number of basis (vector)functions at the adjacent elements (%)", faceBasisFunctionId, nDOFLeft + (isInternal() ? getPtrElementRight()->getNrOfUnknows() * getPtrElementRight()->getNrOfBasisFunctions() : 0));
            return Side::RIGHT;
        }
    }
    
    const std::size_t Face::getElementBasisFunctionId(std::size_t faceBasisFunctionId) const
    {
        std::size_t nDOFLeft = getPtrElementLeft()->getNrOfUnknows() * getPtrElementLeft()->getNrOfBasisFunctions();
        if(faceBasisFunctionId < nDOFLeft)
        {
            return faceBasisFunctionId;
        }
        else
        {
            logger.assert(faceBasisFunctionId < nDOFLeft + (isInternal() ? getPtrElementRight()->getNrOfUnknows() * getPtrElementRight()->getNrOfBasisFunctions() : 0), "The index for the face basis (vector)function (%) is larger than the number of basis (vector)functions at the adjacent elements (%)", faceBasisFunctionId, nDOFLeft + (isInternal() ? getPtrElementRight()->getNrOfUnknows() * getPtrElementRight()->getNrOfBasisFunctions() : 0));
            return faceBasisFunctionId - nDOFLeft;
        }
    }
}
;

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
    
    /// \details The user does not need to worry about the contruction of faces. This is done by mesh-generators. For example the interface HpgemAPIBase can be used to create meshes.
    Face::Face(Element* ptrElemL, const LocalFaceNrTypeT& localFaceNumL, Element* ptrElemR, const LocalFaceNrTypeT& localFaceNumR, std::size_t faceID, std::size_t numberOfFaceMatrixes, std::size_t numberOfFaceVectors)
            : FaceGeometryT(ptrElemL, localFaceNumL, ptrElemR, localFaceNumR),
            FaceData(ptrElemL->getNrOfBasisFunctions() * ptrElemL->getNrOfUnknows() + ptrElemR->getNrOfBasisFunctions() * ptrElemR->getNrOfUnknows(), numberOfFaceMatrixes, numberOfFaceVectors), 
            elementLeft_(ptrElemL), elementRight_(ptrElemR), nrOfConformingDOFOnTheFace_(0), faceID_(faceID)
    {
        logger.assert(ptrElemL != nullptr, "Invalid element passed");
        logger.assert(ptrElemR != nullptr, "Error: passing a boundary face to the constructor for internal faces!");
        createQuadratureRules();
        ptrElemL->setFace(localFaceNumL, this);
        ptrElemR->setFace(localFaceNumR, this);
        
        std::vector<std::size_t> leftNodes, rightNodes;
        std::vector<std::size_t> localLeftNodes = ptrElemL->getPhysicalGeometry()->getLocalFaceNodeIndices(localFaceNumL);
        std::vector<std::size_t> localRightNodes = ptrElemR->getPhysicalGeometry()->getLocalFaceNodeIndices(localFaceNumR);
        for (std::size_t i = 0; i < getReferenceGeometry()->getNumberOfNodes(); ++i)
        {
            leftNodes.push_back(ptrElemL->getNode(localLeftNodes[i])->getID());
            rightNodes.push_back(ptrElemR->getNode(localRightNodes[i])->getID());
        }
        initialiseFaceToFaceMapIndex(leftNodes, rightNodes);
    }
    
    Face::Face(Element* ptrElemL, const LocalFaceNrTypeT& localFaceNumL, const Geometry::FaceType& faceType, std::size_t faceID, std::size_t numberOfFaceMatrixes, std::size_t numberOfFaceVectors)
            : FaceGeometryT(ptrElemL, localFaceNumL, faceType), FaceData(ptrElemL->getNrOfBasisFunctions() * ptrElemL->getNrOfUnknows(), numberOfFaceMatrixes, numberOfFaceVectors), elementLeft_(ptrElemL), elementRight_(nullptr), nrOfConformingDOFOnTheFace_(0), faceID_(faceID)
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
        std::size_t numBasisFuncs = getPtrElementLeft()->getNrOfBasisFunctions();
        if (i < numBasisFuncs)
        {
            return getPtrElementLeft()->basisFunction(i, mapRefFaceToRefElemL(p));
        }
        else
        {
            return getPtrElementRight()->basisFunction(i - numBasisFuncs, mapRefFaceToRefElemR(p));
        }
    }
    
    void Face::basisFunction(std::size_t i, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret) const
    {
        logger.assert(i<getNrOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNrOfBasisFunctions());
        std::size_t n(getPtrElementLeft()->getNrOfBasisFunctions());
        if (i < n)
        {
            getPtrElementLeft()->basisFunction(i, mapRefFaceToRefElemL(p), ret);
        }
        else
        {
            getPtrElementRight()->basisFunction(i - n, mapRefFaceToRefElemR(p), ret);
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
        std::size_t n = getPtrElementLeft()->getNrOfBasisFunctions();
        if (i < n)
        {
            ret = normal;
            ret *= getPtrElementLeft()->basisFunction(i, mapRefFaceToRefElemL(p)) / Base::L2Norm(normal);
        }
        else
        {
            ret = normal;
            ret *= -getPtrElementRight()->basisFunction(i - n, mapRefFaceToRefElemR(p)) / Base::L2Norm(normal);
        }
        return ret;
    }
    
    /// \param[in] iSide The index corresponding to the side of the face.
    /// \param[in] iBasisFunction The index corresponding to the basis function.
    /// \param[in] normal The normal vector (pointing outwards with respect to the element on the left side).
    /// \param[in] p The reference point on the reference element.
    LinearAlgebra::NumericalVector Face::basisFunctionNormal(Side iSide, std::size_t iBasisFunction, const LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p) const
    {
        if (iSide == Side::LEFT)
        {
            logger.assert(iBasisFunction < getPtrElementLeft()->getNrOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", iBasisFunction, getPtrElementLeft()->getNrOfBasisFunctions());
            return getPtrElementLeft()->basisFunction(iBasisFunction, mapRefFaceToRefElemL(p)) * normal / Base::L2Norm(normal);
        }
        else
        {
            logger.assert(isInternal(), "boundary faces only have a \"left\" element");
            logger.assert(iBasisFunction < getPtrElementRight()->getNrOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", iBasisFunction, getPtrElementRight()->getNrOfBasisFunctions());
            return -getPtrElementRight()->basisFunction(iBasisFunction, mapRefFaceToRefElemR(p)) * normal / Base::L2Norm(normal);
        }
    }
    
    double Face::basisFunctionDeriv(std::size_t i, std::size_t jDir, const Geometry::PointReference& p) const
    {
        logger.assert(i<getNrOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNrOfBasisFunctions());
        std::size_t n = getPtrElementLeft()->getNrOfBasisFunctions();
        if (i < n)
        {
            return getPtrElementLeft()->basisFunctionDeriv(i, jDir, mapRefFaceToRefElemL(p));
        }
        else
        {
            return getPtrElementRight()->basisFunctionDeriv(i - n, jDir, mapRefFaceToRefElemR(p));
        }
    }
    
    LinearAlgebra::NumericalVector Face::basisFunctionDeriv(std::size_t i, const Geometry::PointReference& p) const
    {
        logger.assert(i<getNrOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNrOfBasisFunctions());
        std::size_t n = getPtrElementLeft()->getNrOfBasisFunctions();
        if (i < n)
        {
            return getPtrElementLeft()->basisFunctionDeriv(i, mapRefFaceToRefElemL(p));
        }
        else
        {
            return getPtrElementRight()->basisFunctionDeriv(i - n, mapRefFaceToRefElemR(p));
        }
    }
    
    /// \param[in] iSide The index corresponding to the side of the face.
    /// \param[in] iBasisFunction The index corresponding to the basis function.
    /// \param[in] p The reference point on the reference element.
    LinearAlgebra::NumericalVector Face::basisFunctionDeriv(Side iSide, std::size_t iBasisFunction, const Geometry::PointReference& p) const
    {
        if (iSide == Side::LEFT)
        {
            logger.assert(iBasisFunction < getPtrElementLeft()->getNrOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", iBasisFunction, getPtrElementLeft()->getNrOfBasisFunctions());
            return getPtrElementLeft()->basisFunctionDeriv(iBasisFunction, mapRefFaceToRefElemL(p));
        }
        else
        {
            logger.assert(isInternal(), "boundary faces only have a \"left\" element");
            logger.assert(iBasisFunction < getPtrElementRight()->getNrOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", iBasisFunction, getPtrElementRight()->getNrOfBasisFunctions());
            return getPtrElementRight()->basisFunctionDeriv(iBasisFunction, mapRefFaceToRefElemR(p));
        }
    }
    
    LinearAlgebra::NumericalVector Face::basisFunctionCurl(std::size_t i, const Geometry::PointReference& p) const
    {
        logger.assert(i<getNrOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNrOfBasisFunctions());
        std::size_t numBasisFuncsLeft = getPtrElementLeft()->getNrOfBasisFunctions();
        if (i < numBasisFuncsLeft)
        {
            return getPtrElementLeft()->basisFunctionCurl(i, mapRefFaceToRefElemL(p));
        }
        else
        {
            return getPtrElementRight()->basisFunctionCurl(i - numBasisFuncsLeft, mapRefFaceToRefElemR(p));
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
    
    /// \param[in] side The side of the face.
    /// \param[in] varId The index corresponding to the variable.
    /// \param[in] scalarBasisFunctionId The index corresponding to the basisfunction.
    std::size_t Face::convertToSingleIndex(Side side, std::size_t scalarBasisFunctionId, std::size_t varId) const
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
    
    Side Face::getSide(std::size_t faceBasisFunctionId) const
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
    
    std::size_t Face::getElementBasisFunctionId(std::size_t faceBasisFunctionId) const
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

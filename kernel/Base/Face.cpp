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
    Face::Face(ElementT* ptrElemL, const LocalFaceNrTypeT& localFaceNumL, ElementT* ptrElemR, const LocalFaceNrTypeT& localFaceNumR, std::size_t faceID, std::size_t numberOfFaceMatrixes, std::size_t numberOfFaceVectors)
            : FaceGeometryT((ElementGeometryT*) ptrElemL, localFaceNumL, (ElementGeometryT*) ptrElemR, localFaceNumR), elementLeft_(ptrElemL), elementRight_(ptrElemR), faceID_(faceID), nrOfConformingDOFOnTheFace_(0), FaceData(ptrElemL->getNrOfBasisFunctions() * ptrElemL->getNrOfUnknows() + ptrElemR->getNrOfBasisFunctions() * ptrElemR->getNrOfUnknows(), numberOfFaceMatrixes, numberOfFaceVectors)
    {
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
    Face::Face(ElementT* ptrElemL, const LocalFaceNrTypeT& localFaceNumL, const Geometry::FaceType& faceType, std::size_t faceID, std::size_t numberOfFaceMatrixes, std::size_t numberOfFaceVectors)
            : FaceGeometryT((ElementGeometryT*) ptrElemL, localFaceNumL, faceType), elementLeft_(ptrElemL), elementRight_(nullptr), faceID_(faceID), nrOfConformingDOFOnTheFace_(0), FaceData(ptrElemL->getNrOfBasisFunctions() * ptrElemL->getNrOfUnknows(), numberOfFaceMatrixes, numberOfFaceVectors)
    {
        createQuadratureRules();
        ptrElemL->setFace(localFaceNumL, this);
    }
    
    void Face::createQuadratureRules()
    {
        std::size_t rightOrder = (elementRight_ == nullptr ? 0 : elementRight_->getGaussQuadratureRule()->order());
        std::size_t leftOrder = elementLeft_->getGaussQuadratureRule()->order();
        if (leftOrder >= rightOrder)
        {
            quadratureRule_ = elementLeft_->getReferenceGeometry()->getCodim1ReferenceGeometry(FaceGeometryT::localFaceNumberLeft_)->getGaussQuadratureRule(leftOrder);
        }
        else
        {
            std::cout << "again..... Face<DIM>::createQuadratureRules(): " << leftOrder << " " << rightOrder << std::endl;
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
            return getPtrElementLeft()->basisFunction(iBasisFunction, mapRefFaceToRefElemL(p));
        }
        else
        {
            return getPtrElementRight()->basisFunction(iBasisFunction, mapRefFaceToRefElemR(p));
        }
    }
    
    LinearAlgebra::NumericalVector Face::basisFunctionNormal(std::size_t i, const LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p) const
    {
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
    /// \param[in] normal The normal vector (pointing outwards with respect to the elemeent on the left side).
    /// \param[in] p The reference point on the reference element.
    LinearAlgebra::NumericalVector Face::basisFunctionNormal(Side iSide, std::size_t iBasisFunction, const LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p) const
    {
        Geometry::PointReference pElement(p.size() + 1);
        if (iSide == Side::LEFT)
        {
            pElement = mapRefFaceToRefElemL(p);
            return getPtrElementLeft()->basisFunction(iBasisFunction, pElement) * normal / Base::L2Norm(normal);
        }
        else
        {
            pElement = mapRefFaceToRefElemR(p);
            return -getPtrElementRight()->basisFunction(iBasisFunction, pElement) * normal / Base::L2Norm(normal);
        }
    }
    
    double Face::basisFunctionDeriv(std::size_t i, std::size_t jDir, const Geometry::PointReference& p) const
    {
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
            pElement = mapRefFaceToRefElemL(p);
            return getPtrElementLeft()->basisFunctionDeriv(iBasisFunction, pElement);
        }
        else
        {
            pElement = mapRefFaceToRefElemR(p);
            return getPtrElementRight()->basisFunctionDeriv(iBasisFunction, pElement);
        }
    }
    
    LinearAlgebra::NumericalVector Face::basisFunctionCurl(std::size_t i, const Geometry::PointReference& p) const
    {
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
    LinearAlgebra::NumericalVector Face::getTimeLevelData(std::size_t timeLevel)
    {
        LinearAlgebra::NumericalVector resLeft = getPtrElementLeft()->getTimeLevelData(timeLevel);
        if (isInternal())
        {
            std::size_t numBasisFuncs = getNrOfBasisFunctions();
            std::size_t numBasisFuncsLeft = getPtrElementLeft()->getNrOfBasisFunctions();
            resLeft.resize(numBasisFuncs);
            LinearAlgebra::NumericalVector resRight = getPtrElementRight()->getTimeLevelData(timeLevel);
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
    
    /// \param[in] iSide The index corresponding to the side of the face.
    /// \param[in] iVar The index corresponding to the variable.
    /// \param[in] iBasisFunction The index corresponding to the basisfunction.
    const std::size_t Face::convertToSingleIndex(Side iSide, std::size_t iBasisFunction, std::size_t iVar) const
    {
        if (iSide == Side::LEFT)
        {
            return iVar * getPtrElementLeft()->getNrOfBasisFunctions() + iBasisFunction;
        }
        else
        {
            std::size_t nDOFLeft = getPtrElementLeft()->getNrOfUnknows() * getPtrElementLeft()->getNrOfBasisFunctions();
            return nDOFLeft + iVar * getPtrElementRight()->getNrOfBasisFunctions() + iBasisFunction;
        }
    }
}
;

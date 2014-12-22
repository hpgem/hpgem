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
#include "Face.hpp"
#include "Element.hpp"

#include "TestErrorDebug.hpp"
#include <iostream>
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"
#include "Geometry/ReferenceGeometry.hpp"
#include "Geometry/PointReference.hpp"
#include "LinearAlgebra/NumericalVector.hpp"
#include "L2Norm.hpp"
#include "FaceCacheData.hpp"
#include "ElementCacheData.hpp"
#include "Node.hpp"
#include "Geometry/PhysicalGeometry.hpp"
#include "Geometry/PointPhysical.hpp"

namespace Base
{

    class Face;
    Face::Face(ElementT* ptrElemL, const LocalFaceNrTypeT& localFaceNumL,
               ElementT* ptrElemR, const LocalFaceNrTypeT& localFaceNumR,
               size_t faceID, size_t numberOfFaceMatrixes, size_t numberOfFaceVectors) :
    FaceGeometryT((ElementGeometryT*) ptrElemL, localFaceNumL, (ElementGeometryT*) ptrElemR, localFaceNumR),
    elementLeft_(ptrElemL),
    elementRight_(ptrElemR),
    faceID_(faceID),
    nrOfConformingDOFOnTheFace_(0),
    FaceData(ptrElemL->getNrOfBasisFunctions() * ptrElemL->getNrOfUnknows() +
             ptrElemR->getNrOfBasisFunctions() * ptrElemR->getNrOfUnknows(), numberOfFaceMatrixes, numberOfFaceVectors)
    {
        TestErrorDebug(ptrElemR != nullptr, "Error: passing a boundary face to the constructor for internal faces!");
        createQuadratureRules();
        ptrElemL->setFace(localFaceNumL, this);
        ptrElemR->setFace(localFaceNumR, this);

        std::vector<unsigned int> leftVertices, rightVertices;
        std::vector<unsigned int> localLeftVertices, localRightVertices;
        ptrElemL->getPhysicalGeometry()->getLocalFaceNodeIndices(localFaceNumL, localLeftVertices);
        ptrElemR->getPhysicalGeometry()->getLocalFaceNodeIndices(localFaceNumR, localRightVertices);
        for (size_t i = 0; i < getReferenceGeometry()->getNumberOfNodes(); ++i)
        {
            leftVertices.push_back(ptrElemL->getNode(localLeftVertices[i])->getID());
            rightVertices.push_back(ptrElemR->getNode(localRightVertices[i])->getID());
        }
        initialiseFaceToFaceMapIndex(leftVertices, rightVertices);
    }
    Face::Face(ElementT* ptrElemL, const LocalFaceNrTypeT& localFaceNumL, const Geometry::FaceType& faceType,
               int faceID, size_t numberOfFaceMatrixes, size_t numberOfFaceVectors) :
    FaceGeometryT((ElementGeometryT*) ptrElemL, localFaceNumL, faceType),
    elementLeft_(ptrElemL),
    elementRight_(nullptr),
    faceID_(faceID),
    nrOfConformingDOFOnTheFace_(0),
    FaceData(ptrElemL->getNrOfBasisFunctions() * ptrElemL->getNrOfUnknows(), numberOfFaceMatrixes, numberOfFaceVectors)
    {
        createQuadratureRules();
        ptrElemL->setFace(localFaceNumL, this);
    }
    
    void Face::createQuadratureRules()
    {
        unsigned int rightOrder = (elementRight_ == NULL ? 0 : elementRight_->getGaussQuadratureRule()->order());
        unsigned int leftOrder = elementLeft_->getGaussQuadratureRule()->order();
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
    
    int Face::getNrOfBasisFunctions() const
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
    
    double Face::basisFunction(size_t i, const Geometry::PointReference& p) const
    {
        Geometry::PointReference pElement(p.size() + 1);
        int n(getPtrElementLeft()->getNrOfBasisFunctions());
        if (i < n)
        {
            mapRefFaceToRefElemL(p, pElement);
            return getPtrElementLeft()->basisFunction(i, pElement);
        }
        else
        {
            mapRefFaceToRefElemR(p, pElement);
            return getPtrElementRight()->basisFunction(i - n, pElement);
        }
    }
    
    void Face::basisFunction(size_t i, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret) const
    {
        Geometry::PointReference pElement(p.size() + 1);
        int n(getPtrElementLeft()->getNrOfBasisFunctions());
        if (i < n)
        {
            mapRefFaceToRefElemL(p, pElement);
            getPtrElementLeft()->basisFunction(i, pElement, ret);
        }
        else
        {
            mapRefFaceToRefElemR(p, pElement);
            getPtrElementRight()->basisFunction(i - n, pElement, ret);
        }
    }
    
    void Face::basisFunctionNormal(size_t i, const LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret) const
    {
        Geometry::PointReference pElement(p.size() + 1);
        size_t n = getPtrElementLeft()->getNrOfBasisFunctions();
        if (i < n)
        {
            mapRefFaceToRefElemL(p, pElement);
            ret = normal;
            ret *= getPtrElementLeft()->basisFunction(i, pElement) / Base::L2Norm(normal);
        }
        else
        {
            mapRefFaceToRefElemR(p, pElement);
            ret = normal;
            ret *= -getPtrElementRight()->basisFunction(i - n, pElement) / Base::L2Norm(normal);
        }
    }
    
    double Face::basisFunctionDeriv(size_t i, size_t jDir, const Geometry::PointReference& p) const
    {
        Geometry::PointReference pElement(p.size() + 1);
        size_t n = getPtrElementLeft()->getNrOfBasisFunctions();
        if (i < n)
        {
            mapRefFaceToRefElemL(p, pElement);
            return getPtrElementLeft()->basisFunctionDeriv(i, jDir, pElement);
        }
        else
        {
            mapRefFaceToRefElemR(p, pElement);
            return getPtrElementRight()->basisFunctionDeriv(i - n, jDir, pElement);
        }
    }
    
    void Face::basisFunctionDeriv(size_t i, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret) const
    {
        Geometry::PointReference pElement(p.size() + 1);
        size_t n = getPtrElementLeft()->getNrOfBasisFunctions();
        if (i < n)
        {
            mapRefFaceToRefElemL(p, pElement);
            getPtrElementLeft()->basisFunctionDeriv(i, pElement, ret);
        }
        else
        {
            mapRefFaceToRefElemR(p, pElement);
            getPtrElementRight()->basisFunctionDeriv(i - n, pElement, ret);
        }
    }
    
    void Face::basisFunctionCurl(size_t i, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret) const
    {
        Geometry::PointReference pElement(p.size() + 1);
        size_t numBasisFuncsLeft = getPtrElementLeft()->getNrOfBasisFunctions();
        if (i < numBasisFuncsLeft)
        {
            mapRefFaceToRefElemL(p, pElement);
            getPtrElementLeft()->basisFunctionCurl(i, pElement, ret);
        }
        else
        {
            mapRefFaceToRefElemR(p, pElement);
            getPtrElementRight()->basisFunctionCurl(i - numBasisFuncsLeft, pElement, ret);
        }
    }
    
    ///Get the time level data from both elements and concatenate them. 
    ///Note that we assume that the data is stored as column "vectors".
    LinearAlgebra::NumericalVector Face::getTimeLevelData(size_t timeLevel)
    {
        LinearAlgebra::NumericalVector resLeft = getPtrElementLeft()->getTimeLevelData(timeLevel);
        if (isInternal())
        {
            size_t numBasisFuncs = getNrOfBasisFunctions();
            size_t numBasisFuncsLeft = getPtrElementLeft()->getNrOfBasisFunctions();
            resLeft.resize(numBasisFuncs);
            LinearAlgebra::NumericalVector resRight = getPtrElementRight()->getTimeLevelData(timeLevel);
            for (size_t i = numBasisFuncsLeft; i < numBasisFuncs; ++i)
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
            size_t numBasisFuncs = getNrOfBasisFunctions();
            size_t numBasisFuncsLeft = getPtrElementLeft()->getNrOfBasisFunctions();
            dataLeft.resize(numBasisFuncs);
            LinearAlgebra::NumericalVector dataRight = getPtrElementRight()->getCurrentData();
            for (size_t i = numBasisFuncsLeft; i < numBasisFuncs; ++i)
            {
                dataLeft[i] = dataRight[i - numBasisFuncsLeft];
            }
        }
        return dataLeft;
    }
};

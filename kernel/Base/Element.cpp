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

#include "Element.h"
#include "PhysGradientOfBasisFunction.h"
#include "Edge.h"
#include "Face.h"

#include "BasisFunctionSet.h"
#include "Geometry/PhysicalGeometry.h"
#include "LinearAlgebra/NumericalVector.h"
#include "Geometry/PointPhysical.h"
#include "FaceCacheData.h"
#include "ElementCacheData.h"
#include "Geometry/ReferenceGeometry.h"
#include "Geometry/PointReference.h"
#include "Logger.h"
#include "Node.h"
#include "Integration/QuadratureRules/GaussQuadratureRule.h"
#include "Geometry/Jacobian.h"

#include <limits>
#include <algorithm>

namespace Base
{
    
    Element::Element(const VectorOfPointIndexesT& globalNodeIndexes, 
                     const CollectionOfBasisFunctionSets *basisFunctionSet, 
                     VectorOfPhysicalPointsT& allNodes, 
                     std::size_t nrOfUnkowns, 
                     std::size_t nrOfTimeLevels, 
                     std::size_t nrOfBasisFunc, 
                     std::size_t id, 
                     std::size_t numberOfElementMatrixes, 
                     std::size_t numberOfElementVectors, 
                     const std::vector<int>& basisFunctionSetPositions)
            : ElementGeometryT(globalNodeIndexes, allNodes), 
        ElementData(nrOfTimeLevels, nrOfUnkowns, nrOfBasisFunc, numberOfElementMatrixes, numberOfElementVectors), 
        quadratureRule_(nullptr), basisFunctionSet_(basisFunctionSet), 
        vecCacheData_(), id_(id), basisFunctionSetPositions_(basisFunctionSetPositions)
    {
        logger.assert(basisFunctionSet!=nullptr, "Invalid basis function set passed");
        logger.assert(basisFunctionSet->size()>0, "Not enough basis function sets passed");
        logger(VERBOSE, "numberOfElementMatrixes: %", numberOfElementMatrixes);
        logger(VERBOSE, "numberOfElementVectors: %", numberOfElementVectors);
        
        orderCoeff_ = 2; // for safety
        std::size_t numberOfBasisFunctions = 0;
        for (std::size_t i = 0; i < basisFunctionSetPositions_.size(); ++i)
        {
            //basisFunctionSetPositions_ may be set to the special value -1 for the empty set, so this must be an integer comparison
            logger.assert(basisFunctionSetPositions_[i]<static_cast<int>(basisFunctionSet->size()), "Not enough basis function sets passed");
            logger.assert(basisFunctionSetPositions_[i] == -1 || basisFunctionSet->at(basisFunctionSetPositions_[i])!=nullptr, "Invalid basis function set passed");
            if(basisFunctionSetPositions_[i] != -1)
            {
                numberOfBasisFunctions += basisFunctionSet_->at(basisFunctionSetPositions_[i])->size();
            }
        }
        logger.assert(nrOfBasisFunc==numberOfBasisFunctions, "Redundant argument set to the wrong value");
        setNumberOfBasisFunctions(numberOfBasisFunctions);
        setQuadratureRulesWithOrder(orderCoeff_ * basisFunctionSet_->at(basisFunctionSetPositions_[0])->getOrder() + 1);
        nrOfDOFinTheElement_ = basisFunctionSet_->at(basisFunctionSetPositions_[0])->size();
        facesList_.assign(getReferenceGeometry()->getNrOfCodim1Entities(), nullptr);
        if (getReferenceGeometry()->getNrOfCodim3Entities() > 0)
        {
            edgesList_.assign(getReferenceGeometry()->getNrOfCodim2Entities(), nullptr);
        }
        nodesList_.assign(getReferenceGeometry()->getNumberOfNodes(), nullptr);
    }
        
    Element::Element(const ElementData& otherData, const ElementGeometry& otherGeometry)
        : ElementData(otherData),
        ElementGeometry(otherGeometry)
    {        
    }
    
    ///Very ugly default constructor that's only here because it is needed in
    ///ShortTermStorageElementBase.
    Element::Element()
            : ElementData(0, 0, 0, 0, 0)
    {
    }
    
    Element::~Element()
    {        
    }
    
    Element* Element::copyWithoutFacesEdgesNodes()
    {
        //Make a new element with the data and geometry of this element
        Element* newElement = new Element(*this, *this);
        
        //copy the pointers to singletons
        newElement->quadratureRule_ = quadratureRule_;
        newElement->basisFunctionSet_ = basisFunctionSet_;
        
        //copy other data
        newElement->basisFunctionSetPositions_ = basisFunctionSetPositions_;
        newElement->id_ = id_;
        newElement->nrOfDOFinTheElement_ = nrOfDOFinTheElement_;
        newElement->orderCoeff_ = orderCoeff_;
        newElement->vecCacheData_ = vecCacheData_;
        
        //allocate memory for nodesList, facesList and edgesList
        newElement->nodesList_.resize(getNrOfNodes());
        newElement->edgesList_.resize(getNrOfEdges());
        newElement->facesList_.resize(getNrOfFaces());
        
        return newElement;
    }
    
    void Element::setDefaultBasisFunctionSet(std::size_t position)
    {
        logger.assert(position < basisFunctionSet_->size(), "Not enough basis function sets passed");
        basisFunctionSetPositions_.resize(1, -1);
        basisFunctionSetPositions_[0] = position;
        std::size_t numberOfBasisFunctions(0);
        for (int i : basisFunctionSetPositions_)
        {
            if (i != -1)
                numberOfBasisFunctions += basisFunctionSet_->at(i)->size();
        }
        setNumberOfBasisFunctions(numberOfBasisFunctions);
        setQuadratureRulesWithOrder(orderCoeff_ * basisFunctionSet_->at(position)->getOrder() + 1);
        nrOfDOFinTheElement_ = basisFunctionSet_->at(position)->size();
    }
    
    void Element::setFaceBasisFunctionSet(std::size_t position, std::size_t localFaceIndex)
    {
        logger.assert(position < basisFunctionSet_->size(), "Not enough basis function sets passed");
        logger.assert(localFaceIndex < getNrOfFaces(), "Asked for face %, but there are only % faces", localFaceIndex, getNrOfFaces());
        if (basisFunctionSetPositions_.size() < 1 + getNrOfFaces())
        {
            basisFunctionSetPositions_.resize(1 + getNrOfFaces(), -1);
        }
        basisFunctionSetPositions_[1 + localFaceIndex] = position;
        std::size_t numberOfBasisFunctions(0);
        for (int i : basisFunctionSetPositions_)
        {
            if (i != -1)
                numberOfBasisFunctions += basisFunctionSet_->at(i)->size();
        }
        setNumberOfBasisFunctions(numberOfBasisFunctions);
    }
    
    void Element::setEdgeBasisFunctionSet(std::size_t position, std::size_t localEdgeIndex)
    {
        logger.assert(position < basisFunctionSet_->size(), "Not enough basis function sets passed");
        logger.assert(localEdgeIndex < getNrOfEdges(), "Asked for edge %, but there are only % edges", localEdgeIndex, getNrOfEdges());
        if (basisFunctionSetPositions_.size() < 1 + getNrOfFaces() + getNrOfEdges())
        {
            basisFunctionSetPositions_.resize(1 + getNrOfFaces() + getNrOfEdges(), -1);
        }
        basisFunctionSetPositions_[1 + getNrOfFaces() + localEdgeIndex] = position;
        std::size_t numberOfBasisFunctions(0);
        for (int i : basisFunctionSetPositions_)
        {
            if (i != -1)
                numberOfBasisFunctions += basisFunctionSet_->at(i)->size();
        }
        setNumberOfBasisFunctions(numberOfBasisFunctions);
    }
    
    void Element::setVertexBasisFunctionSet(std::size_t position, std::size_t localVertexIndex)
    {
        logger.assert(position < basisFunctionSet_->size(), "Not enough basis function sets passed");
        logger.assert(localVertexIndex < getNrOfNodes(), "Asked for node %, but there are only % nodes", localVertexIndex, getNrOfNodes());
        if (basisFunctionSetPositions_.size() < 1 + getNrOfFaces() + getNrOfEdges() + getNrOfNodes())
        {
            basisFunctionSetPositions_.resize(1 + getNrOfFaces() + getNrOfEdges() + getNrOfNodes(), -1);
        }
        basisFunctionSetPositions_[1 + getNrOfFaces() + getNrOfEdges() + localVertexIndex] = position;
        std::size_t numberOfBasisFunctions(0);
        for (int i : basisFunctionSetPositions_)
        {
            if (i != -1)
                numberOfBasisFunctions += basisFunctionSet_->at(i)->size();
        }
        setNumberOfBasisFunctions(numberOfBasisFunctions);
    }
    
    double Element::basisFunctionDeriv(std::size_t i, std::size_t jDir, const PointReferenceT& p) const
    {
        logger.assert(i<getNrOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNrOfBasisFunctions());
        logger.assert((jDir < p.size()), "Error in BasisFunctionSet.EvalDeriv: invalid derivative direction!");
        
        int basePosition(0);
        for (int j : basisFunctionSetPositions_)
        {
            if (j != -1)
            {
                std::size_t n = basisFunctionSet_->at(j)->size();
                if (i - basePosition < n)
                {
                    return basisFunctionSet_->at(j)->evalDeriv(i - basePosition, jDir, p);
                }
                else
                {
                    basePosition += n;
                }
            }
        }
        logger(ERROR, "This is not supposed to happen, please try again with assertions turned on");
        return std::numeric_limits<double>::quiet_NaN();
    }
    
    double Element::basisFunction(std::size_t i, const PointReferenceT& p) const
    {
        logger.assert(i<getNrOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNrOfBasisFunctions());
        const Base::BaseBasisFunction* function = nullptr;
        int basePosition(0);
        for (int j : basisFunctionSetPositions_)
        {
            if (j != -1)
            {
                std::size_t n = basisFunctionSet_->at(j)->size();
                if (i - basePosition < n)
                {
                    function = basisFunctionSet_->at(j)->operator[](i - basePosition);
                    basePosition += n;
                }
                else
                {
                    basePosition += n;
                }
            }
        }
        return getReferenceGeometry()->getBasisFunctionValue(function, p);
    }
    
    std::size_t Element::getID() const
    {
        return id_;
    }
    
    std::size_t Element::getID()
    {
        return id_;
    }
    
    void Element::setQuadratureRulesWithOrder(std::size_t quadrROrder)
    {
        quadratureRule_ = Geometry::ElementGeometry::referenceGeometry_->getGaussQuadratureRule(quadrROrder);
    }
    
    void Element::setGaussQuadratureRule(GaussQuadratureRuleT* const quadR)
    {
        logger.assert(quadR!=nullptr, "Invalid quadrature rule passed");
        quadratureRule_ = quadR;
    }
    
    const QuadratureRules::GaussQuadratureRule* Element::getGaussQuadratureRule() const
    {
        return quadratureRule_;
    }
    
    std::vector<Base::ElementCacheData>& Element::getVecCacheData()
    {
        return vecCacheData_;
    }
    
    Element::SolutionVector Element::getSolution(std::size_t timeLevel, const PointReferenceT& p) const
    {
        std::size_t numberOfUnknows = ElementData::getNrOfUnknows();
        std::size_t numberOfBasisFunctions = ElementData::getNrOfBasisFunctions();
        SolutionVector solution(numberOfUnknows);
        
        LinearAlgebra::NumericalVector data(numberOfBasisFunctions * numberOfUnknows);
        data = ElementData::getTimeLevelDataVector(timeLevel);
        
        std::size_t iVB = 0;
        for (std::size_t iV = 0; iV < numberOfUnknows; ++iV)
        {
            for (std::size_t iB = 0; iB < numberOfBasisFunctions; ++iB)
            {
                iVB = convertToSingleIndex(iB, iV);
                solution[iV] += data(iVB) * basisFunction(iB, p);
            }
        }
        return solution;
    }
    
    std::vector<LinearAlgebra::NumericalVector> Element::getSolutionGradient(std::size_t timeLevel, const PointReferenceT& p) const
    {
        std::size_t numberOfUnknows = ElementData::getNrOfUnknows();
        std::size_t numberOfBasisFunctions = ElementData::getNrOfBasisFunctions();
        std::vector<LinearAlgebra::NumericalVector> solution(numberOfUnknows, LinearAlgebra::NumericalVector(p.size()));

        LinearAlgebra::NumericalVector data(numberOfBasisFunctions * numberOfUnknows);
        data = ElementData::getTimeLevelDataVector(timeLevel);

        std::size_t iVB = 0;
        for (std::size_t iV = 0; iV < numberOfUnknows; ++iV)
        {
            for (std::size_t iB = 0; iB < numberOfBasisFunctions; ++iB)
            {
                iVB = convertToSingleIndex(iB, iV);
                solution[iV] += data(iVB) * basisFunctionDeriv(iB, p);
            }
        }
        return solution;
    }

    void Element::basisFunction(std::size_t i, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret) const
    {
        logger.assert(i<getNrOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNrOfBasisFunctions());
        int basePosition(0);
        for (int j : basisFunctionSetPositions_)
        {
            if (j != -1)
            {
                int n = basisFunctionSet_->at(j)->size();
                if (i - basePosition < n)
                {
                    basisFunctionSet_->at(j)->eval(i - basePosition, p, ret);
                    return;
                }
                else
                {
                    basePosition += n;
                }
            }
        }
    }
    
    LinearAlgebra::NumericalVector Element::basisFunctionCurl(std::size_t i, const PointReferenceT& p) const
    {
        logger.assert(i<getNrOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNrOfBasisFunctions());
        int basePosition(0);
        for (int j : basisFunctionSetPositions_)
        {
            if (j != -1)
            {
                int n = basisFunctionSet_->at(j)->size();
                if (i - basePosition < n)
                {
                    return basisFunctionSet_->at(j)->evalCurl(i - basePosition, p);
                }
                else
                {
                    basePosition += n;
                }
            }
        }
        logger(ERROR, "This is not supposed to happen, please try again with assertions turned on");
        return LinearAlgebra::NumericalVector(0);
    }
    
    LinearAlgebra::NumericalVector Element::basisFunctionDeriv(std::size_t i, const PointReferenceT& p, const Element* wrapper) const
    {
        logger.assert(i<getNrOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNrOfBasisFunctions());
        if (wrapper == nullptr)
        {
            wrapper = this; 
        }
        int basePosition(0);
        for (int j : basisFunctionSetPositions_)
        {
            if (j != -1)
            {
                std::size_t n = basisFunctionSet_->at(j)->size();
                if (i - basePosition < n)
                {
                    Utilities::PhysGradientOfBasisFunction functionGradient(wrapper, basisFunctionSet_->at(j)->operator[](i - basePosition));
                    return functionGradient(p);
                }
                else
                {
                    basePosition += n;
                }
            }
        }
        logger(ERROR, "It should not be possible to reach this line");
        return LinearAlgebra::NumericalVector(0);
    }
    
#ifndef NDEBUG
    const Base::BaseBasisFunction* Element::getBasisFunction(std::size_t i) const
    {
        logger.assert(i<getNrOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNrOfBasisFunctions());
        int basePosition = 0;
        for (int j : basisFunctionSetPositions_)
        {
            if (j != -1)
            {
                if (i - basePosition < basisFunctionSet_->at(j)->size())
                {
                    return basisFunctionSet_->at(j)->operator[](i - basePosition);
                }
                else
                {
                    basePosition += basisFunctionSet_->at(j)->size();
                }
            }
        }
        logger(ERROR, "This is not supposed to happen, please try again with assertions turned on");
        return nullptr;
    }
#endif
    
    void Element::setFace(std::size_t localFaceNr, const Face* face)
    {
        logger.assert(localFaceNr < getNrOfFaces(), "Asked for face %, but there are only % faces", localFaceNr, getNrOfFaces());
        logger.assert(face!=nullptr, "Invalid face passed");
        logger.assert((face->getPtrElementLeft() == this && face->localFaceNumberLeft() == localFaceNr) 
                || (face->getPtrElementRight() == this && face->localFaceNumberRight() == localFaceNr),
                      "You are only allowed to set a face to a local face index that matches");
        if (facesList_.size() < localFaceNr + 1)
        {
            logger(WARN, "Resizing the facesList, since it's smaller(%) than to localFaceNr + 1(%)", facesList_.size(), localFaceNr + 1);
            facesList_.resize(localFaceNr + 1);
        }
        facesList_[localFaceNr] = face;
    }
    
    void Element::setEdge(std::size_t localEdgeNr, const Edge* edge)
    {
        logger.assert(localEdgeNr < getNrOfEdges(), "Asked for edge %, but there are only % edges", localEdgeNr, getNrOfEdges());
        logger.assert(edge!=nullptr, "Invalid edge passed");
        //This if-statement is needed, since it could happen in 4D
        if (edgesList_.size() < localEdgeNr + 1)
        {
            edgesList_.resize(localEdgeNr + 1);
        }
        edgesList_[localEdgeNr] = edge;
    }
    
    void Element::setNode(std::size_t localNodeNr, const Node* node)
    {
        logger.assert(node!=nullptr, "Invalid node passed");
        logger.assert(std::find(nodesList_.begin(), nodesList_.end(), node) == nodesList_.end(), "Trying to add node %, but it was already added", node->getID());
        logger.assert(localNodeNr < getNrOfNodes(), "Asked for node %, but there are only % nodes", localNodeNr, getNrOfNodes());
        nodesList_[localNodeNr] = node;
    }
    
    ///Function that computes the mass matrix. First resize the mass matrix to 
    ///the correct size, then for all quadrature points, compute the values of 
    ///all the products of basisfunctions and add this with the appropriate weight
    ///to the mass matrix.
    void Element::computeMassMatrix()
    {
        //get the number of basisfunctions, dimension and number of quadrature 
        //points on this element.
        std::size_t numBasisFuncs = getNrOfBasisFunctions();
        std::size_t dim = quadratureRule_->dimension();
        std::size_t numQuadPoints = quadratureRule_->nrOfPoints();
        
        //make the mass matrix of the correct size and set all entries to zero.
        massMatrix_.resize(numBasisFuncs, numBasisFuncs);
        massMatrix_ *= 0;
        
        //declare the relevant auxiliary variables
        LinearAlgebra::Matrix tempMatrix(numBasisFuncs, numBasisFuncs);
        Geometry::Jacobian jac(dim, dim);
        
        //for each quadrature point, compute the value of the product of the 
        //basisfunctions, then add it with the correct weight to massMatrix_
        for (std::size_t pIndex = 0; pIndex < numQuadPoints; ++pIndex)
        {
            Geometry::PointReference p = quadratureRule_->getPoint(pIndex);
            jac = calcJacobian(p);
            for (std::size_t i = 0; i < numBasisFuncs; ++i)
            {
                for (std::size_t j = 0; j < numBasisFuncs; ++j)
                {
                    tempMatrix(i, j) = basisFunction(i, p) * basisFunction(j, p);
                }
            }
            massMatrix_.axpy((quadratureRule_->weight(pIndex)) * std::abs(jac.determinant()), tempMatrix);
        }
    }

    std::ostream& operator<<(std::ostream& os, const Element& element)
    {
        os << '(';
        const Geometry::ElementGeometry& elemG = static_cast<const Geometry::ElementGeometry&>(element);
        operator<<(os, elemG);
        os << std::endl;
        return os;
    }

}

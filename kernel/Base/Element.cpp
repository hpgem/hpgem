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
#include "LinearAlgebra/MiddleSizeVector.h"
#include "Geometry/PointPhysical.h"
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
        
    Element::Element(std::size_t owner, bool owned, const ElementData& otherData, const Geometry::ElementGeometry& otherGeometry)
        : ElementGeometry(otherGeometry), ElementData(otherData), owned_ (owned), owner_ (owner)
    {        
    }
    
    Element* Element::copyWithoutFacesEdgesNodes()
    {
        //Make a new element with the data and geometry of this element
        Element* newElement = new Element(this->owner_, this->owned_, *this, *this);
        
        //copy the pointers to singletons
        newElement->quadratureRule_ = quadratureRule_;
        newElement->basisFunctions_ = basisFunctions_;
        
        //copy other data
        newElement->id_ = id_;
        newElement->orderCoeff_ = orderCoeff_;
        
        //allocate memory for nodesList, facesList and edgesList
        newElement->nodesList_.resize(getNumberOfNodes());
        newElement->edgesList_.resize(getNumberOfEdges());
        newElement->facesList_.resize(getNumberOfFaces());
        
        return newElement;
    }

    void Element::setDefaultBasisFunctionSet(std::size_t position, std::size_t unknown)
    {
        basisFunctions_.clearBasisFunctionPosition(unknown);
        basisFunctions_.registerBasisFunctionPosition(unknown, 0, position);
        std::size_t numberOfBasisFunctions = basisFunctions_.getNumberOfBasisFunctions(unknown);
        setNumberOfBasisFunctions(numberOfBasisFunctions, unknown);
        setQuadratureRulesWithOrder(orderCoeff_ * basisFunctions_.getMaximumOrder() + 1);
    }
    
    void Element::setFaceBasisFunctionSet(std::size_t position, std::size_t localFaceIndex, std::size_t unknown)
    {
        logger.assert_debug(localFaceIndex < getNumberOfFaces(), "Asked for face %, but there are only % faces", localFaceIndex, getNumberOfFaces());
        basisFunctions_.registerBasisFunctionPosition(unknown, 1 + localFaceIndex, position);
        std::size_t numberOfBasisFunctions = basisFunctions_.getNumberOfBasisFunctions(unknown);
        setNumberOfBasisFunctions(numberOfBasisFunctions, unknown);
    }
    
    void Element::setEdgeBasisFunctionSet(std::size_t position, std::size_t localEdgeIndex, std::size_t unknown)
    {
        logger.assert_debug(localEdgeIndex < getNumberOfEdges(), "Asked for edge %, but there are only % edges", localEdgeIndex, getNumberOfEdges());
        basisFunctions_.registerBasisFunctionPosition(unknown, 1 + getNumberOfFaces() + localEdgeIndex, position);
        std::size_t numberOfBasisFunctions = basisFunctions_.getNumberOfBasisFunctions(unknown);
        setNumberOfBasisFunctions(numberOfBasisFunctions, unknown);
    }
    
    void Element::setVertexBasisFunctionSet(std::size_t position, std::size_t localNodeIndex, std::size_t unknown)
    {
        logger.assert_debug(localNodeIndex < getNumberOfNodes(), "Asked for node %, but there are only % nodes", localNodeIndex, getNumberOfNodes());
        basisFunctions_.registerBasisFunctionPosition(unknown,
                1 + getNumberOfFaces() + getNumberOfEdges() + localNodeIndex, position);
        std::size_t numberOfBasisFunctions = basisFunctions_.getNumberOfBasisFunctions(unknown);
        setNumberOfBasisFunctions(numberOfBasisFunctions, unknown);
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
        bool changed = false;
        if (quadratureRule_ == nullptr)
        {
            quadratureRule_ = Geometry::ElementGeometry::referenceGeometry_->getGaussQuadratureRule(quadrROrder);
            changed = true;
        }
        else if(quadrROrder > quadratureRule_->order())
        {
            quadratureRule_ = Geometry::ElementGeometry::referenceGeometry_->getGaussQuadratureRule(quadrROrder);
            changed = true;
        }
        // The quadrature orders of the faces are assigned when they are created
        // based on the current quadrature order on its neighbouring elements.
        // If the quadrature order changes on this element it probably also has
        // to change on the face.
        if(changed)
        {
            for(std::size_t i = 0; i < facesList_.size(); ++i)
            {
                Face* face = facesList_[i];
                if(face != nullptr)
                {
                    QuadratureRules::GaussQuadratureRule* rule = face->getGaussQuadratureRule();
                    std::size_t currentOrder = rule == nullptr ? 0 : rule->order();
                    if(quadrROrder > currentOrder)
                    {
                        face->setGaussQuadratureRule(
                                this->getReferenceGeometry()
                                    ->getCodim1ReferenceGeometry(i)
                                    ->getGaussQuadratureRule(quadrROrder)
                        );
                    }
                }
            }
        }
    }
    
    void Element::setGaussQuadratureRule(QuadratureRules::GaussQuadratureRule* const quadR)
    {
        logger.assert_debug(quadR != nullptr, "Invalid quadrature rule passed");
        quadratureRule_ = quadR;
    }
    
    QuadratureRules::GaussQuadratureRule* Element::getGaussQuadratureRule() const
    {
        return quadratureRule_;
    }
    
#ifndef NDEBUG
    const Base::BaseBasisFunction* Element::getBasisFunction(std::size_t i) const
    {
        logger.assert_debug(i < getNumberOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNumberOfBasisFunctions());
        const BasisFunctionSet *subSet;
        std::size_t subIndex;
        std::tie(subSet, subIndex) = basisFunctions_.getBasisFunctionSetAndIndex(i);
        return (*subSet)[subIndex];
    }
#endif
    
    void Element::setFace(std::size_t localFaceNumber, Face* face)
    {
        logger.assert_debug(localFaceNumber < getNumberOfFaces(), "Asked for face %, but there are only % faces", localFaceNumber, getNumberOfFaces());
        logger.assert_debug(face != nullptr, "Invalid face passed");
        logger.assert_debug((face->getPtrElementLeft() == this && face->localFaceNumberLeft() == localFaceNumber)
                            || (face->getPtrElementRight() == this && face->localFaceNumberRight() == localFaceNumber),
                            "You are only allowed to set a face to a local face index that matches");
        if (facesList_.size() < localFaceNumber + 1)
        {
            logger(WARN, "Resizing the facesList, since it's smaller(%) than to localFaceNumber + 1(%)", facesList_.size(), localFaceNumber + 1);
            facesList_.resize(localFaceNumber + 1);
        }
        facesList_[localFaceNumber] = face;
    }
    
    void Element::setEdge(std::size_t localEdgeNumber, Edge* edge)
    {
        logger.assert_debug(localEdgeNumber < getNumberOfEdges(), "Asked for edge %, but there are only % edges", localEdgeNumber, getNumberOfEdges());
        logger.assert_debug(edge != nullptr, "Invalid edge passed");
        //This if-statement is needed, since it could happen in 4D
        if (edgesList_.size() < localEdgeNumber + 1)
        {
            edgesList_.resize(localEdgeNumber + 1);
        }
        edgesList_[localEdgeNumber] = edge;
    }
    
    void Element::setNode(std::size_t localNodeNumber, Node* node)
    {
        logger.assert_debug(node != nullptr, "Invalid node passed");
        logger.assert_debug(std::find(nodesList_.begin(), nodesList_.end(), node) == nodesList_.end(), "Trying to add node %, but it was already added",
                            node->getID());
        logger.assert_debug(localNodeNumber < getNumberOfNodes(), "Asked for node %, but there are only % nodes", localNodeNumber, getNumberOfNodes());
        nodesList_[localNodeNumber] = node;
    }

    std::ostream& operator<<(std::ostream& os, const Element& element)
    {
        os << '(';
        const Geometry::ElementGeometry& elemG = static_cast<const Geometry::ElementGeometry&>(element);
        operator<<(os, elemG);
        os << std::endl;
        return os;
    }

    double Element::basisFunction(std::size_t i, QuadratureRules::GaussQuadratureRule *quadratureRule,
                                  std::size_t quadraturePointIndex) const
    {
        logger.assert_debug(i < getNumberOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNumberOfBasisFunctions());
        const BasisFunctionSet *subSet;
        std::size_t subIndex;
        std::tie(subSet, subIndex) = basisFunctions_.getBasisFunctionSetAndIndex(i);
        return quadratureRule->eval(subSet, subIndex, quadraturePointIndex);
    }
    
    double Element::basisFunction(std::size_t i, QuadratureRules::GaussQuadratureRule *quadratureRule,
                                  std::size_t quadraturePointIndex, std::size_t unknown) const
    {
        logger.assert_debug(i < getNumberOfBasisFunctions(unknown), "Asked for basis function %, but there are only % basis functions", i,
                            getNumberOfBasisFunctions(unknown));
        const BasisFunctionSet *subSet;
        std::size_t subIndex;
        std::tie(subSet, subIndex) = basisFunctions_.getBasisFunctionSetAndIndex(i, unknown);
        return quadratureRule->eval(subSet, subIndex, quadraturePointIndex);
    }

    double Element::basisFunctionDiv(std::size_t i, QuadratureRules::GaussQuadratureRule *quadratureRule,
                                     std::size_t quadraturePointIndex) const
    {
        logger.assert_debug(i < getNumberOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNumberOfBasisFunctions());
        const BasisFunctionSet *subSet;
        std::size_t subIndex;
        std::tie(subSet, subIndex) = basisFunctions_.getBasisFunctionSetAndIndex(i);
        return quadratureRule->evalDiv(subSet, subIndex, quadraturePointIndex);
    }
    
    double Element::basisFunctionDiv(std::size_t i, QuadratureRules::GaussQuadratureRule *quadratureRule,
                                     std::size_t quadraturePointIndex, std::size_t unknown) const
    {
        logger.assert_debug(i < getNumberOfBasisFunctions(unknown), "Asked for basis function %, but there are only % basis functions", i,
                            getNumberOfBasisFunctions(unknown));
        const BasisFunctionSet *subSet;
        std::size_t subIndex;
        std::tie(subSet, subIndex) = basisFunctions_.getBasisFunctionSetAndIndex(i, unknown);
        return quadratureRule->evalDiv(subSet, subIndex, quadraturePointIndex);
    }

    double Element::basisFunction(std::size_t i, QuadratureRules::GaussQuadratureRule *quadratureRule,
                                  std::size_t quadraturePointIndex, const Geometry::MappingReferenceToReference<1> *map) const
    {
        logger.assert_debug(i < getNumberOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNumberOfBasisFunctions());
        const BasisFunctionSet *subSet;
        std::size_t subIndex;
        std::tie(subSet, subIndex) = basisFunctions_.getBasisFunctionSetAndIndex(i);
        return quadratureRule->eval(subSet, subIndex, quadraturePointIndex, map);
    }
    
    double Element::basisFunction(std::size_t i, QuadratureRules::GaussQuadratureRule *quadratureRule,
                                  std::size_t quadraturePointIndex, const Geometry::MappingReferenceToReference<1> *map, std::size_t unknown) const
    {
        logger.assert_debug(i < getNumberOfBasisFunctions(unknown), "Asked for basis function %, but there are only % basis functions", i,
                            getNumberOfBasisFunctions(unknown));
        const BasisFunctionSet *subSet;
        std::size_t subIndex;
        std::tie(subSet, subIndex) = basisFunctions_.getBasisFunctionSetAndIndex(i, unknown);
        return quadratureRule->eval(subSet, subIndex, quadraturePointIndex, map);
    }

    double Element::basisFunctionDiv(std::size_t i, QuadratureRules::GaussQuadratureRule *quadratureRule,
                                     std::size_t quadraturePointIndex, const Geometry::MappingReferenceToReference<1> *map) const
    {
        logger.assert_debug(i < getNumberOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNumberOfBasisFunctions());
        const BasisFunctionSet *subSet;
        std::size_t subIndex;
        std::tie(subSet, subIndex) = basisFunctions_.getBasisFunctionSetAndIndex(i);
        return quadratureRule->evalDiv(subSet, subIndex, quadraturePointIndex, map);
    }
    
    double Element::basisFunctionDiv(std::size_t i, QuadratureRules::GaussQuadratureRule *quadratureRule,
                                     std::size_t quadraturePointIndex, const Geometry::MappingReferenceToReference<1> *map, std::size_t unknown) const
    {
        logger.assert_debug(i < getNumberOfBasisFunctions(unknown), "Asked for basis function %, but there are only % basis functions", i,
                            getNumberOfBasisFunctions(unknown));
        const BasisFunctionSet *subSet;
        std::size_t subIndex;
        std::tie(subSet, subIndex) = basisFunctions_.getBasisFunctionSetAndIndex(i, unknown);
        return quadratureRule->evalDiv(subSet, subIndex, quadraturePointIndex, map);
    }

    void Element::setOwnedByCurrentProcessor(std::size_t owner, bool owned)
    {
        owner_ = owner;
        owned_ = owned;
    }

    bool Element::isOwnedByCurrentProcessor() const
    {
        return owned_;
    }

    std::size_t Element::getOwner() const
    {
        return owner_;
    }
}

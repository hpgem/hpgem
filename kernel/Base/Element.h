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

//----------------------------------------------------------------
#ifndef Element_h
#define Element_h
//----------------------------------------------------------------
#include "Base/ElementData.h"
#include "Geometry/ElementGeometry.h"
#include "LinearAlgebra/MiddleSizeVector.h"
#include "BasisFunctionSet.h"
#include "TreeEntry.h"

//#include "PhysGradientOfBasisFunction.h"//included after class definition due to cross dependencies (needed for templated function definition)
//#include "PhysicalElement.h"//included after class definition due to cross dependencies (needed for templated function definition)

#include <vector>
#include <iostream>
#include <memory>
#include <Integration/QuadratureRules/GaussQuadratureRule.h>

namespace QuadratureRules
{
    class GaussQuadratureRule;
}

namespace Base
{
    class Node;
    class Edge;
    class Face;
    class BasisFunctionSet;
    class ElementCacheData;
    class BaseBasisFunction;
    template<std::size_t DIM>
    class PhysicalElement;

    //class is final as a reminder that the virtual default destructor should be added once something inherits from this class
    class Element final: public Geometry::ElementGeometry, public ElementData
    {
    public:
        using SolutionVector = LinearAlgebra::MiddleSizeVector;
        using CollectionOfBasisFunctionSets = std::vector<std::shared_ptr<const BasisFunctionSet>>;
        
        template<std::size_t DIM>
        Element(const std::vector<std::size_t>& globalNodeIndexes, 
            const CollectionOfBasisFunctionSets *basisFunctionSet, 
            std::vector<Geometry::PointPhysical<DIM> >& allNodes, 
            std::size_t numberOfUnkowns, 
            std::size_t numberOfTimeLevels, 
            std::size_t numberOfBasisFunctions, 
            std::size_t id, 
            std::size_t numberOfElementMatrices = 0, 
            std::size_t numberOfElementVectors = 0, 
            const std::vector<int>& basisFunctionSetPositions = std::vector<int>(1, 0));

        Element(const Element &other) = delete;
        Element& operator=(const Element &other) = delete;
        
        Element* copyWithoutFacesEdgesNodes();

        std::size_t getID() const;

        std::size_t getID();

        void setQuadratureRulesWithOrder(std::size_t quadrROrder);

        void setGaussQuadratureRule(QuadratureRules::GaussQuadratureRule * const quadR);

        void setDefaultBasisFunctionSet(std::size_t position);

        void setVertexBasisFunctionSet(std::size_t position, std::size_t localIndex);
        void setEdgeBasisFunctionSet(std::size_t position, std::size_t localIndex);
        void setFaceBasisFunctionSet(std::size_t position, std::size_t localIndex);

        /// \brief Get a pointer to the quadrature rule used to do integration on this element.
        QuadratureRules::GaussQuadratureRule* getGaussQuadratureRule() const;

        std::vector<Base::ElementCacheData>& getVecCacheData();

        /// \brief Get the value of the basis function (corresponding to index i) at the physical point corresponding to reference point p.
        template<std::size_t DIM>
        double basisFunction(std::size_t i, const Geometry::PointReference<DIM>& p) const;

        double basisFunction(std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule, std::size_t quadraturePointIndex) const;

        double basisFunction(std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule, std::size_t quadraturePointIndex, const Geometry::MappingReferenceToReference<1> *map) const;

        /// \brief returns the value of the i-th basisfunction at point p in ret.
        template<std::size_t DIM>
        void basisFunction(std::size_t i, const Geometry::PointReference<DIM>& p, LinearAlgebra::SmallVector<DIM>& ret) const;

        template<std::size_t DIM>
        void basisFunction(std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule, std::size_t quadraturePointIndex, LinearAlgebra::SmallVector<DIM>& ret) const;

        template<std::size_t DIM>
        void basisFunction(std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule, std::size_t quadraturePointIndex, const Geometry::MappingReferenceToReference<1> *map, LinearAlgebra::SmallVector<DIM>& ret) const;

        /// \brief Get the value of the derivative of the physical basisfunction (corresponding to index i) in the direction jDir at the physical point corresponding to reference point p.
        template<std::size_t DIM>
        double basisFunctionDeriv(std::size_t i, std::size_t jDir, const Geometry::PointReference<DIM>& p) const;
        
        /// \brief Get the gradient of the physical basis function (corresponding to index i) at the physical point corresponding to reference point p.
        /// \details WARNING: contrary to some previous versions of hpGEM, this will NOT do any coordinate transformations!!
        template<std::size_t DIM>
        LinearAlgebra::SmallVector<DIM> basisFunctionDeriv(std::size_t i, const Geometry::PointReference<DIM>& p) const;

        template<std::size_t DIM>
        LinearAlgebra::SmallVector<DIM> basisFunctionDeriv(std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule, std::size_t quadraturePointIndex) const;

        template<std::size_t DIM>
        LinearAlgebra::SmallVector<DIM> basisFunctionDeriv(std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule, std::size_t quadraturePointIndex, const Geometry::MappingReferenceToReference<1> *map) const;

        /// \brief Returns the curl of the i-th basisfunction at point p in ret.
        template<std::size_t DIM>
        LinearAlgebra::SmallVector<DIM> basisFunctionCurl(std::size_t i, const Geometry::PointReference<DIM>& p) const;

        template<std::size_t DIM>
        LinearAlgebra::SmallVector<DIM> basisFunctionCurl(std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule, std::size_t quadraturePointIndex) const;

        template<std::size_t DIM>
        LinearAlgebra::SmallVector<DIM> basisFunctionCurl(std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule, std::size_t quadraturePointIndex, const Geometry::MappingReferenceToReference<1> *map) const;

        /// \brief Returns the divergence of the i-th basisfunction at point p in ret.
        template<std::size_t DIM>
        double basisFunctionDiv(std::size_t i, const Geometry::PointReference<DIM>& p) const;

        double basisFunctionDiv(std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule, std::size_t quadraturePointIndex) const;

        double basisFunctionDiv(std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule, std::size_t quadraturePointIndex, const Geometry::MappingReferenceToReference<1> *map) const;

        /// \brief Get the solution at the given timeLevel at the physical point corresponding to reference point p.
        /// \details This routine assumes the result of evaluating a basis function has to be transformed using the identity transformation
        template<std::size_t DIM>
        SolutionVector getSolution(std::size_t timeLevel, const Geometry::PointReference<DIM>& p) const;

        /// \brief Get the gradient of the solution at the given timeLevel at the physical point corresponding to reference point p.
        /// \details returns a vector of gradients. This routine assumes the result of evaluating a gradient of a basis function has to be transformed using the identity transformation
        template<std::size_t DIM>
        std::vector<LinearAlgebra::SmallVector<DIM> > getSolutionGradient(std::size_t timeLevel, const Geometry::PointReference<DIM>& p) const;

        /// \brief Get the solution at the given timeLevel at the physical point corresponding to reference point p.
        /// \details uses the physical element for evaluation and transformation of the basis functions
        template<std::size_t DIM>
        SolutionVector getSolution(std::size_t timeLevel, PhysicalElement<DIM>& element) const;

        /// \brief Get the gradient of the solution at the given timeLevel at the physical point corresponding to reference point p.
        /// \details returns a vector of gradients. Uses the physical element for evaluation and transformation of gradients of the basis functions
        template<std::size_t DIM>
        std::vector<LinearAlgebra::SmallVector<DIM> > getSolutionGradient(std::size_t timeLevel, PhysicalElement<DIM>& element) const;

        ///\todo not implemented
        void initialiseSolution(std::size_t timeLevel, std::size_t solutionId, const SolutionVector& solution);
                
        void setFace(std::size_t localFaceNumber, const Face* face);

        void setEdge(std::size_t localEdgeNumber, const Edge* edge);

        void setNode(std::size_t localNodeNumber, Node* node);

        ///\deprecated Does not follow naming conventions, use getLocalNumberOfBasisFunctions instead
        std::size_t getLocalNrOfBasisFunctions() const
        {
            return getLocalNumberOfBasisFunctions();
        }
        
        ///return the number of basis functions that are associated with this element only. This always includes functions with compact support on the interior of the element and DG basis function, but never include conforming basis functions that are nonzero on a face, edge or node
        std::size_t getLocalNumberOfBasisFunctions() const
        {
            return numberOfDOFinTheElement_;
        }
        
        const Face* getFace(std::size_t localFaceNumber) const
        {
            logger.assert(localFaceNumber<getNumberOfFaces(), "Asked for face %, but there are only % faces", localFaceNumber, getNumberOfFaces());
            return facesList_[localFaceNumber];
        }
        
        const std::vector<const Face*> getFacesList() const
        {
            return facesList_;
        }

        const Edge* getEdge(std::size_t localEdgeNumber) const
        {
            logger.assert(localEdgeNumber<getNumberOfEdges(), "Asked for edge %, but there are only % edges", localEdgeNumber, getNumberOfEdges());
            return edgesList_[localEdgeNumber];
        }
        
        const std::vector<const Edge*> getEdgesList() const
        {
            return edgesList_;
        }

        const Node* getNode(std::size_t localNodeNumber) const
        {
            logger.assert(localNodeNumber<getNumberOfNodes(), "Asked for node %, but there are only % nodes", localNodeNumber, getNumberOfNodes());
            return nodesList_[localNodeNumber];
        }

        Node* getNode(std::size_t localNodeNumber)
        {
            logger.assert(localNodeNumber<getNumberOfNodes(), "Asked for node %, but there are only % nodes", localNodeNumber, getNumberOfNodes());
            return nodesList_[localNodeNumber];
        }
        
        const std::vector<Node*> getNodesList() const
        {
            return nodesList_;
        }

        ///\deprecated Does not follow naming conventions, use getNumberOfFaces instead
        std::size_t getNrOfFaces() const
        {
            return getNumberOfFaces();
        }
        
        ///\deprecated Does not follow naming conventions, use getNumberOfEdges instead
        std::size_t getNrOfEdges() const
        {
            return getNumberOfEdges();
        }
        
        ///\deprecated Does not follow naming conventions, use getNumberOfNodes instead
        std::size_t getNrOfNodes() const
        {
            return getNumberOfNodes();
        }
        
        std::size_t getNumberOfFaces() const
        {
            return facesList_.size();
        }
        
        std::size_t getNumberOfEdges() const
        {
            return edgesList_.size();
        }

        ///\todo This function overwrites the non-virtual function in ElementGeometry. Either remove this function, rename one of the two functions or make the one in ElementGeometry virtual
        std::size_t getNumberOfNodes() const
        {
            return nodesList_.size();
        }
        
#ifndef NDEBUG
        const Base::BaseBasisFunction* getBasisFunction(std::size_t i) const;
#endif

        void setPositionInTree(const TreeEntry<Element*>* position) {
            logger.assert(position->getData() == this, "Trying to set the position of another element as this element");
            positionInTheTree_ = position;
        }

        const TreeEntry<Element*>* getPositionInTree() const {
            return positionInTheTree_;
        }

        /// Output operator.        
        friend std::ostream& operator<<(std::ostream& os, const Element& element);
        
    private:
        ///Constructor that copies the data and geometry of the given ElementData and ElementGeometry.
        Element(const ElementData& otherData, const Geometry::ElementGeometry& otherGeometry);

        std::tuple<const BasisFunctionSet*, std::size_t> getBasisFunctionSetAndIndex(std::size_t index) const;

        /// Quadrature rule used to do the integration on this element.
        QuadratureRules::GaussQuadratureRule *quadratureRule_;
        
        /// Vector of basis function sets. Usually you only need one basis function set.
        const CollectionOfBasisFunctionSets *basisFunctionSet_;
        
        const TreeEntry<Element*>* positionInTheTree_;

        /// Identifier (index) of the element.
        std::size_t id_;
        
        /// Constant that describes a relation between the polynomial order of the basis function set and the accuracy of the quadrature rule.
        std::size_t orderCoeff_;
        
        /// Indices of the basis function sets that are used.
        std::vector<int> basisFunctionSetPositions_;
        
        std::vector<const Face*> facesList_;
        std::vector<const Edge*> edgesList_;
        std::vector<Node*> nodesList_;

        /// Degrees of freedom corresponding to this element. When using conforming basis functions only the basis functions with support on only this element are counted.
        std::size_t numberOfDOFinTheElement_;
        
        /// Vector of data which the user might want to store. For example determinants of the Jacobian for each quadrature point.
        std::vector<Base::ElementCacheData> vecCacheData_;
        
    };
}
#include "PhysGradientOfBasisFunction.h"
#include "PhysicalElement.h"
namespace Base
{
    /// \details The user does not need to worry about the construction of elements. This is done by mesh-generators. For example the interface HpgemAPIBase can be used to create meshes.
    template<std::size_t DIM>
    Element::Element(const std::vector<std::size_t>& globalNodeIndexes,
                     const CollectionOfBasisFunctionSets *basisFunctionSet,
                     std::vector<Geometry::PointPhysical<DIM> >& allNodes,
                     std::size_t numberOfUnkowns,
                     std::size_t numberOfTimeLevels,
                     std::size_t numberOfBasisFuncs,
                     std::size_t id,
                     std::size_t numberOfElementMatrices,
                     std::size_t numberOfElementVectors,
                     const std::vector<int>& basisFunctionSetPositions)
            : ElementGeometry(globalNodeIndexes, allNodes),
        ElementData(numberOfTimeLevels, numberOfUnkowns, numberOfBasisFuncs, numberOfElementMatrices, numberOfElementVectors),
        quadratureRule_(nullptr), basisFunctionSet_(basisFunctionSet),
        id_(id), basisFunctionSetPositions_(basisFunctionSetPositions), vecCacheData_()
    {
        logger.assert(basisFunctionSet!=nullptr, "Invalid basis function set passed");
        logger.assert(basisFunctionSet->size()>0, "Not enough basis function sets passed");
        logger(VERBOSE, "numberOfElementMatrices: %", numberOfElementMatrices);
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
        logger.assert(numberOfBasisFunctions==numberOfBasisFuncs, "Redundant argument set to the wrong value");
        setNumberOfBasisFunctions(numberOfBasisFunctions);
        setQuadratureRulesWithOrder(orderCoeff_ * basisFunctionSet_->at(basisFunctionSetPositions_[0])->getOrder() + 1);
        numberOfDOFinTheElement_ = basisFunctionSet_->at(basisFunctionSetPositions_[0])->size();
        facesList_.assign(getReferenceGeometry()->getNumberOfCodim1Entities(), nullptr);
        if (getReferenceGeometry()->getNumberOfCodim3Entities() > 0)
        {
            edgesList_.assign(getReferenceGeometry()->getNumberOfCodim2Entities(), nullptr);
        }
        nodesList_.assign(getReferenceGeometry()->getNumberOfNodes(), nullptr);
    }

    template<std::size_t DIM>
    double Element::basisFunctionDeriv(std::size_t i, std::size_t jDir, const Geometry::PointReference<DIM>& p) const
    {
        logger.assert(i<getNumberOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNumberOfBasisFunctions());
        logger.assert((jDir < p.size()), "Error in BasisFunctionSet.EvalDeriv: invalid derivative direction!");

        const BasisFunctionSet *subSet;
        std::size_t subIndex;
        std::tie(subSet, subIndex) = getBasisFunctionSetAndIndex(i);
        return subSet->evalDeriv(subIndex, jDir, p);
    }

    template<std::size_t DIM>
    double Element::basisFunction(std::size_t i, const Geometry::PointReference<DIM>& p) const
    {
        logger.assert(i<getNumberOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNumberOfBasisFunctions());
        const BasisFunctionSet *subSet;
        std::size_t subIndex;
        std::tie(subSet, subIndex) = getBasisFunctionSetAndIndex(i);
        return subSet->eval(subIndex, p);
    }

    template<std::size_t DIM>
    void Element::basisFunction(std::size_t i, const Geometry::PointReference<DIM>& p, LinearAlgebra::SmallVector<DIM>& ret) const
    {
        logger.assert(i<getNumberOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNumberOfBasisFunctions());
        const BasisFunctionSet *subSet;
        std::size_t subIndex;
        std::tie(subSet, subIndex) = getBasisFunctionSetAndIndex(i);
        subSet->eval(subIndex, p, ret);
    }

    template<std::size_t DIM>
    LinearAlgebra::SmallVector<DIM> Element::basisFunctionCurl(std::size_t i, const Geometry::PointReference<DIM>& p) const
    {
        logger.assert(i<getNumberOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNumberOfBasisFunctions());
        const BasisFunctionSet *subSet;
        std::size_t subIndex;
        std::tie(subSet, subIndex) = getBasisFunctionSetAndIndex(i);
        return subSet->evalCurl(subIndex, p);
    }

    template<std::size_t DIM>
    double Element::basisFunctionDiv(std::size_t i, const Geometry::PointReference<DIM>& p) const
    {
        logger.assert(i<getNumberOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNumberOfBasisFunctions());
        const BasisFunctionSet *subSet;
        std::size_t subIndex;
        std::tie(subSet, subIndex) = getBasisFunctionSetAndIndex(i);
        return subSet->evalDiv(subIndex, p);
    }

    template<std::size_t DIM>
    LinearAlgebra::SmallVector<DIM> Element::basisFunctionDeriv(std::size_t i, const Geometry::PointReference<DIM>& p) const
    {
        logger.assert(i<getNumberOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNumberOfBasisFunctions());
        const BasisFunctionSet *subSet;
        std::size_t subIndex;
        std::tie(subSet, subIndex) = getBasisFunctionSetAndIndex(i);
        return subSet->evalDeriv(subIndex, p);
    }

    template<std::size_t DIM>
    Element::SolutionVector Element::getSolution(std::size_t timeIntegrationVectorId, const Geometry::PointReference<DIM>& p) const
    {
        std::size_t numberOfUnknows = ElementData::getNumberOfUnknowns();
        std::size_t numberOfBasisFunctions = ElementData::getNumberOfBasisFunctions();
        SolutionVector solution(numberOfUnknows);

        LinearAlgebra::MiddleSizeVector data = ElementData::getTimeIntegrationVector(timeIntegrationVectorId);

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

    template<std::size_t DIM>
    std::vector<LinearAlgebra::SmallVector<DIM> > Element::getSolutionGradient(std::size_t timeIntegrationVectorId, const Geometry::PointReference<DIM>& p) const
    {
        std::size_t numberOfUnknows = ElementData::getNumberOfUnknowns();
        std::size_t numberOfBasisFunctions = ElementData::getNumberOfBasisFunctions();
        std::vector<LinearAlgebra::SmallVector<DIM> > solution(numberOfUnknows);

        LinearAlgebra::MiddleSizeVector data = ElementData::getTimeIntegrationVector(timeIntegrationVectorId);

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

    template<std::size_t DIM>
    Element::SolutionVector Element::getSolution(std::size_t timeIntegrationVectorId, PhysicalElement<DIM>& element) const
    {
        logger.assert(element.getElement() == this, "Cannot find the solution in a different element!");
        std::size_t numberOfUnknowns = ElementData::getNumberOfUnknowns();
        std::size_t numberOfBasisFunctions = ElementData::getNumberOfBasisFunctions();
        SolutionVector solution(numberOfUnknowns);

        const LinearAlgebra::MiddleSizeVector& data = ElementData::getTimeIntegrationVector(timeIntegrationVectorId);

        std::size_t iVB = 0;
        for (std::size_t iV = 0; iV < numberOfUnknowns; ++iV)
        {
            for (std::size_t iB = 0; iB < numberOfBasisFunctions; ++iB)
            {
                iVB = convertToSingleIndex(iB, iV);
                solution[iV] += data(iVB) * element.basisFunction(iB);
            }
        }
        return solution;
    }

    template<std::size_t DIM>
    std::vector<LinearAlgebra::SmallVector<DIM> > Element::getSolutionGradient(std::size_t timeIntegrationVectorId, PhysicalElement<DIM>& element) const
    {
        logger.assert(element.getElement() == this, "Cannot find the gradient of the solution in a different element!");
        std::size_t numberOfUnknowns = ElementData::getNumberOfUnknowns();
        std::size_t numberOfBasisFunctions = ElementData::getNumberOfBasisFunctions();
        std::vector<LinearAlgebra::SmallVector<DIM> > solution(numberOfUnknowns);

        LinearAlgebra::MiddleSizeVector data = ElementData::getTimeIntegrationVector(timeIntegrationVectorId);

        std::size_t iVB = 0;
        for (std::size_t iV = 0; iV < numberOfUnknowns; ++iV)
        {
            for (std::size_t iB = 0; iB < numberOfBasisFunctions; ++iB)
            {
                iVB = convertToSingleIndex(iB, iV);
                solution[iV] += data(iVB) * element.basisFunctionDeriv(iB);
            }
        }
        return solution;
    }

    template<std::size_t DIM>
    void Element::basisFunction(std::size_t i, QuadratureRules::GaussQuadratureRule *quadratureRule,
                                std::size_t quadraturePointIndex, LinearAlgebra::SmallVector<DIM> &ret) const
    {
        logger.assert(i<getNumberOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNumberOfBasisFunctions());
        const BasisFunctionSet *subSet;
        std::size_t subIndex;
        std::tie(subSet, subIndex) = getBasisFunctionSetAndIndex(i);
        quadratureRule->eval(subSet, subIndex, quadraturePointIndex, ret);
    }

    template<std::size_t DIM>
    LinearAlgebra::SmallVector<DIM> Element::basisFunctionDeriv(std::size_t i,
                                                                QuadratureRules::GaussQuadratureRule *quadratureRule,
                                                                std::size_t quadraturePointIndex) const
    {
        logger.assert(i<getNumberOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNumberOfBasisFunctions());
        const BasisFunctionSet *subSet;
        std::size_t subIndex;
        std::tie(subSet, subIndex) = getBasisFunctionSetAndIndex(i);
        return quadratureRule->evalGrad(subSet, subIndex, quadraturePointIndex);
    }

    template<std::size_t DIM>
    LinearAlgebra::SmallVector<DIM> Element::basisFunctionCurl(std::size_t i,
                                                               QuadratureRules::GaussQuadratureRule *quadratureRule,
                                                               std::size_t quadraturePointIndex) const
    {
        logger.assert(i<getNumberOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNumberOfBasisFunctions());
        const BasisFunctionSet *subSet;
        std::size_t subIndex;
        std::tie(subSet, subIndex) = getBasisFunctionSetAndIndex(i);
        return subSet->evalCurl<DIM>(subIndex, quadratureRule, quadraturePointIndex);
    }

    template<std::size_t DIM>
    void Element::basisFunction(std::size_t i, QuadratureRules::GaussQuadratureRule *quadratureRule,
                                std::size_t quadraturePointIndex, const Geometry::MappingReferenceToReference<1> *map, LinearAlgebra::SmallVector<DIM> &ret) const
    {
        logger.assert(i<getNumberOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNumberOfBasisFunctions());
        const BasisFunctionSet *subSet;
        std::size_t subIndex;
        std::tie(subSet, subIndex) = getBasisFunctionSetAndIndex(i);
        quadratureRule->eval(subSet, subIndex, quadraturePointIndex, map, ret);
    }

    template<std::size_t DIM>
    LinearAlgebra::SmallVector<DIM> Element::basisFunctionDeriv(std::size_t i,
                                                                QuadratureRules::GaussQuadratureRule *quadratureRule,
                                                                std::size_t quadraturePointIndex, const Geometry::MappingReferenceToReference<1> *map) const
    {
        logger.assert(i<getNumberOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNumberOfBasisFunctions());
        const BasisFunctionSet *subSet;
        std::size_t subIndex;
        std::tie(subSet, subIndex) = getBasisFunctionSetAndIndex(i);
        return quadratureRule->evalGrad(subSet, subIndex, quadraturePointIndex, map);
    }

    template<std::size_t DIM>
    LinearAlgebra::SmallVector<DIM> Element::basisFunctionCurl(std::size_t i,
                                                               QuadratureRules::GaussQuadratureRule *quadratureRule,
                                                               std::size_t quadraturePointIndex, const Geometry::MappingReferenceToReference<1> *map) const
    {
        logger.assert(i<getNumberOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNumberOfBasisFunctions());
        const BasisFunctionSet *subSet;
        std::size_t subIndex;
        std::tie(subSet, subIndex) = getBasisFunctionSetAndIndex(i);
        return subSet->evalCurl<DIM>(subIndex, quadratureRule, quadraturePointIndex, map);
    }
}

#endif

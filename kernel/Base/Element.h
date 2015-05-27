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
//#include "PhysGradientOfBasisFunction.h"//included after class definition due to cross dependencies (needed for templated function definition)

#include <vector>
#include <iostream>
#include <memory>

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

    //Programmer note: functions that are needed during integration should be overridden in ShortTermStorageFaceBase
    class Element : public Geometry::ElementGeometry, public ElementData
    {
    public:
        using ReferenceGeometryT = Geometry::ReferenceGeometry;
        using MappingReferenceToPhysicalT = Geometry::MappingReferenceToPhysical;
        using ElementGeometryT = Geometry::ElementGeometry;
        using CacheT = Base::ElementCacheData;
        using GaussQuadratureRuleT = QuadratureRules::GaussQuadratureRule;
        using VecCacheT = std::vector<CacheT>;
        using SolutionVector = LinearAlgebra::MiddleSizeVector;
        using CollectionOfBasisFunctionSets = std::vector<std::shared_ptr<const BasisFunctionSet>>;

    public:
        
        template<std::size_t DIM>
        Element(const std::vector<std::size_t>& globalNodeIndexes, const CollectionOfBasisFunctionSets *basisFunctionSet, std::vector<Geometry::PointPhysical<DIM> >& allNodes, std::size_t nrOfUnkowns, std::size_t nrOfTimeLevels, std::size_t nrOfBasisFunc, std::size_t id, std::size_t numberOfElementMatrices = 0, std::size_t numberOfElementVectors = 0, const std::vector<int>& basisFunctionSetPositions = std::vector<int>(1, 0));

        Element(const Element &other) = delete;
        Element& operator=(const Element &other) = delete;        
        
        virtual ~ Element();
        
        Element* copyWithoutFacesEdgesNodes();

        virtual std::size_t getID() const;

        virtual std::size_t getID();

        void setQuadratureRulesWithOrder(std::size_t quadrROrder);

        void setGaussQuadratureRule(GaussQuadratureRuleT * const quadR);

        void setDefaultBasisFunctionSet(std::size_t position);

        void setVertexBasisFunctionSet(std::size_t position, std::size_t localIndex);
        void setEdgeBasisFunctionSet(std::size_t position, std::size_t localIndex);
        void setFaceBasisFunctionSet(std::size_t position, std::size_t localIndex);

        /// \brief Get a pointer to the quadrature rule used to do integration on this element.
        virtual const GaussQuadratureRuleT* getGaussQuadratureRule() const;

        virtual VecCacheT& getVecCacheData();

        /// \brief Get the value of the basis function (corresponding to index i) at the physical point corresponding to reference point p.
        template<std::size_t DIM>
        double basisFunction(std::size_t i, const Geometry::PointReference<DIM>& p) const;

        /// \brief returns the value of the i-th basisfunction at point p in ret.
        template<std::size_t DIM>
        void basisFunction(std::size_t i, const Geometry::PointReference<DIM>& p, LinearAlgebra::SmallVector<DIM>& ret) const;

        /// \brief Get the value of the derivative of the physical basisfunction (corresponding to index i) in the direction jDir at the physical point corresponding to reference point p.
        template<std::size_t DIM>
        double basisFunctionDeriv(std::size_t i, std::size_t jDir, const Geometry::PointReference<DIM>& p) const;
        
        /// \brief Get the gradient of the physical basis function (corresponding to index i) at the physical point corresponding to reference point p.
        /// \details If some of the data needed for the reference to physical gradient mapping is already stored on a wrapper class, you can pass the class to this function for more efficient computation
        template<std::size_t DIM>
        LinearAlgebra::SmallVector<DIM> basisFunctionDeriv(std::size_t i, const Geometry::PointReference<DIM>& p) const;

        /// \brief Returns the curl of the i-th basisfunction at point p in ret.
        template<std::size_t DIM>
        LinearAlgebra::SmallVector<DIM> basisFunctionCurl(std::size_t i, const Geometry::PointReference<DIM>& p) const;

        /// \brief Get the solution at the given timeLevel at the physical point corresponding to reference point p.
        template<std::size_t DIM>
        SolutionVector getSolution(std::size_t timeLevel, const Geometry::PointReference<DIM>& p) const;

        /// \brief Get the gradient of the solution at the given timeLevel at the physical point corresponding to reference point p.
        /// \details returns a vector of gradients
        template<std::size_t DIM>
        std::vector<LinearAlgebra::SmallVector<DIM> > getSolutionGradient(std::size_t timeLevel, const Geometry::PointReference<DIM>& p) const;

        void initialiseSolution(std::size_t timeLevel, std::size_t solutionId, const SolutionVector& solution); ///\todo not implemented
                
        void setFace(std::size_t localFaceNr, const Face* face);

        void setEdge(std::size_t localEdgeNr, const Edge* edge);

        void setNode(std::size_t localNodeNr, const Node* node);

        
        virtual std::size_t getLocalNrOfBasisFunctions() const
        {
            return nrOfDOFinTheElement_;
        }
        
        virtual const Face* getFace(std::size_t localFaceNr) const
        {
            logger.assert(localFaceNr<getNrOfFaces(), "Asked for face %, but there are only % faces", localFaceNr, getNrOfFaces());
            return facesList_[localFaceNr];
        }
        
        virtual const std::vector<const Face*> getFacesList() const
        {
            return facesList_;
        }

        virtual const Edge* getEdge(std::size_t localEdgeNr) const
        {
            logger.assert(localEdgeNr<getNrOfEdges(), "Asked for edge %, but there are only % edges", localEdgeNr, getNrOfEdges());
            return edgesList_[localEdgeNr];
        }
        
        virtual const std::vector<const Edge*> getEdgesList() const
        {
            return edgesList_;
        }

        virtual const Node* getNode(std::size_t localNodeNr) const
        {
            logger.assert(localNodeNr<getNrOfNodes(), "Asked for node %, but there are only % nodes", localNodeNr, getNrOfNodes());
            return nodesList_[localNodeNr];
        }
        
        virtual const std::vector<const Node*> getNodesList() const
        {
            return nodesList_;
        }

        virtual std::size_t getNrOfFaces() const
        {
            return facesList_.size();
        }
        
        virtual std::size_t getNrOfEdges() const
        {
            return edgesList_.size();
        }
        
        virtual std::size_t getNrOfNodes() const
        {
            return nodesList_.size();
        }
        
#ifndef NDEBUG
        virtual const Base::BaseBasisFunction* getBasisFunction(std::size_t i) const;
#endif
        
    protected:
        
        ///\brief default constructor - for use with wrapper classes (that can delegate functionality of Element in another way)
        Element();

    public:
        /// Output operator.        
        friend std::ostream& operator<<(std::ostream& os, const Element& element);
        
    private:
        ///Constructor that copies the data and geometry of the given ElementData and ElementGeometry.
        Element(const ElementData& otherData, const ElementGeometry& otherGeometry);

        /// Quadrature rule used to do the integration on this element.
        const GaussQuadratureRuleT *quadratureRule_;
        
        /// Vector of basis function sets. Usually you only need one basis function set.
        const CollectionOfBasisFunctionSets *basisFunctionSet_;
        
        /// Identifier (index) of the element.
        std::size_t id_;
        
        /// Constant that describes a relation between the polynomial order of the basis function set and the accuracy of the quadrature rule.
        std::size_t orderCoeff_;
        
        /// Indices of the basis function sets that are used.
        std::vector<int> basisFunctionSetPositions_;
        
        std::vector<const Face*> facesList_;
        std::vector<const Edge*> edgesList_;
        std::vector<const Node*> nodesList_;

        /// Degrees of freedom corresponding to this element. When using conforming basis functions only the basis functions with support on only this element are counted.
        std::size_t nrOfDOFinTheElement_;
        
        /// Vector of data which the user might want to store. For example determinants of the Jacobian for each quadrature point.
        VecCacheT vecCacheData_;
        
    };
}
#include "PhysGradientOfBasisFunction.h"
namespace Base
{
    /// \details The user does not need to worry about the contruction of elements. This is done by mesh-generators. For example the interface HpgemAPIBase can be used to create meshes.
    template<std::size_t DIM>
    Element::Element(const VectorOfPointIndexesT& globalNodeIndexes,
                     const CollectionOfBasisFunctionSets *basisFunctionSet,
                     std::vector<Geometry::PointPhysical<DIM> >& allNodes,
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
        id_(id), basisFunctionSetPositions_(basisFunctionSetPositions), vecCacheData_()
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

    template<std::size_t DIM>
    double Element::basisFunctionDeriv(std::size_t i, std::size_t jDir, const Geometry::PointReference<DIM>& p) const
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

    template<std::size_t DIM>
    double Element::basisFunction(std::size_t i, const Geometry::PointReference<DIM>& p) const
    {
        logger.assert(i<getNrOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNrOfBasisFunctions());
        const Base::BaseBasisFunction* function = nullptr;
        for (int j : basisFunctionSetPositions_)
        {
            if (j != -1)
            {
                std::size_t n = basisFunctionSet_->at(j)->size();
                if (i < n)
                {
                    function = basisFunctionSet_->at(j)->operator[](i);
                    return p.getBasisFunctionValue(function);
                }
                else
                {
                    i -= n;
                }
            }
        }
        return p.getBasisFunctionValue(function);
    }

    template<std::size_t DIM>
    void Element::basisFunction(std::size_t i, const Geometry::PointReference<DIM>& p, LinearAlgebra::SmallVector<DIM>& ret) const
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

    template<std::size_t DIM>
    LinearAlgebra::SmallVector<DIM> Element::basisFunctionCurl(std::size_t i, const Geometry::PointReference<DIM>& p) const
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
        return LinearAlgebra::SmallVector<DIM>();
    }

    template<std::size_t DIM>
    LinearAlgebra::SmallVector<DIM> Element::basisFunctionDeriv(std::size_t i, const Geometry::PointReference<DIM>& p) const
    {
        logger.assert(i<getNrOfBasisFunctions(), "Asked for basis function %, but there are only % basis functions", i, getNrOfBasisFunctions());
        int basePosition(0);
        for (int j : basisFunctionSetPositions_)
        {
            if (j != -1)
            {
                std::size_t n = basisFunctionSet_->at(j)->size();
                if (i - basePosition < n)
                {
                    Utilities::PhysGradientOfBasisFunction functionGradient(this, basisFunctionSet_->at(j)->operator[](i - basePosition));
                    return functionGradient(p);
                }
                else
                {
                    basePosition += n;
                }
            }
        }
        logger(ERROR, "It should not be possible to reach this line");
        return LinearAlgebra::SmallVector<DIM>();
    }

    template<std::size_t DIM>
    Element::SolutionVector Element::getSolution(std::size_t timeLevel, const Geometry::PointReference<DIM>& p) const
    {
        std::size_t numberOfUnknows = ElementData::getNrOfUnknows();
        std::size_t numberOfBasisFunctions = ElementData::getNrOfBasisFunctions();
        SolutionVector solution(numberOfUnknows);

        LinearAlgebra::MiddleSizeVector data = ElementData::getTimeLevelDataVector(timeLevel);

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
    std::vector<LinearAlgebra::SmallVector<DIM> > Element::getSolutionGradient(std::size_t timeLevel, const Geometry::PointReference<DIM>& p) const
    {
        std::size_t numberOfUnknows = ElementData::getNrOfUnknows();
        std::size_t numberOfBasisFunctions = ElementData::getNrOfBasisFunctions();
        std::vector<LinearAlgebra::SmallVector<DIM> > solution(numberOfUnknows);

        LinearAlgebra::MiddleSizeVector data = ElementData::getTimeLevelDataVector(timeLevel);

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
}

#endif

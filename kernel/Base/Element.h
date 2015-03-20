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
#include "LinearAlgebra/NumericalVector.h"

#include <vector>
#include <iostream>

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

    class Element : public Geometry::ElementGeometry, public ElementData
    {
    public:
        using PointPhysicalT = Geometry::PointPhysical;
        using PointReferenceT = Geometry::PointReference;
        using ReferenceGeometryT = Geometry::ReferenceGeometry;
        using MappingReferenceToPhysicalT = Geometry::MappingReferenceToPhysical;
        using ElementGeometryT = Geometry::ElementGeometry;
        using CacheT = Base::ElementCacheData;
        using GaussQuadratureRuleT = QuadratureRules::GaussQuadratureRule;
        using VecCacheT = std::vector<CacheT>;
        using SolutionVector = LinearAlgebra::NumericalVector;

    public:
        
        Element(const std::vector<std::size_t>& globalNodeIndexes, const std::vector<const BasisFunctionSet*>* basisFunctionSet, std::vector<Geometry::PointPhysical>& allNodes, std::size_t nrOfUnkowns, std::size_t nrOfTimeLevels, std::size_t nrOfBasisFunc, std::size_t id, std::size_t numberOfElementMatrices = 0, std::size_t numberOfElementVectors = 0, const std::vector<int>& basisFunctionSetPositions = std::vector<int>(1, 0));

        Element(const Element& other) = delete;

        ~ Element();

        virtual std::size_t getID() const;

        virtual std::size_t getID();

        void setQuadratureRulesWithOrder(std::size_t quadrROrder);

        void setGaussQuadratureRule(GaussQuadratureRuleT * const quadR);

        void setDefaultBasisFunctionSet(std::size_t position);

        void setVertexBasisFunctionSet(std::size_t position, std::size_t localIndex);
        void setEdgeBasisFunctionSet(std::size_t position, std::size_t localIndex);
        void setFaceBasisFunctionSet(std::size_t position, std::size_t localIndex);

        virtual const GaussQuadratureRuleT* getGaussQuadratureRule() const;

        virtual VecCacheT& getVecCacheData();

        virtual double basisFunction(std::size_t i, const PointReferenceT& p) const;

        ///\brief returns the value of the i-th basisfunction at point p in ret
        virtual void basisFunction(std::size_t i, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret) const;

        /// \param[in] jDir Direction of the derivative, jDir=0 means x, and etc.
        virtual double basisFunctionDeriv(std::size_t i, std::size_t jDir, const PointReferenceT& p) const;
        
        ///\brief the all directions in one go edition of basisFunctionDeriv. Also applies the scaling gained from transforming to the reference element.
        ///\details if some of the data needed for this mapping is already stored on a wrapper class, you can pass the class to this function for more efficient computation
        virtual LinearAlgebra::NumericalVector basisFunctionDeriv(std::size_t i, const PointReferenceT& p, const Element* wrapper = nullptr) const;

        ///\brief returns the curl of the i-th basisfunction at point p in ret
        virtual LinearAlgebra::NumericalVector basisFunctionCurl(std::size_t i, const PointReferenceT& p) const;

        virtual SolutionVector getSolution(std::size_t timeLevel, const PointReferenceT& p) const;

        void initialiseSolution(std::size_t timeLevel, std::size_t solutionId, const SolutionVector& solution); ///\TODO not implemented  
                
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
        
        virtual const Edge* getEdge(std::size_t localEdgeNr) const
        {
            logger.assert(localEdgeNr<getNrOfEdges(), "Asked for edge %, but there are only % edges", localEdgeNr, getNrOfEdges());
            return edgesList_[localEdgeNr];
        }
        
        virtual const Node* getNode(std::size_t localNodeNr) const
        {
            logger.assert(localNodeNr<getNrOfNodes(), "Asked for node %, but there are only % nodes", localNodeNr, getNrOfNodes());
            return nodesList_[localNodeNr];
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
        
        ///\brief Return the mass int(phi_i phi_j) matrix of this element.
        ///\details If the mass matrix is computed earlier, we just return it. 
        ///Otherwise, the mass matrix is computed, stored and returned. 
        ///\bug does not work for moving meshes
        const LinearAlgebra::Matrix& getMassMatrix()
        {
            if (massMatrix_.size() == 0)
            {
                computeMassMatrix();
            }
            return massMatrix_;
        }
        
#ifndef NDEBUG
        virtual const Base::BaseBasisFunction* getBasisFunction(std::size_t i) const;
#endif
        
    protected:
        
        ///\brief default constructor - for use with wrapper classes (that can delegate functionality of Element in another way)
        Element();

    public:
        /// Output operator.        
        friend std::ostream& operator<<(std::ostream& os, const Element& element)
        {
            os << '(';
            const Geometry::ElementGeometry& elemG = static_cast<const Geometry::ElementGeometry&>(element);
            operator<<(os, elemG);
            os << std::endl;
            return os;
        }
        
    private:
        
        ///Compute the mass matrix of this element.
        void computeMassMatrix();

        const GaussQuadratureRuleT* quadratureRule_;
        const std::vector<const BasisFunctionSet*>* basisFunctionSet_;
        VecCacheT vecCacheData_;
        std::size_t id_;
        double orderCoeff_;
        std::vector<int> basisFunctionSetPositions_;
        std::vector<const Face*> facesList_;
        std::vector<const Edge*> edgesList_;
        std::vector<const Node*> nodesList_;

        //IN the element, so don't count conforming DOF from faces/...
        std::size_t nrOfDOFinTheElement_;

        ///Stores that mass matrix for this element
        LinearAlgebra::Matrix massMatrix_;
    };
}

#endif

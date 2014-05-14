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
#ifndef Element_hpp
#define Element_hpp
//----------------------------------------------------------------
#include "Base/ElementData.hpp"
#include "Base/BasisFunctionSet.hpp"
#include "Geometry/ElementGeometry.hpp"
#include "Geometry/PointPhysical.hpp"
#include "Geometry/ReferenceSquare.hpp"
#include "Geometry/PointPhysical.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"
#include "Base/ElementCacheData.hpp"


#include <vector>
#include <iostream>

namespace Base
{
	class Edge;
	class Face;

    class Element: public Geometry::ElementGeometry,
                   public ElementData
    {
    public:
        typedef Geometry::PointPhysical                PointPhysicalT;
        typedef Geometry::PointReference               PointReferenceT;
        typedef Geometry::ReferenceGeometry            ReferenceGeometryT;
        typedef Geometry::MappingReferenceToPhysical   MappingReferenceToPhysicalT;
        typedef Geometry::ElementGeometry              ElementGeometryT;
        typedef unsigned int                                PointIndexT;
        typedef unsigned int                                UId;
        typedef std::vector<PointPhysicalT>                 VectorOfPhysicalPointsT;
        typedef std::vector<PointIndexT>                    VectorOfPointIndexesT;
        typedef Base::ElementCacheData                 CacheT;
        typedef Base::BasisFunctionSet                 BasisFunctionSetT;
        typedef QuadratureRules::GaussQuadratureRule   GaussQuadratureRuleT;
        typedef Base::ElementData                      ElementDataT;
        typedef std::vector<CacheT>                         VecCacheT;
        typedef LinearAlgebra::NumericalVector              SolutionVector;
        
    public:
        
        Element(const VectorOfPointIndexesT& globalNodeIndexes,
                const std::vector<const BasisFunctionSetT*>* basisFunctionSet,
                const VectorOfPhysicalPointsT& allNodes,
                unsigned int nrOfUnkowns,
                unsigned int nrOfTimeLevels,
                unsigned int nrOfBasisFunc,
                unsigned int id,
                unsigned int numberOfElementMatrices=0,
                unsigned int numberOfElementVectors=0,
                const std::vector< int>& basisFunctionSetPositions=std::vector< int>(1,0));
    
        Element(const Element& other);
        
        ~Element();
	

        virtual unsigned int                    getID()const;

        virtual unsigned int                    getID();

        void                            setQuadratureRulesWithOrder(unsigned int quadrROrder);

        void                            setGaussQuadratureRule(GaussQuadratureRuleT* const quadR);

        void							setDefaultBasisFunctionSet(unsigned int position);

        void                            setVertexBasisFunctionSet(unsigned int position, int localIndex);
        void                            setEdgeBasisFunctionSet(unsigned int position,int localIndex);
        void                            setFaceBasisFunctionSet(unsigned int position,int localIndex);

        virtual const GaussQuadratureRuleT*     getGaussQuadratureRule() const;

        virtual VecCacheT&                      getVecCacheData();
        
        virtual double                          basisFunction(unsigned int i, const PointReferenceT& p) const;

		///\brief returns the value of the i-th basisfunction at point p in ret
        virtual void                            basisFunction(unsigned int i, const PointReferenceT& p, NumericalVector& ret) const;
 
            /// jDir=0 means x, and etc.
        virtual double                          basisFunctionDeriv(unsigned int i, unsigned int jDir, const PointReferenceT& p) const;
            //unsigned int                    getNumberOfDegreesOfFreedom()const;
            //unsigned int                    getNumberOfDegreesOfFreedom();
        
        ///\brief the all directions in one go edition of basisFunctionDeriv. Also applies the scaling gained from transforming to the reference element.
        ///if some of the data needed for this mapping is already stored on a wrapper class, you can pass the class to this function for more efficient computation
        virtual void                            basisFunctionDeriv(unsigned int i,const PointReferenceT& p, NumericalVector& ret,const Element* wrapper=NULL) const;
	
		///\brief returns the curl of the i-th basisfunction at point p in ret
        virtual void                            basisFunctionCurl(unsigned int i, const PointReferenceT& p, NumericalVector& ret) const;
	
        virtual void                            getSolution(unsigned int timeLevel, const PointReferenceT& p, SolutionVector& solution) const;
        
        void                            initialiseSolution(unsigned int timeLevel, unsigned int solutionId, const SolutionVector& solution);///\TODO not implemented
        
        void                            setFace(int localFaceNr, const Face* face);

        void                            setEdge(int localEdgeNr, const Edge* edge);

        virtual int                             getLocalNrOfBasisFunctions() const{return nrOfDOFinTheElement_;}

        virtual int                             getLocalNrOfBasisFunctionsVertex() const{return nrOfDOFperVertex_;}

        virtual const Face*                     getFace(int localFaceNr)const {return facesList_[localFaceNr];}

        virtual const Edge*                     getEdge(int localEdgeNr)const {return edgesList_[localEdgeNr];}

        virtual int                             getNrOfEdges() const{return edgesList_.size();}

#ifndef NDEBUG
        virtual const Base::BaseBasisFunction*  getBasisFunction(int i)const;
#endif

    protected:

        ///\brief default constructor - for use with wrapper classes (that can delegate functionality of Element in another way)
        Element():ElementData(0,0,0,0,0),ElementGeometry(),quadratureRule_(NULL),basisFunctionSet_(NULL){}

    public:
            /// Output operator.
        friend ostream& operator<<(ostream& os, const Element& element)
        {
            os << '(' ;
            const Geometry::ElementGeometry& elemG = static_cast<const Geometry::ElementGeometry&>(element);
            operator<<(os, elemG);
            os<<endl;
            return os;
        }
    
    
    private:
        GaussQuadratureRuleT*                        quadratureRule_;
        const std::vector<const BasisFunctionSetT*>* basisFunctionSet_;
        VecCacheT                                    vecCacheData_;
        UId                                          id_;
        double                                       orderCoeff_;
        std::vector<int>                    basisFunctionSetPositions_;
        std::vector<const Face*>					 facesList_;
        std::vector<const Edge*>      				 edgesList_;

        //IN the element, so don't count conforming DOF from faces/...
        unsigned int                                 nrOfDOFinTheElement_;

        //the number of DOF per vertex is just a number so if you want to do p-adaptation
        //use a set of basisfunctions where the DOF per vertex are constant in p (this is usually the case)
        unsigned int                                 nrOfDOFperVertex_;
    };
}

#endif

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
#ifndef _Element_Impl_hpp
#define _Element_Impl_hpp

#include "Element.hpp"
#include "PhysGradientOfBasisFunction.hpp"
#include "Edge.hpp"
#include "Face.hpp"

#include "BasisFunctionSet.hpp"
#include "Geometry/PhysicalGeometry.hpp"
#include "LinearAlgebra/NumericalVector.hpp"
#include "Geometry/PointPhysical.hpp"
#include "FaceCacheData.hpp"
#include "ElementCacheData.hpp"
#include "Geometry/ReferenceGeometry.hpp"
#include "Geometry/PointReference.hpp"
#include "TestErrorDebug.hpp"
#include "Node.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"
#include "Geometry/Jacobian.hpp"

namespace Base
{
    //class Element;
    
    Element::Element(const VectorOfPointIndexesT& globalNodeIndexes,
                          const std::vector<const BasisFunctionSetT*>* basisFunctionSet,
                          const VectorOfPhysicalPointsT& allNodes,
                          size_t nrOfUnkowns,
                          size_t nrOfTimeLevels,
                          size_t nrOfBasisFunc,
                          size_t id,
                          size_t numberOfElementMatrixes,
                          size_t numberOfElementVectors,
                          const std::vector< int>& basisFunctionSetPositions):
        ElementGeometryT(globalNodeIndexes, allNodes),
        ElementDataT(nrOfTimeLevels, nrOfUnkowns, nrOfBasisFunc,numberOfElementMatrixes,numberOfElementVectors),
        basisFunctionSet_(basisFunctionSet),
        quadratureRule_(nullptr),
        vecCacheData_(),
        id_(id),
        basisFunctionSetPositions_(basisFunctionSetPositions)
    {
        orderCoeff_ = 2;// for safety
    	size_t numberOfBasisFunctions = 0;
    	for (size_t i = 0; i < basisFunctionSetPositions_.size(); ++i)
    	{
    		numberOfBasisFunctions += basisFunctionSet_->at(basisFunctionSetPositions_[i])->size();
    	}
    	setNumberOfBasisFunctions(numberOfBasisFunctions);
        setQuadratureRulesWithOrder(orderCoeff_ * basisFunctionSet_->at(basisFunctionSetPositions_[0])->getOrder()+1);
        nrOfDOFinTheElement_=basisFunctionSet_->at(basisFunctionSetPositions_[0])->size();
        facesList_.assign(getReferenceGeometry()->getNrOfCodim1Entities(),nullptr);
        if (getReferenceGeometry()->getNrOfCodim3Entities() > 0)
        {
            edgesList_.assign(getReferenceGeometry()->getNrOfCodim2Entities(),nullptr);
        }
        nodesList_.assign(getReferenceGeometry()->getNumberOfNodes(),nullptr);
    }

    Element::Element(const Element& other):
        ElementGeometryT(other),
        ElementDataT(other),
        basisFunctionSet_(other.basisFunctionSet_),
        quadratureRule_(other.quadratureRule_),
        vecCacheData_(other.vecCacheData_),
        id_(other.id_),
        orderCoeff_(other.orderCoeff_),
        basisFunctionSetPositions_(other.basisFunctionSetPositions_),
        nrOfDOFinTheElement_(other.nrOfDOFinTheElement_),
        facesList_(other.facesList_),
        edgesList_(other.edgesList_),
        nodesList_(other.nodesList_)
            
    {
        std::cout << "In the copy constructor of Element " << std::endl;
    }
    
    ///Very ugly default constructor that's only here because it is needed in
    ///ShortTermStorageElementBase. 
    ///\todo Remove the default constructor everywhere.
    Element::Element():
    ElementDataT(0, 0, 0, 0, 0)
    {
    }

    Element::~Element()
    {
        
    }
    
    void Element::setDefaultBasisFunctionSet(unsigned int position)
    {
    	basisFunctionSetPositions_.resize(1,-1);
    	basisFunctionSetPositions_[0]=position;
    	int numberOfBasisFunctions(0);
    	for(int i:basisFunctionSetPositions_)
    	{
    		if(i!=-1) numberOfBasisFunctions+=basisFunctionSet_->at(i)->size();
    	}
    	setNumberOfBasisFunctions(numberOfBasisFunctions);
        setQuadratureRulesWithOrder(orderCoeff_ * basisFunctionSet_->at(position)->getOrder()+1);
        nrOfDOFinTheElement_=basisFunctionSet_->at(position)->size();
    }

    void Element::setFaceBasisFunctionSet(unsigned int position,int localFaceIndex)
    {
    	if(basisFunctionSetPositions_.size()<1+getNrOfFaces()){
			basisFunctionSetPositions_.resize(1+getNrOfFaces(),-1);
		}
    	basisFunctionSetPositions_[1+localFaceIndex]=position;
    	int numberOfBasisFunctions(0);
    	for(int i:basisFunctionSetPositions_)
    	{
    		if(i!=-1) numberOfBasisFunctions+=basisFunctionSet_->at(i)->size();
    	}
    	setNumberOfBasisFunctions(numberOfBasisFunctions);
    }

    void Element::setEdgeBasisFunctionSet(unsigned int position, int localEdgeIndex)
    {
    	if(basisFunctionSetPositions_.size()<1+getNrOfFaces()+getNrOfEdges()){
			basisFunctionSetPositions_.resize(1+getNrOfFaces()+getNrOfEdges(),-1);
		}
    	basisFunctionSetPositions_[1+getNrOfFaces()+localEdgeIndex]=position;
    	int numberOfBasisFunctions(0);
    	for(int i:basisFunctionSetPositions_)
    	{
    		if(i!=-1) numberOfBasisFunctions+=basisFunctionSet_->at(i)->size();
    	}
    	setNumberOfBasisFunctions(numberOfBasisFunctions);
    }

    void Element::setVertexBasisFunctionSet(unsigned int position, int localVertexIndex)
    {
    	if(basisFunctionSetPositions_.size()<1+getNrOfFaces()+getNrOfEdges()+getNrOfNodes()){
			basisFunctionSetPositions_.resize(1+getNrOfFaces()+getNrOfEdges()+getNrOfNodes(),-1);
		}
    	basisFunctionSetPositions_[1+getNrOfFaces()+getNrOfEdges()+localVertexIndex]=position;
    	int numberOfBasisFunctions(0);
    	for(int i:basisFunctionSetPositions_)
    	{
    		if(i!=-1) numberOfBasisFunctions+=basisFunctionSet_->at(i)->size();
    	}
    	setNumberOfBasisFunctions(numberOfBasisFunctions);
    }

    double Element::basisFunctionDeriv(unsigned int i, unsigned int jDir, const PointReferenceT& p)const
    {
        TestErrorDebug((jDir<p.size()),"Error in BasisFunctionSet.EvalDeriv: invalid derivative direction!");

        /*if (jDir>= DIM)
            return -1.e50;
        else*/
        int basePosition(0);
        for(int j:basisFunctionSetPositions_){
            if(j!=-1){
                int n=basisFunctionSet_->at(j)->size();
                if(i-basePosition<n){
                    return basisFunctionSet_->at(j)->evalDeriv(i-basePosition, jDir, p);
                }else{
                    basePosition+=n;
                }
            }
        }
        throw "in basisFunctionDeriv(jdir): asked for a basisFunction that doesn't exist!";
    }
    
    double Element::basisFunction(unsigned int i, const PointReferenceT& p)const
    {
    	const Base::BaseBasisFunction* function;
        int basePosition(0);
        for(int j:basisFunctionSetPositions_)
        {
        	if(j!=-1)
            {
        		int n=basisFunctionSet_->at(j)->size();
				if(i-basePosition<n)
                {
					function=basisFunctionSet_->at(j)->operator[](i-basePosition);
					basePosition+=n;
				}
                else
                {
					basePosition+=n;
				}
        	}
        }
    	return getReferenceGeometry()->getBasisFunctionValue(function,p);
        //return basisFunctionSet_->eval(i,p);
    }
    
    unsigned int Element::getID()const
    {
        return id_;
    }
    
    unsigned int Element::getID()
    {
        return id_;
    }
    
    void Element::setQuadratureRulesWithOrder(unsigned int quadrROrder)
    {
        quadratureRule_ =  Geometry::ElementGeometry::referenceGeometry_->getGaussQuadratureRule(quadrROrder);
    }
    
    void Element::setGaussQuadratureRule(GaussQuadratureRuleT* const quadR)
    {
        quadratureRule_ = quadR;
    }
    
    const QuadratureRules::GaussQuadratureRule* Element::getGaussQuadratureRule() const
    {
        return quadratureRule_;
    }
    
    std::vector<Base::ElementCacheData >& Element::getVecCacheData()
    {
        return vecCacheData_;
    }
    
    
    void Element::getSolution(unsigned int timeLevel, const PointReferenceT& p, SolutionVector& solution) const
    {
        size_t numberOfUnknows = ElementData::getNrOfUnknows();
        solution.resize(numberOfUnknows);
        
        for(size_t k = 0; k < numberOfUnknows; ++k)
        {
            solution[k] = 0;
        }
        
        const LinearAlgebra::NumericalVector& data = ElementData::getTimeLevelData(0);
        for (size_t i=0; i < ElementData::getNrOfBasisFunctions(); ++i)
        {
            for (size_t k=0; k < numberOfUnknows; ++k)
            {
                solution[k] += data(i) * basisFunction(i, p);
            }
        }
    }
    
    void Element::basisFunction(unsigned int i, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret) const
    {
        int basePosition(0);
        for(int j:basisFunctionSetPositions_){
        	if(j!=-1){
        		int n=basisFunctionSet_->at(j)->size();
				if(i-basePosition<n){
					basisFunctionSet_->at(j)->eval(i-basePosition,p,ret);
					return;
				}else{
					basePosition+=n;
				}
        	}
        }
        throw "in basisFunction: asked for a basisFunction that doesn't exist!";
    }
    
    void Element::basisFunctionCurl(unsigned int i, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret) const
    {
        int basePosition(0);
        for(int j:basisFunctionSetPositions_)
        {
        	if(j!=-1){
        		int n=basisFunctionSet_->at(j)->size();
				if(i-basePosition<n)
                {
					basisFunctionSet_->at(j)->evalCurl(i-basePosition,p,ret);
					return;
				}
                else
                {
					basePosition+=n;
				}
        	}
        }
        throw "in basisFunctionCurl: asked for a basisFunction that doesn't exist!";
    }

    void Element::basisFunctionDeriv(unsigned int i, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret,const Element* wrapper) const
    {
    	if(wrapper==nullptr)
        {
    		wrapper=this;//Apparently you can't default to this
    	}
    	const Base::BaseBasisFunction* function;
        int basePosition(0);
        for(int j:basisFunctionSetPositions_){
        	if(j!=-1){
        		unsigned int n=basisFunctionSet_->at(j)->size();
				if(i-basePosition<n){
					function=basisFunctionSet_->at(j)->operator[](i-basePosition);
					basePosition+=n;
				}else{
					basePosition+=n;
				}
        	}
        }
    	Utilities::PhysGradientOfBasisFunction functionGradient(wrapper,function);
    	functionGradient(p,ret);
    }

#ifndef NDEBUG
    const Base::BaseBasisFunction* Element::getBasisFunction(int i) const
    {
    	int basePosition(0);
    	for(int j:basisFunctionSetPositions_){
    		if(j!=-1){
    			if(i-basePosition<basisFunctionSet_->at(j)->size()){
    				return basisFunctionSet_->at(j)->operator[](i-basePosition);
    			}else{
					basePosition+=basisFunctionSet_->at(j)->size();
    			}
    		}
    	}
        throw "in getBasisFunction(): asked for a basisFunction that doesn't exist!";
    }
#endif

    void Element::setFace(int localFaceNr, const Face* face)
    {
    	TestErrorDebug((face->getPtrElementLeft()==this&&face->localFaceNumberLeft()==localFaceNr)||
    			       (face->getPtrElementRight()==this&&face->localFaceNumberRight()==localFaceNr),"You are only allowed to set a face to a local face index that matches");
    	if(facesList_.size()<localFaceNr+1)
    	{
            //should not happen
    		facesList_.resize(localFaceNr+1);
    	}
    	facesList_[localFaceNr]=face;
    }

    void Element::setEdge(int localEdgeNr, const Edge* edge)
    {
    	if(edgesList_.size()<localEdgeNr+1)
    	{
            //could happen in 4D
    		edgesList_.resize(localEdgeNr+1);
    	}
    	edgesList_[localEdgeNr]=edge;
    }

    void Element::setNode(int localNodeNr, const Node* node)
    {
    	if(nodesList_.size()<localNodeNr+1)
    	{
            //should not happen
    		nodesList_.resize(localNodeNr+1);
    	}
    	nodesList_[localNodeNr]=node;
    }
    
    ///Function that computes the mass matrix. First resize the mass matrix to 
    ///the correct size, then for all quadrature points, compute the values of 
    ///all the products of basisfunctions and add this with the appropriate weight
    ///to the mass matrix.
    void Element::computeMassMatrix()
    {
        //get the number of basisfunctions, dimension and number of quadrature 
        //points on this element.
        size_t numBasisFuncs = getNrOfBasisFunctions();
        size_t dim = quadratureRule_->dimension();
        size_t numQuadPoints = quadratureRule_->nrOfPoints();
        
        //make the mass matrix of the correct size and set all entries to zero.
        massMatrix_.resize(numBasisFuncs,numBasisFuncs);
        massMatrix_ *= 0;
        
        //declare the relevant auxiliary variables
        LinearAlgebra::Matrix tempMatrix(numBasisFuncs,numBasisFuncs);        
        PointReferenceT p(dim);
        Geometry::Jacobian jac(dim,dim);        
        
        //for each quadrature point, compute the value of the product of the 
        //basisfunctions, then add it with the correct weight to massMatrix_
        for(size_t pIndex = 0; pIndex < numQuadPoints; ++pIndex)
        {
            quadratureRule_->getPoint(pIndex,p);
            calcJacobian(p, jac);
            for (size_t i = 0; i < numBasisFuncs; ++i)
            {
                for (size_t j = 0; j < numBasisFuncs; ++j)
                {
                    tempMatrix(i,j) = basisFunction(i,p) * basisFunction(j,p);
                }
            }
            massMatrix_.axpy((quadratureRule_->weight(pIndex))*std::abs(jac.determinant()),tempMatrix);
        }
    }
}
#endif

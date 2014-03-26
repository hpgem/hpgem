//
//  Element_Impl.hpp
//  
//
//  Created by Shavarsh Nurijanyan on 2/25/13.
//
//

#ifndef _Element_Impl_hpp
#define _Element_Impl_hpp

#include "Base/Element.hpp"
#include "PhysGradientOfBasisFunction.hpp"
#include "edge.hpp"
#include "Face.hpp"

namespace Base
{
    //class Element;
    
    Element::Element(const VectorOfPointIndexesT& globalNodeIndexes,
                          const std::vector<const BasisFunctionSetT*>* basisFunctionSet,
                          const VectorOfPhysicalPointsT& allNodes,
                          unsigned int nrOfUnkowns,
                          unsigned int nrOfTimeLevels,
                          unsigned int nrOfBasisFunc,
                          unsigned int id,
                          unsigned int numberOfElementMatrixes,
                          unsigned int numberOfElementVectors,
                          const std::vector< int>& basisFunctionSetPositions):
        ElementGeometryT(globalNodeIndexes, allNodes),
        ElementDataT(nrOfTimeLevels, nrOfUnkowns, nrOfBasisFunc,numberOfElementMatrixes,numberOfElementVectors),
        basisFunctionSet_(basisFunctionSet),
        quadratureRule_(NULL),
        vecCacheData_(),
        id_(id),
        basisFunctionSetPositions_(basisFunctionSetPositions)
    {
        orderCoeff_ = 2;// for safety
    	int numberOfBasisFunctions(0);
    	for(int i=0;i<basisFunctionSetPositions_.size();++i)
    	{
    		numberOfBasisFunctions+=basisFunctionSet_->at(basisFunctionSetPositions_[i])->size();
    	}
    	setNumberOfBasisFunctions(numberOfBasisFunctions);
        setQuadratureRulesWithOrder(orderCoeff_ * basisFunctionSet_->at(basisFunctionSetPositions_[0])->getOrder()+1);
        nrOfDOFinTheElement_=basisFunctionSet_->at(basisFunctionSetPositions_[0])->size();
        nrOfDOFperVertex_=0;
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
        nrOfDOFperVertex_(other.nrOfDOFperVertex_)
    {
        std::cout << "In the copy constructor of Element " << std::endl;
    }

    Element::~Element()
    {
        
    }
//    template<unsigned int DIM>
//    unsigned int
//    Element<DIM>::getNumberOfDegreesOfFreedom()
//    {
//        return basisFunctionSet_.size();
//    }
//    
//    template<unsigned int DIM>
//    unsigned int 
//    Element<DIM>::getNumberOfDegreesOfFreedom()const
//    {
//        return basisFunctionSet_.size();
//    }
    
    void
    Element::setDefaultBasisFunctionSet(unsigned int position)
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
        nrOfDOFperVertex_=0;
    }

    void
    Element::setFaceBasisFunctionSet(unsigned int position,int localFaceIndex)
    {
    	if(basisFunctionSetPositions_.size()<1+physicalGeometry_->getNrOfFaces()){
			basisFunctionSetPositions_.resize(1+physicalGeometry_->getNrOfFaces(),-1);
		}
    	basisFunctionSetPositions_[1+localFaceIndex]=position;
    	int numberOfBasisFunctions(0);
    	for(int i:basisFunctionSetPositions_)
    	{
    		if(i!=-1) numberOfBasisFunctions+=basisFunctionSet_->at(i)->size();
    	}
    	setNumberOfBasisFunctions(numberOfBasisFunctions);
    }

    void
    Element::setEdgeBasisFunctionSet(unsigned int position, int localEdgeIndex)
    {
    	if(basisFunctionSetPositions_.size()<1+physicalGeometry_->getNrOfFaces()+getNrOfEdges()){
			basisFunctionSetPositions_.resize(1+physicalGeometry_->getNrOfFaces()+getNrOfEdges(),-1);
		}
    	basisFunctionSetPositions_[1+physicalGeometry_->getNrOfFaces()+localEdgeIndex]=position;
    	int numberOfBasisFunctions(0);
    	for(int i:basisFunctionSetPositions_)
    	{
    		if(i!=-1) numberOfBasisFunctions+=basisFunctionSet_->at(i)->size();
    	}
    	setNumberOfBasisFunctions(numberOfBasisFunctions);
    }

    void
    Element::setVertexBasisFunctionSet(unsigned int position, int localVertexIndex)
    {
    	if(basisFunctionSetPositions_.size()<1+physicalGeometry_->getNrOfFaces()+getNrOfEdges()+getNrOfNodes()){
			basisFunctionSetPositions_.resize(1+physicalGeometry_->getNrOfFaces()+getNrOfEdges()+getNrOfNodes(),-1);
		}
    	basisFunctionSetPositions_[1+physicalGeometry_->getNrOfFaces()+getNrOfEdges()+localVertexIndex]=position;
    	int numberOfBasisFunctions(0);
    	for(int i:basisFunctionSetPositions_)
    	{
    		if(i!=-1) numberOfBasisFunctions+=basisFunctionSet_->at(i)->size();
    	}
    	setNumberOfBasisFunctions(numberOfBasisFunctions);
    	nrOfDOFperVertex_=basisFunctionSet_->at(position)->size();
    }

    double
    Element::basisFunctionDeriv(unsigned int i, unsigned int jDir, const PointReferenceT& p)const
    {
        TestErrorDebug((jDir<p.size()),"Error in BasisFunctionSet.EvalDeriv: invalid derivative direction!");

        /*if (jDir>= DIM)
            return -1.e50;
        else*/
        int basePosition(0);
        for(int j:basisFunctionSetPositions_){
        	if(j!=-1){
				if(i-basePosition<basisFunctionSet_->at(j)->size()){
					return basisFunctionSet_->at(j)->evalDeriv(i-basePosition, jDir, p);
				}else{
					basePosition+=basisFunctionSet_->at(j)->size();
				}
        	}
        }
        throw "in basisFunctionDeriv(jdir): asked for a basisFunction that doesn't exist!";
    }
    
    double
    Element::basisFunction(unsigned int i, const PointReferenceT& p)const
    {
    	const Base::BaseBasisFunction* function;
        int basePosition(0);
        for(int j:basisFunctionSetPositions_){
        	if(j!=-1){
				if(i-basePosition<basisFunctionSet_->at(j)->size()){
					function=basisFunctionSet_->at(j)->operator[](i-basePosition);
					basePosition+=basisFunctionSet_->at(j)->size();
				}else{
					basePosition+=basisFunctionSet_->at(j)->size();
				}
        	}
        }
    	return getReferenceGeometry()->getBasisFunctionValue(function,p);
        //return basisFunctionSet_->eval(i,p);
    }
    
    unsigned int
    Element::getID()const
    {
        return id_;
    }
    
    unsigned int
    Element::getID()
    {
        return id_;
    }
    
    void
    Element::setQuadratureRulesWithOrder(unsigned int quadrROrder)
    {
        quadratureRule_ =  Geometry::ElementGeometry::referenceGeometry_->getGaussQuadratureRule(quadrROrder);
    }
    
    void
    Element::setGaussQuadratureRule(GaussQuadratureRuleT* const quadR)
    {
        quadratureRule_ = quadR;
    }
    
    const QuadratureRules::GaussQuadratureRule*
    Element::getGaussQuadratureRule() const
    {
        return quadratureRule_;
    }
    
    std::vector<Base::ElementCacheData >&
    Element::getVecCacheData()
    {
        return vecCacheData_;
    }
    
    
    void
    Element::getSolution(unsigned int timeLevel, const PointReferenceT& p, SolutionVector& solution) const
    {
        unsigned int numberOfUnknows = ElementData::getNrOfUnknows();
        solution.resize(numberOfUnknows);
        
        const LinearAlgebra::Matrix& data = ElementData::getTimeLevelData(0);
        for (int i=0; i < ElementData::getNrOfBasisFunctions(); ++i)
        {
            for (int k=0; k < numberOfUnknows; ++k)
            {
                solution[k] += data(k, i) * basisFunction(i, p);
            }
        }
    }
    
    void
    Element::basisFunction(unsigned int i, const PointReferenceT& p, NumericalVector& ret) const
    {
        int basePosition(0);
        for(int j:basisFunctionSetPositions_){
        	if(j!=-1){
				if(i-basePosition<basisFunctionSet_->at(j)->size()){
					basisFunctionSet_->at(j)->eval(i-basePosition,p,ret);
					return;
				}else{
					basePosition+=basisFunctionSet_->at(j)->size();
				}
        	}
        }
        throw "in basisFunctionDeriv(jdir): asked for a basisFunction that doesn't exist!";
    }
    
    void
    Element::basisFunctionCurl(unsigned int i, const PointReferenceT& p, NumericalVector& ret) const
    {
        int basePosition(0);
        for(int j:basisFunctionSetPositions_){
        	if(j!=-1){
				if(i-basePosition<basisFunctionSet_->at(j)->size()){
					basisFunctionSet_->at(j)->evalCurl(i-basePosition,p,ret);
					return;
				}else{
					basePosition+=basisFunctionSet_->at(j)->size();
				}
        	}
        }
        throw "in basisFunctionDeriv(jdir): asked for a basisFunction that doesn't exist!";
    }

    void
    Element::basisFunctionDeriv(unsigned int i, const PointReferenceT& p, NumericalVector& ret) const
    {
    	const Base::BaseBasisFunction* function;
        int basePosition(0);
        for(int j:basisFunctionSetPositions_){
        	if(j!=-1){
				if(i-basePosition<basisFunctionSet_->at(j)->size()){
					function=basisFunctionSet_->at(j)->operator[](i-basePosition);
					basePosition+=basisFunctionSet_->at(j)->size();
				}else{
					basePosition+=basisFunctionSet_->at(j)->size();
				}
        	}
        }
    	Utilities::PhysGradientOfBasisFunction functionGradient(this,function);
    	functionGradient(p,ret);
    }

#ifndef NDEBUG
    const Base::BaseBasisFunction*
    Element::getBasisFunction(int i)
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

    void
    Element::setFace(int localFaceNr, const Face* face)
    {
    	TestErrorDebug((face->getPtrElementLeft()==this&&face->localFaceNumberLeft()==localFaceNr)||
    			       (face->getPtrElementRight()==this&&face->localFaceNumberRight()==localFaceNr),"You are only allowed to set a face to a local face index that matches");
    	if(facesList_.size()<localFaceNr+1)
    	{
    		facesList_.resize(localFaceNr+1);
    	}
    	facesList_[localFaceNr]=face;
    }

    void
    Element::setEdge(int localEdgeNr, const Edge* edge)
    {
    	if(edgesList_.size()<localEdgeNr+1)
    	{
    		edgesList_.resize(localEdgeNr+1);
    	}
    	edgesList_[localEdgeNr]=edge;
    }
}
#endif

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

namespace Base
{
    
    class Face;
    
    Face::Face(ElementT*  ptrElemL, const LocalFaceNrTypeT& localFaceNumL, ElementT* ptrElemR, const LocalFaceNrTypeT& localFaceNumR,int faceID,unsigned int numberOfFaceMatrixes, unsigned int numberOfFaceVectors) :
        FaceGeometryT((ElementGeometryT*)ptrElemL, localFaceNumL,(ElementGeometryT*)ptrElemR,localFaceNumR),
        elementLeft_(ptrElemL),
        elementRight_(ptrElemR),
        faceID_(faceID),
    	FaceData(ptrElemL->getNrOfBasisFunctions()*ptrElemL->getNrOfUnknows()+
    			 ptrElemR->getNrOfBasisFunctions()*ptrElemR->getNrOfUnknows(),numberOfFaceMatrixes,numberOfFaceVectors)
    {
    	TestErrorDebug(ptrElemR!=NULL,"Error: passing a boundary face to the constructor for internal faces!");
        createQuadratureRules();
        ptrElemL->setFace(localFaceNumL,this);
        ptrElemR->setFace(localFaceNumR,this);
    }
    
    Face::Face(ElementT* ptrElemL, const LocalFaceNrTypeT& localFaceNumL, const Geometry::FaceType&  faceType,int faceID,unsigned int numberOfFaceMatrixes, unsigned int numberOfFaceVectors):
        FaceGeometryT((ElementGeometryT*)ptrElemL, localFaceNumL, faceType),
        elementLeft_(ptrElemL),
        elementRight_(NULL),
        faceID_(faceID),
    	FaceData(ptrElemL->getNrOfBasisFunctions()*ptrElemL->getNrOfUnknows(),numberOfFaceMatrixes,numberOfFaceVectors)
    {
        createQuadratureRules();
        ptrElemL->setFace(localFaceNumL,this);
    }

    
   /* void
    Face::setPtrElementLeft( ElementT* left)
    {
        elementLeft_ = left;
        FaceGeometry::setPtrElementLeft(left);
        //setElementGLeft(left);
    }
        
    void
    Face::setPtrElementRight( ElementT* right)
    {
        elementRight_ = right;
        FaceGeometry::setPtrElementRight(right);
        //setElementGRight(right);
    }*/
    
    bool
    Face::isInternal()const
    {
        bool internal;
        if((FaceGeometryT::faceType_==Geometry::INTERNAL && elementRight_==NULL)||(FaceGeometryT::faceType_!=Geometry::INTERNAL && elementRight_!=NULL))
            std::cout << "Something wrong with boundaries";
            
            
        return (FaceGeometryT::faceType_==Geometry::INTERNAL? true: false);
    }
    
    void
    Face::createQuadratureRules()
    {
        unsigned int rightOrder = (elementRight_==NULL? 0 : elementRight_->getGaussQuadratureRule()->order());
        unsigned int leftOrder  = elementLeft_->getGaussQuadratureRule()->order();
        if (leftOrder >=rightOrder)
        {
             quadratureRule_ = elementLeft_->getReferenceGeometry()->getCodim1ReferenceGeometry(FaceGeometryT::localFaceNumberLeft_)->getGaussQuadratureRule(leftOrder);
        }
         else
        {
            std::cout << "again..... Face<DIM>::createQuadratureRules(): " << leftOrder << " " << rightOrder << std::endl;
            quadratureRule_ = elementRight_->getReferenceGeometry()->getCodim1ReferenceGeometry(FaceGeometryT::localFaceNumberRight_)->getGaussQuadratureRule(rightOrder);
        }
        
    }

    int
    Face::getNrOfBasisFunctions() const
    {
    	if(isInternal())
    	{
    		return getPtrElementLeft()->getNrOfBasisFunctions()+getPtrElementRight()->getNrOfBasisFunctions();
    	}
    	else
    	{
    		return getPtrElementLeft()->getNrOfBasisFunctions();
    	}
    }

    double
    Face::basisFunction(unsigned int i, const Geometry::PointReference& p) const
    {
    	Geometry::PointReference pElement(p.size()+1);
    	int n(getPtrElementLeft()->getNrOfBasisFunctions());
    	if(i<n){
    		mapRefFaceToRefElemL(p,pElement);
    		return getPtrElementLeft()->basisFunction(i,pElement);
    	}else{
    		mapRefFaceToRefElemR(p,pElement);
    		return getPtrElementRight()->basisFunction(i-n,pElement);
    	}
    }

    void
    Face::basisFunction(unsigned int i, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret) const
    {
    	Geometry::PointReference pElement(p.size()+1);
    	int n(getPtrElementLeft()->getNrOfBasisFunctions());
    	if(i<n){
    		mapRefFaceToRefElemL(p,pElement);
    		getPtrElementLeft()->basisFunction(i,pElement,ret);
    	}else{
    		mapRefFaceToRefElemR(p,pElement);
    		getPtrElementRight()->basisFunction(i-n,pElement,ret);
    	}
    }

    void
    Face::basisFunctionNormal(unsigned int i, const LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret) const
    {
    	Geometry::PointReference pElement(p.size()+1);
    	int n(getPtrElementLeft()->getNrOfBasisFunctions());
    	if(i<n){
    		mapRefFaceToRefElemL(p,pElement);
    		ret=normal;
    		ret*=getPtrElementLeft()->basisFunction(i,pElement)/Base::L2Norm(normal);
    	}else{
    		mapRefFaceToRefElemR(p,pElement);
    		ret=normal;
    		ret*=-getPtrElementRight()->basisFunction(i-n,pElement)/Base::L2Norm(normal);
    	}
    }

    double
    Face::basisFunctionDeriv(unsigned int i, unsigned int jDir, const Geometry::PointReference& p) const
    {
    	Geometry::PointReference pElement(p.size()+1);
    	int n(getPtrElementLeft()->getNrOfBasisFunctions());
    	if(i<n){
    		mapRefFaceToRefElemL(p,pElement);
    		return getPtrElementLeft()->basisFunctionDeriv(i,jDir,pElement);
    	}else{
    		mapRefFaceToRefElemR(p,pElement);
    		return getPtrElementRight()->basisFunctionDeriv(i-n,jDir,pElement);
    	}
    }

    void
    Face::basisFunctionDeriv(unsigned int i, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret) const
    {
    	Geometry::PointReference pElement(p.size()+1);
    	int n(getPtrElementLeft()->getNrOfBasisFunctions());
    	if(i<n){
    		mapRefFaceToRefElemL(p,pElement);
    		getPtrElementLeft()->basisFunctionDeriv(i,pElement,ret);
    	}else{
    		mapRefFaceToRefElemR(p,pElement);
    		getPtrElementRight()->basisFunctionDeriv(i-n,pElement,ret);
    	}
    }

    void
    Face::basisFunctionCurl(unsigned int i, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret) const
    {
    	Geometry::PointReference pElement(p.size()+1);
    	int n(getPtrElementLeft()->getNrOfBasisFunctions());
    	if(i<n){
    		mapRefFaceToRefElemL(p,pElement);
    		getPtrElementLeft()->basisFunctionCurl(i,pElement,ret);
    	}else{
    		mapRefFaceToRefElemR(p,pElement);
    		getPtrElementRight()->basisFunctionCurl(i-n,pElement,ret);
    	}
    }
};

//
//  Face.cpp
//  
//
//  Created by Shavarsh Nurijanyan on 2/9/13.
//
//

#include "Face.hpp"
#include "Element.hpp"

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
            cout << "Something wrong with boundaries";
            
            
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
    Face::basisFunction(unsigned int i, const Geometry::PointReference& p, NumericalVector& ret) const
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
    Face::basisFunctionNormal(unsigned int i, const NumericalVector& normal, const Geometry::PointReference& p, NumericalVector& ret) const
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
    Face::basisFunctionDeriv(unsigned int i, const Geometry::PointReference& p, NumericalVector& ret) const
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
    Face::basisFunctionCurl(unsigned int i, const Geometry::PointReference& p, NumericalVector& ret) const
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

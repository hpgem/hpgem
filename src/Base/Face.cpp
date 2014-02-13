//
//  Face.cpp
//  
//
//  Created by Shavarsh Nurijanyan on 2/9/13.
//
//

#include "Face.hpp"

namespace Base
{
    
    class Face;
    
    Face::Face(ElementT*  ptrElemL, const LocalFaceNrTypeT& localFaceNumL, ElementT* ptrElemR, const LocalFaceNrTypeT& localFaceNumR,unsigned int numberOfFaceMatrixes, unsigned int numberOfFaceVectors) :
        FaceGeometryT((ElementGeometryT*)ptrElemL, localFaceNumL,(ElementGeometryT*)ptrElemR,localFaceNumR),
        elementLeft_(ptrElemL),
        elementRight_(ptrElemR),
    	FaceData(ptrElemL->getNrOfBasisFunctions()*ptrElemL->getNrOfUnknows()+
    			 ptrElemR->getNrOfBasisFunctions()*ptrElemR->getNrOfUnknows(),numberOfFaceMatrixes,numberOfFaceVectors)
    {
    	//TestErrorDebug(ptrElemR!=NULL,"Error: passing a boundary face to the constructor for internal faces!");
        createQuadratureRules();
	
    }
    
    Face::Face(ElementT* ptrElemL, const LocalFaceNrTypeT& localFaceNumL, const Geometry::FaceType&  faceType,unsigned int numberOfFaceMatrixes, unsigned int numberOfFaceVectors):
        FaceGeometryT((ElementGeometryT*)ptrElemL, localFaceNumL, faceType),
        elementLeft_(ptrElemL),
        elementRight_(NULL),
    	FaceData(ptrElemL->getNrOfBasisFunctions()*ptrElemL->getNrOfUnknows(),numberOfFaceMatrixes,numberOfFaceVectors)
    {
        createQuadratureRules();
    }

    
    void
    Face::setPtrElementLeft(const ElementT* left)
    {
        elementLeft_ = left;
        FaceGeometry::setPtrElementLeft(left);
        //setElementGLeft(left);
    }
        
    void
    Face::setPtrElementRight(const ElementT* right)
    {
        elementRight_ = right;
        FaceGeometry::setPtrElementRight(right);
        //setElementGLeft(right);
    }
    
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
};

//
//  Face.cpp
//  
//
//  Created by Shavarsh Nurijanyan on 2/9/13.
//
//

namespace Base
{
    
    template<unsigned int DIM>
    class Face;
    
    template<unsigned int DIM>
    Face<DIM>::Face(ElementT*  ptrElemL, const LocalFaceNrTypeT& localFaceNumL, ElementT* ptrElemRight, const LocalFaceNrTypeT& localFaceNumR) :
        Geometry::FaceGeometry<DIM>((ElementGeometryT*)ptrElemL, localFaceNumL,(ElementGeometryT*)ptrElemRight,localFaceNumR),
        elementLeft_(ptrElemL),
        elementRight_(ptrElemL)
    {
       createQuadratureRules();
    }
    
    template<unsigned int DIM>
    Face<DIM>::Face(ElementT* ptrElemL, const LocalFaceNrTypeT& localFaceNumL, const Geometry::FaceType&  faceType):
        Geometry::FaceGeometry<DIM>((ElementGeometryT*)ptrElemL, localFaceNumL, faceType),
        elementLeft_(ptrElemL),
        elementRight_(ptrElemL)
    {
       createQuadratureRules();
    }

    
    template<unsigned int DIM>
    void
    Face<DIM>::setPtrElementLeft(const ElementT* left)
    {
        elementLeft_ = left;
        setElementGLeft(left);
    }
        
    template<unsigned int DIM>
    void
    Face<DIM>::setPtrElementRight(const ElementT* right)
    {
        elementRight_ = right;
        setElementGLeft(right);
    }
    
    template<unsigned int DIM>
    void
    Face<DIM>::createQuadratureRules()
    {
        unsigned int rightOrder = (elementRight_==NULL? 0:elementRight_->getGaussQuadratureRule()->order());
        unsigned int leftOrder  = elementLeft_->getGaussQuadratureRule()->order();
        std::cout << "Face<DIM>::createQuadratureRules(): " << leftOrder << " " << rightOrder << std::endl;
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
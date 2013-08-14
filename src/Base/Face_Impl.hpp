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
    Face<DIM>::Face(ElementT*  ptrElemL, const LocalFaceNrTypeT& localFaceNumL, ElementT* ptrElemR, const LocalFaceNrTypeT& localFaceNumR) :
        FaceGeometryT((ElementGeometryT*)ptrElemL, localFaceNumL,(ElementGeometryT*)ptrElemR,localFaceNumR),
        elementLeft_(ptrElemL),
        elementRight_(ptrElemR)
    {
        createQuadratureRules();
    }
    
    template<unsigned int DIM>
    Face<DIM>::Face(ElementT* ptrElemL, const LocalFaceNrTypeT& localFaceNumL, const Geometry::FaceType&  faceType):
        FaceGeometryT((ElementGeometryT*)ptrElemL, localFaceNumL, faceType),
        elementLeft_(ptrElemL),
        elementRight_(NULL)
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
    bool
    Face<DIM>::isInternal()const
    {
        bool internal;
        if((FaceGeometryT::faceType_==Geometry::INTERNAL && elementRight_==NULL)||(FaceGeometryT::faceType_!=Geometry::INTERNAL && elementRight_!=NULL))
            cout << "Something wrong with boundaries";
            
            
        return (FaceGeometryT::faceType_==Geometry::INTERNAL? true: false);
    }
    
    template<unsigned int DIM>
    void
    Face<DIM>::createQuadratureRules()
    {
        if (elementRight_->getGaussQuadratureRule()==NULL)
            cout <<"MOTHERFUCKERRRR!"<<endl;
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
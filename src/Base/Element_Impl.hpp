//
//  Element_Impl.hpp
//  
//
//  Created by Shavarsh Nurijanyan on 2/25/13.
//
//

#ifndef _Element_Impl_hpp
#define _Element_Impl_hpp

namespace Base
{
    template<unsigned int DIM>
    class Element;
    
    template<unsigned int DIM>
    Element<DIM>::Element(const VectorOfPointIndexesT& globalNodeIndexes,
                          const BasisFunctionSetT* const& basisFunctionSet,
                          const VectorOfPhysicalPointsT& allNodes,
                          unsigned int nrOfUnkowns,
                          unsigned int nrOfTimeLevels,
                          unsigned int counter):
        ElementGeometryT(globalNodeIndexes, allNodes),
        ElementDataT(nrOfTimeLevels, nrOfUnkowns, basisFunctionSet->size()),
        basisFunctionSet_(basisFunctionSet),
        quadratureRule_(NULL),
        vecCacheData_(),
        id_(counter)
    {
        orderCoeff_ = 2;
        setQuadratureRulesWithOrder(orderCoeff_ * basisFunctionSet_->getOrder());
    }

    template<unsigned int DIM>
    Element<DIM>::Element(const Element& other):
        ElementGeometryT(other),
        ElementDataT(other),
        basisFunctionSet_(other.basisFunctionSet_),
        quadratureRule_(other.quadratureRule_),
        vecCacheData_(other.vecCacheData_),
        id_(other.id_)
    {
        std::cout << "In the copy constructor of Elemenet " << std::endl;
    }

    template<unsigned int DIM>
    Element<DIM>::~Element()
    {
        
    }
    
    template<unsigned int DIM>
    unsigned int
    Element<DIM>::getID()const
    {
        return id_;
    }
    
    template<unsigned int DIM>
    unsigned int
    Element<DIM>::getID()
    {
        return id_;
    }
    
    template<unsigned int DIM>
    void
    Element<DIM>::setQuadratureRulesWithOrder(unsigned int quadrROrder)
    {
        quadratureRule_ =  Geometry::ElementGeometry<DIM>::referenceGeometry_->getGaussQuadratureRule(quadrROrder);
    }
    
    template<unsigned int DIM>
    void
    Element<DIM>::setGaussQuadratureRule(GaussQuadratureRuleT* const quadR)
    {
        quadratureRule_ = quadR;
    }
    
    template<unsigned int DIM>
    const QuadratureRules::GaussQuadratureRule<DIM>*
    Element<DIM>::getGaussQuadratureRule() const
    {
        return quadratureRule_;
    }
    
    template<unsigned int DIM>
    std::vector<Base::ElementCacheData<DIM> >&
    Element<DIM>::getVecCacheData()
    {
        return vecCacheData_;
    }
}
#endif

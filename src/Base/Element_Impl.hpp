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
                          const BasisFunctionSetT* const basisFunctionSet,
                          const VectorOfPhysicalPointsT& allNodes,
                          unsigned int nrOfUnkowns,
                          unsigned int nrOfTimeLevels,
                          unsigned int id):
        ElementGeometryT(globalNodeIndexes, allNodes),
        ElementDataT(nrOfTimeLevels, nrOfUnkowns, basisFunctionSet->size()),
        basisFunctionSet_(basisFunctionSet),
        quadratureRule_(NULL),
        vecCacheData_(),
        id_(id)
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
    
    template<unsigned int DIM>
    double
    Element<DIM>::basisFunctionDeriv(unsigned int i, unsigned int jDir, const PointReferenceT& p)const
    {
        TestErrorDebug((jDir<DIM),"Error in BasisFunctionSet.EvalDeriv: invalid derivative direction!");

        if (jDir>= DIM)
            return -1.e50;
        else
            return basisFunctionSet_->evalDeriv(i, jDir, p);
    }
    
    template<unsigned int DIM>
    double
    Element<DIM>::basisFunction(unsigned int i, const PointReferenceT& p)const
    {
        return basisFunctionSet_->eval(i,p);
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
    
    
    template<unsigned int DIM>
    void
    Element<DIM>::getSolution(unsigned int timeLevel, const PointReferenceT& p, SolutionVector& solution) const
    {
        unsigned int numberOfUnknows = ElementData<DIM>::getNrOfUnknows();
        solution.resize(numberOfUnknows);
        
        const LinearAlgebra::Matrix& data = ElementData<DIM>::getTimeLevelData(0);
        for (int i=0; i < ElementData<DIM>::getNrOfBasisFunctions(); ++i)
        {
            for (int k=0; k < numberOfUnknows; ++k)
            {
                solution[k] += data(k, i) * basisFunction(i, p);
            }
        }
    }
}
#endif

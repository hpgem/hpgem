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

namespace Base
{
    //class Element;
    
    Element::Element(const VectorOfPointIndexesT& globalNodeIndexes,
                          const BasisFunctionSetT* const basisFunctionSet,
                          const VectorOfPhysicalPointsT& allNodes,
                          unsigned int nrOfUnkowns,
                          unsigned int nrOfTimeLevels,
                          unsigned int nrOfBasisFunc,
                          unsigned int id,
                          unsigned int numberOfElementMatrixes,
                          unsigned int numberOfElementVectors):
        ElementGeometryT(globalNodeIndexes, allNodes),
        ElementDataT(nrOfTimeLevels, nrOfUnkowns, nrOfBasisFunc,numberOfElementMatrixes,numberOfElementVectors),
        basisFunctionSet_(basisFunctionSet),
        quadratureRule_(NULL),
        vecCacheData_(),
        id_(id)
    {
        orderCoeff_ = 2;// for safety
        setQuadratureRulesWithOrder(orderCoeff_ * basisFunctionSet_->getOrder()+1);
    }

    Element::Element(const Element& other):
        ElementGeometryT(other),
        ElementDataT(other),
        basisFunctionSet_(other.basisFunctionSet_),
        quadratureRule_(other.quadratureRule_),
        vecCacheData_(other.vecCacheData_),
        id_(other.id_),
        orderCoeff_(other.orderCoeff_)
    {
        std::cout << "In the copy constructor of Elemenet " << std::endl;
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
    Element::setDefaultBasisFunctionSet(BasisFunctionSetT* bFSet)
    {
    	basisFunctionSet_=bFSet;
    	setNumberOfBasisFunctions(bFSet->size());
    }


    double
    Element::basisFunctionDeriv(unsigned int i, unsigned int jDir, const PointReferenceT& p)const
    {
        TestErrorDebug((jDir<p.size()),"Error in BasisFunctionSet.EvalDeriv: invalid derivative direction!");

        /*if (jDir>= DIM)
            return -1.e50;
        else*/
            return basisFunctionSet_->evalDeriv(i, jDir, p);
    }
    
    double
    Element::basisFunction(unsigned int i, const PointReferenceT& p)const
    {
        return basisFunctionSet_->eval(i,p);
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
        basisFunctionSet_->eval(i,p,ret);
    }
    
    void
    Element::basisFunctionCurl(unsigned int i, const PointReferenceT& p, NumericalVector& ret) const
    {
        basisFunctionSet_->evalCurl(i,p,ret);
    }

    void
    Element::basisFunctionDeriv(unsigned int i, const PointReferenceT& p, NumericalVector& ret) const
    {
    	Utilities::PhysGradientOfBasisFunction function(this,i);
    	function(p,ret);
    }
}
#endif

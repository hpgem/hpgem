//
//  BasisFunctionSet.cpp
//  
//
//  Created by Shavarsh Nurijanyan on 6/24/13.
//
//

#include "Base/BasisFunctionSet.hpp"

namespace Base {

    //class BasisFunctionSet;

    BasisFunctionSet::BasisFunctionSet(unsigned int order):
    order_(order)
    {
    }
            
    BasisFunctionSet::~BasisFunctionSet()
    {
        while(!vecOfBasisFcn_.empty())
        {
            delete vecOfBasisFcn_.back();
            vecOfBasisFcn_.pop_back();
        }
    }

    unsigned int
    BasisFunctionSet::size() const
    {
        return vecOfBasisFcn_.size();
    }

    unsigned int
    BasisFunctionSet::getOrder() const
    {
        return order_;
    }
    void
    BasisFunctionSet::addBasisFunction(BaseBasisFunctionT* bf)
    {
        vecOfBasisFcn_.push_back(bf);
    }

    double
    BasisFunctionSet::eval(unsigned int i, const PointReferenceT& p) const
    {
        return vecOfBasisFcn_[i]->eval(p);
    }
    
    /*double
    BasisFunctionSet::evalDeriv(unsigned int i, unsigned int jDir, const PointReferenceT& p) const
    {
        TestErrorDebug((jDir<1),"Error in BasisFunctionSet.EvalDeriv: invalid derivative direction!");
                   
        switch (jDir)
        {
            case 0:
                return vecOfBasisFcn_[i]->evalDeriv0(p);
                
                break;
            default: -1.e50;
        }
    }
    template<>
    inline double
    BasisFunctionSet<1>::evalDeriv(unsigned int i, unsigned int jDir, const PointReferenceT& p) const
    {
        TestErrorDebug((jDir<1),"Error in BasisFunctionSet.EvalDeriv: invalid derivative direction!");
        
        switch (jDir)
        {
            case 0:
                return vecOfBasisFcn_[i]->evalDeriv0(p);
                
                break;
            default: -1.e50;
        }
        
    }
    template<>
    inline double
    BasisFunctionSet<2>::evalDeriv(unsigned int i, unsigned int jDir, const PointReferenceT& p) const
    {
        TestErrorDebug((jDir<2),"Error in BasisFunctionSet.EvalDeriv: invalid derivative direction!");
        
        switch (jDir)
        {
            case 0:
                return vecOfBasisFcn_[i]->evalDeriv0(p);
                break;
            case 1:
                return vecOfBasisFcn_[i]->evalDeriv1(p);
                break;
            default: -1.e50;
        }
        
    }
    template<>
    inline double
    BasisFunctionSet<3>::evalDeriv(unsigned int i, unsigned int jDir, const PointReferenceT& p) const
    {
        TestErrorDebug((jDir<3),"Error in BasisFunctionSet.EvalDeriv: invalid derivative direction!");
        
        switch (jDir)
        {
            case 0:
                return vecOfBasisFcn_[i]->evalDeriv0(p);
                break;
            case 1:
                return vecOfBasisFcn_[i]->evalDeriv1(p);
                break;
            case 2:
                return vecOfBasisFcn_[i]->evalDeriv2(p);
                break;
            default: -1.e50;
        }
        
    }*/
     double
    BasisFunctionSet::evalDeriv(unsigned int i, unsigned int jDir, const PointReferenceT& p) const
    {
        TestErrorDebug((jDir<4),"Error in BasisFunctionSet.EvalDeriv: invalid derivative direction!");
        
        switch (jDir)
        {
            case 0:
                return vecOfBasisFcn_[i]->evalDeriv0(p);
                break;
            case 1:
                return vecOfBasisFcn_[i]->evalDeriv1(p);
                break;
            case 2:
                return vecOfBasisFcn_[i]->evalDeriv2(p);
                break;
            case 3:
                return vecOfBasisFcn_[i]->evalDeriv3(p);;
                break;
            default: return -1.e50;
        }
        
    }   
    
     void
    BasisFunctionSet::eval(unsigned int i, const PointReferenceT& p, NumericalVector& ret) const
    {
        vecOfBasisFcn_[i]->eval(p,ret);
    }

     void
    BasisFunctionSet::evalCurl(unsigned int i, const PointReferenceT& p, NumericalVector& ret) const
    {
        vecOfBasisFcn_[i]->evalCurl(p,ret);
    }
}


//
//  BasisFunctionSet.cpp
//  
//
//  Created by Shavarsh Nurijanyan on 6/24/13.
//
//
namespace Base {

    template<unsigned int DIM>
    class BasisFunctionSet;

    template<unsigned int DIM>
    BasisFunctionSet<DIM>::BasisFunctionSet(unsigned int order):
    order_(order)
    {
    }
            
    template<unsigned int DIM>
    BasisFunctionSet<DIM>::~BasisFunctionSet()
    {
        while(!vecOfBasisFcn_.empty())
        {
            delete vecOfBasisFcn_.back();
            vecOfBasisFcn_.pop_back();
        }
    }

    template<unsigned int DIM>
    unsigned int
    BasisFunctionSet<DIM>::size() const
    {
        return vecOfBasisFcn_.size();
    }

    template<unsigned int DIM>
    unsigned int
    BasisFunctionSet<DIM>::getOrder() const
    {
        return order_;
    }
    template<unsigned int DIM>
    void
    BasisFunctionSet<DIM>::addBasisFunction(BaseBasisFunctionT* bf)
    {
        vecOfBasisFcn_.push_back(bf);
    }

    template<unsigned int DIM>
    double
    BasisFunctionSet<DIM>::eval(unsigned int i, const PointReferenceT& p) const
    {
        return vecOfBasisFcn_[i]->eval(p);
    }
    
    template<unsigned int DIM>
    double
    BasisFunctionSet<DIM>::evalDeriv(unsigned int i, unsigned int jDir, const PointReferenceT& p) const
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
        
    }
    template<>
    inline double
    BasisFunctionSet<4>::evalDeriv(unsigned int i, unsigned int jDir, const PointReferenceT& p) const
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
            default: -1.e50;
        }
        
    }   
    
    template<unsigned int DIM>
    inline void
    BasisFunctionSet<DIM>::eval(unsigned int i, const PointReferenceT& p, NumericalVector& ret) const
    {
        vecOfBasisFcn_[i]->eval(p,ret);
    }

    template<unsigned int DIM>
    inline void
    BasisFunctionSet<DIM>::evalCurl(unsigned int i, const PointReferenceT& p, NumericalVector& ret) const
    {
        vecOfBasisFcn_[i]->evalCurl(p,ret);
    }
}


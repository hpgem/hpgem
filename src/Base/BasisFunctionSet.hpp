#ifndef BasisFunctionSet_hpp
#define BasisFunctionSet_hpp

#include <vector>
#include "Base/BaseBasisFunction.hpp"
#include "Base/TestErrorDebug.hpp"
#include "Geometry/PointReference.hpp"

namespace Base
{
    template <unsigned int DIM>
    class BasisFunctionSet
    {
    public:
        BasisFunctionSet(unsigned int order) : order_(order) {}
        
        ~BasisFunctionSet()
        {
            while(!vecOfBasisFcn_.empty()) 
            {
                delete vecOfBasisFcn_.back();
                vecOfBasisFcn_.pop_back();
            }
        }
        
        unsigned int size() const
        {
            return vecOfBasisFcn_.size();
        }
        
        unsigned int getOrder() const
        {
            return order_;
        }
        
        void AddBasisFunction(BaseBasisFunction<DIM>* bf)
        {
            vecOfBasisFcn_.push_back(bf);
        }
        
        double Eval(unsigned int i, const Geometry::PointReference<DIM>& p) const
        {
            return vecOfBasisFcn_[i]->Eval(p);
        }
        
        double EvalDeriv(unsigned int i, unsigned int jDir, const Geometry::PointReference<DIM>& p) const
        {
            TestErrorDebug((jDir<DIM),"Error in BasisFunctionSet.EvalDeriv: invalid derivative direction!");
          
            switch (jDir)
            {
            case 0:
              return vecOfBasisFcn_[i]->EvalDeriv0(p);
              break;
            case 1:
              return vecOfBasisFcn_[i]->EvalDeriv1(p);
              break;
            case 2:
              return vecOfBasisFcn_[i]->EvalDeriv2(p);
              break;
            case 3:
              return vecOfBasisFcn_[i]->EvalDeriv3(p);
              break;
            }

            return -1.e50;
        }

    private:
        BasisFunctionSet();
        BasisFunctionSet(BasisFunctionSet& other);

    private:
        unsigned int                          order_;
        std::vector<BaseBasisFunction<DIM>*>  vecOfBasisFcn_;
    };

};

#endif

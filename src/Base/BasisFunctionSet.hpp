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
        typedef BaseBasisFunction<DIM>              BaseBasisFunctionT;
        typedef std::vector<BaseBasisFunctionT*>    BaseBasisFunctions;
        typedef Geometry::PointReference<DIM>       PointReferenceT;
        
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
        
        void addBasisFunction(BaseBasisFunctionT* bf)
        {
            vecOfBasisFcn_.push_back(bf);
        }
        
        double eval(unsigned int i, const PointReferenceT& p) const
        {
            return vecOfBasisFcn_[i]->eval(p);
        }
        
        double evalDeriv(unsigned int i, unsigned int jDir, const PointReferenceT& p) const
        {
            TestErrorDebug((jDir<DIM),"Error in BasisFunctionSet.EvalDeriv: invalid derivative direction!");
          
            double (BaseBasisFunction<DIM>::*evalDeriv)(const PointReferenceT&)const;
            
            switch (jDir)
            {
                case 0:
                    evalDeriv=&BaseBasisFunction<DIM>::evalDeriv0;
                  break;
                case 1:
                    evalDeriv=&BaseBasisFunction<DIM>::evalDeriv1;
                    break;
                case 2:
                    evalDeriv=&BaseBasisFunction<DIM>::evalDeriv2;
                  break;
                case 3:
                        //evalDeriv=&BaseBasisFunction<DIM>::evalDeriv3;
                  break;
                default: evalDeriv=NULL;
            }
            return (vecOfBasisFcn_[i]->*evalDeriv)(p);
        }

    private:
        BasisFunctionSet();
        BasisFunctionSet(const BasisFunctionSet& other);

    private:
        unsigned int                          order_;
        BaseBasisFunctions                    vecOfBasisFcn_;
    };

};

#endif

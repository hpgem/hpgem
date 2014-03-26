#ifndef BasisFunctionSet_hpp
#define BasisFunctionSet_hpp

#include <vector>
#include "Base/BaseBasisFunction.hpp"
#include "Base/TestErrorDebug.hpp"
#include "Geometry/PointReference.hpp"

namespace Base
{
    class BasisFunctionSet
    {
    public:
        typedef BaseBasisFunction              BaseBasisFunctionT;
        typedef std::vector<BaseBasisFunctionT*>    BaseBasisFunctions;
        typedef Geometry::PointReference       PointReferenceT;
        
    public:
        BasisFunctionSet(unsigned int order);
        
        virtual ~BasisFunctionSet();
        
         virtual unsigned int size() const;
        
         virtual unsigned int getOrder() const;
        
         virtual void         addBasisFunction(BaseBasisFunctionT* bf);
        
         virtual double       eval(unsigned int i, const PointReferenceT& p) const;
	
	///\brief returns the value of the i-th basisfunction at point p in ret
	 virtual void         eval(unsigned int i, const PointReferenceT& p, NumericalVector& ret) const;
        
         virtual double       evalDeriv(unsigned int i, unsigned int jDir, const PointReferenceT& p) const;
	
	///\brief returns the curl of the i-th basisfunction at point p in ret
	 virtual void         evalCurl(unsigned int i, const PointReferenceT& p, NumericalVector& ret) const;
	
	    virtual const BaseBasisFunction* operator[](int i) const {return vecOfBasisFcn_[i];}

    private:
        BasisFunctionSet();
        BasisFunctionSet(const BasisFunctionSet& other);

    private:
        unsigned int                          order_;
        BaseBasisFunctions                    vecOfBasisFcn_;
    };
}

#endif

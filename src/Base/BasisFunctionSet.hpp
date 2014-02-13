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
        
        ~BasisFunctionSet();
        
         unsigned int size() const;
        
         unsigned int getOrder() const;
        
         void         addBasisFunction(BaseBasisFunctionT* bf);
        
         double       eval(unsigned int i, const PointReferenceT& p) const;
	
	///\brief returns the value of the i-th basisfunction at point p in ret
	 void         eval(unsigned int i, const PointReferenceT& p, NumericalVector& ret) const;
        
         double       evalDeriv(unsigned int i, unsigned int jDir, const PointReferenceT& p) const;
	
	///\brief returns the curl of the i-th basisfunction at point p in ret
	 void         evalCurl(unsigned int i, const PointReferenceT& p, NumericalVector& ret) const;
	

    private:
        BasisFunctionSet();
        BasisFunctionSet(const BasisFunctionSet& other);

    private:
        unsigned int                          order_;
        BaseBasisFunctions                    vecOfBasisFcn_;
    };
}

#endif

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
        BasisFunctionSet(unsigned int order);
        
        ~BasisFunctionSet();
        
        inline unsigned int size() const;
        
        inline unsigned int getOrder() const;
        
        inline void         addBasisFunction(BaseBasisFunctionT* bf);
        
        inline double       eval(unsigned int i, const PointReferenceT& p) const;
        
        inline double       evalDeriv(unsigned int i, unsigned int jDir, const PointReferenceT& p) const;

    private:
        BasisFunctionSet();
        BasisFunctionSet(const BasisFunctionSet& other);

    private:
        unsigned int                          order_;
        BaseBasisFunctions                    vecOfBasisFcn_;
    };
}
#include "BasisFunctionSet_Impl.hpp"

#endif

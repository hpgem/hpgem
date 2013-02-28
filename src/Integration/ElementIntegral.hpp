//------------
// Modified by MT Julianto, Fri Feb 15 17:55:34 WET 2013
//------------------------------------------------------------------------------
#ifndef ElementIntegral_hpp
#define ElementIntegral_hpp
//------------------------------------------------------------------------------
// Package configuration (namespace Integration):
#include "Integration/GlobalNamespaceIntegration.hpp"
//------------------------------------------------------------------------------
// Package includes:
#include "Base/Element.hpp"
#include "Base/TestErrorDebug.hpp"
#include "Geometry/Jacobian.hpp"
#include "Geometry/PointReference.hpp"
#include "Integration/ReturnTrait1.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"
//------------------------------------------------------------------------------

namespace Integration 
{

    template <unsigned int DIM>
    class ElementIntegral
    {
     public:
        typedef typename Base::Element<DIM>::CacheT                 CacheT;
        typedef typename Base::Element<DIM>::VecCacheT              VecCacheT;
        typedef typename QuadratureRules::GaussQuadratureRule<DIM>  QuadratureRulesT;
        typedef typename Base::Element<DIM>                         ElementT;

    public:
        
        //! \brief Construct an ElementIntegral, either with or without cache.
        ElementIntegral(bool useCache=false);
        //! \brief Class destructor
        ~ElementIntegral();

            //! \brief Start caching (geometry) information now.
        void    cacheOn();
        
            //! \brief Stop using cache.
        void    cacheOff();
        
            //! \brief Set recompute the cache ON.
        void    recomputeCacheOn();
        
            //! \brief Set recompute the cache OFF.
        void    recomputeCacheOff();

            //! \brief Directly integrate the inegrand and return ReturnTraits1.
        template <typename IntegrandT>
        void    integrate(ElementT& el, IntegrandT& integrand, typename ReturnTrait1<IntegrandT>::ReturnType& result,
                          const QuadratureRulesT* const qdrRule = NULL);

    
    private:
       
        bool useCache_;
        bool recomputeCache_;
    };

} // close namespace Integration

#include "ElementIntegral_Impl.hpp"

#endif

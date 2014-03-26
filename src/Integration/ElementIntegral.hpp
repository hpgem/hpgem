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
//#include "Integration/ReturnTrait1.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"
#include "ElementIntegrandBase.hpp"
//------------------------------------------------------------------------------

namespace Integration 
{

    class ElementIntegral
    {
     public:
        typedef typename Base::Element::CacheT                 CacheT;
        typedef typename Base::Element::VecCacheT              VecCacheT;
        typedef typename QuadratureRules::GaussQuadratureRule  QuadratureRulesT;
        typedef typename Base::Element                         ElementT;

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

            //! \brief Directly integrate the integrand and return ReturnTraits1.
            //! ReturnTrait1 needs to have the function axpy() implemented
        template<class ReturnTrait1>
        void    integrate(ElementT* el, ElementIntegrandBase<ReturnTrait1>* integrand, ReturnTrait1& result,
                          const QuadratureRulesT* const qdrRule = NULL);
             // \brief Directly integrate the inegrand and return ReturnTraits1, member function version.
        /*template <typename OBJ, typename IntegrandT>
        void    integrate(ElementT* el, IntegrandT& integrand, typename ReturnTrait1<IntegrandT>::ReturnType& result, OBJ* objPtr, const QuadratureRulesT* const qdrRule = NULL);*/
        
        
            /// Probably not needed, this is for classes which are template!
//        template <class OBJ, typename IntegrandT>
//        void    integrate(ElementT* el, IntegrandT& integrand, typename ReturnTrait1<IntegrandT>::ReturnType& result, OBJ* objPtr, const QuadratureRulesT* const qdrRule = NULL);

    
    private:
       
        bool useCache_;
        bool recomputeCache_;
    };

} // close namespace Integration

#include "ElementIntegral_Impl.hpp"

#endif

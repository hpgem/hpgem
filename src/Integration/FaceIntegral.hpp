//------------------------------------------------------------------------------
// File: FaceIntegral.hpp
// Class with members that allow to integrate on faces, 
// with reused cached data.
//------------------------------------------------------------------------------
#ifndef FaceIntegral_hpp
#define FaceIntegral_hpp
//------------------------------------------------------------------------------
// Package configuration (namespace Integration):
#include "GlobalNamespaceIntegration.hpp"
//------------------------------------------------------------------------------
// System includes and names imported from them:
//------------------------------------------------------------------------------
// Package includes:
#include "Base/Face.hpp"
#include "Base/L2Norm.hpp"
#include "Base/TestErrorDebug.hpp"
#include "Geometry/Jacobian.hpp"
#include "Geometry/PointReference.hpp"
#include "Integration/ReturnTrait1.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"

#include "LinearAlgebra/NumericalVector.hpp"

//------------------------------------------------------------------------------

namespace Integration 
{
    template <unsigned int DIM>
    class FaceIntegral
    {
    public:
        typedef typename Base::Face<DIM>::CacheT                        CacheT;
        typedef typename Base::Face<DIM>::VecCacheT                     VecCacheT;
            
        typedef typename QuadratureRules::GaussQuadratureRule<DIM-1>    QuadratureRulesT;
        typedef typename Base::Face<DIM>                                FaceT;

    public:
            //! \brief Construct an FaceIntegral with cache on.
        FaceIntegral(bool useCache = false);
       
            //! \brief Free the memory used for the data storage.
        ~FaceIntegral() ;

        //! \brief Start caching (geometry) information now.
        void cacheOn();

        //! \brief Stop using cache.
        void cacheOff();
        //! \brief Stop using cache.
        void recomputeCacheOn();

        //! \brief Stop using cache.
        void recomputeCacheOff();

        //! \brief Do the face integration using given Gauss integration rule.
        template <typename IntegrandT>
        void integrate(Base::Face<DIM>& fa,
                       IntegrandT& integrand,
                       typename ReturnTrait1<IntegrandT>::ReturnType& result,
                       const QuadratureRules::GaussQuadratureRule<DIM-1>* qdrRule = NULL);
    private:
        
        bool        useCache_;
        bool        recomputeCache_;

    };  // class FaceIntegral
   
} // close namespace Integration

#include "FaceIntegral_Impl.hpp"
#endif

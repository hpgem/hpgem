//
//  FaceIntegral_Impl.hpp
//  
//
//  Created by Shavarsh Nurijanyan on 2/27/13.
//
//

#ifndef _FaceIntegral_Impl_hpp
#define _FaceIntegral_Impl_hpp

namespace Integration
{
    
    template<unsigned int DIM>
    class FaceIntegral;
    
    
        //! \brief Construct an FaceIntegral with cache on.
    template<unsigned int DIM>
    FaceIntegral<DIM>::FaceIntegral(bool useCache):
        useCache_(useCache),
        recomputeCache_(false)
    {}

        //! \brief Free the memory used for the data storage.
    template<unsigned int DIM>
    FaceIntegral<DIM>::~FaceIntegral()
    {
    }

        //! \brief Start caching (geometry) information now.
    template<unsigned int DIM>
    void
    FaceIntegral<DIM>::cacheOn()
    {
        useCache_ = true;
        recomputeCache_ = false;
    }

        //! \brief Stop using cache.
    template<unsigned int DIM>
    void
    FaceIntegral<DIM>::cacheOff()
    {
        useCache_ = true;
        recomputeCache_ = false;
    }

        //! \brief Stop using cache.
    template<unsigned int DIM>
    void
    FaceIntegral<DIM>::recomputeCacheOn()
    {
        recomputeCache_ = true;
    }

        //! \brief Stop using cache.
    template<unsigned int DIM>
    void
    FaceIntegral<DIM>::recomputeCacheOff()
    {
        recomputeCache_ = false;
    }

        //! \brief Do the face integration using given Gauss integration rule.
    template<unsigned int DIM>
    template <typename IntegrandT>
    void
    FaceIntegral<DIM>::integrate(FaceT&                                           fa,
                                 IntegrandT&                                      integrand,
                                 typename ReturnTrait1<IntegrandT>::ReturnType&   result,
                                 const QuadratureRulesT* const                    qdrRule )
    {
        const QuadratureRulesT* const qdrRuleLoc = (qdrRule==NULL? fa.getGaussQuadratureRule(): qdrRule);
        
            // check whether the GaussIntegrationRule is actually for the
            // Element's ReferenceGeometry
        TestErrorDebug((qdrRuleLoc->forReferenceGeometry() == fa.getReferenceGeometry()),
                       "FaceIntegral: " + qdrRuleLoc->getName() + " rule is not for THIS ReferenceGeometry!");
        
            // value returned by the integrand
        typename ReturnTrait1<IntegrandT>::ReturnType value(result);
        
            // number of Gauss quadrature points
        unsigned int nrOfPoints = qdrRuleLoc->nrOfPoints();
        
            // Gauss quadrature point
        Geometry::PointReference<DIM-1> p;
        
        if (!useCache_)
        {
            Geometry::PointPhysical<DIM> Normal;
            
                // first Gauss point
            qdrRuleLoc->getPoint(0, p);
            fa.getNormalVector(p, Normal);
            integrand(fa, Normal, p, result);
            result *= (qdrRuleLoc->weight(0) * Base::L2Norm<DIM>(Normal));
            
                // next Gauss points
            for (unsigned int i = 1; i < nrOfPoints; ++i)
            {
                qdrRuleLoc->getPoint(i, p);
                fa.getNormalVector(p, Normal);
                integrand(fa, Normal, p, value);
                
                 //Y = alpha * X + Y
                result.axpy(qdrRuleLoc->weight(i) * Base::L2Norm<DIM>(Normal),value);
                
            }
        }
        else  // useCache_
        {
                // get vector of cache data
            VecCacheT vecCache = fa.getVecCacheData();
            
                // Calculate the cache
            if ((vecCache.size()!=nrOfPoints) || recomputeCache_)
            {
                std::cout << qdrRuleLoc->getName() << " ";
                std::cout << "FaceIntegral: filling up the cache (" << nrOfPoints << "points)!\n";
                
                vecCache.resize(nrOfPoints);
                for (unsigned int i=0; i<nrOfPoints; ++i)
                {
                    qdrRuleLoc->getPoint(i, p);
                    vecCache[i](fa,p);
                }
            }
            
                // first Gauss point
            qdrRuleLoc->getPoint(0, p);
            integrand(fa, vecCache[0].Normal, p, result);
            result *= (qdrRuleLoc->weight(0) * vecCache[0].L2Normal);
            
                // next Gauss point(s)
            for (unsigned int i = 1; i < nrOfPoints; ++i)
            {
                qdrRuleLoc->getPoint(i, p);
                integrand(fa, vecCache[i].Normal, p, value);
                
                    //Y = alpha * X + Y
                result.axpy(qdrRuleLoc->weight(i) * vecCache[i].L2Normal,value);
                
            } // for integration points
        } // if cached data (else)
    } // function

    
    /*! \brief Integration class for face integrals, a specialization for 1D case
     
     *  Call to the integrand: must pass the face reference, since that allows
     *  access to the elements' data.  We use no cache for 1D case. */
        //-----------------------------------------
        //----------------------------------------- Specialization for 1D case
        //-----------------------------------------
    template <>
    class FaceIntegral<1>
    {
    public:
            //! Construct a FaceIntegral, without geometry cache.
        FaceIntegral(bool) {}
        
        ~FaceIntegral()
        {}
        
        template <class IntegrandT>
        void integrate(Base::Face<1>& fa,
                       IntegrandT& integrand,
                       typename ReturnTrait1<IntegrandT>::ReturnType& result,
                       const QuadratureRules::GaussQuadratureRule<0>* qdrRule = NULL)
        {
            Geometry::PointReference<0> p;
            Geometry::PointPhysical<1> Normal;
            
            fa.getNormalVector(p, Normal);
            integrand(fa, Normal, p, result);
        }
    };
    
}
#endif

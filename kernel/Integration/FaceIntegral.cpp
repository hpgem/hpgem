//
//  FaceIntegral_Impl.hpp
//  
//
//  Created by Shavarsh Nurijanyan on 2/27/13.
//
//

#include "FaceIntegral.hpp"

namespace Integration
{
    
    class FaceIntegral;
    
    
        //! \brief Construct an FaceIntegral with cache on.
    FaceIntegral::FaceIntegral(bool useCache):
        useCache_(useCache),
        recomputeCache_(false)
    {}

        //! \brief Free the memory used for the data storage.
    FaceIntegral::~FaceIntegral()
    {
    }

        //! \brief Start caching (geometry) information now.
    void
    FaceIntegral::cacheOn()
    {
        useCache_ = true;
        recomputeCache_ = false;
    }

        //! \brief Stop using cache.
    void
    FaceIntegral::cacheOff()
    {
        useCache_ = true;
        recomputeCache_ = false;
    }

        //! \brief Stop using cache.
    void
    FaceIntegral::recomputeCacheOn()
    {
        recomputeCache_ = true;
    }

        //! \brief Stop using cache.
    void
    FaceIntegral::recomputeCacheOff()
    {
        recomputeCache_ = false;
    }

        //! \brief Do the face integration using given Gauss integration rule.
    


    
    /*! \brief Integration class for face integrals, a specialization for 1D case
     *  Call to the integrand: must pass the face reference, since that allows
     *  access to the elements' data.  We use no cache for 1D case. */
        //-----------------------------------------
        //----------------------------------------- Specialization for 1D case
        //-----------------------------------------
    /*template <>
    class FaceIntegral<1>
    {
    public:
            //! Construct a FaceIntegral, without geometry cache.
        FaceIntegral(bool) {}
        
        ~FaceIntegral()
        {}
        
        template <class IntegrandT>
        void integrate(Base::Face<1>* fa,
                       IntegrandT& integrand,
                       typename ReturnTrait1<IntegrandT>::ReturnType& result,
                       const QuadratureRules::GaussQuadratureRule<0>* qdrRule = NULL)
        {
            Geometry::PointReference<0> p;
            Geometry::PointPhysical<1> Normal;
            
            fa->getNormalVector(p, Normal);
            integrand(fa, Normal, p, result);
        }

        template <class IntegrandT,class OBJ>
		void integrate(Base::Face<1>* fa,
					   IntegrandT& integrand,
					   typename ReturnTrait1<IntegrandT>::ReturnType& result,OBJ* objPtr,
					   const QuadratureRules::GaussQuadratureRule<0>* qdrRule = NULL)
		{
			Geometry::PointReference<0> p;
			Geometry::PointPhysical<1> Normal;

			fa->getNormalVector(p, Normal);
			(objPtr->*integrand)(fa, Normal, p, result);
		}
    };*/
    
}

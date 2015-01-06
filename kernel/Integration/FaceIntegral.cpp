/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, University of Twente
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#include "FaceIntegral.hpp"

#include "FaceCacheData.hpp"

namespace Integration
{
    
    class FaceIntegral;
    
    
        //! \brief Construct an FaceIntegral with cache on.
    FaceIntegral::FaceIntegral(bool useCache):
        useCache_(useCache),
    localFace_(nullptr)
    {}

        //! \brief Free the memory used for the data storage.
    FaceIntegral::~FaceIntegral()
    {
    	delete localFace_;
    }

        //! \brief Start caching (geometry) information now.
    void
    FaceIntegral::cacheOn()
    {
        useCache_ = true;
        if (localFace_ != nullptr)
        {
        	localFace_->cacheOn();
        }
    }

        //! \brief Stop using cache.
    void
    FaceIntegral::cacheOff()
    {
        useCache_ = false;
        if (localFace_ != nullptr)
        {
        	localFace_->cacheOff();
        }
    }

        //! \brief Force a recomputation of the cache, the next time it is needed
    void
    FaceIntegral::recomputeCacheOn()
    {
        if (localFace_ != nullptr)
        {
        	localFace_->recomputeCacheOn();
        }
    }

        //! \brief Stop forcing a recomputation of the cache
    void
    FaceIntegral::recomputeCacheOff()
    {
        if (localFace_ != nullptr)
        {
        	localFace_->recomputeCacheOff();
        }
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
                       const QuadratureRules::GaussQuadratureRule<0>* qdrRule = nullptr)
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
                       const QuadratureRules::GaussQuadratureRule<0>* qdrRule = nullptr)
		{
			Geometry::PointReference<0> p;
			Geometry::PointPhysical<1> Normal;

			fa->getNormalVector(p, Normal);
			(objPtr->*integrand)(fa, Normal, p, result);
		}
    };*/
    

}

void Integration::FaceIntegral::setStorageWrapper(Base::ShortTermStorageFaceBase* transform) {
	delete localFace_;
	localFace_=transform;
	if(useCache_){
		localFace_->cacheOn();
	}else{
		localFace_->cacheOff();
	}
}

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
#ifndef FaceIntegral_hpp
#define FaceIntegral_hpp

#include <functional>
//------------------------------------------------------------------------------
// Package configuration (namespace Integration):
//------------------------------------------------------------------------------
// System includes and names imported from them:
//------------------------------------------------------------------------------
// Package includes:

//------------------------------------------------------------------------------

namespace LinearAlgebra {
    class NumericalVector;
}

namespace Base {
	class Face;
	class ShortTermStorageFaceBase;
}

namespace Geometry {
    class PointReference;
}

namespace QuadratureRules {
	class GaussQuadratureRule;
}

namespace Integration 
{
	template<class ReturnType1>
	class FaceIntegrandBase;

    class FaceIntegral
    {
    public:
        //typedef typename Base::Face::CacheT                        CacheT;
        //typedef typename Base::Face::VecCacheT                     VecCacheT;
            
        using QuadratureRulesT = QuadratureRules::GaussQuadratureRule;
        using FaceT = Base::Face;

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

        ///\brief provide a Face wrapper that can be used to store transformed function data
        ///this wrapper is responsible for transforming the functions to the reference coordinates
        ///the default is suitable for 2D H1 conforming bases (no transformation for values and multiply with Jac^-T for derivatives)
        ///this class will take over responsibility for the data management
        void setStorageWrapper(Base::ShortTermStorageFaceBase *transform);

        //! \brief Do the face integration using given Gauss integration rule.
        template <class ReturnTrait1>
        ReturnTrait1 integrate(FaceT* fa, FaceIntegrandBase<ReturnTrait1>* integrand, const QuadratureRulesT* qdrRule = nullptr);
        
        //! \brief Nice version accepting an appropriate std::function
        template <class ReturnTrait1>
        ReturnTrait1 integrate(FaceT* fa, std::function<ReturnTrait1(const Base::Face*, const LinearAlgebra::NumericalVector&, const Geometry::PointReference&)> integrandFunc, const QuadratureRulesT* qdrRule = nullptr);
        
        /*template <typename OBJ, typename IntegrandT>
        void    integrate(FaceT* el, IntegrandT& integrand, typename ReturnTrait1<IntegrandT>::ReturnType& result, OBJ* objPtr, const QuadratureRulesT* const qdrRule = nullptr);*/
	
    private:
        
        bool        useCache_;

        Base::ShortTermStorageFaceBase *localFace_;

    };  // class FaceIntegral
   
} // close namespace Integration

#include "FaceIntegral_Impl.hpp"
#endif

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
#ifndef ElementIntegral_hpp
#define ElementIntegral_hpp
//------------------------------------------------------------------------------
// Package configuration (namespace Integration):
//------------------------------------------------------------------------------
// Package includes:
//------------------------------------------------------------------------------

namespace Base{
	class Element;
	class ShortTermStorageElementBase;
}

namespace QuadratureRules {
	class GaussQuadratureRule;
}

namespace Integration 
{
	template<class returntrait1>
	class ElementIntegrandBase;

    class ElementIntegral
    {
     public:
        //typedef typename Base::Element::CacheT                 CacheT;
        //typedef typename Base::Element::VecCacheT              VecCacheT;
        typedef typename QuadratureRules::GaussQuadratureRule  QuadratureRulesT;
        typedef typename Base::Element                         ElementT;

    public:
        
        //! \brief Construct an ElementIntegral, either with or without cache.
        ElementIntegral(bool useCache=false);
        //! \brief Class destructor
        ~ElementIntegral();

            //! \brief Start caching (geometry) information now.
        void    cacheOn();
        
            //! \brief Stop using cache. This routine is not required to delete any stored data.
        void    cacheOff();
        
            //! \brief Set recompute the cache ON.
        void    recomputeCacheOn();
        
            //! \brief Set recompute the cache OFF.
        void    recomputeCacheOff();

        ///\brief provide an Element wrapper that can be used to store transformed function data
        ///this wrapper is responsible for transforming the functions to the reference coordinates
        ///the default is suitable for 2D H1 conforming bases (no transformation for values and multiply with Jac^-T for derivatives)
        ///this class will take over responsibility for the data management
        void setStorageWrapper(Base::ShortTermStorageElementBase *transform);

            //! \brief Directly integrate the integrand and return ReturnTraits1.
            //! ReturnTrait1 needs to have the function axpy() implemented
        template<class ReturnTrait1>
        void    integrate(ElementT* el, ElementIntegrandBase<ReturnTrait1>* integrand, ReturnTrait1& result,
                          const QuadratureRulesT* const qdrRule = nullptr);
             // \brief Directly integrate the inegrand and return ReturnTraits1, member function version.
        /*template <typename OBJ, typename IntegrandT>
        void    integrate(ElementT* el, IntegrandT& integrand, typename ReturnTrait1<IntegrandT>::ReturnType& result, OBJ* objPtr, const QuadratureRulesT* const qdrRule = NULL);*/
        
        
            /// Probably not needed, this is for classes which are template!
//        template <class OBJ, typename IntegrandT>
//        void    integrate(ElementT* el, IntegrandT& integrand, typename ReturnTrait1<IntegrandT>::ReturnType& result, OBJ* objPtr, const QuadratureRulesT* const qdrRule = NULL);

    
    private:
       
        bool useCache_;

        Base::ShortTermStorageElementBase* localElement_;
    };

} // close namespace Integration

#include "ElementIntegral_Impl.hpp"

#endif

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
#ifndef ElementIntegral_h
#define ElementIntegral_h
//------------------------------------------------------------------------------
// Package configuration (namespace Integration):
//------------------------------------------------------------------------------
// Package includes:
#include <functional>
#include "Base/SerializationInclude.h"
//------------------------------------------------------------------------------

namespace Base
{
    class Element;
}

namespace QuadratureRules
{
    class GaussQuadratureRule;
}

namespace Geometry
{
    template<std::size_t DIM>
    class PointReference;
}

namespace LinearAlgebra
{
    class MiddleSizeVector;
}

namespace Integration
{
    template<class returntrait1, std::size_t DIM>
    class ElementIntegrandBase;
    
    template<std::size_t DIM>
    class ElementIntegral
    {
    public:
        
        //! \brief Construct an ElementIntegral, either with or without cache.
        ElementIntegral(bool useCache = false);

        //! \brief Class destructor
        ~ ElementIntegral();
        
        ElementIntegral(const ElementIntegral &other) = delete;

        //! \brief Start caching (geometry) information now.
        void cacheOn();

        //! \brief Stop using cache. This routine is not required to delete any stored data.
        void cacheOff();

        //! \brief Set recompute the cache ON.
        void recomputeCacheOn();

        //! \brief Set recompute the cache OFF.
        void recomputeCacheOff();

        /// you can set the coordinate transformation that is to be used here, before calling integrate
        void setTransformation(std::shared_ptr<Base::CoordinateTransformation<DIM> > transform);

        Base::CoordinateTransformation<DIM>& getTransformation();

        Base::PhysicalElement<DIM>& getPhysicalElement();

        //! \brief Directly integrate the integrand and return ReturnTrait1.
        //! ReturnTrait1 needs to have the function LinearAlgebra::axpy() implemented
        template<class ReturnTrait1>
        ReturnTrait1 integrate(const Base::Element* el, ElementIntegrandBase<ReturnTrait1, DIM>* integrand, QuadratureRules::GaussQuadratureRule * qdrRule = nullptr);

        template<class ReturnType>
        ReturnType integrate(const Base::Element* el, std::function<ReturnType(Base::PhysicalElement<DIM>&)> integrand, QuadratureRules::GaussQuadratureRule * qdrRule = nullptr);

        /// \brief Compute the integral on a reference element. IntegrandType needs to have the function LinearAlgebra::axpy() implemented.
        template<typename IntegrandType>
        IntegrandType referenceElementIntegral(const QuadratureRules::GaussQuadratureRule *ptrQdrRule, std::function<IntegrandType()> integrandFunction);

        template <class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & useCache_;
            ar & element_;
        }

    private:
        
        bool useCache_;

        Base::PhysicalElement<DIM> element_;
    };

} // close namespace Integration

#include "ElementIntegral_Impl.h"

#endif

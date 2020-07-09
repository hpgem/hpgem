/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2014, University of Twente
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors
 may be used to endorse or promote products derived from this software without
 specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef HPGEM_KERNEL_FACEINTEGRAL_H
#define HPGEM_KERNEL_FACEINTEGRAL_H

#include <functional>

#include "Base/PhysicalElement.h"
#include "Base/Face.h"

namespace LinearAlgebra {
class MiddleSizeVector;
}

namespace Base {
class Face;
}

namespace Geometry {
template <std::size_t DIM>
class PointReference;
}

namespace QuadratureRules {
class GaussQuadratureRule;
}

namespace Integration {
template <class ReturnType1, std::size_t DIM>
class FaceIntegrandBase;

template <std::size_t DIM>
class FaceIntegral {
   public:
    //! \brief Construct an FaceIntegral with cache on.
    FaceIntegral();

    ///! \brief Copy constructor should not be used and is therefore deleted.
    FaceIntegral(const FaceIntegral& other) = delete;

    //! \brief Free the memory used for the data storage.
    ~FaceIntegral();

    /// you can set the coordinate transformation that is to be used here,
    /// before calling integrate
    void setTransformation(
        std::shared_ptr<Base::CoordinateTransformation<DIM> > transform,
        std::size_t unknown = 0);

    Base::CoordinateTransformation<DIM>& getTransformation(
        std::size_t unknown = 0);

    //! \brief Do the face integration using given Gauss integration rule.
    template <typename ReturnTrait1>
    ReturnTrait1 integrate(
        const Base::Face* fa, FaceIntegrandBase<ReturnTrait1, DIM>* integrand,
        QuadratureRules::GaussQuadratureRule* qdrRule = nullptr);

    //! \brief Nice version accepting an appropriate std::function
    template <typename FunctionType>
    std::result_of_t<FunctionType(Base::PhysicalFace<DIM>&)> integrate(
        const Base::Face* fa, FunctionType integrandFunc,
        QuadratureRules::GaussQuadratureRule* qdrRule = nullptr);

    //! \brief Nice version accepting an appropriate std::function
    template <typename FunctionType>
    std::result_of_t<FunctionType(Base::PhysicalFace<DIM>&)> integratePair(
        const Base::Face* fa, FunctionType integrandFunc,
        QuadratureRules::GaussQuadratureRule* qdrRule = nullptr);

    /// \brief Compute the integral on a reference element. IntegrandType needs
    /// to have the function LinearAlgebra::axpy() implemented.
    // need to know information about the face to perform the integration
    // template<typename IntegrandType>
    // IntegrandType referenceFaceIntegral(const
    // QuadratureRules::GaussQuadratureRule *ptrQdrRule,
    // std::function<IntegrandType()> integrandFunction) const;

   private:
    Base::PhysicalFace<DIM> internalFace_;
    Base::PhysicalFace<DIM> boundaryFace_;

};  // class FaceIntegral

}  // namespace Integration

#include "FaceIntegral_Impl.h"
#endif // HPGEM_KERNEL_FACEINTEGRAL_H

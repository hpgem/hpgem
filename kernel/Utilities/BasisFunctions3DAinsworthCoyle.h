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

#ifndef HPGEM_KERNEL_BASISFUNCTIONS3DAINSWORTHCOYLE_H
#define HPGEM_KERNEL_BASISFUNCTIONS3DAINSWORTHCOYLE_H

#include "Base/BaseBasisFunction.h"
#include <vector>

namespace hpgem {

namespace Base {
class BasisFunctionSet;
}

namespace Geometry {
template <std::size_t DIM>
class PointReference;
}

namespace Utilities {
//! Curl conforming edge functions.
class BasisCurlEdgeAinsworthCoyle : public Base::BaseBasisFunction {
   public:
    BasisCurlEdgeAinsworthCoyle(std::size_t degree,
                                std::size_t localFirstVertex,
                                std::size_t localSecondVertex);

    void eval(const Geometry::PointReference<3>& p,
              LinearAlgebra::SmallVector<3>& ret) const override;

    LinearAlgebra::SmallVector<3> evalCurl(
        const Geometry::PointReference<3>& p) const override;

   private:
    const std::size_t deg, o, i;
};

//! Curl conforming edge based face functions.
class BasisCurlEdgeFaceAinsworthCoyle : public Base::BaseBasisFunction {
   public:
    BasisCurlEdgeFaceAinsworthCoyle(std::size_t degree,
                                    std::size_t localOpposingVertex,
                                    std::size_t localSpecialVertex);

    void eval(const Geometry::PointReference<3>& p,
              LinearAlgebra::SmallVector<3>& ret) const override;

    LinearAlgebra::SmallVector<3> evalCurl(
        const Geometry::PointReference<3>& p) const override;

   private:
    std::size_t deg, a, b, c;
};

//! Curl conforming face functions.
class BasisCurlFaceAinsworthCoyle : public Base::BaseBasisFunction {
   public:
    BasisCurlFaceAinsworthCoyle(std::size_t degree1, std::size_t degree2,
                                std::size_t localOpposingVertex,
                                std::size_t direction);

    void eval(const Geometry::PointReference<3>& p,
              LinearAlgebra::SmallVector<3>& ret) const override;

    LinearAlgebra::SmallVector<3> evalCurl(
        const Geometry::PointReference<3>& p) const override;

   private:
    std::size_t deg1, deg2, a, b, c;
};

//! Curl conforming face based interior functions.
class BasisCurlFaceinteriorAinsworthCoyle : public Base::BaseBasisFunction {
   public:
    BasisCurlFaceinteriorAinsworthCoyle(std::size_t degree1,
                                        std::size_t degree2,
                                        std::size_t localOpposingVertex);

    void eval(const Geometry::PointReference<3>& p,
              LinearAlgebra::SmallVector<3>& ret) const override;

    LinearAlgebra::SmallVector<3> evalCurl(
        const Geometry::PointReference<3>& p) const override;

   private:
    std::size_t deg1, deg2, a, b, c, d;
};

//! curl conforming interior functions
class BasisCurlinteriorAinsworthCoyle : public Base::BaseBasisFunction {
   public:
    BasisCurlinteriorAinsworthCoyle(std::size_t degree1, std::size_t degree2,
                                    std::size_t degree3, std::size_t direction);

    void eval(const Geometry::PointReference<3>& p,
              LinearAlgebra::SmallVector<3>& ret) const override;

    LinearAlgebra::SmallVector<3> evalCurl(
        const Geometry::PointReference<3>& p) const override;

   private:
    const std::size_t deg1, deg2, deg3, direction;
};

Base::BasisFunctionSet* createDGBasisFunctionSet3DAinsworthCoyle(
    std::size_t order);

}  // namespace Utilities

}  // namespace hpgem

#endif  // HPGEM_KERNEL_BASISFUNCTIONS3DAINSWORTHCOYLE_H

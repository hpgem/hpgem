/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2019, University of Twente
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

#include "BasisFunctionsMonomials.h"
#include "BaseBasisFunction.h"
#include "BasisFunctionSet.h"

#include <array>

namespace hpgem {

namespace Utilities {
double ipow(double value, std::size_t power) {
    double result = 1.0;
    for (std::size_t i = 0; i < power; ++i) result *= value;
    return result;
}

template <std::size_t DIM>
class BasisFunctionsMonomials {
   public:
    explicit BasisFunctionsMonomials(std::array<std::size_t, DIM> powers)
        : powers_(powers) {}

    double evalDIM(const Geometry::PointReference<DIM>& p) const {
        double result = 1;
        for (std::size_t i = 0; i < DIM; ++i) {
            result *= ipow(p[i], powers_[i]);
        }
        return result;
    }

    double evalDerivDIM(const Geometry::PointReference<DIM>& p,
                        std::size_t derivCoord) const {
        double result = powers_[derivCoord];
        for (std::size_t i = 0; i < DIM; ++i) {
            if (i != derivCoord) {
                result *= ipow(p[i], powers_[i]);
            } else if (powers_[i] != 0) {
                // The prefactor is already at the start, multiply by the rest
                result *= ipow(p[i], powers_[i] - 1);
            } else {
                result = 0;
            }
        }
        return result;
    }

   private:
    std::array<std::size_t, DIM> powers_;
};

class BasisFunctionsMonomials1 : public Base::BaseBasisFunction,
                                 private BasisFunctionsMonomials<1> {
   public:
    explicit BasisFunctionsMonomials1(std::size_t power)
        : BasisFunctionsMonomials<1>({power}) {}

    double eval(const Geometry::PointReference<1>& p) const final {
        return evalDIM(p);
    }
    double evalDeriv0(const Geometry::PointReference<1>& p) const final {
        return evalDerivDIM(p, 0);
    }
};

class BasisFunctionsMonomials2 : public Base::BaseBasisFunction,
                                 private BasisFunctionsMonomials<2> {
   public:
    BasisFunctionsMonomials2(std::size_t power1, std::size_t power2)
        : BasisFunctionsMonomials<2>({power1, power2}) {}

    double eval(const Geometry::PointReference<2>& p) const final {
        return evalDIM(p);
    }

    double evalDeriv0(const Geometry::PointReference<2>& p) const final {
        return evalDerivDIM(p, 0);
    }

    double evalDeriv1(const Geometry::PointReference<2>& p) const final {
        return evalDerivDIM(p, 1);
    }
};

class BasisFunctionsMonomials3 : public Base::BaseBasisFunction,
                                 private BasisFunctionsMonomials<3> {
   public:
    BasisFunctionsMonomials3(std::size_t p1, std::size_t p2, std::size_t p3)
        : BasisFunctionsMonomials<3>({p1, p2, p3}) {}

    double eval(const Geometry::PointReference<3>& p) const final {
        return evalDIM(p);
    }

    double evalDeriv0(const Geometry::PointReference<3>& p) const final {
        return evalDerivDIM(p, 0);
    }

    double evalDeriv1(const Geometry::PointReference<3>& p) const final {
        return evalDerivDIM(p, 1);
    }

    double evalDeriv2(const Geometry::PointReference<3>& p) const final {
        return evalDerivDIM(p, 2);
    }
};

class BasisFunctionsMonomials4 : public Base::BaseBasisFunction,
                                 private BasisFunctionsMonomials<4> {
   public:
    BasisFunctionsMonomials4(std::size_t p1, std::size_t p2, std::size_t p3,
                             std::size_t p4)
        : BasisFunctionsMonomials<4>({p1, p2, p3, p4}) {}

    double eval(const Geometry::PointReference<4>& p) const final {
        return evalDIM(p);
    }

    double evalDeriv0(const Geometry::PointReference<4>& p) const final {
        return evalDerivDIM(p, 0);
    }

    double evalDeriv1(const Geometry::PointReference<4>& p) const final {
        return evalDerivDIM(p, 1);
    }

    double evalDeriv2(const Geometry::PointReference<4>& p) const final {
        return evalDerivDIM(p, 2);
    }

    double evalDeriv3(const Geometry::PointReference<4>& p) const final {
        return evalDerivDIM(p, 3);
    }
};

// Actual construction methods
void assembleMonomialBasisFunctions1D(Base::BasisFunctionSet& set,
                                      std::size_t maxPower) {
    for (std::size_t i = 0; i <= maxPower; ++i) {
        set.addBasisFunction(new BasisFunctionsMonomials1(i));
    }
}
void assembleMonomialBasisFunctions2D(Base::BasisFunctionSet& set,
                                      std::size_t maxPower) {
    for (std::size_t i = 0; i <= maxPower; ++i) {
        for (std::size_t j = 0; j + i <= maxPower; ++j) {
            set.addBasisFunction(new BasisFunctionsMonomials2(i, j));
        }
    }
}
void assembleMonomialBasisFunctions3D(Base::BasisFunctionSet& set,
                                      std::size_t maxPower) {
    for (std::size_t i = 0; i <= maxPower; ++i) {
        for (std::size_t j = 0; i + j <= maxPower; ++j) {
            for (std::size_t k = 0; i + j + k <= maxPower; ++k) {
                set.addBasisFunction(new BasisFunctionsMonomials3(i, j, k));
            }
        }
    }
}
void assembleMonomialBasisFunctions4D(Base::BasisFunctionSet& set,
                                      std::size_t maxPower) {
    for (std::size_t i = 0; i <= maxPower; ++i) {
        for (std::size_t j = 0; i + j <= maxPower; ++j) {
            for (std::size_t k = 0; i + j + k <= maxPower; ++k) {
                for (std::size_t l = 0; i + j + k + l <= maxPower; ++l)
                    set.addBasisFunction(
                        new BasisFunctionsMonomials4(i, j, k, l));
            }
        }
    }
}
}  // namespace Utilities

}  // namespace hpgem

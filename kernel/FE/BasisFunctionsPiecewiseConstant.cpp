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

#include "BasisFunctionsPiecewiseConstant.h"

#include <memory>
#include "BaseBasisFunction.h"
#include "BasisFunctionSet.h"

namespace hpgem {

namespace FE {
class PiecewiseConstant1 : public BaseBasisFunction {
   public:
    double eval(const Geometry::PointReference<1>& p) const final {
        return 1.0;
    }

    double evalDeriv0(const Geometry::PointReference<1>& p) const final {
        return 0.0;
    }
};

void addPiecewiseConstantBasisFunction1D(BasisFunctionSet& set) {
    set.addBasisFunction(createPiecewiseConstant1D());
}

BaseBasisFunction* createPiecewiseConstant1D() {
    return new PiecewiseConstant1();
}

class PiecewiseConstant2 : public BaseBasisFunction {
   public:
    double eval(const Geometry::PointReference<2>& p) const final {
        return 1.0;
    }

    double evalDeriv0(const Geometry::PointReference<2>& p) const final {
        return 0.0;
    }

    double evalDeriv1(const Geometry::PointReference<2>& p) const final {
        return 0.0;
    }
};

BaseBasisFunction* createPiecewiseConstant2D() {
    return new PiecewiseConstant2();
}

void addPiecewiseConstantBasisFunction2D(BasisFunctionSet& set) {
    set.addBasisFunction(createPiecewiseConstant2D());
}

class PiecewiseConstant3 : public BaseBasisFunction {
   public:
    double eval(const Geometry::PointReference<3>& p) const final {
        return 1.0;
    }

    double evalDeriv0(const Geometry::PointReference<3>& p) const final {
        return 0.0;
    }

    double evalDeriv1(const Geometry::PointReference<3>& p) const final {
        return 0.0;
    }

    double evalDeriv2(const Geometry::PointReference<3>& p) const final {
        return 0.0;
    }
};

void addPiecewiseConstantBasisFunction3D(BasisFunctionSet& set) {
    set.addBasisFunction(new PiecewiseConstant3());
}

class PiecewiseConstant4 : public BaseBasisFunction {
   public:
    double eval(const Geometry::PointReference<4>& p) const final {
        return 1.0;
    }

    double evalDeriv0(const Geometry::PointReference<4>& p) const final {
        return 0.0;
    }

    double evalDeriv1(const Geometry::PointReference<4>& p) const final {
        return 0.0;
    }

    double evalDeriv2(const Geometry::PointReference<4>& p) const final {
        return 0.0;
    }

    double evalDeriv3(const Geometry::PointReference<4>& p) const final {
        return 0.0;
    }
};

void addPiecewiseConstantBasisFunction4D(BasisFunctionSet& set) {
    set.addBasisFunction(new PiecewiseConstant4());
}
}  // namespace FE
}  // namespace hpgem

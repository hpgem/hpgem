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

#include "BasisFunctions2DNedelec.h"
#include "helperFunctions.h"
#include "Base/BasisFunctionSet.h"
#include "Geometry/ReferenceTriangle.h"
#include "Geometry/PointReference.h"

#include "BasisFunctions1DH1ConformingLine.h"
#include "BasisFunctions2DH1ConformingTriangle.h"

namespace hpgem {

namespace Utilities {
namespace {
LinearAlgebra::SmallVector<2> baricentricDeriv(std::size_t node) {
    LinearAlgebra::SmallVector<2> ret;
    if (node == 0) {
        ret[0] = -1;
        ret[1] = -1;
    } else {
        // clear the return vector so we don't return trash
        ret[0] = 0;
        ret[1] = 0;
        ret[node - 1] = 1;
    }
    return ret;
}
}  // namespace

BasisCurlEdgeNedelec2D::BasisCurlEdgeNedelec2D(
    std::size_t side, std::shared_ptr<const Base::BaseBasisFunction> modulator)
    : side_(side), modulator_(std::move(modulator)) {
    logger.assert_debug(side <= 2, "Triangle has only three sides");
}

void BasisCurlEdgeNedelec2D::eval(const Geometry::PointReference<2>& p,
                                  LinearAlgebra::SmallVector<2>& ret) const {

    evalSideFunction(p, ret);
    Geometry::PointReference<1> coord;
    coord[0] = evalModulatorCoord(p);
    ret *= modulator_->eval(coord);
}

LinearAlgebra::SmallVector<2> BasisCurlEdgeNedelec2D::evalCurl(
    const Geometry::PointReference<2>& p) const {

    double result;

    Geometry::PointReference<1> modulatorCoord;
    modulatorCoord[0] = evalModulatorCoord(p);
    double modulatorVal = modulator_->eval(modulatorCoord);
    LinearAlgebra::SmallVector<2> sideVal;
    evalSideFunction(p, sideVal);

    // Compute rotation Fy/dx - Fx/dy from the side function
    result = -2.0;
    // By chain rule, multiply by modulator (q) value
    result *= modulatorVal;
    // Now second contribution, applying the derivatives to the modulator.
    // Note that the modulator q(t) depends on only one coordinate, where
    // t depends on side_.
    if (side_ == 0) {
        // t=x, so add dq/dt Fy
        result += modulator_->evalDeriv0(modulatorCoord) * sideVal[1];
    } else if (side_ == 1) {
        // t=y, so add -dq/dt Fx
        result -= modulator_->evalDeriv0(modulatorCoord) * sideVal[0];
    } else {
        // t = 1-y, so add +dq/dt Fx
        result += modulator_->evalDeriv0(modulatorCoord) * sideVal[0];
    }
    return LinearAlgebra::SmallVector<2>({result, 0});
}

void BasisCurlEdgeNedelec2D::evalSideFunction(
    const Geometry::PointReference<2>& p,
    LinearAlgebra::SmallVector<2>& ret) const {

    if (side_ == 0) {
        ret[0] = p[1] - 1.0;
        ret[1] = -p[0];
    } else if (side_ == 1) {
        ret[0] = p[1];
        ret[1] = -p[0];
    } else {
        ret[0] = p[1];
        ret[1] = 1.0 - p[0];
    }
}

double BasisCurlEdgeNedelec2D::evalModulatorCoord(
    const Geometry::PointReference<2>& p) const {
    if (side_ == 0) {
        return p[0];
    } else if (side_ == 1) {
        return p[1];
    } else {
        return 1 - p[1];
    }
}

BasisFunctionCurlInteriorNedelec2D::BasisFunctionCurlInteriorNedelec2D(
    std::size_t type, std::shared_ptr<const Base::BaseBasisFunction> modulator)
    : type_(type), modulator_(std::move(modulator)) {
    logger.assert_always(type == 0 || type == 1,
                         "Only two interior types supported");
}

void BasisFunctionCurlInteriorNedelec2D::eval(
    const Geometry::PointReference<2>& p,
    LinearAlgebra::SmallVector<2>& ret) const {

    // Type 0: x (y; 1-x) q(x,y), Type 1: y(1-y;x) q(x,y)
    // where q is the modulating function

    evalTypePart(p, ret);
    double pval = modulator_->eval(p);
    ret *= pval;
}

LinearAlgebra::SmallVector<2> BasisFunctionCurlInteriorNedelec2D::evalCurl(
    const Geometry::PointReference<2>& p) const {

    // Curl = rotation, thus F_y/dx - F_x/dy evaluated via the chain rule
    double result;

    LinearAlgebra::SmallVector<2> typeValue;
    double typeFydx, typeFxdy;
    double modulatorValue;
    LinearAlgebra::SmallVector<2> modulatorGrad;

    evalTypePart(p, typeValue);
    if (type_ == 0) {
        typeFydx = 1 - 2.0 * p[0];
        typeFxdy = p[0];
    } else {
        typeFydx = p[1];
        typeFxdy = 1 - 2.0 * p[1];
    }
    modulatorValue = modulator_->eval(p);
    modulatorGrad = modulator_->evalDeriv(p);

    // Chain rule part 1, F_y/dx
    result = typeFydx * modulatorValue + typeValue[1] * modulatorGrad[0];
    // Chain rule part 2, -F_x/dy
    result -= typeFxdy * modulatorValue + typeValue[0] * modulatorGrad[1];

    return LinearAlgebra::SmallVector<2>({result, 0.0});
}

void BasisFunctionCurlInteriorNedelec2D::evalTypePart(
    const Geometry::PointReference<2>& p,
    LinearAlgebra::SmallVector<2>& ret) const {
    if (type_ == 0) {
        ret[0] = p[1];
        ret[1] = 1.0 - p[0];
        ret *= p[0];
    } else {
        ret[0] = 1.0 - p[1];
        ret[1] = p[0];
        ret *= p[1];
    }
}

std::vector<Base::BaseBasisFunction*> createDGBasisFunctions2DNedelec(
    std::size_t order) {
    //    logger.assert_debug(
    //        order < 2, "2D Nedelec basis functions only implemented for p =
    //        1");

    std::vector<Base::BaseBasisFunction*> result;

    std::vector<Base::BaseBasisFunction*> lineFunctions =
        createDGBasisFunctions1DH1Line(order - 1);
    for (auto& i : lineFunctions) {
        auto lineFunction = std::shared_ptr<Base::BaseBasisFunction>(i);
        result.emplace_back(new BasisCurlEdgeNedelec2D(0, lineFunction));
        result.emplace_back(new BasisCurlEdgeNedelec2D(1, lineFunction));
        result.emplace_back(new BasisCurlEdgeNedelec2D(2, lineFunction));
    }

    if (order > 1) {
        std::vector<Base::BaseBasisFunction*> triangleFunctions =
            createDGBasisFunctions2DH1Triangle(order - 2);
        for (auto& i : triangleFunctions) {
            auto triangleFunction = std::shared_ptr<Base::BaseBasisFunction>(i);
            result.emplace_back(
                new BasisFunctionCurlInteriorNedelec2D(0, triangleFunction));
            result.emplace_back(
                new BasisFunctionCurlInteriorNedelec2D(1, triangleFunction));
        }
    }
    return result;
}
Base::BasisFunctionSet* createDGBasisFunctionSet2DNedelec(std::size_t order) {
    Base::BasisFunctionSet* set = new Base::BasisFunctionSet(order);
    for (auto* bf : createDGBasisFunctions2DNedelec(order)) {
        set->addBasisFunction(bf);
    }
    return set;
}
}  // namespace Utilities

}  // namespace hpgem

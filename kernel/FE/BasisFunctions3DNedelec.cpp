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

#include "BasisFunctions3DNedelec.h"
#include "helperFunctions.h"
#include "BasisFunctionSet.h"
#include "Geometry/ReferenceTetrahedron.h"
#include "Geometry/PointReference.h"

namespace hpgem {

namespace Utilities {
namespace {
void OuterProduct(const LinearAlgebra::SmallVector<3>& a,
                  const LinearAlgebra::SmallVector<3>& b,
                  LinearAlgebra::SmallVector<3>& ret) {
    // ret.resize(3);
    // direct computation using the definition
    ret[0] = a[1] * b[2] - a[2] * b[1];
    ret[1] = a[2] * b[0] - a[0] * b[2];
    ret[2] = a[0] * b[1] - a[1] * b[0];
}

LinearAlgebra::SmallVector<3> baricentricDeriv(std::size_t node) {
    LinearAlgebra::SmallVector<3> ret;
    // ret.resize(3);
    if (node == 0) {
        ret[0] = -1;
        ret[1] = -1;
        ret[2] = -1;
    } else {
        // clear the return vector so we don't return trash
        ret[0] = 0;
        ret[1] = 0;
        ret[2] = 0;
        ret[node - 1] = 1;
    }
    return ret;
}
}  // namespace

BasisCurlEdgeNedelec::BasisCurlEdgeNedelec(std::size_t degree1,
                                           std::size_t degree2,
                                           std::size_t localFirstVertex,
                                           std::size_t localSecondVertex)
    : deg1(degree1), deg2(degree2), i(localFirstVertex), j(localSecondVertex) {
    logger.assert_debug(i < 4 && j < 4, "A tetrahedron only has 4 nodes");
}

void BasisCurlEdgeNedelec::eval(const Geometry::PointReference<3>& p,
                                LinearAlgebra::SmallVector<3>& ret) const {
    LinearAlgebra::SmallVector<3> dummy;

    ret = baricentricDeriv(i);
    dummy = baricentricDeriv(j);

    double valI(baricentric_3D(i, p)), valJ(baricentric_3D(j, p));

    ret *= valJ;
    dummy *= valI;
    ret -= dummy;

    ret *= std::pow(valI, deg1) * std::pow(valJ, deg2);
}

LinearAlgebra::SmallVector<3> BasisCurlEdgeNedelec::evalCurl(
    const Geometry::PointReference<3>& p) const {
    LinearAlgebra::SmallVector<3> dummy, dummy2, ret;

    dummy = baricentricDeriv(i);
    dummy2 = baricentricDeriv(j);
    OuterProduct(dummy2, dummy, ret);

    double valI(baricentric_3D(i, p)), valJ(baricentric_3D(j, p));

    ret *=
        double(deg1 + deg2 + 2) * std::pow(valI, deg1) * std::pow(valJ, deg2);
    return ret;
}

BasisCurlFace1Nedelec::BasisCurlFace1Nedelec(std::size_t degree1,
                                             std::size_t degree2,
                                             std::size_t degree3,
                                             std::size_t localOpposingVertex)
    : deg1(degree1), deg2(degree2), deg3(degree3), d(localOpposingVertex) {
    logger.assert_debug(d < 4, "A tetrahedron only has 4 nodes");
    a = 0;
    if (d == a) {
        a++;
    }
    b = a + 1;
    if (d == b) {
        b++;
    }
    c = b + 1;
    if (d == c) {
        c++;
    }
}

void BasisCurlFace1Nedelec::eval(const Geometry::PointReference<3>& p,
                                 LinearAlgebra::SmallVector<3>& ret) const {
    LinearAlgebra::SmallVector<3> dummy;

    double valI(baricentric_3D(a, p)), valJ(baricentric_3D(b, p)),
        valK(baricentric_3D(c, p));

    ret = baricentricDeriv(a);
    dummy = baricentricDeriv(b);

    ret *= valJ;
    dummy *= valI;
    ret -= dummy;

    ret *=
        std::pow(valI, deg1) * std::pow(valJ, deg2) * std::pow(valK, deg3 + 1);
}

LinearAlgebra::SmallVector<3> BasisCurlFace1Nedelec::evalCurl(
    const Geometry::PointReference<3>& p) const {
    LinearAlgebra::SmallVector<3> dummy, dummy2, dummy3, ret;

    double valI(baricentric_3D(a, p)), valJ(baricentric_3D(b, p)),
        valK(baricentric_3D(c, p));
    dummy = baricentricDeriv(a);
    dummy2 = baricentricDeriv(b);
    dummy3 = baricentricDeriv(c);

    OuterProduct(dummy2, dummy, ret);
    ret *= double(deg1 + deg2 + 2) * std::pow(valI, deg1) *
           std::pow(valJ, deg2) * std::pow(valK, deg3 + 1);

    dummy *= valJ;
    dummy2 *= valI;
    dummy -= dummy2;

    OuterProduct(dummy3, dummy, dummy2);
    dummy2 *= double(1 + deg3) * std::pow(valI, deg1) * std::pow(valJ, deg2) *
              std::pow(valK, deg3);

    ret += dummy2;
    return ret;
}

BasisCurlFace2Nedelec::BasisCurlFace2Nedelec(std::size_t degree1,
                                             std::size_t degree2,
                                             std::size_t degree3,
                                             std::size_t localOpposingVertex)
    : deg1(degree1), deg2(degree2), deg3(degree3), d(localOpposingVertex) {
    logger.assert_debug(d < 4, "A tetrahedron only has 4 nodes");
    a = 0;
    if (d == a) {
        a++;
    }
    b = a + 1;
    if (d == b) {
        b++;
    }
    c = b + 1;
    if (d == c) {
        c++;
    }
}

void BasisCurlFace2Nedelec::eval(const Geometry::PointReference<3>& p,
                                 LinearAlgebra::SmallVector<3>& ret) const {
    LinearAlgebra::SmallVector<3> dummy;

    double valI(baricentric_3D(a, p)), valJ(baricentric_3D(b, p)),
        valK(baricentric_3D(c, p));

    ret = baricentricDeriv(b);
    dummy = baricentricDeriv(c);

    ret *= valK;
    dummy *= valJ;
    ret -= dummy;

    ret *=
        std::pow(valI, deg1 + 1) * std::pow(valJ, deg2) * std::pow(valK, deg3);
}

LinearAlgebra::SmallVector<3> BasisCurlFace2Nedelec::evalCurl(
    const Geometry::PointReference<3>& p) const {
    LinearAlgebra::SmallVector<3> dummy, dummy2, ret, dummy3;

    double valI(baricentric_3D(a, p)), valJ(baricentric_3D(b, p)),
        valK(baricentric_3D(c, p));
    dummy3 = baricentricDeriv(a);
    dummy = baricentricDeriv(b);
    dummy2 = baricentricDeriv(c);

    OuterProduct(dummy2, dummy, ret);
    ret *= double(deg2 + deg3 + 2) * std::pow(valI, deg1 + 1) *
           std::pow(valJ, deg2) * std::pow(valK, deg3);

    dummy *= valK;
    dummy2 *= valJ;
    dummy -= dummy2;

    OuterProduct(dummy3, dummy, dummy2);
    dummy2 *= double(1 + deg1) * std::pow(valI, deg1) * std::pow(valJ, deg2) *
              std::pow(valK, deg3);

    ret += dummy2;
    return ret;
}

BasisCurlinterior1Nedelec::BasisCurlinterior1Nedelec(std::size_t degree1,
                                                     std::size_t degree2,
                                                     std::size_t degree3,
                                                     std::size_t degree4)
    : deg1(degree1), deg2(degree2), deg3(degree3), deg4(degree4) {}

void BasisCurlinterior1Nedelec::eval(const Geometry::PointReference<3>& p,
                                     LinearAlgebra::SmallVector<3>& ret) const {
    LinearAlgebra::SmallVector<3> dummy;

    double val0(baricentric_3D(0, p)), val1(baricentric_3D(1, p)),
        val2(baricentric_3D(2, p)), val3(baricentric_3D(3, p));

    ret = baricentricDeriv(0);
    dummy = baricentricDeriv(1);

    ret *= val1;
    dummy *= val0;
    ret -= dummy;

    ret *= std::pow(val0, deg1) * std::pow(val1, deg2) *
           std::pow(val2, deg3 + 1) * std::pow(val3, deg4 + 1);
}

LinearAlgebra::SmallVector<3> BasisCurlinterior1Nedelec::evalCurl(
    const Geometry::PointReference<3>& p) const {
    LinearAlgebra::SmallVector<3> dummy, dummy2, dummy3, dummy4, ret;

    double val0(baricentric_3D(0, p)), val1(baricentric_3D(1, p)),
        val2(baricentric_3D(2, p)), val3(baricentric_3D(3, p));

    dummy = baricentricDeriv(0);
    dummy2 = baricentricDeriv(1);
    dummy3 = baricentricDeriv(2);
    dummy4 = baricentricDeriv(3);

    OuterProduct(dummy2, dummy, ret);

    dummy *= val1;
    dummy2 *= val0;
    dummy -= dummy2;

    OuterProduct(dummy3, dummy, dummy2);
    OuterProduct(dummy4, dummy, dummy3);

    dummy2 *= double(1 + deg3) * std::pow(val0, deg1) * std::pow(val1, deg2) *
              std::pow(val2, deg3) * std::pow(val3, deg4 + 1);
    dummy3 *= double(1 + deg4) * std::pow(val0, deg1) * std::pow(val1, deg2) *
              std::pow(val2, deg3 + 1) * std::pow(val3, deg4);

    ret *= double(deg1 + deg2 + 2) * std::pow(val0, deg1) *
           std::pow(val1, deg2) * std::pow(val2, deg3 + 1) *
           std::pow(val3, deg4 + 1);
    ret += dummy2 + dummy3;
    return ret;
}

BasisCurlinterior2Nedelec::BasisCurlinterior2Nedelec(std::size_t degree1,
                                                     std::size_t degree2,
                                                     std::size_t degree3,
                                                     std::size_t degree4)
    : deg1(degree1), deg2(degree2), deg3(degree3), deg4(degree4) {}

void BasisCurlinterior2Nedelec::eval(const Geometry::PointReference<3>& p,
                                     LinearAlgebra::SmallVector<3>& ret) const {
    LinearAlgebra::SmallVector<3> dummy;

    double val0(baricentric_3D(0, p)), val1(baricentric_3D(1, p)),
        val2(baricentric_3D(2, p)), val3(baricentric_3D(3, p));

    ret = baricentricDeriv(1);
    dummy = baricentricDeriv(2);

    ret *= val2;
    dummy *= val1;
    ret -= dummy;

    ret *= std::pow(val0, deg1 + 1) * std::pow(val1, deg2) *
           std::pow(val2, deg3) * std::pow(val3, deg4 + 1);
}

LinearAlgebra::SmallVector<3> BasisCurlinterior2Nedelec::evalCurl(
    const Geometry::PointReference<3>& p) const {
    LinearAlgebra::SmallVector<3> dummy, dummy2, dummy3, dummy4, ret;

    double val0(baricentric_3D(0, p)), val1(baricentric_3D(1, p)),
        val2(baricentric_3D(2, p)), val3(baricentric_3D(3, p));

    dummy = baricentricDeriv(1);
    dummy2 = baricentricDeriv(2);
    dummy3 = baricentricDeriv(0);
    dummy4 = baricentricDeriv(3);

    OuterProduct(dummy2, dummy, ret);

    dummy *= val2;
    dummy2 *= val1;
    dummy -= dummy2;

    OuterProduct(dummy3, dummy, dummy2);
    OuterProduct(dummy4, dummy, dummy3);

    dummy2 *= double(1 + deg1) * std::pow(val0, deg1) * std::pow(val1, deg2) *
              std::pow(val2, deg3) * std::pow(val3, deg4 + 1);
    dummy3 *= double(1 + deg4) * std::pow(val0, deg1 + 1) *
              std::pow(val1, deg2) * std::pow(val2, deg3) *
              std::pow(val3, deg4);

    ret *= double(deg2 + deg3 + 2) * std::pow(val0, deg1 + 1) *
           std::pow(val1, deg2) * std::pow(val2, deg3) *
           std::pow(val3, deg4 + 1);
    ret += dummy2 + dummy3;
    return ret;
}

BasisCurlinterior3Nedelec::BasisCurlinterior3Nedelec(std::size_t degree1,
                                                     std::size_t degree2,
                                                     std::size_t degree3,
                                                     std::size_t degree4)
    : deg1(degree1), deg2(degree2), deg3(degree3), deg4(degree4) {}

void BasisCurlinterior3Nedelec::eval(const Geometry::PointReference<3>& p,
                                     LinearAlgebra::SmallVector<3>& ret) const {
    LinearAlgebra::SmallVector<3> dummy;

    double val0(baricentric_3D(0, p)), val1(baricentric_3D(1, p)),
        val2(baricentric_3D(2, p)), val3(baricentric_3D(3, p));

    ret = baricentricDeriv(2);
    dummy = baricentricDeriv(3);

    ret *= val3;
    dummy *= val2;
    ret -= dummy;

    ret *= std::pow(val0, deg1 + 1) * std::pow(val1, deg2 + 1) *
           std::pow(val2, deg3) * std::pow(val3, deg4);
}

LinearAlgebra::SmallVector<3> BasisCurlinterior3Nedelec::evalCurl(
    const Geometry::PointReference<3>& p) const {
    LinearAlgebra::SmallVector<3> dummy, dummy2, dummy3, dummy4, ret;

    double val0(baricentric_3D(0, p)), val1(baricentric_3D(1, p)),
        val2(baricentric_3D(2, p)), val3(baricentric_3D(3, p));

    dummy = baricentricDeriv(2);
    dummy2 = baricentricDeriv(3);
    dummy3 = baricentricDeriv(0);
    dummy4 = baricentricDeriv(1);

    OuterProduct(dummy2, dummy, ret);

    dummy *= val3;
    dummy2 *= val2;
    dummy -= dummy2;

    OuterProduct(dummy3, dummy, dummy2);
    OuterProduct(dummy4, dummy, dummy3);

    dummy2 *= double(1 + deg1) * std::pow(val0, deg1) *
              std::pow(val1, deg2 + 1) * std::pow(val2, deg3) *
              std::pow(val3, deg4);
    dummy3 *= double(1 + deg2) * std::pow(val0, deg1 + 1) *
              std::pow(val1, deg2) * std::pow(val2, deg3) *
              std::pow(val3, deg4);

    ret *= double(deg3 + deg4 + 2) * std::pow(val0, deg1 + 1) *
           std::pow(val1, deg2 + 1) * std::pow(val2, deg3) *
           std::pow(val3, deg4);
    ret += dummy2 + dummy3;
    return ret;
}

Base::BasisFunctionSet* createDGBasisFunctionSet3DNedelec(std::size_t order) {
    Base::BasisFunctionSet* bFset = new Base::BasisFunctionSet(order);
    for (std::size_t l = 0; l + 1 <= order; ++l) {
        std::size_t m((order - 1) - l);
        bFset->addBasisFunction(new BasisCurlEdgeNedelec(l, m, 0, 1));
        bFset->addBasisFunction(new BasisCurlEdgeNedelec(l, m, 0, 2));
        bFset->addBasisFunction(new BasisCurlEdgeNedelec(l, m, 0, 3));
        bFset->addBasisFunction(new BasisCurlEdgeNedelec(l, m, 2, 3));
        bFset->addBasisFunction(new BasisCurlEdgeNedelec(l, m, 1, 3));
        bFset->addBasisFunction(new BasisCurlEdgeNedelec(l, m, 1, 2));
        //	    std::cout<<"constructed edge functions with
        // p="<<p<<std::endl;
    }
    if (order > 1) {
        for (std::size_t l = 0; l + 2 <= order; ++l) {
            for (std::size_t m = 0; m + l + 2 <= order; ++m) {
                std::size_t n((order - 2) - (l + m));
                bFset->addBasisFunction(new BasisCurlFace1Nedelec(l, m, n, 0));
                bFset->addBasisFunction(new BasisCurlFace1Nedelec(l, m, n, 1));
                bFset->addBasisFunction(new BasisCurlFace1Nedelec(l, m, n, 2));
                bFset->addBasisFunction(new BasisCurlFace1Nedelec(l, m, n, 3));
                bFset->addBasisFunction(new BasisCurlFace2Nedelec(l, m, n, 0));
                bFset->addBasisFunction(new BasisCurlFace2Nedelec(l, m, n, 1));
                bFset->addBasisFunction(new BasisCurlFace2Nedelec(l, m, n, 2));
                bFset->addBasisFunction(new BasisCurlFace2Nedelec(l, m, n, 3));
                //	    std::cout<<"constructed face functions and face
                // based interior functions with l="<<l<<" and
                // m="<<m<<std::endl;
            }
        }
    }
    if (order > 2) {
        for (std::size_t l = 0; l + 3 <= order; ++l) {
            for (std::size_t m = 0; m + l + 3 <= order; ++m) {
                for (std::size_t n = 0; n + m + l + 3 <= order; ++n) {
                    std::size_t o((order - 3) - (l + m + n));
                    bFset->addBasisFunction(
                        new BasisCurlinterior1Nedelec(l, m, n, o));
                    bFset->addBasisFunction(
                        new BasisCurlinterior2Nedelec(l, m, n, o));
                    bFset->addBasisFunction(
                        new BasisCurlinterior3Nedelec(l, m, n, o));
                    //	    std::cout<<"constructed interior functions with
                    // l="<<l<<" and m="<<m<<" and n="<<n<<std::endl;
                }
            }
        }
    }
    return bFset;
}
}  // namespace Utilities

}  // namespace hpgem

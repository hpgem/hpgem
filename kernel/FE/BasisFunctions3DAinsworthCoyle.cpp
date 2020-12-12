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

#include "BasisFunctions3DAinsworthCoyle.h"
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

BasisCurlEdgeAinsworthCoyle::BasisCurlEdgeAinsworthCoyle(
    std::size_t degree, std::size_t localFirstVertex,
    std::size_t localSecondVertex)
    : deg(degree), o(localFirstVertex), i(localSecondVertex) {
    logger.assert_debug(i < 4 && o < 4, "A tetrahedron only has 4 nodes");
}

void BasisCurlEdgeAinsworthCoyle::eval(
    const Geometry::PointReference<3>& p,
    LinearAlgebra::SmallVector<3>& ret) const {
    LinearAlgebra::SmallVector<3> dummy;  // dummy vectors are used to store
                                          // parial results
    ret = baricentricDeriv(o);
    dummy = baricentricDeriv(i);

    double valI(baricentric_3D(i, p)), valO(baricentric_3D(o, p));

    ret *= valI;  // use only *= and += and so on near numerical vectors or you
                  // will be rapidly creating
    dummy *=
        valO;  // and destroying a gazilion of them, which is a waste of time
    // switch between the special cases in the definition
    switch (deg) {
        case 0:
            ret -= dummy;
            break;
        case 1:
            ret += dummy;
            break;
        default:
            ret -= dummy;  // degree 0 value
            dummy *= 2;
            dummy += ret;  // degree 1 value
            ret *= -double(deg - 1) / double(deg) *
                   LegendrePolynomial(deg - 2, valI - valO);
            dummy *= double(2 * deg - 1) / double(deg) *
                     LegendrePolynomial(deg - 1, valI - valO);
            ret += dummy;
    }
}

LinearAlgebra::SmallVector<3> BasisCurlEdgeAinsworthCoyle::evalCurl(
    const Geometry::PointReference<3>& p) const {
    LinearAlgebra::SmallVector<3> ret, dummy, dummy2;
    // ret.resize(3);
    switch (deg) {
        case 0:
            dummy = baricentricDeriv(o);
            dummy2 = baricentricDeriv(i);
            OuterProduct(dummy2, dummy, ret);
            ret *= 2;
            break;
        case 1:
            ret[0] = 0;
            ret[1] = 0;
            ret[2] = 0;
            break;
        default:
            dummy = baricentricDeriv(o);
            dummy2 = baricentricDeriv(i);
            LinearAlgebra::SmallVector<3> dummy3, dummy4(dummy2);
            dummy4 -= dummy;  // dummy4=nambla(labda_i-labda_o)
            double valI(baricentric_3D(i, p)), valO(baricentric_3D(o, p));
            OuterProduct(dummy2, dummy, ret);
            ret *= -double(deg - 1) / double(deg) *
                   LegendrePolynomial(deg - 2, valI - valO) * 2;
            dummy *= valI;
            dummy2 *= valO;
            dummy += dummy2;  // dummy=phi_1
            OuterProduct(dummy4, dummy, dummy3);
            dummy3 *= double(2 * deg - 1) / double(deg) *
                      LegendrePolynomialDerivative(deg - 1, valI - valO);
            ret += dummy3;
            dummy2 *= 2;
            dummy -= dummy2;  // dummy=phi_0
            OuterProduct(dummy4, dummy, dummy3);
            dummy3 *= double(deg - 1) / double(deg) *
                      LegendrePolynomialDerivative(deg - 2, valI - valO);
            ret -= dummy3;
    }
    return ret;
}

BasisCurlEdgeFaceAinsworthCoyle::BasisCurlEdgeFaceAinsworthCoyle(
    std::size_t degree, std::size_t localOpposingVertex,
    std::size_t localSpecialVertex)
    : deg(degree), c(localSpecialVertex) {
    logger.assert_debug(c < 4, "A tetrahedron only has 4 nodes");
    a = 0;
    // find the edge (a,b) this functions is based on
    while (c == a || localOpposingVertex == a) {
        a++;
    }
    b = a + 1;
    while (c == b || localOpposingVertex == b) {
        b++;
    }
}

void BasisCurlEdgeFaceAinsworthCoyle::eval(
    const Geometry::PointReference<3>& p,
    LinearAlgebra::SmallVector<3>& ret) const {
    ret = baricentricDeriv(c);
    double valA(baricentric_3D(a, p)), valB(baricentric_3D(b, p));
    ret *= valA * valB * LegendrePolynomial(deg - 2, valB - valA);
}

LinearAlgebra::SmallVector<3> BasisCurlEdgeFaceAinsworthCoyle::evalCurl(
    const Geometry::PointReference<3>& p) const {
    LinearAlgebra::SmallVector<3> ret, dummy, dummy2, dummy3;
    // ret.resize(3);
    double valA(baricentric_3D(a, p)), valB(baricentric_3D(b, p));
    dummy = baricentricDeriv(a);
    dummy2 = baricentricDeriv(b);
    dummy3 = baricentricDeriv(c);
    OuterProduct(dummy, dummy3, ret);
    OuterProduct(dummy2, dummy3, dummy);
    ret *= valB * LegendrePolynomial(deg - 2, valB - valA) -
           valA * valB * LegendrePolynomialDerivative(deg - 2, valB - valA);
    dummy *= valA * LegendrePolynomial(deg - 2, valB - valA) +
             valA * valB * LegendrePolynomialDerivative(deg - 2, valB - valA);
    ret += dummy;
    return ret;
}

BasisCurlFaceAinsworthCoyle::BasisCurlFaceAinsworthCoyle(
    std::size_t degree1, std::size_t degree2, std::size_t localOpposingVertex,
    std::size_t direction)
    : deg1(degree1), deg2(degree2) {
    a = 0;
    // construct the face this function lives on
    if (localOpposingVertex == a) {
        a++;
    }
    b = a + 1;
    if (localOpposingVertex == b) {
        b++;
    }
    c = b + 1;
    if (localOpposingVertex == c) {
        c++;
    }
    // lets the tangent always point to vertex b in the rest of this class
    if (direction == 2) {      // this also swaps the use of l and m in the
                               // ainsworth definition
        std::size_t temp = c;  // but that should be inconsequential
        c = b;
        b = temp;
    }
}

void BasisCurlFaceAinsworthCoyle::eval(
    const Geometry::PointReference<3>& p,
    LinearAlgebra::SmallVector<3>& ret) const {
    // ret.resize(3);
    double valA(baricentric_3D(a, p)), valB(baricentric_3D(b, p)),
        valC(baricentric_3D(c, p));
    ret[0] = (b == 1 - a == 1);
    ret[1] = (b == 2);
    ret[2] = (b == 3);
    ret *= valA * valB * valC * LegendrePolynomial(deg1, valB - valA) *
           LegendrePolynomial(deg2, valC - valA);
}

LinearAlgebra::SmallVector<3> BasisCurlFaceAinsworthCoyle::evalCurl(
    const Geometry::PointReference<3>& p) const {
    double valA(baricentric_3D(a, p)), valB(baricentric_3D(b, p)),
        valC(baricentric_3D(c, p));
    LinearAlgebra::SmallVector<3> ret, dummy, dummy2, dummy3;
    // ret.resize(3);
    dummy[0] = (b == 1 - a == 1);
    dummy[1] = (b == 2);
    dummy[2] = (b == 3);
    dummy2 = baricentricDeriv(a);
    OuterProduct(dummy2, dummy, ret);
    ret *= valB * valC * LegendrePolynomial(deg1, valB - valA) *
               LegendrePolynomial(deg2, valC - valA) -
           valA * valB * valC *
               (LegendrePolynomialDerivative(deg1, valB - valA) *
                    LegendrePolynomial(deg2, valC - valA) +
                LegendrePolynomial(deg1, valB - valA) *
                    LegendrePolynomialDerivative(deg2, valC - valA));
    dummy2 = baricentricDeriv(b);
    OuterProduct(dummy2, dummy, dummy3);
    dummy3 *= valA * valC * LegendrePolynomial(deg1, valB - valA) *
                  LegendrePolynomial(deg2, valC - valA) +
              valA * valB * valC *
                  LegendrePolynomialDerivative(deg1, valB - valA) *
                  LegendrePolynomial(deg2, valC - valA);
    ret += dummy3;
    dummy2 = baricentricDeriv(c);
    OuterProduct(dummy2, dummy, dummy3);
    dummy3 *= valA * valB * LegendrePolynomial(deg1, valB - valA) *
                  LegendrePolynomial(deg2, valC - valA) +
              valA * valB * valC * LegendrePolynomial(deg1, valB - valA) *
                  LegendrePolynomialDerivative(deg2, valC - valA);
    ret += dummy3;
    return ret;
}

BasisCurlFaceinteriorAinsworthCoyle::BasisCurlFaceinteriorAinsworthCoyle(
    std::size_t degree1, std::size_t degree2, std::size_t localOpposingVertex)
    : deg1(degree1), deg2(degree2), d(localOpposingVertex) {
    a = 0;
    // find the face...
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

void BasisCurlFaceinteriorAinsworthCoyle::eval(
    const Geometry::PointReference<3>& p,
    LinearAlgebra::SmallVector<3>& ret) const {
    double valA(baricentric_3D(a, p)), valB(baricentric_3D(b, p)),
        valC(baricentric_3D(c, p));
    ret = baricentricDeriv(d);
    ret *= valA * valB * valC * LegendrePolynomial(deg1, valB - valA) *
           LegendrePolynomial(deg2, valC - valA);
}

LinearAlgebra::SmallVector<3> BasisCurlFaceinteriorAinsworthCoyle::evalCurl(
    const Geometry::PointReference<3>& p) const {
    double valA(baricentric_3D(a, p)), valB(baricentric_3D(b, p)),
        valC(baricentric_3D(c, p));
    LinearAlgebra::SmallVector<3> ret, dummy, dummy2, dummy3;
    // ret.resize(3);
    dummy = baricentricDeriv(d);
    dummy2 = baricentricDeriv(a);
    OuterProduct(dummy2, dummy, ret);
    ret *= valB * valC * LegendrePolynomial(deg1, valB - valA) *
               LegendrePolynomial(deg2, valC - valA) -
           valA * valB * valC *
               (LegendrePolynomialDerivative(deg1, valB - valA) *
                    LegendrePolynomial(deg2, valC - valA) +
                LegendrePolynomial(deg1, valB - valA) *
                    LegendrePolynomialDerivative(deg2, valC - valA));
    dummy2 = baricentricDeriv(b);
    OuterProduct(dummy2, dummy, dummy3);
    dummy3 *= valA * valC * LegendrePolynomial(deg1, valB - valA) *
                  LegendrePolynomial(deg2, valC - valA) +
              valA * valB * valC *
                  LegendrePolynomialDerivative(deg1, valB - valA) *
                  LegendrePolynomial(deg2, valC - valA);
    ret += dummy3;
    dummy2 = baricentricDeriv(c);
    OuterProduct(dummy2, dummy, dummy3);
    dummy3 *= valA * valB * LegendrePolynomial(deg1, valB - valA) *
                  LegendrePolynomial(deg2, valC - valA) +
              valA * valB * valC * LegendrePolynomial(deg1, valB - valA) *
                  LegendrePolynomialDerivative(deg2, valC - valA);
    ret += dummy3;
    return ret;
}

BasisCurlinteriorAinsworthCoyle::BasisCurlinteriorAinsworthCoyle(
    std::size_t degree1, std::size_t degree2, std::size_t degree3,
    std::size_t direction)
    : deg1(degree1), deg2(degree2), deg3(degree3), direction(direction) {}

void BasisCurlinteriorAinsworthCoyle::eval(
    const Geometry::PointReference<3>& p,
    LinearAlgebra::SmallVector<3>& ret) const {
    // ret.resize(3);
    ret[0] = 0;
    ret[1] = 0;
    ret[2] = 0;
    ret[direction - 1] = 1;
    double val0(baricentric_3D(0, p)), val1(baricentric_3D(1, p)),
        val2(baricentric_3D(2, p)), val3(baricentric_3D(3, p));
    ret *= val0 * val1 * val2 * val3 * LegendrePolynomial(deg1, val1 - val0) *
           LegendrePolynomial(deg2, val2 - val0) *
           LegendrePolynomial(deg3, val3 - val0);
}

LinearAlgebra::SmallVector<3> BasisCurlinteriorAinsworthCoyle::evalCurl(
    const Geometry::PointReference<3>& p) const {
    double val0(baricentric_3D(0, p)), val1(baricentric_3D(1, p)),
        val2(baricentric_3D(2, p)), val3(baricentric_3D(3, p));
    LinearAlgebra::SmallVector<3> ret, dummy, dummy2, dummy3;
    dummy[0] = 0;
    dummy[1] = 0;
    dummy[2] = 0;
    dummy[direction - 1] = 1;
    dummy2 = baricentricDeriv(0);
    // ret.resize(3);
    OuterProduct(dummy2, dummy, ret);
    ret *= val1 * val2 * val3 * LegendrePolynomial(deg1, val1 - val0) *
               LegendrePolynomial(deg2, val2 - val0) *
               LegendrePolynomial(deg3, val3 - val0) -
           val0 * val1 * val2 * val3 *
               (LegendrePolynomialDerivative(deg1, val1 - val0) *
                    LegendrePolynomial(deg2, val2 - val0) *
                    LegendrePolynomial(deg3, val3 - val0) +
                LegendrePolynomial(deg1, val1 - val0) *
                    LegendrePolynomialDerivative(deg2, val2 - val0) *
                    LegendrePolynomial(deg3, val3 - val0) +
                LegendrePolynomial(deg1, val1 - val0) *
                    LegendrePolynomial(deg2, val2 - val0) *
                    LegendrePolynomialDerivative(deg3, val3 - val0));
    dummy2 = baricentricDeriv(1);
    OuterProduct(dummy2, dummy, dummy3);
    dummy3 *= val0 * val2 * val3 * LegendrePolynomial(deg1, val1 - val0) *
                  LegendrePolynomial(deg2, val2 - val0) *
                  LegendrePolynomial(deg3, val3 - val0) +
              val0 * val1 * val2 * val3 *
                  LegendrePolynomialDerivative(deg1, val1 - val0) *
                  LegendrePolynomial(deg2, val2 - val0) *
                  LegendrePolynomial(deg3, val3 - val0);
    ret += dummy3;
    dummy2 = baricentricDeriv(2);
    OuterProduct(dummy2, dummy, dummy3);
    dummy3 *= val0 * val1 * val3 * LegendrePolynomial(deg1, val1 - val0) *
                  LegendrePolynomial(deg2, val2 - val0) *
                  LegendrePolynomial(deg3, val3 - val0) +
              val0 * val1 * val2 * val3 *
                  LegendrePolynomial(deg1, val1 - val0) *
                  LegendrePolynomialDerivative(deg2, val2 - val0) *
                  LegendrePolynomial(deg3, val3 - val0);
    ret += dummy3;
    dummy2 = baricentricDeriv(3);
    OuterProduct(dummy2, dummy, dummy3);
    dummy3 *= val0 * val1 * val2 * LegendrePolynomial(deg1, val1 - val0) *
                  LegendrePolynomial(deg2, val2 - val0) *
                  LegendrePolynomial(deg3, val3 - val0) +
              val0 * val1 * val2 * val3 *
                  LegendrePolynomial(deg1, val1 - val0) *
                  LegendrePolynomial(deg2, val2 - val0) *
                  LegendrePolynomialDerivative(deg3, val3 - val0);
    ret += dummy3;
    return ret;
}

Base::BasisFunctionSet* createDGBasisFunctionSet3DAinsworthCoyle(
    std::size_t order) {
    Base::BasisFunctionSet* bFset = new Base::BasisFunctionSet(order);
    for (std::size_t p = 0; p <= order; ++p) {
        // constructor takes first the degree of the function then the two
        // vertices on the edge
        bFset->addBasisFunction(new BasisCurlEdgeAinsworthCoyle(p, 0, 1));
        bFset->addBasisFunction(new BasisCurlEdgeAinsworthCoyle(p, 0, 2));
        bFset->addBasisFunction(new BasisCurlEdgeAinsworthCoyle(p, 0, 3));
        bFset->addBasisFunction(new BasisCurlEdgeAinsworthCoyle(p, 2, 3));
        bFset->addBasisFunction(new BasisCurlEdgeAinsworthCoyle(p, 1, 3));
        bFset->addBasisFunction(new BasisCurlEdgeAinsworthCoyle(p, 1, 2));
        //	    std::cout<<"constructed edge functions with
        // p="<<p<<std::endl;
    }
    for (std::size_t p = 2; p <= order; ++p) {
        // constructor takes first the degree of the function the the vertex
        // opposing the face, then the vertex not on the edge
        bFset->addBasisFunction(new BasisCurlEdgeFaceAinsworthCoyle(p, 3, 0));
        bFset->addBasisFunction(new BasisCurlEdgeFaceAinsworthCoyle(p, 3, 1));
        bFset->addBasisFunction(new BasisCurlEdgeFaceAinsworthCoyle(p, 3, 2));
        bFset->addBasisFunction(new BasisCurlEdgeFaceAinsworthCoyle(p, 2, 0));
        bFset->addBasisFunction(new BasisCurlEdgeFaceAinsworthCoyle(p, 2, 1));
        bFset->addBasisFunction(new BasisCurlEdgeFaceAinsworthCoyle(p, 2, 3));
        bFset->addBasisFunction(new BasisCurlEdgeFaceAinsworthCoyle(p, 1, 0));
        bFset->addBasisFunction(new BasisCurlEdgeFaceAinsworthCoyle(p, 1, 2));
        bFset->addBasisFunction(new BasisCurlEdgeFaceAinsworthCoyle(p, 1, 3));
        bFset->addBasisFunction(new BasisCurlEdgeFaceAinsworthCoyle(p, 0, 1));
        bFset->addBasisFunction(new BasisCurlEdgeFaceAinsworthCoyle(p, 0, 2));
        bFset->addBasisFunction(new BasisCurlEdgeFaceAinsworthCoyle(p, 0, 3));
        //	    std::cout<<"constructed edge based face functions with
        // p="<<p<<std::endl;
    }
    for (std::size_t l = 0; l + 3 <= order; ++l) {
        for (std::size_t m = 0; m + l + 3 <= order; ++m) {
            // constructor takes first the degrees of the legendre polinomials,
            // then the vertex opposing the face, then a 1 or 2 to fix the
            // direction of the bubble function
            bFset->addBasisFunction(
                new BasisCurlFaceAinsworthCoyle(l, m, 3, 1));
            bFset->addBasisFunction(
                new BasisCurlFaceAinsworthCoyle(l, m, 3, 2));
            bFset->addBasisFunction(
                new BasisCurlFaceAinsworthCoyle(l, m, 2, 1));
            bFset->addBasisFunction(
                new BasisCurlFaceAinsworthCoyle(l, m, 2, 2));
            bFset->addBasisFunction(
                new BasisCurlFaceAinsworthCoyle(l, m, 1, 1));
            bFset->addBasisFunction(
                new BasisCurlFaceAinsworthCoyle(l, m, 1, 2));
            bFset->addBasisFunction(
                new BasisCurlFaceAinsworthCoyle(l, m, 0, 1));
            bFset->addBasisFunction(
                new BasisCurlFaceAinsworthCoyle(l, m, 0, 2));
            // constructor takes first the degrees of the legendre polinomials,
            // then the vertex opposing the face
            bFset->addBasisFunction(
                new BasisCurlFaceinteriorAinsworthCoyle(l, m, 3));
            bFset->addBasisFunction(
                new BasisCurlFaceinteriorAinsworthCoyle(l, m, 2));
            bFset->addBasisFunction(
                new BasisCurlFaceinteriorAinsworthCoyle(l, m, 1));
            bFset->addBasisFunction(
                new BasisCurlFaceinteriorAinsworthCoyle(l, m, 0));
            //	    std::cout<<"constructed face functions and face based
            // interior functions with l="<<l<<" and m="<<m<<std::endl;
        }
    }
    for (std::size_t l = 0; l + 4 <= order; ++l) {
        for (std::size_t m = 0; m + l + 4 <= order; ++m) {
            for (std::size_t n = 0; n + m + l + 4 <= order; ++n) {
                // constructor takes first the degrees of the legendre
                // polinomials, then a 1, 2 or 3 to fix the direction of the
                // bubble function
                bFset->addBasisFunction(
                    new BasisCurlinteriorAinsworthCoyle(l, m, n, 1));
                bFset->addBasisFunction(
                    new BasisCurlinteriorAinsworthCoyle(l, m, n, 2));
                bFset->addBasisFunction(
                    new BasisCurlinteriorAinsworthCoyle(l, m, n, 3));
                //	    std::cout<<"constructed interior functions with
                // l="<<l<<" and m="<<m<<" and n="<<n<<std::endl;
            }
        }
    }
    return bFset;
}
}  // namespace Utilities

}  // namespace hpgem

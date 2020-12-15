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

#include "BasisFunctions3DH1ConformingPrism.h"
#include "helperFunctions.h"
#include "BasisFunctionSet.h"
#include "OrientedBasisFunctionSet.h"
#include "Geometry/ReferenceTriangularPrism.h"
#include "Geometry/PointReference.h"

// only uses the constant basis functions
#include "BasisFunctionsPiecewiseConstant.h"

namespace hpgem {

namespace FE {

BasisFunction3DVertexPrism::BasisFunction3DVertexPrism(std::size_t node) {
    logger.assert_debug(node < 6, "A triangular prism only has 6 nodes");
    nodePosition_ = (static_cast<int>(node) / 3) * 2 - 1;
    node_ = node % 3;
}

double BasisFunction3DVertexPrism::eval(
    const Geometry::PointReference<3>& p) const {
    return baricentric_2D(node_, {p[0], p[1]}) * (1 + nodePosition_ * p[2]) /
           2.;
}

double BasisFunction3DVertexPrism::evalDeriv0(
    const Geometry::PointReference<3>& p) const {
    return baricentricDeriv(node_, 0) * (1 + nodePosition_ * p[2]) / 2.;
}

double BasisFunction3DVertexPrism::evalDeriv1(
    const Geometry::PointReference<3>& p) const {
    return baricentricDeriv(node_, 1) * (1 + nodePosition_ * p[2]) / 2.;
}

double BasisFunction3DVertexPrism::evalDeriv2(
    const Geometry::PointReference<3>& p) const {
    return baricentric_2D(node_, {p[0], p[1]}) * nodePosition_ / 2.;
}

double BasisFunction3DEdgePrism_0::eval(
    const Geometry::PointReference<3>& p) const {
    return (1 + edgePosition_ * p[2]) * baricentric_2D(node0_, {p[0], p[1]}) *
           baricentric_2D(node1_, {p[0], p[1]}) *
           LobattoPolynomial(polynomialOrder_,
                             baricentric_2D(node0_, {p[0], p[1]}) -
                                 baricentric_2D(node1_, {p[0], p[1]})) /
           2.;
}

double BasisFunction3DEdgePrism_0::evalDeriv0(
    const Geometry::PointReference<3>& p) const {
    const Geometry::PointReference<2>& p2D = {p[0], p[1]};
    return (1 + edgePosition_ * p[2]) * baricentricDeriv(node0_, 0) *
               (baricentric_2D(node1_, p2D) *
                    LobattoPolynomial(polynomialOrder_,
                                      baricentric_2D(node0_, p2D) -
                                          baricentric_2D(node1_, p2D)) +
                baricentric_2D(node0_, p2D) * baricentric_2D(node1_, p2D) *
                    LobattoPolynomialDerivative(
                        polynomialOrder_, baricentric_2D(node0_, p2D) -
                                              baricentric_2D(node1_, p2D))) /
               2. +
           (1 + edgePosition_ * p[2]) * baricentricDeriv(node1_, 0) *
               (baricentric_2D(node0_, p2D) *
                    LobattoPolynomial(polynomialOrder_,
                                      baricentric_2D(node0_, p2D) -
                                          baricentric_2D(node1_, p2D)) -
                baricentric_2D(node0_, p2D) * baricentric_2D(node1_, p2D) *
                    LobattoPolynomialDerivative(
                        polynomialOrder_, baricentric_2D(node0_, p2D) -
                                              baricentric_2D(node1_, p2D))) /
               2.;
}

double BasisFunction3DEdgePrism_0::evalDeriv1(
    const Geometry::PointReference<3>& p) const {
    const Geometry::PointReference<2>& p2D = {p[0], p[1]};
    return (1 + edgePosition_ * p[2]) * baricentricDeriv(node0_, 1) *
               (baricentric_2D(node1_, p2D) *
                    LobattoPolynomial(polynomialOrder_,
                                      baricentric_2D(node0_, p2D) -
                                          baricentric_2D(node1_, p2D)) +
                baricentric_2D(node0_, p2D) * baricentric_2D(node1_, p2D) *
                    LobattoPolynomialDerivative(
                        polynomialOrder_, baricentric_2D(node0_, p2D) -
                                              baricentric_2D(node1_, p2D))) /
               2. +
           (1 + edgePosition_ * p[2]) * baricentricDeriv(node1_, 1) *
               (baricentric_2D(node0_, p2D) *
                    LobattoPolynomial(polynomialOrder_,
                                      baricentric_2D(node0_, p2D) -
                                          baricentric_2D(node1_, p2D)) -
                baricentric_2D(node0_, p2D) * baricentric_2D(node1_, p2D) *
                    LobattoPolynomialDerivative(
                        polynomialOrder_, baricentric_2D(node0_, p2D) -
                                              baricentric_2D(node1_, p2D))) /
               2.;
}

double BasisFunction3DEdgePrism_0::evalDeriv2(
    const Geometry::PointReference<3>& p) const {
    const Geometry::PointReference<2>& p2D = {p[0], p[1]};
    return edgePosition_ * baricentric_2D(node0_, p2D) *
           baricentric_2D(node1_, p2D) *
           LobattoPolynomial(
               polynomialOrder_,
               baricentric_2D(node0_, p2D) - baricentric_2D(node1_, p2D)) /
           2.;
}

BasisFunction3DEdgePrism_1::BasisFunction3DEdgePrism_1(
    std::size_t node0, std::size_t node1, std::size_t polynomialOrder)
    : node_(node0 % 3), polynomialOrder_(polynomialOrder) {
    logger.assert_debug(node0 < 6, "A triangular prism only has 6 nodes");
    logger.assert_debug(node1 < 6, "A triangular prism only has 6 nodes");
    logger.assert_debug(node0 % 3 == node1 % 3,
                        "This edge is aligned next to a triangular face");
    mirroring_ = node0 < node1 ? -1 : 1;
}

double BasisFunction3DEdgePrism_1::eval(
    const Geometry::PointReference<3>& p) const {
    return baricentric_2D(node_, {p[0], p[1]}) * (1 - p[2]) * (1 + p[2]) *
           LobattoPolynomial(polynomialOrder_, mirroring_ * p[2]) / 4.;
}

double BasisFunction3DEdgePrism_1::evalDeriv0(
    const Geometry::PointReference<3>& p) const {
    return baricentricDeriv(node_, 0) * (1 - p[2]) * (1 + p[2]) *
           LobattoPolynomial(polynomialOrder_, mirroring_ * p[2]) / 4.;
}

double BasisFunction3DEdgePrism_1::evalDeriv1(
    const Geometry::PointReference<3>& p) const {
    return baricentricDeriv(node_, 1) * (1 - p[2]) * (1 + p[2]) *
           LobattoPolynomial(polynomialOrder_, mirroring_ * p[2]) / 4.;
}

double BasisFunction3DEdgePrism_1::evalDeriv2(
    const Geometry::PointReference<3>& p) const {
    return baricentric_2D(node_, {p[0], p[1]}) *
           (-p[2] * LobattoPolynomial(polynomialOrder_, mirroring_ * p[2]) /
                2. +
            (1 - p[2]) * (1 + p[2]) *
                LobattoPolynomialDerivative(polynomialOrder_,
                                            mirroring_ * p[2]) *
                mirroring_ / 4.);
}

BasisFunction3DFacePrism_0::BasisFunction3DFacePrism_0(
    std::size_t node0, std::size_t node1, std::size_t node2,
    std::size_t polynomialOrder0, std::size_t polynomialOrder1)
    : polynomialOrder0_(polynomialOrder0),
      polynomialOrder1_(polynomialOrder1),
      node0_(node0 % 3),
      node1_(node1 % 3),
      node2_(node2 % 3) {
    logger.assert_debug(node0 < 6, "A triangular prism only has 6 nodes");
    logger.assert_debug(node1 < 6, "A triangular prism only has 6 nodes");
    logger.assert_debug(node2 < 6, "A triangular prism only has 6 nodes");
    logger.assert_debug(node0 / 3 == node1 / 3,
                        "This is not a triangular face");
    logger.assert_debug(node0 / 3 == node2 / 3,
                        "This is not a triangular face");
    facePosition_ = (static_cast<int>(node0) / 3) * 2 - 1;
}

double BasisFunction3DFacePrism_0::eval(
    const Geometry::PointReference<3>& p) const {
    const Geometry::PointReference<2>& p2D = {p[0], p[1]};
    double x0(baricentric_2D(node0_, p2D) - baricentric_2D(node1_, p2D)),
        x1(baricentric_2D(node1_, p2D) - baricentric_2D(node2_, p2D));
    return (1 + facePosition_ * p[2]) * baricentric_2D(node0_, p2D) *
           baricentric_2D(node1_, p2D) * baricentric_2D(node2_, p2D) *
           LobattoPolynomial(polynomialOrder0_, x0) *
           LobattoPolynomial(polynomialOrder1_, x1) / 2.;
}

double BasisFunction3DFacePrism_0::evalDeriv0(
    const Geometry::PointReference<3>& p) const {
    const Geometry::PointReference<2>& p2D = {p[0], p[1]};
    double x0(baricentric_2D(node0_, p2D) - baricentric_2D(node1_, p2D)),
        x1(baricentric_2D(node1_, p2D) - baricentric_2D(node2_, p2D));
    return (1 + facePosition_ * p[2]) * baricentricDeriv(node0_, 0) *
               (baricentric_2D(node1_, p2D) * baricentric_2D(node2_, p2D) *
                    LobattoPolynomial(polynomialOrder0_, x0) *
                    LobattoPolynomial(polynomialOrder1_, x1) +
                baricentric_2D(node0_, p2D) * baricentric_2D(node1_, p2D) *
                    baricentric_2D(node2_, p2D) *
                    LobattoPolynomialDerivative(polynomialOrder0_, x0) *
                    LobattoPolynomial(polynomialOrder1_, x1)) /
               2. +
           (1 + facePosition_ * p[2]) * baricentricDeriv(node1_, 0) *
               (baricentric_2D(node0_, p2D) * baricentric_2D(node2_, p2D) *
                    LobattoPolynomial(polynomialOrder0_, x0) *
                    LobattoPolynomial(polynomialOrder1_, x1) -
                baricentric_2D(node0_, p2D) * baricentric_2D(node1_, p2D) *
                    baricentric_2D(node2_, p2D) *
                    LobattoPolynomialDerivative(polynomialOrder0_, x0) *
                    LobattoPolynomial(polynomialOrder1_, x1) +
                baricentric_2D(node0_, p2D) * baricentric_2D(node1_, p2D) *
                    baricentric_2D(node2_, p2D) *
                    LobattoPolynomial(polynomialOrder0_, x0) *
                    LobattoPolynomialDerivative(polynomialOrder1_, x1)) /
               2. +
           (1 + facePosition_ * p[2]) * baricentricDeriv(node2_, 0) *
               (baricentric_2D(node0_, p2D) * baricentric_2D(node1_, p2D) *
                    LobattoPolynomial(polynomialOrder0_, x0) *
                    LobattoPolynomial(polynomialOrder1_, x1) -
                baricentric_2D(node0_, p2D) * baricentric_2D(node1_, p2D) *
                    baricentric_2D(node2_, p2D) *
                    LobattoPolynomial(polynomialOrder0_, x0) *
                    LobattoPolynomialDerivative(polynomialOrder1_, x1)) /
               2.;
}

double BasisFunction3DFacePrism_0::evalDeriv1(
    const Geometry::PointReference<3>& p) const {
    const Geometry::PointReference<2>& p2D = {p[0], p[1]};
    double x0(baricentric_2D(node0_, p2D) - baricentric_2D(node1_, p2D)),
        x1(baricentric_2D(node1_, p2D) - baricentric_2D(node2_, p2D));
    return (1 + facePosition_ * p[2]) * baricentricDeriv(node0_, 1) *
               (baricentric_2D(node1_, p2D) * baricentric_2D(node2_, p2D) *
                    LobattoPolynomial(polynomialOrder0_, x0) *
                    LobattoPolynomial(polynomialOrder1_, x1) +
                baricentric_2D(node0_, p2D) * baricentric_2D(node1_, p2D) *
                    baricentric_2D(node2_, p2D) *
                    LobattoPolynomialDerivative(polynomialOrder0_, x0) *
                    LobattoPolynomial(polynomialOrder1_, x1)) /
               2. +
           (1 + facePosition_ * p[2]) * baricentricDeriv(node1_, 1) *
               (baricentric_2D(node0_, p2D) * baricentric_2D(node2_, p2D) *
                    LobattoPolynomial(polynomialOrder0_, x0) *
                    LobattoPolynomial(polynomialOrder1_, x1) -
                baricentric_2D(node0_, p2D) * baricentric_2D(node1_, p2D) *
                    baricentric_2D(node2_, p2D) *
                    LobattoPolynomialDerivative(polynomialOrder0_, x0) *
                    LobattoPolynomial(polynomialOrder1_, x1) +
                baricentric_2D(node0_, p2D) * baricentric_2D(node1_, p2D) *
                    baricentric_2D(node2_, p2D) *
                    LobattoPolynomial(polynomialOrder0_, x0) *
                    LobattoPolynomialDerivative(polynomialOrder1_, x1)) /
               2. +
           (1 + facePosition_ * p[2]) * baricentricDeriv(node2_, 1) *
               (baricentric_2D(node0_, p2D) * baricentric_2D(node1_, p2D) *
                    LobattoPolynomial(polynomialOrder0_, x0) *
                    LobattoPolynomial(polynomialOrder1_, x1) -
                baricentric_2D(node0_, p2D) * baricentric_2D(node1_, p2D) *
                    baricentric_2D(node2_, p2D) *
                    LobattoPolynomial(polynomialOrder0_, x0) *
                    LobattoPolynomialDerivative(polynomialOrder1_, x1)) /
               2.;
}

double BasisFunction3DFacePrism_0::evalDeriv2(
    const Geometry::PointReference<3>& p) const {
    const Geometry::PointReference<2>& p2D = {p[0], p[1]};
    double x0(baricentric_2D(node0_, p2D) - baricentric_2D(node1_, p2D)),
        x1(baricentric_2D(node1_, p2D) - baricentric_2D(node2_, p2D));
    return facePosition_ * baricentric_2D(node0_, p2D) *
           baricentric_2D(node1_, p2D) * baricentric_2D(node2_, p2D) *
           LobattoPolynomial(polynomialOrder0_, x0) *
           LobattoPolynomial(polynomialOrder1_, x1) / 2.;
}

BasisFunction3DFacePrism_1::BasisFunction3DFacePrism_1(
    std::size_t node0, std::size_t node1, std::size_t node2,
    std::size_t polynomialOrder0, std::size_t polynomialOrder1)
    : node0_(node0 % 3) {
    logger.assert_debug(node0 < 6, "A triangular prism only has 6 nodes");
    logger.assert_debug(node1 < 6, "A triangular prism only has 6 nodes");
    logger.assert_debug(node2 < 6, "A triangular prism only has 6 nodes");
    if (node1 % 3 == node0_) {
        mirroring_ = node0 < node1 ? -1 : 1;
        polynomialOrder0_ = polynomialOrder0;
        polynomialOrder1_ = polynomialOrder1;
        node1_ = node2 % 3;
    } else {
        mirroring_ = node0 < node2 ? -1 : 1;
        polynomialOrder0_ = polynomialOrder1;
        polynomialOrder1_ = polynomialOrder0;
        node1_ = node1 % 3;
    }
}

double BasisFunction3DFacePrism_1::eval(
    const Geometry::PointReference<3>& p) const {
    const Geometry::PointReference<2>& p2D = {p[0], p[1]};
    return (1 + p[2]) * (1 - p[2]) * baricentric_2D(node0_, p2D) *
           baricentric_2D(node1_, p2D) *
           LobattoPolynomial(
               polynomialOrder0_,
               baricentric_2D(node0_, p2D) - baricentric_2D(node1_, p2D)) *
           LobattoPolynomial(polynomialOrder1_, mirroring_ * p[2]) / 4.;
}

double BasisFunction3DFacePrism_1::evalDeriv0(
    const Geometry::PointReference<3>& p) const {
    const Geometry::PointReference<2>& p2D = {p[0], p[1]};
    return (1 + p[2]) * (1 - p[2]) *
               LobattoPolynomial(polynomialOrder1_, mirroring_ * p[2]) *
               baricentricDeriv(node0_, 0) *
               (baricentric_2D(node1_, p2D) *
                    LobattoPolynomial(polynomialOrder0_,
                                      baricentric_2D(node0_, p2D) -
                                          baricentric_2D(node1_, p2D)) +
                baricentric_2D(node0_, p2D) * baricentric_2D(node1_, p2D) *
                    LobattoPolynomialDerivative(
                        polynomialOrder0_, baricentric_2D(node0_, p2D) -
                                               baricentric_2D(node1_, p2D))) /
               4. +
           (1 + p[2]) * (1 - p[2]) *
               LobattoPolynomial(polynomialOrder1_, mirroring_ * p[2]) *
               baricentricDeriv(node1_, 0) *
               (baricentric_2D(node0_, p2D) *
                    LobattoPolynomial(polynomialOrder0_,
                                      baricentric_2D(node0_, p2D) -
                                          baricentric_2D(node1_, p2D)) -
                baricentric_2D(node0_, p2D) * baricentric_2D(node1_, p2D) *
                    LobattoPolynomialDerivative(
                        polynomialOrder0_, baricentric_2D(node0_, p2D) -
                                               baricentric_2D(node1_, p2D))) /
               4.;
}

double BasisFunction3DFacePrism_1::evalDeriv1(
    const Geometry::PointReference<3>& p) const {
    const Geometry::PointReference<2>& p2D = {p[0], p[1]};
    return (1 + p[2]) * (1 - p[2]) *
               LobattoPolynomial(polynomialOrder1_, mirroring_ * p[2]) *
               baricentricDeriv(node0_, 1) *
               (baricentric_2D(node1_, p2D) *
                    LobattoPolynomial(polynomialOrder0_,
                                      baricentric_2D(node0_, p2D) -
                                          baricentric_2D(node1_, p2D)) +
                baricentric_2D(node0_, p2D) * baricentric_2D(node1_, p2D) *
                    LobattoPolynomialDerivative(
                        polynomialOrder0_, baricentric_2D(node0_, p2D) -
                                               baricentric_2D(node1_, p2D))) /
               4. +
           (1 + p[2]) * (1 - p[2]) *
               LobattoPolynomial(polynomialOrder1_, mirroring_ * p[2]) *
               baricentricDeriv(node1_, 1) *
               (baricentric_2D(node0_, p2D) *
                    LobattoPolynomial(polynomialOrder0_,
                                      baricentric_2D(node0_, p2D) -
                                          baricentric_2D(node1_, p2D)) -
                baricentric_2D(node0_, p2D) * baricentric_2D(node1_, p2D) *
                    LobattoPolynomialDerivative(
                        polynomialOrder0_, baricentric_2D(node0_, p2D) -
                                               baricentric_2D(node1_, p2D))) /
               4.;
}

double BasisFunction3DFacePrism_1::evalDeriv2(
    const Geometry::PointReference<3>& p) const {
    const Geometry::PointReference<2>& p2D = {p[0], p[1]};
    return baricentric_2D(node0_, p2D) * baricentric_2D(node1_, p2D) *
           LobattoPolynomial(
               polynomialOrder0_,
               baricentric_2D(node0_, p2D) - baricentric_2D(node1_, p2D)) *
           (-p[2] * LobattoPolynomial(polynomialOrder1_, mirroring_ * p[2]) /
                2. +
            (1 + p[2]) * (1 - p[2]) *
                LobattoPolynomialDerivative(polynomialOrder1_,
                                            mirroring_ * p[2]) *
                mirroring_ / 4.);
}

double BasisFunction3DInteriorPrism::eval(
    const Geometry::PointReference<3>& p) const {
    const Geometry::PointReference<2>& p2D = {p[0], p[1]};
    double x0(baricentric_2D(0, p2D) - baricentric_2D(1, p2D)),
        x1(baricentric_2D(1, p2D) - baricentric_2D(2, p2D));
    return baricentric_2D(0, p2D) * baricentric_2D(1, p2D) *
           baricentric_2D(2, p2D) * (1 - p[2]) * (1 + p[2]) *
           LobattoPolynomial(polnomialOrder0_, x0) *
           LobattoPolynomial(polynomialOrder1_, x1) *
           LobattoPolynomial(polynomialOrder2_, p[2]) / 4.;
}

double BasisFunction3DInteriorPrism::evalDeriv0(
    const Geometry::PointReference<3>& p) const {
    const Geometry::PointReference<2>& p2D = {p[0], p[1]};
    double x0(baricentric_2D(0, p2D) - baricentric_2D(1, p2D)),
        x1(baricentric_2D(1, p2D) - baricentric_2D(2, p2D));
    return (1 - p[2]) * (1 + p[2]) *
           LobattoPolynomial(polynomialOrder2_, p[2]) / 4. *
           (baricentricDeriv(0, 0) *
                (baricentric_2D(1, p2D) * baricentric_2D(2, p2D) *
                     LobattoPolynomial(polnomialOrder0_, x0) *
                     LobattoPolynomial(polynomialOrder1_, x1) +
                 baricentric_2D(0, p2D) * baricentric_2D(1, p2D) *
                     baricentric_2D(2, p2D) *
                     LobattoPolynomialDerivative(polnomialOrder0_, x0) *
                     LobattoPolynomial(polynomialOrder1_, x1)) +
            baricentricDeriv(1, 0) *
                (baricentric_2D(0, p2D) * baricentric_2D(2, p2D) *
                     LobattoPolynomial(polnomialOrder0_, x0) *
                     LobattoPolynomial(polynomialOrder1_, x1) -
                 baricentric_2D(0, p2D) * baricentric_2D(1, p2D) *
                     baricentric_2D(2, p2D) *
                     LobattoPolynomialDerivative(polnomialOrder0_, x0) *
                     LobattoPolynomial(polynomialOrder1_, x1) +
                 baricentric_2D(0, p2D) * baricentric_2D(1, p2D) *
                     baricentric_2D(2, p2D) *
                     LobattoPolynomial(polnomialOrder0_, x0) *
                     LobattoPolynomialDerivative(polynomialOrder1_, x1)));
}

double BasisFunction3DInteriorPrism::evalDeriv1(
    const Geometry::PointReference<3>& p) const {
    const Geometry::PointReference<2>& p2D = {p[0], p[1]};
    double x0(baricentric_2D(0, p2D) - baricentric_2D(1, p2D)),
        x1(baricentric_2D(1, p2D) - baricentric_2D(2, p2D));
    return (1 - p[2]) * (1 + p[2]) *
           LobattoPolynomial(polynomialOrder2_, p[2]) / 4. *
           (baricentricDeriv(0, 1) *
                (baricentric_2D(1, p2D) * baricentric_2D(2, p2D) *
                     LobattoPolynomial(polnomialOrder0_, x0) *
                     LobattoPolynomial(polynomialOrder1_, x1) +
                 baricentric_2D(0, p2D) * baricentric_2D(1, p2D) *
                     baricentric_2D(2, p2D) *
                     LobattoPolynomialDerivative(polnomialOrder0_, x0) *
                     LobattoPolynomial(polynomialOrder1_, x1)) +
            baricentricDeriv(2, 1) *
                (baricentric_2D(0, p2D) * baricentric_2D(1, p2D) *
                     LobattoPolynomial(polnomialOrder0_, x0) *
                     LobattoPolynomial(polynomialOrder1_, x1) -
                 baricentric_2D(0, p2D) * baricentric_2D(1, p2D) *
                     baricentric_2D(2, p2D) *
                     LobattoPolynomial(polnomialOrder0_, x0) *
                     LobattoPolynomialDerivative(polynomialOrder1_, x1)));
}

double BasisFunction3DInteriorPrism::evalDeriv2(
    const Geometry::PointReference<3>& p) const {
    const Geometry::PointReference<2>& p2D = {p[0], p[1]};
    double x0(baricentric_2D(0, p2D) - baricentric_2D(1, p2D)),
        x1(baricentric_2D(1, p2D) - baricentric_2D(2, p2D));
    return baricentric_2D(0, p2D) * baricentric_2D(1, p2D) *
           baricentric_2D(2, p2D) * LobattoPolynomial(polnomialOrder0_, x0) *
           LobattoPolynomial(polynomialOrder1_, x1) *
           (-p[2] * LobattoPolynomial(polynomialOrder2_, p[2]) / 2. +
            (1 - p[2]) * (1 + p[2]) *
                LobattoPolynomialDerivative(polynomialOrder2_, p[2]) / 4.);
}

BasisFunctionSet* createDGBasisFunctionSet3DH1ConformingPrism(
    std::size_t order) {
    BasisFunctionSet* result = new BasisFunctionSet(order);
    if (order > 0) {
        Geometry::ReferenceTriangularPrism& prism =
            Geometry::ReferenceTriangularPrism::Instance();
        std::vector<std::size_t> vectorOfPointIndexes(4);
        for (std::size_t i = 0; i < 6; ++i) {
            result->addBasisFunction(new BasisFunction3DVertexPrism(i));
        }
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            for (std::size_t i = 0; i < 6; ++i) {
                vectorOfPointIndexes = prism.getCodim2EntityLocalIndices(i);
                result->addBasisFunction(new BasisFunction3DEdgePrism_0(
                    vectorOfPointIndexes[0], vectorOfPointIndexes[1], j));
            }
            for (std::size_t i = 6; i < 9; ++i) {
                vectorOfPointIndexes = prism.getCodim2EntityLocalIndices(i);
                result->addBasisFunction(new BasisFunction3DEdgePrism_1(
                    vectorOfPointIndexes[0], vectorOfPointIndexes[1], j));
            }
            for (std::size_t i = 0; i < 2; ++i) {
                vectorOfPointIndexes = prism.getCodim1EntityLocalIndices(i);
                if (j > 0) {
                    for (std::size_t k = 0; k < j; ++k) {
                        result->addBasisFunction(new BasisFunction3DFacePrism_0(
                            vectorOfPointIndexes[0], vectorOfPointIndexes[1],
                            vectorOfPointIndexes[2], j - k - 1, k));
                    }
                }
            }
            for (std::size_t i = 2; i < 5; ++i) {
                vectorOfPointIndexes = prism.getCodim1EntityLocalIndices(i);
                result->addBasisFunction(new BasisFunction3DFacePrism_1(
                    vectorOfPointIndexes[0], vectorOfPointIndexes[1],
                    vectorOfPointIndexes[2], j, j));
                for (std::size_t k = 0; k < j; ++k) {
                    result->addBasisFunction(new BasisFunction3DFacePrism_1(
                        vectorOfPointIndexes[0], vectorOfPointIndexes[1],
                        vectorOfPointIndexes[2], j, k));
                    result->addBasisFunction(new BasisFunction3DFacePrism_1(
                        vectorOfPointIndexes[0], vectorOfPointIndexes[1],
                        vectorOfPointIndexes[2], k, j));
                }
            }
            for (std::size_t i = 0; i < j; ++i) {
                for (std::size_t k = 0; i + k < j; ++k) {
                    result->addBasisFunction(
                        new BasisFunction3DInteriorPrism(i, j, k));
                    result->addBasisFunction(
                        new BasisFunction3DInteriorPrism(j, i, k));
                    result->addBasisFunction(
                        new BasisFunction3DInteriorPrism(i, k, j));
                }
            }
        }
    } else {
        addPiecewiseConstantBasisFunction3D(*result);
    }
    return result;
}

BasisFunctionSet* createInteriorBasisFunctionSet3DH1ConformingPrism(
    std::size_t order) {
    logger.assert_debug(order > 0,
                        "Trying to create a conforming, constant basis "
                        "function set, did you mean the constant solution?");
    BasisFunctionSet* result = new BasisFunctionSet(order);
    for (std::size_t i = 0; i + 3 <= order; ++i) {
        for (std::size_t j = 0; i + j + 3 <= order; ++j) {
            for (std::size_t k = 0; k + 2 <= order; ++k) {
                result->addBasisFunction(
                    new BasisFunction3DInteriorPrism(i, j, k));
            }
        }
    }
    return result;
}

std::vector<const BasisFunctionSet*>
    createVertexBasisFunctionSet3DH1ConformingPrism(std::size_t order) {
    logger.assert_debug(order > 0,
                        "Trying to create a conforming, constant basis "
                        "function set, did you mean the constant solution?");
    std::vector<const BasisFunctionSet*> result;
    BasisFunctionSet* set;
    for (std::size_t i = 0; i < 6; ++i) {
        set = new BasisFunctionSet(order);
        set->addBasisFunction(new BasisFunction3DVertexPrism(i));
        result.push_back(set);
    }
    return result;
}

std::vector<const OrientedBasisFunctionSet*>
    createEdgeBasisFunctionSet3DH1ConformingPrism(std::size_t order) {
    logger.assert_debug(order > 0,
                        "Trying to create a conforming, constant basis "
                        "function set, did you mean the constant solution?");
    std::vector<const OrientedBasisFunctionSet*> result;
    OrientedBasisFunctionSet* set;
    Geometry::ReferenceTriangularPrism& prism =
        Geometry::ReferenceTriangularPrism::Instance();
    std::vector<std::size_t> vectorOfPointIndexes(2);
    for (std::size_t i = 0; i < 6; ++i) {
        vectorOfPointIndexes = prism.getCodim2EntityLocalIndices(i);
        set = new OrientedBasisFunctionSet(order, 0, i);
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DEdgePrism_0(
                vectorOfPointIndexes[0], vectorOfPointIndexes[1], j));
        }
        result.push_back(set);
        set = new OrientedBasisFunctionSet(order, 1, i);
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DEdgePrism_0(
                vectorOfPointIndexes[1], vectorOfPointIndexes[0], j));
        }
        result.push_back(set);
    }
    for (std::size_t i = 6; i < 9; ++i) {
        vectorOfPointIndexes = prism.getCodim2EntityLocalIndices(i);
        set = new OrientedBasisFunctionSet(order, 0, i);
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DEdgePrism_1(
                vectorOfPointIndexes[0], vectorOfPointIndexes[1], j));
        }
        result.push_back(set);
        set = new OrientedBasisFunctionSet(order, 1, i);
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DEdgePrism_1(
                vectorOfPointIndexes[1], vectorOfPointIndexes[0], j));
        }
        result.push_back(set);
    }
    return result;
}

std::vector<const OrientedBasisFunctionSet*>
    createFaceBasisFunctionSet3DH1ConformingPrism(std::size_t order) {
    logger.assert_debug(order > 0,
                        "Trying to create a conforming, constant basis "
                        "function set, did you mean the constant solution?");
    std::vector<const OrientedBasisFunctionSet*> result;
    OrientedBasisFunctionSet* set;
    Geometry::ReferenceTriangularPrism& prism =
        Geometry::ReferenceTriangularPrism::Instance();
    std::vector<std::size_t> vectorOfPointIndexes(4);
    for (std::size_t i = 0; i < 2; ++i) {
        vectorOfPointIndexes = prism.getCodim1EntityLocalIndices(i);
        set = new OrientedBasisFunctionSet(order, 0, i);
        for (std::size_t j = 0; j + 3 <= order; ++j) {
            for (std::size_t k = 0; j + k + 3 <= order; ++k) {
                set->addBasisFunction(new BasisFunction3DFacePrism_0(
                    vectorOfPointIndexes[0], vectorOfPointIndexes[1],
                    vectorOfPointIndexes[2], j, k));
            }
        }
        result.push_back(set);
        set = new OrientedBasisFunctionSet(order, 1, i);
        for (std::size_t j = 0; j + 3 <= order; ++j) {
            for (std::size_t k = 0; j + k + 3 <= order; ++k) {
                set->addBasisFunction(new BasisFunction3DFacePrism_0(
                    vectorOfPointIndexes[0], vectorOfPointIndexes[2],
                    vectorOfPointIndexes[1], j, k));
            }
        }
        result.push_back(set);
        set = new OrientedBasisFunctionSet(order, 2, i);
        for (std::size_t j = 0; j + 3 <= order; ++j) {
            for (std::size_t k = 0; j + k + 3 <= order; ++k) {
                set->addBasisFunction(new BasisFunction3DFacePrism_0(
                    vectorOfPointIndexes[1], vectorOfPointIndexes[2],
                    vectorOfPointIndexes[0], j, k));
            }
        }
        result.push_back(set);
        set = new OrientedBasisFunctionSet(order, 3, i);
        for (std::size_t j = 0; j + 3 <= order; ++j) {
            for (std::size_t k = 0; j + k + 3 <= order; ++k) {
                set->addBasisFunction(new BasisFunction3DFacePrism_0(
                    vectorOfPointIndexes[1], vectorOfPointIndexes[0],
                    vectorOfPointIndexes[2], j, k));
            }
        }
        result.push_back(set);
        set = new OrientedBasisFunctionSet(order, 4, i);
        for (std::size_t j = 0; j + 3 <= order; ++j) {
            for (std::size_t k = 0; j + k + 3 <= order; ++k) {
                set->addBasisFunction(new BasisFunction3DFacePrism_0(
                    vectorOfPointIndexes[2], vectorOfPointIndexes[1],
                    vectorOfPointIndexes[0], j, k));
            }
        }
        result.push_back(set);
        set = new OrientedBasisFunctionSet(order, 5, i);
        for (std::size_t j = 0; j + 3 <= order; ++j) {
            for (std::size_t k = 0; j + k + 3 <= order; ++k) {
                set->addBasisFunction(new BasisFunction3DFacePrism_0(
                    vectorOfPointIndexes[2], vectorOfPointIndexes[0],
                    vectorOfPointIndexes[1], j, k));
            }
        }
        result.push_back(set);
    }
    for (std::size_t i = 2; i < 5; ++i) {
        vectorOfPointIndexes = prism.getCodim1EntityLocalIndices(i);
        set = new OrientedBasisFunctionSet(order, 0, i);
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            for (std::size_t k = 0; k + 2 <= order; ++k) {
                set->addBasisFunction(new BasisFunction3DFacePrism_0(
                    vectorOfPointIndexes[0], vectorOfPointIndexes[1],
                    vectorOfPointIndexes[2], j, k));
            }
        }
        result.push_back(set);
        set = new OrientedBasisFunctionSet(order, 1, i);
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            for (std::size_t k = 0; k + 2 <= order; ++k) {
                set->addBasisFunction(new BasisFunction3DFacePrism_0(
                    vectorOfPointIndexes[1], vectorOfPointIndexes[2],
                    vectorOfPointIndexes[3], j, k));
            }
        }
        result.push_back(set);
        set = new OrientedBasisFunctionSet(order, 2, i);
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            for (std::size_t k = 0; k + 2 <= order; ++k) {
                set->addBasisFunction(new BasisFunction3DFacePrism_0(
                    vectorOfPointIndexes[2], vectorOfPointIndexes[3],
                    vectorOfPointIndexes[0], j, k));
            }
        }
        result.push_back(set);
        set = new OrientedBasisFunctionSet(order, 3, i);
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            for (std::size_t k = 0; k + 2 <= order; ++k) {
                set->addBasisFunction(new BasisFunction3DFacePrism_0(
                    vectorOfPointIndexes[3], vectorOfPointIndexes[0],
                    vectorOfPointIndexes[1], j, k));
            }
        }
        result.push_back(set);
        set = new OrientedBasisFunctionSet(order, 4, i);
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            for (std::size_t k = 0; k + 2 <= order; ++k) {
                set->addBasisFunction(new BasisFunction3DFacePrism_0(
                    vectorOfPointIndexes[3], vectorOfPointIndexes[2],
                    vectorOfPointIndexes[1], j, k));
            }
        }
        result.push_back(set);
        set = new OrientedBasisFunctionSet(order, 5, i);
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            for (std::size_t k = 0; k + 2 <= order; ++k) {
                set->addBasisFunction(new BasisFunction3DFacePrism_0(
                    vectorOfPointIndexes[1], vectorOfPointIndexes[0],
                    vectorOfPointIndexes[3], j, k));
            }
        }
        result.push_back(set);
        set = new OrientedBasisFunctionSet(order, 6, i);
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            for (std::size_t k = 0; k + 2 <= order; ++k) {
                set->addBasisFunction(new BasisFunction3DFacePrism_0(
                    vectorOfPointIndexes[2], vectorOfPointIndexes[1],
                    vectorOfPointIndexes[0], j, k));
            }
        }
        result.push_back(set);
        set = new OrientedBasisFunctionSet(order, 7, i);
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            for (std::size_t k = 0; k + 2 <= order; ++k) {
                set->addBasisFunction(new BasisFunction3DFacePrism_0(
                    vectorOfPointIndexes[0], vectorOfPointIndexes[3],
                    vectorOfPointIndexes[2], j, k));
            }
        }
        result.push_back(set);
    }
    return result;
}

}  // namespace FE

}  // namespace hpgem

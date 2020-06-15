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

#include "BasisFunctions3DH1ConformingCube.h"
#include "helperFunctions.h"
#include "Base/BasisFunctionSet.h"
#include "Base/OrientedBasisFunctionSet.h"
#include "Geometry/ReferenceCube.h"
#include "Geometry/PointReference.h"

// only uses the constant basis functions
#include "BasisFunctionsPiecewiseConstant.h"

namespace Utilities {

BasisFunction3DVertexCube::BasisFunction3DVertexCube(std::size_t node) {
    logger.assert_debug(node < 8, "A cube only has 8 nodes");
    nodePosition0_ = (static_cast<int>(node) % 2) * 2 - 1;
    nodePosition1_ = ((static_cast<int>(node) / 2) % 2) * 2 - 1;
    nodePosition2_ = (static_cast<int>(node) / 4) * 2 - 1;
}

double BasisFunction3DVertexCube::eval(
    const Geometry::PointReference<3>& p) const {
    return (1 + nodePosition0_ * p[0]) * (1 + nodePosition1_ * p[1]) *
           (1 + nodePosition2_ * p[2]) / 8.;
}

double BasisFunction3DVertexCube::evalDeriv0(
    const Geometry::PointReference<3>& p) const {
    return nodePosition0_ * (1 + nodePosition1_ * p[1]) *
           (1 + nodePosition2_ * p[2]) / 8.;
}

double BasisFunction3DVertexCube::evalDeriv1(
    const Geometry::PointReference<3>& p) const {
    return nodePosition1_ * (1 + nodePosition0_ * p[0]) *
           (1 + nodePosition2_ * p[2]) / 8.;
}

double BasisFunction3DVertexCube::evalDeriv2(
    const Geometry::PointReference<3>& p) const {
    return nodePosition2_ * (1 + nodePosition0_ * p[0]) *
           (1 + nodePosition1_ * p[1]) / 8.;
}

BasisFunction3DEdgeCube_0::BasisFunction3DEdgeCube_0(
    std::size_t node0, std::size_t node1, std::size_t polynomialOrder)
    : polynomialOrder_(polynomialOrder) {
    logger.assert_debug(node0 < 8, "A cube only has 8 nodes");
    logger.assert_debug(node1 < 8, "A cube only has 8 nodes");
    logger.assert_debug(
        node0 / 4 == node1 / 4,
        "This class is intended for edges aligned in the x-direction");
    logger.assert_debug(
        (node0 / 2) % 2 == (node1 / 2) % 2,
        "This class is intended for edges aligned in the x-direction");
    edgePosition1_ = ((static_cast<int>(node0) / 2) % 2) * 2 - 1;
    edgePosition2_ = (static_cast<int>(node0) / 4) * 2 - 1;
    mirroring_ = node0 < node1 ? 1 : -1;
}

double BasisFunction3DEdgeCube_0::eval(
    const Geometry::PointReference<3>& p) const {
    return (1 - p[0]) * (1 + p[0]) * (1 + edgePosition1_ * p[1]) *
           (1 + edgePosition2_ * p[2]) *
           LobattoPolynomial(polynomialOrder_, mirroring_ * p[0]) / 16.;
}

double BasisFunction3DEdgeCube_0::evalDeriv0(
    const Geometry::PointReference<3>& p) const {
    return -p[0] * (1 + edgePosition1_ * p[1]) * (1 + edgePosition2_ * p[2]) *
               LobattoPolynomial(polynomialOrder_, mirroring_ * p[0]) / 8. +
           mirroring_ * (1 - p[0]) * (1 + p[0]) * (1 + edgePosition1_ * p[1]) *
               (1 + edgePosition2_ * p[2]) *
               LobattoPolynomialDerivative(polynomialOrder_,
                                           mirroring_ * p[0]) /
               16.;
}

double BasisFunction3DEdgeCube_0::evalDeriv1(
    const Geometry::PointReference<3>& p) const {
    return (1 - p[0]) * (1 + p[0]) * edgePosition1_ *
           (1 + edgePosition2_ * p[2]) *
           LobattoPolynomial(polynomialOrder_, mirroring_ * p[0]) / 16.;
}

double BasisFunction3DEdgeCube_0::evalDeriv2(
    const Geometry::PointReference<3>& p) const {
    return (1 - p[0]) * (1 + p[0]) * edgePosition2_ *
           (1 + edgePosition1_ * p[1]) *
           LobattoPolynomial(polynomialOrder_, mirroring_ * p[0]) / 16.;
}

BasisFunction3DEdgeCube_1::BasisFunction3DEdgeCube_1(
    std::size_t node0, std::size_t node1, std::size_t polynomialOrder)
    : polynomialOrder_(polynomialOrder) {
    logger.assert_debug(node0 < 8, "A cube only has 8 nodes");
    logger.assert_debug(node1 < 8, "A cube only has 8 nodes");
    logger.assert_debug(
        node0 / 4 == node1 / 4,
        "This class is intended for edges aligned in the y-direction");
    logger.assert_debug(
        node0 % 2 == node1 % 2,
        "This class is intended for edges aligned in the y-direction");
    edgePosition0_ = (static_cast<int>(node0) % 2) * 2 - 1;
    edgePosition2_ = (static_cast<int>(node0) / 4) * 2 - 1;
    mirroring_ = node0 < node1 ? 1 : -1;
}

double BasisFunction3DEdgeCube_1::eval(
    const Geometry::PointReference<3>& p) const {
    return (1 - p[1]) * (1 + p[1]) * (1 + edgePosition0_ * p[0]) *
           (1 + edgePosition2_ * p[2]) *
           LobattoPolynomial(polynomialOrder_, mirroring_ * p[1]) / 16.;
}

double BasisFunction3DEdgeCube_1::evalDeriv0(
    const Geometry::PointReference<3>& p) const {
    return (1 - p[1]) * (1 + p[1]) * edgePosition0_ *
           (1 + edgePosition2_ * p[2]) *
           LobattoPolynomial(polynomialOrder_, mirroring_ * p[1]) / 16.;
}

double BasisFunction3DEdgeCube_1::evalDeriv1(
    const Geometry::PointReference<3>& p) const {
    return -p[1] * (1 + edgePosition0_ * p[0]) * (1 + edgePosition2_ * p[2]) *
               LobattoPolynomial(polynomialOrder_, mirroring_ * p[1]) / 8. +
           mirroring_ * (1 - p[1]) * (1 + p[1]) * (1 + edgePosition0_ * p[0]) *
               (1 + edgePosition2_ * p[2]) *
               LobattoPolynomialDerivative(polynomialOrder_,
                                           mirroring_ * p[1]) /
               16.;
}

double BasisFunction3DEdgeCube_1::evalDeriv2(
    const Geometry::PointReference<3>& p) const {
    return (1 - p[1]) * (1 + p[1]) * edgePosition2_ *
           (1 + edgePosition0_ * p[0]) *
           LobattoPolynomial(polynomialOrder_, mirroring_ * p[1]) / 16.;
}

BasisFunction3DEdgeCube_2::BasisFunction3DEdgeCube_2(
    std::size_t node0, std::size_t node1, std::size_t polynomialOrder)
    : polynomialOrder_(polynomialOrder) {
    logger.assert_debug(node0 < 8, "A cube only has 8 nodes");
    logger.assert_debug(node1 < 8, "A cube only has 8 nodes");
    logger.assert_debug(
        node0 % 2 == node1 % 2,
        "This class is intended for edges aligned in the z-direction");
    logger.assert_debug(
        (node0 / 2) % 2 == (node1 / 2) % 2,
        "This class is intended for edges aligned in the z-direction");
    edgePosition0_ = (static_cast<int>(node0) % 2) * 2 - 1;
    edgePosition1_ = ((static_cast<int>(node0) / 2) % 2) * 2 - 1;
    mirroring_ = node0 < node1 ? 1 : -1;
}

double BasisFunction3DEdgeCube_2::evalDeriv0(
    const Geometry::PointReference<3>& p) const {
    return (1 - p[2]) * (1 + p[2]) * edgePosition0_ *
           (1 + edgePosition1_ * p[1]) *
           LobattoPolynomial(polynomialOrder_, mirroring_ * p[2]) / 16.;
}

double BasisFunction3DEdgeCube_2::eval(
    const Geometry::PointReference<3>& p) const {
    return (1 - p[2]) * (1 + p[2]) * (1 + edgePosition0_ * p[0]) *
           (1 + edgePosition1_ * p[1]) *
           LobattoPolynomial(polynomialOrder_, mirroring_ * p[2]) / 16.;
}

double BasisFunction3DEdgeCube_2::evalDeriv1(
    const Geometry::PointReference<3>& p) const {
    return (1 - p[2]) * (1 + p[2]) * edgePosition1_ *
           (1 + edgePosition0_ * p[0]) *
           LobattoPolynomial(polynomialOrder_, mirroring_ * p[2]) / 16.;
}

double BasisFunction3DEdgeCube_2::evalDeriv2(
    const Geometry::PointReference<3>& p) const {
    return -p[2] * (1 + edgePosition0_ * p[0]) * (1 + edgePosition1_ * p[1]) *
               LobattoPolynomial(polynomialOrder_, mirroring_ * p[2]) / 8. +
           mirroring_ * (1 - p[2]) * (1 + p[2]) * (1 + edgePosition0_ * p[0]) *
               (1 + edgePosition1_ * p[1]) *
               LobattoPolynomialDerivative(polynomialOrder_,
                                           mirroring_ * p[2]) /
               16.;
}

BasisFunction3DFaceCube_0::BasisFunction3DFaceCube_0(
    std::size_t node0, std::size_t node1, std::size_t node2,
    std::size_t polynomialOrder1, std::size_t polynomialOrder2)
    : polynomialOrder1_(polynomialOrder1), polynomialOrder2_(polynomialOrder2) {
    logger.assert_debug(node0 < 8, "A cube only has 8 nodes");
    logger.assert_debug(node1 < 8, "A cube only has 8 nodes");
    logger.assert_debug(node2 < 8, "A cube only has 8 nodes");
    logger.assert_debug(
        node0 % 2 == node1 % 2,
        "This class is intended for faces aligned in the x-direction");
    logger.assert_debug(
        node0 % 2 == node2 % 2,
        "This class is intended for faces aligned in the x-direction");
    mirroring1_ = node0 < node1 ? 1 : -1;  // choices about mirroring need only
                                           // be consistent for one face
    mirroring2_ = node0 < node2 ? 1 : -1;
    facePosition_ = (static_cast<int>(node0) % 2) * 2 - 1;
}

double BasisFunction3DFaceCube_0::eval(
    const Geometry::PointReference<3>& p) const {
    return (1 + facePosition_ * p[0]) * (1 - p[1]) * (1 + p[1]) * (1 - p[2]) *
           (1 + p[2]) *
           LobattoPolynomial(polynomialOrder1_, mirroring1_ * p[1]) *
           LobattoPolynomial(polynomialOrder2_, mirroring2_ * p[2]) / 32.;
}

double BasisFunction3DFaceCube_0::evalDeriv1(
    const Geometry::PointReference<3>& p) const {
    return -p[1] * (1 + facePosition_ * p[0]) * (1 - p[2]) * (1 + p[2]) *
               LobattoPolynomial(polynomialOrder1_, mirroring1_ * p[1]) *
               LobattoPolynomial(polynomialOrder2_, mirroring2_ * p[2]) / 16. +
           mirroring1_ * (1 + facePosition_ * p[0]) * (1 - p[1]) * (1 + p[1]) *
               (1 - p[2]) * (1 + p[2]) *
               LobattoPolynomialDerivative(polynomialOrder1_,
                                           mirroring1_ * p[1]) *
               LobattoPolynomial(polynomialOrder2_, mirroring2_ * p[2]) / 32.;
}

double BasisFunction3DFaceCube_0::evalDeriv2(
    const Geometry::PointReference<3>& p) const {
    return -p[2] * (1 + facePosition_ * p[0]) * (1 - p[1]) * (1 + p[1]) *
               LobattoPolynomial(polynomialOrder1_, mirroring1_ * p[1]) *
               LobattoPolynomial(polynomialOrder2_, mirroring2_ * p[2]) / 16. +
           mirroring2_ * (1 + facePosition_ * p[0]) * (1 - p[1]) * (1 + p[1]) *
               (1 - p[2]) * (1 + p[2]) *
               LobattoPolynomial(polynomialOrder1_, mirroring1_ * p[1]) *
               LobattoPolynomialDerivative(polynomialOrder2_,
                                           mirroring2_ * p[2]) /
               32.;
}

BasisFunction3DFaceCube_1::BasisFunction3DFaceCube_1(
    std::size_t node0, std::size_t node1, std::size_t node2,
    std::size_t polynomialOrder0, std::size_t polynomialOrder2)
    : polynomialOrder0_(polynomialOrder0), polynomialOrder2_(polynomialOrder2) {
    logger.assert_debug(node0 < 8, "A cube only has 8 nodes");
    logger.assert_debug(node1 < 8, "A cube only has 8 nodes");
    logger.assert_debug(node2 < 8, "A cube only has 8 nodes");
    logger.assert_debug(
        (node0 / 2) % 2 == (node1 / 2) % 2,
        "This class is intended for faces aligned in the y-direction");
    logger.assert_debug(
        (node0 / 2) % 2 == (node2 / 2) % 2,
        "This class is intended for faces aligned in the y-direction");
    mirroring0_ = node0 < node1 ? 1 : -1;
    mirroring2_ = node0 < node2 ? 1 : -1;
    facePosition_ = ((static_cast<int>(node0) / 2) % 2) * 2 - 1;
}

double BasisFunction3DFaceCube_1::eval(
    const Geometry::PointReference<3>& p) const {
    return (1 + facePosition_ * p[1]) * (1 - p[0]) * (1 + p[0]) * (1 - p[2]) *
           (1 + p[2]) *
           LobattoPolynomial(polynomialOrder0_, mirroring0_ * p[0]) *
           LobattoPolynomial(polynomialOrder2_, mirroring2_ * p[2]) / 32.;
}

double BasisFunction3DFaceCube_1::evalDeriv0(
    const Geometry::PointReference<3>& p) const {
    return -p[0] * (1 + facePosition_ * p[1]) * (1 - p[2]) * (1 + p[2]) *
               LobattoPolynomial(polynomialOrder0_, mirroring0_ * p[0]) *
               LobattoPolynomial(polynomialOrder2_, mirroring2_ * p[2]) / 16. +
           mirroring0_ * (1 + facePosition_ * p[1]) * (1 - p[0]) * (1 + p[0]) *
               (1 - p[2]) * (1 + p[2]) *
               LobattoPolynomialDerivative(polynomialOrder0_,
                                           mirroring0_ * p[0]) *
               LobattoPolynomial(polynomialOrder2_, mirroring2_ * p[2]) / 32.;
}

double BasisFunction3DFaceCube_1::evalDeriv1(
    const Geometry::PointReference<3>& p) const {
    return facePosition_ * (1 - p[0]) * (1 + p[0]) * (1 - p[2]) * (1 + p[2]) *
           LobattoPolynomial(polynomialOrder0_, mirroring0_ * p[0]) *
           LobattoPolynomial(polynomialOrder2_, mirroring2_ * p[2]) / 32.;
}

double BasisFunction3DFaceCube_1::evalDeriv2(
    const Geometry::PointReference<3>& p) const {
    return -p[2] * (1 + facePosition_ * p[1]) * (1 - p[0]) * (1 + p[0]) *
               LobattoPolynomial(polynomialOrder0_, mirroring0_ * p[0]) *
               LobattoPolynomial(polynomialOrder2_, mirroring2_ * p[2]) / 16. +
           mirroring2_ * (1 + facePosition_ * p[1]) * (1 - p[0]) * (1 + p[0]) *
               (1 - p[2]) * (1 + p[2]) *
               LobattoPolynomial(polynomialOrder0_, mirroring0_ * p[0]) *
               LobattoPolynomialDerivative(polynomialOrder2_,
                                           mirroring2_ * p[2]) /
               32.;
}

BasisFunction3DFaceCube_2::BasisFunction3DFaceCube_2(
    std::size_t node0, std::size_t node1, std::size_t node2,
    std::size_t polynomialOrder0, std::size_t polynomialOrder1)
    : polynomialOrder0_(polynomialOrder0), polynomialOrder1_(polynomialOrder1) {
    logger.assert_debug(node0 < 8, "A cube only has 8 nodes");
    logger.assert_debug(node1 < 8, "A cube only has 8 nodes");
    logger.assert_debug(node2 < 8, "A cube only has 8 nodes");
    logger.assert_debug(
        node0 / 4 == node1 / 4,
        "This class is intended for faces aligned in the x-direction");
    logger.assert_debug(
        node0 / 4 == node2 / 4,
        "This class is intended for faces aligned in the x-direction");
    mirroring0_ = node0 < node1 ? 1 : -1;
    mirroring1_ = node0 < node2 ? 1 : -1;
    facePosition_ = (static_cast<int>(node0) / 4) * 2 - 1;
}

double BasisFunction3DFaceCube_2::eval(
    const Geometry::PointReference<3>& p) const {
    return (1 + facePosition_ * p[2]) * (1 - p[0]) * (1 + p[0]) * (1 - p[1]) *
           (1 + p[1]) *
           LobattoPolynomial(polynomialOrder0_, mirroring0_ * p[0]) *
           LobattoPolynomial(polynomialOrder1_, mirroring1_ * p[1]) / 32.;
}

double BasisFunction3DFaceCube_2::evalDeriv0(
    const Geometry::PointReference<3>& p) const {
    return -p[0] * (1 + facePosition_ * p[2]) * (1 - p[1]) * (1 + p[1]) *
               LobattoPolynomial(polynomialOrder1_, mirroring1_ * p[1]) *
               LobattoPolynomial(polynomialOrder0_, mirroring0_ * p[0]) / 16. +
           mirroring0_ * (1 + facePosition_ * p[2]) * (1 - p[0]) * (1 + p[0]) *
               (1 - p[1]) * (1 + p[1]) *
               LobattoPolynomial(polynomialOrder1_, mirroring1_ * p[1]) *
               LobattoPolynomialDerivative(polynomialOrder0_,
                                           mirroring0_ * p[0]) /
               32.;
}

double BasisFunction3DFaceCube_2::evalDeriv1(
    const Geometry::PointReference<3>& p) const {
    return -p[1] * (1 + facePosition_ * p[2]) * (1 - p[0]) * (1 + p[0]) *
               LobattoPolynomial(polynomialOrder1_, mirroring1_ * p[1]) *
               LobattoPolynomial(polynomialOrder0_, mirroring0_ * p[0]) / 16. +
           mirroring1_ * (1 + facePosition_ * p[2]) * (1 - p[0]) * (1 + p[0]) *
               (1 - p[1]) * (1 + p[1]) *
               LobattoPolynomialDerivative(polynomialOrder1_,
                                           mirroring1_ * p[1]) *
               LobattoPolynomial(polynomialOrder0_, mirroring0_ * p[0]) / 32.;
}

double BasisFunction3DFaceCube_2::evalDeriv2(
    const Geometry::PointReference<3>& p) const {
    return facePosition_ * (1 - p[0]) * (1 + p[0]) * (1 - p[1]) * (1 + p[1]) *
           LobattoPolynomial(polynomialOrder1_, mirroring1_ * p[1]) *
           LobattoPolynomial(polynomialOrder0_, mirroring0_ * p[0]) / 32.;
}

double BasisFunction3DInteriorCube::eval(
    const Geometry::PointReference<3>& p) const {
    return (1 - p[0]) * (1 + p[0]) * (1 - p[1]) * (1 + p[1]) * (1 - p[2]) *
           (1 + p[2]) * LobattoPolynomial(polynomialOrder0_, p[0]) *
           LobattoPolynomial(polynomialOrder1_, p[1]) *
           LobattoPolynomial(polynomialOrder2_, p[2]) / 64.;
}

double BasisFunction3DInteriorCube::evalDeriv0(
    const Geometry::PointReference<3>& p) const {
    return (1 - p[1]) * (1 + p[1]) * (1 - p[2]) * (1 + p[2]) *
           LobattoPolynomial(polynomialOrder1_, p[1]) *
           LobattoPolynomial(polynomialOrder2_, p[2]) *
           (-p[0] * LobattoPolynomial(polynomialOrder0_, p[0]) / 32. +
            (1 - p[0]) * (1 + p[0]) *
                LobattoPolynomialDerivative(polynomialOrder0_, p[0]) / 64.);
}

double BasisFunction3DInteriorCube::evalDeriv1(
    const Geometry::PointReference<3>& p) const {
    return (1 - p[0]) * (1 + p[0]) * (1 - p[2]) * (1 + p[2]) *
           LobattoPolynomial(polynomialOrder0_, p[0]) *
           LobattoPolynomial(polynomialOrder2_, p[2]) *
           (-p[1] * LobattoPolynomial(polynomialOrder1_, p[1]) / 32. +
            (1 - p[1]) * (1 + p[1]) *
                LobattoPolynomialDerivative(polynomialOrder1_, p[1]) / 64.);
}

double BasisFunction3DInteriorCube::evalDeriv2(
    const Geometry::PointReference<3>& p) const {
    return (1 - p[0]) * (1 + p[0]) * (1 - p[1]) * (1 + p[1]) *
           LobattoPolynomial(polynomialOrder0_, p[0]) *
           LobattoPolynomial(polynomialOrder1_, p[1]) *
           (-p[2] * LobattoPolynomial(polynomialOrder2_, p[2]) / 32. +
            (1 - p[2]) * (1 + p[2]) *
                LobattoPolynomialDerivative(polynomialOrder2_, p[2]) / 64.);
}

Base::BasisFunctionSet* createDGBasisFunctionSet3DH1Cube(std::size_t order) {
    Base::BasisFunctionSet* result(new Base::BasisFunctionSet(order));
    if (order > 0) {
        Geometry::ReferenceCube& cube = Geometry::ReferenceCube::Instance();
        std::vector<std::size_t> vectorOfPointIndices(4);
        for (std::size_t i = 0; i < cube.getNumberOfCodim3Entities(); ++i) {
            result->addBasisFunction(new BasisFunction3DVertexCube(i));
        }
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            for (std::size_t i = 0; i < 4; ++i) {
                vectorOfPointIndices = cube.getCodim2EntityLocalIndices(i);
                result->addBasisFunction(new BasisFunction3DEdgeCube_0(
                    vectorOfPointIndices[0], vectorOfPointIndices[1], j));
            }
            for (std::size_t i = 4; i < 8; ++i) {
                vectorOfPointIndices = cube.getCodim2EntityLocalIndices(i);
                result->addBasisFunction(new BasisFunction3DEdgeCube_1(
                    vectorOfPointIndices[0], vectorOfPointIndices[1], j));
            }
            for (std::size_t i = 8; i < 12; ++i) {
                vectorOfPointIndices = cube.getCodim2EntityLocalIndices(i);
                result->addBasisFunction(new BasisFunction3DEdgeCube_2(
                    vectorOfPointIndices[0], vectorOfPointIndices[1], j));
            }
            result->addBasisFunction(
                new BasisFunction3DFaceCube_2(0, 1, 2, j, j));
            result->addBasisFunction(
                new BasisFunction3DFaceCube_1(0, 1, 4, j, j));
            result->addBasisFunction(
                new BasisFunction3DFaceCube_0(0, 2, 4, j, j));
            result->addBasisFunction(
                new BasisFunction3DFaceCube_0(1, 3, 5, j, j));
            result->addBasisFunction(
                new BasisFunction3DFaceCube_1(2, 3, 6, j, j));
            result->addBasisFunction(
                new BasisFunction3DFaceCube_2(4, 5, 6, j, j));
            result->addBasisFunction(new BasisFunction3DInteriorCube(j, j, j));
            for (std::size_t i = 0; i < j; ++i) {
                result->addBasisFunction(
                    new BasisFunction3DFaceCube_2(0, 1, 2, i, j));
                result->addBasisFunction(
                    new BasisFunction3DFaceCube_1(0, 1, 4, i, j));
                result->addBasisFunction(
                    new BasisFunction3DFaceCube_0(0, 2, 4, i, j));
                result->addBasisFunction(
                    new BasisFunction3DFaceCube_0(1, 3, 5, i, j));
                result->addBasisFunction(
                    new BasisFunction3DFaceCube_1(2, 3, 6, i, j));
                result->addBasisFunction(
                    new BasisFunction3DFaceCube_2(4, 5, 6, i, j));
                result->addBasisFunction(
                    new BasisFunction3DFaceCube_2(0, 1, 2, j, i));
                result->addBasisFunction(
                    new BasisFunction3DFaceCube_1(0, 1, 4, j, i));
                result->addBasisFunction(
                    new BasisFunction3DFaceCube_0(0, 2, 4, j, i));
                result->addBasisFunction(
                    new BasisFunction3DFaceCube_0(1, 3, 5, j, i));
                result->addBasisFunction(
                    new BasisFunction3DFaceCube_1(2, 3, 6, j, i));
                result->addBasisFunction(
                    new BasisFunction3DFaceCube_2(4, 5, 6, j, i));
                result->addBasisFunction(
                    new BasisFunction3DInteriorCube(i, i, j));
                result->addBasisFunction(
                    new BasisFunction3DInteriorCube(i, j, i));
                result->addBasisFunction(
                    new BasisFunction3DInteriorCube(j, i, i));
                result->addBasisFunction(
                    new BasisFunction3DInteriorCube(j, j, i));
                result->addBasisFunction(
                    new BasisFunction3DInteriorCube(j, i, j));
                result->addBasisFunction(
                    new BasisFunction3DInteriorCube(i, j, j));
                for (std::size_t k = 0; k < i; ++k) {
                    result->addBasisFunction(
                        new BasisFunction3DInteriorCube(i, j, k));
                    result->addBasisFunction(
                        new BasisFunction3DInteriorCube(i, k, j));
                    result->addBasisFunction(
                        new BasisFunction3DInteriorCube(j, i, k));
                    result->addBasisFunction(
                        new BasisFunction3DInteriorCube(j, k, i));
                    result->addBasisFunction(
                        new BasisFunction3DInteriorCube(k, i, j));
                    result->addBasisFunction(
                        new BasisFunction3DInteriorCube(k, j, i));
                }
            }
        }
    } else {
        addPiecewiseConstantBasisFunction3D(*result);
    }
    return result;
}

Base::BasisFunctionSet* createInteriorBasisFunctionSet3DH1Cube(
    std::size_t order) {
    logger.assert_debug(order > 0,
                        "Trying to create a conforming, constant basis "
                        "function set, did you mean the constant solution?");
    Base::BasisFunctionSet* result(new Base::BasisFunctionSet(order));
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            for (std::size_t k = 0; k + 2 <= order; ++k) {
                result->addBasisFunction(
                    new BasisFunction3DInteriorCube(i, j, k));
            }
        }
    }
    return result;
}

std::vector<const Base::BasisFunctionSet*> createVertexBasisFunctionSet3DH1Cube(
    std::size_t order) {
    logger.assert_debug(order > 0,
                        "Trying to create a conforming, constant basis "
                        "function set, did you mean the constant solution?");
    std::vector<const Base::BasisFunctionSet*> result;
    Base::BasisFunctionSet* set;
    Geometry::ReferenceCube& cube = Geometry::ReferenceCube::Instance();
    for (std::size_t i = 0; i < cube.getNumberOfCodim3Entities(); ++i) {
        set = new Base::BasisFunctionSet(order);
        set->addBasisFunction(new BasisFunction3DVertexCube(i));
        result.push_back(set);
    }
    return result;
}

std::vector<const Base::OrientedBasisFunctionSet*>
    createEdgeBasisFunctionSet3DH1Cube(std::size_t order) {
    logger.assert_debug(order > 0,
                        "Trying to create a conforming, constant basis "
                        "function set, did you mean the constant solution?");
    std::vector<const Base::OrientedBasisFunctionSet*> result;
    Base::OrientedBasisFunctionSet* set;
    Geometry::ReferenceCube& cube = Geometry::ReferenceCube::Instance();
    std::vector<std::size_t> vectorOfPointIndices(2);
    for (std::size_t i = 0; i < 4; ++i) {
        vectorOfPointIndices = cube.getCodim2EntityLocalIndices(i);
        set = new Base::OrientedBasisFunctionSet(order, 0, i);
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DEdgeCube_0(
                vectorOfPointIndices[0], vectorOfPointIndices[1], j));
        }
        result.push_back(set);
        set = new Base::OrientedBasisFunctionSet(order, 1, i);
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DEdgeCube_0(
                vectorOfPointIndices[1], vectorOfPointIndices[0], j));
        }
        result.push_back(set);
    }
    for (std::size_t i = 4; i < 8; ++i) {
        vectorOfPointIndices = cube.getCodim2EntityLocalIndices(i);
        set = new Base::OrientedBasisFunctionSet(order, 0, i);
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DEdgeCube_1(
                vectorOfPointIndices[0], vectorOfPointIndices[1], j));
        }
        result.push_back(set);
        set = new Base::OrientedBasisFunctionSet(order, 1, i);
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DEdgeCube_1(
                vectorOfPointIndices[1], vectorOfPointIndices[0], j));
        }
        result.push_back(set);
    }
    for (std::size_t i = 8; i < 12; ++i) {
        vectorOfPointIndices = cube.getCodim2EntityLocalIndices(i);
        set = new Base::OrientedBasisFunctionSet(order, 0, i);
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DEdgeCube_2(
                vectorOfPointIndices[0], vectorOfPointIndices[1], j));
        }
        result.push_back(set);
        set = new Base::OrientedBasisFunctionSet(order, 1, i);
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DEdgeCube_2(
                vectorOfPointIndices[1], vectorOfPointIndices[0], j));
        }
        result.push_back(set);
    }
    return result;
}

std::vector<const Base::OrientedBasisFunctionSet*>
    createFaceBasisFunctionSet3DH1Cube(std::size_t order) {
    logger.assert_debug(order > 0,
                        "Trying to create a conforming, constant basis "
                        "function set, did you mean the constant solution?");
    std::vector<const Base::OrientedBasisFunctionSet*> result;
    Base::OrientedBasisFunctionSet* set;
    Geometry::ReferenceCube& cube = Geometry::ReferenceCube::Instance();
    std::vector<std::size_t> vectorOfPointIndices(4);
    vectorOfPointIndices = cube.getCodim1EntityLocalIndices(0);
    set = new Base::OrientedBasisFunctionSet(order, 0, 0);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_2(
                vectorOfPointIndices[0], vectorOfPointIndices[1],
                vectorOfPointIndices[2], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 1, 0);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_2(
                vectorOfPointIndices[1], vectorOfPointIndices[3],
                vectorOfPointIndices[0], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 2, 0);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_2(
                vectorOfPointIndices[3], vectorOfPointIndices[2],
                vectorOfPointIndices[1], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 3, 0);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_2(
                vectorOfPointIndices[2], vectorOfPointIndices[0],
                vectorOfPointIndices[3], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 4, 0);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_2(
                vectorOfPointIndices[2], vectorOfPointIndices[3],
                vectorOfPointIndices[0], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 5, 0);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_2(
                vectorOfPointIndices[1], vectorOfPointIndices[0],
                vectorOfPointIndices[3], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 6, 0);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_2(
                vectorOfPointIndices[3], vectorOfPointIndices[1],
                vectorOfPointIndices[2], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 7, 0);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_2(
                vectorOfPointIndices[0], vectorOfPointIndices[2],
                vectorOfPointIndices[1], i, j));
        }
    }
    result.push_back(set);
    vectorOfPointIndices = cube.getCodim1EntityLocalIndices(1);
    set = new Base::OrientedBasisFunctionSet(order, 0, 1);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_1(
                vectorOfPointIndices[0], vectorOfPointIndices[1],
                vectorOfPointIndices[2], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 1, 1);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_1(
                vectorOfPointIndices[1], vectorOfPointIndices[3],
                vectorOfPointIndices[0], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 2, 1);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_1(
                vectorOfPointIndices[3], vectorOfPointIndices[2],
                vectorOfPointIndices[1], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 3, 1);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_1(
                vectorOfPointIndices[2], vectorOfPointIndices[0],
                vectorOfPointIndices[3], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 4, 1);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_1(
                vectorOfPointIndices[2], vectorOfPointIndices[3],
                vectorOfPointIndices[0], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 5, 1);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_1(
                vectorOfPointIndices[1], vectorOfPointIndices[0],
                vectorOfPointIndices[3], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 6, 1);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_1(
                vectorOfPointIndices[3], vectorOfPointIndices[1],
                vectorOfPointIndices[2], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 7, 1);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_1(
                vectorOfPointIndices[0], vectorOfPointIndices[2],
                vectorOfPointIndices[1], i, j));
        }
    }
    result.push_back(set);
    vectorOfPointIndices = cube.getCodim1EntityLocalIndices(2);
    set = new Base::OrientedBasisFunctionSet(order, 0, 2);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_0(
                vectorOfPointIndices[0], vectorOfPointIndices[1],
                vectorOfPointIndices[2], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 1, 2);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_0(
                vectorOfPointIndices[1], vectorOfPointIndices[3],
                vectorOfPointIndices[0], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 2, 2);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_0(
                vectorOfPointIndices[3], vectorOfPointIndices[2],
                vectorOfPointIndices[1], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 3, 2);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_0(
                vectorOfPointIndices[2], vectorOfPointIndices[0],
                vectorOfPointIndices[3], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 4, 2);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_0(
                vectorOfPointIndices[2], vectorOfPointIndices[3],
                vectorOfPointIndices[0], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 5, 2);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_0(
                vectorOfPointIndices[1], vectorOfPointIndices[0],
                vectorOfPointIndices[3], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 6, 2);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_0(
                vectorOfPointIndices[3], vectorOfPointIndices[1],
                vectorOfPointIndices[2], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 7, 2);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_0(
                vectorOfPointIndices[0], vectorOfPointIndices[2],
                vectorOfPointIndices[1], i, j));
        }
    }
    result.push_back(set);
    vectorOfPointIndices = cube.getCodim1EntityLocalIndices(3);
    set = new Base::OrientedBasisFunctionSet(order, 0, 3);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_0(
                vectorOfPointIndices[0], vectorOfPointIndices[1],
                vectorOfPointIndices[2], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 1, 3);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_0(
                vectorOfPointIndices[1], vectorOfPointIndices[3],
                vectorOfPointIndices[0], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 2, 3);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_0(
                vectorOfPointIndices[3], vectorOfPointIndices[2],
                vectorOfPointIndices[1], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 3, 3);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_0(
                vectorOfPointIndices[2], vectorOfPointIndices[0],
                vectorOfPointIndices[3], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 4, 3);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_0(
                vectorOfPointIndices[2], vectorOfPointIndices[3],
                vectorOfPointIndices[0], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 5, 3);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_0(
                vectorOfPointIndices[1], vectorOfPointIndices[0],
                vectorOfPointIndices[3], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 6, 3);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_0(
                vectorOfPointIndices[3], vectorOfPointIndices[1],
                vectorOfPointIndices[2], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 7, 3);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_0(
                vectorOfPointIndices[0], vectorOfPointIndices[2],
                vectorOfPointIndices[1], i, j));
        }
    }
    result.push_back(set);
    vectorOfPointIndices = cube.getCodim1EntityLocalIndices(4);
    set = new Base::OrientedBasisFunctionSet(order, 0, 4);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_1(
                vectorOfPointIndices[0], vectorOfPointIndices[1],
                vectorOfPointIndices[2], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 1, 4);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_1(
                vectorOfPointIndices[1], vectorOfPointIndices[3],
                vectorOfPointIndices[0], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 2, 4);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_1(
                vectorOfPointIndices[3], vectorOfPointIndices[2],
                vectorOfPointIndices[1], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 3, 4);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_1(
                vectorOfPointIndices[2], vectorOfPointIndices[0],
                vectorOfPointIndices[3], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 4, 4);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_1(
                vectorOfPointIndices[2], vectorOfPointIndices[3],
                vectorOfPointIndices[0], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 5, 4);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_1(
                vectorOfPointIndices[1], vectorOfPointIndices[0],
                vectorOfPointIndices[3], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 6, 4);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_1(
                vectorOfPointIndices[3], vectorOfPointIndices[1],
                vectorOfPointIndices[2], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 7, 4);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_1(
                vectorOfPointIndices[0], vectorOfPointIndices[2],
                vectorOfPointIndices[1], i, j));
        }
    }
    result.push_back(set);
    vectorOfPointIndices = cube.getCodim1EntityLocalIndices(5);
    set = new Base::OrientedBasisFunctionSet(order, 0, 5);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_2(
                vectorOfPointIndices[0], vectorOfPointIndices[1],
                vectorOfPointIndices[2], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 1, 5);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_2(
                vectorOfPointIndices[1], vectorOfPointIndices[3],
                vectorOfPointIndices[0], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 2, 5);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_2(
                vectorOfPointIndices[3], vectorOfPointIndices[2],
                vectorOfPointIndices[1], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 3, 5);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_2(
                vectorOfPointIndices[2], vectorOfPointIndices[0],
                vectorOfPointIndices[3], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 4, 5);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_2(
                vectorOfPointIndices[2], vectorOfPointIndices[3],
                vectorOfPointIndices[0], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 5, 5);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_2(
                vectorOfPointIndices[1], vectorOfPointIndices[0],
                vectorOfPointIndices[3], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 6, 5);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_2(
                vectorOfPointIndices[3], vectorOfPointIndices[1],
                vectorOfPointIndices[2], i, j));
        }
    }
    result.push_back(set);
    set = new Base::OrientedBasisFunctionSet(order, 7, 5);
    for (std::size_t i = 0; i + 2 <= order; ++i) {
        for (std::size_t j = 0; j + 2 <= order; ++j) {
            set->addBasisFunction(new BasisFunction3DFaceCube_2(
                vectorOfPointIndices[0], vectorOfPointIndices[2],
                vectorOfPointIndices[1], i, j));
        }
    }
    result.push_back(set);
    return result;
}

double BasisFunction3DFaceCube_0::evalDeriv0(
    const Geometry::PointReference<3>& p) const {
    return facePosition_ * (1 - p[1]) * (1 + p[1]) * (1 - p[2]) * (1 + p[2]) *
           LobattoPolynomial(polynomialOrder1_, mirroring1_ * p[1]) *
           LobattoPolynomial(polynomialOrder2_, mirroring2_ * p[2]) / 32.;
}

}  // namespace Utilities

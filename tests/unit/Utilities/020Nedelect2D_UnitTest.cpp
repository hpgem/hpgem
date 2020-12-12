/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2020, University of Twente
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

#include "Utilities/BasisFunctions2DNedelec.h"
#include "Integration/QuadratureRules/AllGaussQuadratureRules.h"
#include "Integration/QuadratureRules/GaussQuadratureRule.h"
#include "Geometry/ReferenceTriangle.h"
#include "../catch.hpp"

using namespace hpgem;

TEST_CASE("Nedelec 2D: Basic properties", "[Nedelec2D]") {
    auto p = GENERATE(1, 2, 3, 4);
    auto* bfSet = Utilities::createDGBasisFunctionSet2DNedelec(p);

    // Basic self check
    REQUIRE(bfSet->getOrder() == p);
    // There should be p*(p+2) basis functions
    REQUIRE(bfSet->size() == p * (p + 2));
}

TEST_CASE("Nedelec 2D: trace values", "[Nedelec2D]") {
    LinearAlgebra::SmallVector<2> tangents[3];
    Geometry::PointReference<2> sidePoints[3];
    // Nedelec elements are for H(curl), for conforming elements this requires
    // continuous tangential part. Hence, they are constructed to have zero
    // tangential trace on all sides but (at most) 1. For order p there are
    // exactly p basis functions with non zero trace.

    // -y
    tangents[0][0] = 1.0;
    tangents[0][1] = 0.0;
    // x,y
    tangents[1][0] = 1.0;
    tangents[1][1] = -1.0;
    //-x
    tangents[2][0] = 0.0;
    tangents[2][1] = -1.0;

    // Random points on the side to prevent hitting an accidental zero
    sidePoints[0][0] = 1 / M_PI;
    sidePoints[1][0] = 1 / M_PI;
    sidePoints[1][1] = 1 - sidePoints[1][0];
    sidePoints[2][1] = 1 / M_PI;

    // There should be exactly p basis functions that have non zero trace
    auto p = GENERATE(1, 2, 3, 4);
    auto* bfSet = Utilities::createDGBasisFunctionSet2DNedelec(p);

    for (std::size_t side = 0; side < 3; ++side) {
        std::size_t nonZeroTangents = 0;
        LinearAlgebra::SmallVector<2> val;
        for (auto& bf : *bfSet) {
            // Intentional reuse of val to possibly catch when it is not zeroed
            bf->eval(sidePoints[side], val);
            if (std::abs(val * tangents[side]) > 1e-5) {
                nonZeroTangents++;
            }
        }
        REQUIRE(nonZeroTangents == p);
    }
}

TEST_CASE("Nedelec 2D: curl values", "[Nedelec2D]") {
    // Test whether the curl of the basis functions is computed correctly. For
    // this we compute the curl through Finite differences and compare it with
    // the value from the basisfunctions.

    auto p = GENERATE(1, 2, 3, 4);
    auto* bfSet = Utilities::createDGBasisFunctionSet2DNedelec(p);

    const double eps = 0.0001;
    const double margin = eps * 100;

    // Test points so that there is not just an accidental correct point
    Geometry::PointReference<2> testPoints[3];
    testPoints[0] = {0.3, 0.1};
    testPoints[1] = {0.9, 0.01};
    testPoints[2] = {0.4, 0.4};

    for (auto& testPoint : testPoints) {
        LinearAlgebra::SmallVector<2> curl;
        LinearAlgebra::SmallVector<2> val;
        LinearAlgebra::SmallVector<2> val2;
        double fydx, fxdy;
        for (std::size_t i = 0; i < bfSet->size(); ++i) {
            bfSet->eval(i, testPoint, val);
            curl = bfSet->evalCurl(i, testPoint);

            CHECK(curl[1] == 0.0);  // By contract

            // Use finite differences to compute the curl
            LinearAlgebra::SmallVector<2> temp({eps, 0});
            Geometry::PointReference<2> pdx(testPoint.getCoordinates() + temp);
            bfSet->eval(i, pdx, val2);
            fydx = (val2[1] - val[1]) / eps;

            temp[0] = 0.0;
            temp[1] = eps;
            Geometry::PointReference<2> pdy(testPoint.getCoordinates() + temp);
            bfSet->eval(i, pdy, val2);
            fxdy = (val2[0] - val[0]) / eps;
            double fdCurl = fydx - fxdy;
            // Actual check
            Approx target = Approx(fdCurl).margin(margin);
            INFO("Curl of bf " << i << " p=" << p);
            REQUIRE(curl[0] == target);
        }
    }
}

TEST_CASE("Nedelec 2D: Non singular mass matrix") {
    // Test to see if the mass matrix is non singular
    auto p = GENERATE(1, 2, 3, 4);
    auto* bfSet = Utilities::createDGBasisFunctionSet2DNedelec(p);

    auto* quadRule =
        QuadratureRules::AllGaussQuadratureRules::instance().getRule(
            &Geometry::ReferenceTriangle::Instance(), p);

    // Mass matrix
    LinearAlgebra::MiddleSizeMatrix mat(bfSet->size(), bfSet->size());
    // Manually perform the quadrature.
    for (std::size_t qPoint = 0; qPoint < quadRule->nrOfPoints(); ++qPoint) {
        Geometry::PointReference<2> qp = quadRule->getPoint(qPoint);
        double weight = quadRule->weight(qPoint);

        for (std::size_t i = 0; i < bfSet->size(); ++i) {
            LinearAlgebra::SmallVector<2> phiI;
            // Diagonal entry
            bfSet->eval(i, qp, phiI);
            mat(i, i) += phiI * phiI * weight;
            // Off diagonals
            for (std::size_t j = i + 1; j < bfSet->size(); ++j) {
                LinearAlgebra::SmallVector<2> phiJ;
                double val = phiI * phiJ * weight;
                mat(i, j) += val;
                mat(j, i) += val;
            }
        }
    }
    // Check for non singularity by doing a cholesky decomposition
    INFO("Checking for singular mass matrix for p=" << p)
    mat.cholesky();
}
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

#include "SavageHutter2DInteriorTest.h"
#include "HeightLimiters/PositiveLayerLimiter.h"

using namespace hpgem;


SavageHutter2DInteriorTest::SavageHutter2DInteriorTest(std::size_t polyOrder,
                                                       std::string meshName)
    : SavageHutter2DBase(3, polyOrder) {
    chuteAngle_ = M_PI / 180 * 29;
    epsilon_ = .1;
    // const PointPhysicalT &pPhys = createMeshDescription(1).bottomLeft_;
    const PointPhysicalT pPhys = {0., 0.};
    inflowBC_ = getInitialSolution(pPhys, 0);

    std::vector<std::string> variableNames = {"h", "hu", "hv"};
    setOutputNames("output2D", "SavageHutter", "SavageHutter", variableNames);

    readMesh(meshName);
}

LinearAlgebra::MiddleSizeVector SavageHutter2DInteriorTest::getInitialSolution(
    const PointPhysicalT &pPhys, const double &startTime,
    const std::size_t orderTimeDerivative) {
    return getExactSolution(pPhys, 0, orderTimeDerivative);
}

LinearAlgebra::MiddleSizeVector SavageHutter2DInteriorTest::getExactSolution(
    const PointPhysicalT &pPhys, const double &time,
    const std::size_t orderTimeDerivative) {
    const double x = pPhys[0];
    const double y = pPhys[1];
    const double h = 1 + .1 * std::sin(2 * M_PI * x) * std::sin(2 * M_PI * y) *
                             std::cos(2 * M_PI * time);
    const double u = 1;
    const double v = 1;
    return MiddleSizeVector({h, h * u, h * v});
}

LinearAlgebra::MiddleSizeVector SavageHutter2DInteriorTest::computeSourceTerm(
    const LinearAlgebra::MiddleSizeVector &numericalSolution,
    const PointPhysicalT &pPhys, const double time) {
    logger.assert_debug(chuteAngle_ < M_PI / 2,
                        "Angle must be in radians, not degrees!");
    const double x = pPhys[0];
    const double y = pPhys[1];
    // abbreviations: first letter for sine(s) or cosine (c), second letter for
    // direction
    const double sx = std::sin(2 * M_PI * x);
    const double sy = std::sin(2 * M_PI * y);
    const double st = std::sin(2 * M_PI * time);
    const double cx = std::cos(2 * M_PI * x);
    const double cy = std::cos(2 * M_PI * y);
    const double ct = std::cos(2 * M_PI * time);
    LinearAlgebra::MiddleSizeVector source(3);
    source[0] = -2 * M_PI / 10 * (sx * sy * st - cx * sy * ct - sx * cy * ct);
    source[1] = -2 * M_PI / 10 *
                (sx * sy * st - cx * sy * ct -
                 epsilon_ * std::cos(chuteAngle_) * (1 + .1 * sx * sy * ct) *
                     cx * sy * ct -
                 sx * cy * ct);
    source[2] = -2 * M_PI / 10 *
                (sx * sy * st - sx * cy * ct -
                 epsilon_ * std::cos(chuteAngle_) * (1 + .1 * sx * sy * ct) *
                     sx * cy * ct -
                 cx * sy * ct);
    return source;
}

LinearAlgebra::MiddleSizeVector SavageHutter2DInteriorTest::computePhysicalFlux(
    const LinearAlgebra::MiddleSizeVector &numericalSolution) {
    const double h = numericalSolution(0);
    logger.assert_debug(h > -1e-16, "Negative height (%)", h);
    double hu = numericalSolution(1);
    double hv = numericalSolution(2);
    double u = 0;
    double v = 0;
    if (h > dryLimit_) {
        u = hu / h;
        v = hv / h;
    }
    MiddleSizeVector flux(6);
    flux(0) = hu;
    flux(1) = hv;
    flux(2) = hu * u + epsilon_ / 2 * std::cos(chuteAngle_) * h * h;
    flux(3) = hu * v;
    flux(4) = hu * v;
    flux(5) = hv * v + epsilon_ / 2 * std::cos(chuteAngle_) * h * h;
    logger(DEBUG, "flux values: %, ", flux);
    return flux;
}

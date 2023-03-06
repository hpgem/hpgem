/*
This file forms part of hpGEM. This package has been developed over a number of
years by various people at the University of Twente and a full list of
contributors can be found at http://hpgem.org/about-the-code/team

This code is distributed using BSD 3-Clause License. A copy of which can found
below.


Copyright (c) 2023, University of Twente
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <iomanip>

#include <Base/CommandLineOptions.h>
#include <DGMaxLogger.h>
#include <DGMaxProgramUtils.h>
#include <Utils/KSpacePath.h>
#include <Utils/BandStructure.h>
#include <Utils/HomogeneousBandStructure.h>
#include <Utils/BraggStackBandstructure.h>

using namespace hpgem;

auto& dimensionArg = Base::register_argument<std::size_t>(
    'd', "dimension", "The dimension", false, 3);
auto& structureArg = Base::register_argument<std::size_t>(
    's', "structure", "The structure (0): Homogenous, (1) Bragg stack", true);
auto& backgroundEpsArg = Base::register_argument<double>(
    'b', "epsilon1", "Background epsilon", false, 1.0);
auto& foregroundEpsArg = Base::register_argument<double>(
    'f', "epsilon2", "Foreground epsilon (Bragg stack)", false, 12.1);

// Points in reciprocal space that are visted
// Usual format of (2@)[point1]:[point2]: .... :[pointn]
// Specified in units of 1/a = k/2pi
auto& pointsArg = Base::register_argument<std::string>(
    '\0', "points",
    "Compute a single point in or a set of lines through k-space", false);

auto& maxFreqArg =
    Base::register_argument<double>('\0', "freqMax", "Maximum frequency", true);

template <std::size_t DIM>
std::unique_ptr<BandStructure<DIM>> getStructure();

template <std::size_t DIM>
void runWithDimension() {
    using namespace DGMax;

    PointPath<DIM> path = DGMax::parsePath<DIM>(pointsArg.getValue());
    for (auto& point : path.points_) {
        // Convert from regular reduced units to actual k-vectors
        point *= (2 * M_PI);
    }
    if (path.steps_ < 0) {
        path.steps_ = 1;
    }
    KSpacePath<DIM> kpath(path.points_, path.steps_);

    auto structure = getStructure<DIM>();
    double omegaMax = maxFreqArg.getValue();

    for (std::size_t kindex = 0; kindex < kpath.totalNumberOfSteps();
         ++kindex) {
        auto kpoint = kpath.k(kindex);
        auto freqs = structure->computeSpectrum(kpoint, omegaMax * (2 * M_PI));
        std::cout << (kindex + 1);
        for (std::size_t i = 0; i < DIM; ++i) {
            std::cout << "," << kpoint[i] / (2 * M_PI);
        }
        std::cout << kpoint.l2Norm() / (2 * M_PI);
        for (auto& freq : freqs) {
            for (std::size_t i = 0; i < freq.second; ++i) {
                std::cout << "," << std::setprecision(8)
                          << (freq.first / (2 * M_PI));
            }
        }
        std::cout << std::endl;
    }
}

int main(int argc, char** argv) {
    Base::parse_options(argc, argv);
    initDGMaxLogging();
    DGMax::printArguments(argc, argv);
    switch (dimensionArg.getValue()) {
        case 2:
            runWithDimension<2>();
            break;
        case 3:
            runWithDimension<3>();
            break;
        default:
            logger(ERROR, "Dimension % is not supported",
                   dimensionArg.getValue());
    }
}

template <>
std::unique_ptr<BandStructure<3>> getStructure() {
    std::array<LinearAlgebra::SmallVector<3>, 3> reciprocals;
    for (std::size_t i = 0; i < 3; ++i) {
        reciprocals[i][i] = 2.0 * M_PI;
    }
    switch (structureArg.getValue()) {
        case 0:
            return std::make_unique<HomogeneousBandStructure<3>>(
                reciprocals, backgroundEpsArg.getValue());
        case 1:
            return std::make_unique<BraggStackBandstructure<3>>(
                backgroundEpsArg.getValue(), foregroundEpsArg.getValue());
        default:
            logger.fail("Unknown structure %", structureArg.getValue());
    }
}

template <>
std::unique_ptr<BandStructure<2>> getStructure() {
    std::array<LinearAlgebra::SmallVector<2>, 2> reciprocals;
    for (std::size_t i = 0; i < 2; ++i) {
        reciprocals[i][i] = 1.0;
    }
    switch (structureArg.getValue()) {
        case 0:
            return std::make_unique<HomogeneousBandStructure<2>>(
                reciprocals, backgroundEpsArg.getValue());
        case 1:
            logger(ERROR, "Braggstack not supported in 2D");
        default:
            logger.fail("Unknown structure %", structureArg.getValue());
    }
}
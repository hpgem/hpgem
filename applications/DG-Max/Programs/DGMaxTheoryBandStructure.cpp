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
#include <Utils/BandstructureGNUPlot.h>
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
auto& modeTypes = Base::register_argument<std::string>(
    '\0', "modes", "Compute \"TE\", \"TM\" or \"both\" (Braggstack only)",
    false, "both");

auto& maxModesArgs = Base::register_argument<std::size_t>(
    '\0', "maxmodes", "Maximum number of modes printed", false);

auto& transverseLattice =Base::register_argument<double>(
    '\0', "braggtransversesize", "Relative size of the transverse lattice ", false, 1.0);

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

    // Header
    std::cout << "k-index,kx,ky,kz,kmag/2pi";
    if (maxModesArgs.isUsed()) {
        for (std::size_t i = 0; i < maxModesArgs.getValue(); ++i) {
            std::cout << ",mode" << (i+1);
        }
    }
    std::cout << std::endl;

    // Compute modes
    for (std::size_t kindex = 0; kindex < kpath.totalNumberOfSteps();
         ++kindex) {
        auto kpoint = kpath.k(kindex);
        auto freqs =
            structure->computeLinearSpectrum(kpoint, omegaMax * (2 * M_PI));
        std::cout << (kindex + 1);
        for (std::size_t i = 0; i < 3; ++i) {
            if (i < DIM) {
                std::cout << "," << kpoint[i] / (2 * M_PI);
            } else {
                std::cout << ",0";
            }
        }
        std::cout << "," << kpoint.l2Norm() / (2 * M_PI);
        std::size_t modeCount = 0;
        for (auto& freq : freqs) {
            modeCount++;
            std::cout << "," << std::setprecision(8) << (freq / (2 * M_PI));
            if (maxModesArgs.isUsed() && maxModesArgs.getValue() == modeCount) {
                break;
            }
        }
        std::cout << std::endl;
    }

    BandstructureGNUPlot<DIM> plotter (kpath, {}, *structure);
    plotter.plot("theory.p");
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

template <std::size_t DIM>
std::unique_ptr<HomogeneousBandStructure<DIM>> getHomogenousStructure() {
    std::array<LinearAlgebra::SmallVector<DIM>, DIM> reciprocals;
    for (std::size_t i = 0; i < DIM; ++i) {
        reciprocals[i][i] = 2.0 * M_PI;
    }
    if (modeTypes.hasArgument()) {
        DGMaxLogger(WARN, "Mode types argument is ignored");
    }
    if (transverseLattice.hasArgument()) {
        DGMaxLogger(WARN, "Argument % unused", transverseLattice.getLongTag());
    }
    return std::make_unique<HomogeneousBandStructure<DIM>>(
        reciprocals, backgroundEpsArg.getValue());
}

template <std::size_t DIM>
std::unique_ptr<BraggStackBandstructure<DIM>> getBraggStack() {
    auto structure = std::make_unique<BraggStackBandstructure<DIM>>(
        backgroundEpsArg.getValue(), foregroundEpsArg.getValue(),
        0.5, transverseLattice.getValue());
    if (modeTypes.getValue() == "TE") {
        structure->setComputeTM(false);
    } else if (modeTypes.getValue() == "TM") {
        structure->setComputeTE(false);
    } else if (modeTypes.getValue() != "both") {
        logger.fail("Unknown mode type \"%\"", modeTypes.getValue());
    }
    return structure;
}

template <std::size_t DIM>
std::unique_ptr<BandStructure<DIM>> getStructure() {
    switch (structureArg.getValue()) {
        case 0:
            return getHomogenousStructure<DIM>();
        case 1:
            return getBraggStack<DIM>();
        default:
            logger.fail("Unknown structure %", structureArg.getValue());
    }
}
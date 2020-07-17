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

// Simple program that can compute and output the two analytic spectra

#include "Base/CommandLineOptions.h"

#include "DGMaxLogger.h"
#include "DGMaxProgramUtils.h"
#include "Utils/BandStructure.h"
#include "Utils/BraggStackBandstructure.h"
#include "Utils/HomogeneousBandStructure.h"
#include "Utils/KSpacePath.h"

auto& structureArg = Base::register_argument<std::size_t>(
    '\0', "structure",
    "The structure to use for the spectrum vacuum (0) or Bragg stack (1)",
    true);

auto& points = Base::register_argument<std::string>(
    '\0', "points", "The points in k-space for which to output the frequencies",
    true);

auto& maxFrequency = Base::register_argument<double>(
    '\0', "maxFrequency", "The maximum frequency to output", true);

auto& dimensionArg = Base::register_argument<std::size_t>(
    'd', "dimension", "The dimension to consider", true);

template <std::size_t DIM>
void runWithDimension();

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
            logger.assert_always(
                false, "Only supports 2 and 3 dimensional structures");
    }
}

template <std::size_t DIM>
std::unique_ptr<BandStructure<DIM>> createBraggStack() {
    logger.assert_always(false,
                         "Bragg stack is not implemented for dimension %", DIM);
    // Return is never reached
    return nullptr;
}

std::unique_ptr<BandStructure<3>> createBraggStack() {
    return std::make_unique<BraggStackBandstructure>(13, 1);
}

template <std::size_t DIM>
void runWithDimension() {
    // Construct reciprocal lattice for a unit square
    std::array<LinearAlgebra::SmallVector<DIM>, DIM> reciprocalVectors;
    for (std::size_t i = 0; i < DIM; ++i) {
        reciprocalVectors[i].set(0);
        reciprocalVectors[i][i] = 2.0 * M_PI;
    }
    std::unique_ptr<BandStructure<DIM>> structure;
    switch (structureArg.getValue()) {
        case 0:
            structure = std::make_unique<HomogeneousBandStructure<DIM>>(
                reciprocalVectors);
            break;
        case 1:
            structure = createBraggStack<DIM>();
            break;
        default:
            logger.assert_always(false, "No analytic spectrum for structure %",
                                 structureArg.getValue());
    }
    DGMax::PointPath<DIM> path = DGMax::parsePath<DIM>(points.getValue());

    std::vector<typename KSpacePath<DIM>::KPoint> rawKPoints = path.points_;
    // Rescale to go from relative points to absolute ones.
    // Note: Works only with unit cube lattice.
    for (auto& point : rawKPoints) {
        point *= M_PI;
    }
    KSpacePath<DIM> kSpacePath(rawKPoints, path.steps_);
    for (std::size_t i = 0; i < kSpacePath.totalNumberOfSteps(); ++i) {
        std::cout << i;
        std::vector<double> spectrum = structure->computeLinearSpectrum(
            kSpacePath.k(i), maxFrequency.getValue());
        for (double& frequency : spectrum) {
            std::cout << "\t" << frequency;
        }
        std::cout << std::endl;
    }
}
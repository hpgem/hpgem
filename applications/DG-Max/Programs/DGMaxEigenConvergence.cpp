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

// Program to run convergence tests

#include <iomanip>
#include <utility>
#include "Base/CommandLineOptions.h"

#include "DGMaxLogger.h"
#include "DGMaxProgramUtils.h"
#include "Algorithms/DivDGMaxEigenvalue.h"
#include "Algorithms/DGMaxEigenvalue.h"
#include "Utils/HomogeneousBandStructure.h"
#include "Utils/Verification/AbstractEVConvergenceTest.h"
#include "Utils/Verification/DGMaxEVConvergenceTest.h"
#include "Utils/Verification/DivDGMaxEVConvergenceTest.h"
#include "Utils/Verification/EVConvergenceResult.h"

auto &dimensionArg = Base::register_argument<std::size_t>(
    'd', "dimension", "The dimension of the problem", true);

auto &structureArg = Base::register_argument<std::size_t>(
    '\0', "structure", "The structure to use", true, 0);

auto &meshFiles = Base::register_argument<std::string>(
    'm', "meshes", "The mesh files to use, comma separated", true);

auto &method = Base::register_argument<std::string>(
    '\0', "method",
    "The method to be used, either 'DGMAX' or 'DIVDGMAX' (default)", false,
    "DIVDGMAX");

/// Create a reference bandstructure based on the structure index
template <std::size_t DIM>
std::unique_ptr<BandStructure<DIM>> createStructure(
    std::size_t structureIndex) {
    // Using same structure indices as for the Jelmer structure
    if (structureIndex == 0) {
        std::array<LinearAlgebra::SmallVector<DIM>, DIM> reciprocalVectors;
        for (std::size_t i = 0; i < DIM; ++i) {
            for (std::size_t j = 0; j < DIM; ++j) {
                reciprocalVectors[i][j] = 0;
            }
            // Assuming unit cube as unit cell
            reciprocalVectors[i][i] = 2 * M_PI;
        }
        return std::make_unique<HomogeneousBandStructure<DIM>>(
            reciprocalVectors);
    } else {
        return std::unique_ptr<BandStructure<DIM>>();
    }
}

/// Parse the input argument with mesh files
/// \return The list of filenames with meshes
std::vector<std::string> parseMeshFiles();

/// Run the actual program with fixed dimension template.
template <std::size_t DIM>
void runWithDimension();

int main(int argc, char **argv) {
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
            logger.assert_always(false, "Can only run with dimension 2 or 3");
    }
    return 0;
}

template <std::size_t DIM>
void runWithDimension() {
    const std::size_t numFrequencies = 10;

    std::vector<std::string> meshFiles = parseMeshFiles();
    std::size_t structureId = structureArg.getValue();

    std::unique_ptr<DGMax::AbstractEVConvergenceTest<DIM>> convergenceTest;
    DGMax::EVTestPoint<DIM> testPoint({0.5, 0.8}, structureId, numFrequencies);
    if (method.getValue() == "DGMAX") {
        DGMaxEigenvalueBase::SolverConfig config;
        config.useHermitian_ = true;
        config.stab_ = 100;
        config.shiftFactor_ = 0.0;
        config.useProjector_ = false;
        convergenceTest = std::make_unique<DGMax::DGMaxEVConvergenceTest<DIM>>(
            testPoint, meshFiles,
            0.0,  // No expecations
            1, config, nullptr);
    } else if (method.getValue() == "DIVDGMAX") {
        // Some default stabilization parameters
        typename DivDGMaxDiscretization<DIM>::Stab stab;
        stab.stab1 = 5;
        stab.stab2 = 0;
        stab.stab3 = 5;
        stab.setAllFluxeTypes(DivDGMaxDiscretization<DIM>::FluxType::BREZZI);

        convergenceTest =
            std::make_unique<DGMax::DivDGMaxEVConvergenceTest<DIM>>(
                testPoint, meshFiles,
                0.0,  // No expectations
                1, stab, nullptr);
    } else {
        DGMaxLogger(ERROR, "Unknown method %", method.getValue());
        return;
    }

    DGMax::EVConvergenceResult refinementResult = convergenceTest->run(false);

    // Compute Errors //
    ////////////////////
    std::vector<double> spectrum;

    std::unique_ptr<BandStructure<DIM>> structure =
        createStructure<DIM>(testPoint.getStructureId());
    if (structure) {
        // TODO: It would be nice to have a way to compute the N lowest instead
        // of up to a frequency
        spectrum = structure->computeLinearSpectrum(testPoint.getKPoint(), 20);
        logger.assert_always(spectrum.size() >= numFrequencies,
                             "Not enough analytical frequencies");
    }

    // Print results //
    ///////////////////

    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    if (rank == 0) {
        refinementResult.filterResults(0.01, false);
        refinementResult.printFrequencyTable(spectrum);
        if (!spectrum.empty()) {
            refinementResult.printErrorTable(spectrum);
        }
    }
}

std::vector<std::string> parseMeshFiles() {
    std::string files = meshFiles.getValue();

    logger.assert_always(!files.empty(), "Empty list of meshes");

    std::vector<std::string> result;

    std::size_t pos = 0;
    while (true) {
        std::size_t next = files.find_first_of(',', pos);
        if (next == pos) {
            logger(WARN, "Empty mesh string at position %", pos);
            pos++;
        }
        if (next == std::string::npos) {
            // Last entry
            result.push_back(files.substr(pos));
            break;
        } else {
            result.push_back(files.substr(pos, next - pos));
            pos = next + 1;  // Skip the comma itself.
        }
    }
    return result;
}
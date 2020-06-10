

#include <iomanip>
#include <utility>
#include "Base/CommandLineOptions.h"

#include "DGMaxLogger.h"
#include "DGMaxProgramUtils.h"
#include "Algorithms/DivDGMaxEigenvalue.h"
#include "Algorithms/DGMaxEigenvalue.h"
#include "Utils/HomogeneousBandStructure.h"
#include "Utils/Verification/RunnableEVTestCase.h"
#include "Utils/Verification/DGMaxEVTestCase.h"
#include "Utils/Verification/DivDGMaxEVTestCase.h"
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
        return std::unique_ptr<BandStructure<DIM>>(
            new HomogeneousBandStructure<DIM>(reciprocalVectors));
    } else {
        return std::unique_ptr<BandStructure<DIM>>();
    }
}

std::vector<std::string> parseMeshFiles();

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

    std::unique_ptr<DGMax::RunnableEVTestCase<DIM>> testCase2;
    DGMax::EVTestCase<DIM> rawTestCase({0.5, 0.8}, structureId, numFrequencies);
    if (method.getValue() == "DGMAX") {
        testCase2 = std::unique_ptr<DGMax::RunnableEVTestCase<DIM>>(
            new DGMax::DGMaxEVTestCase<DIM>(rawTestCase, meshFiles,
                                            0.0,  // No expecations
                                            1, 100, nullptr));
    } else if (method.getValue() == "DIVDGMAX") {
        typename DivDGMaxDiscretization<DIM>::Stab stab;
        stab.stab1 = 5;
        stab.stab2 = 0;
        stab.stab3 = 5;
        stab.setAllFluxeTypes(DivDGMaxDiscretization<DIM>::FluxType::BREZZI);

        testCase2 = std::unique_ptr<DGMax::RunnableEVTestCase<DIM>>(
            new DGMax::DivDGMaxEVTestCase<DIM>(rawTestCase, meshFiles,
                                               0.0,  // No expectations
                                               1, stab, nullptr));
    } else {
        DGMaxLogger(ERROR, "Unknown method %", method.getValue());
        return;
    }

    DGMax::EVConvergenceResult refinementResult = testCase2->run(false);

    // Compute Errors //
    ////////////////////
    std::vector<double> spectrum;

    std::unique_ptr<BandStructure<DIM>> structure =
        createStructure<DIM>(rawTestCase.getStructureId());
    if (structure) {
        // TODO: It would be nice to have a way to compute the N lowest instead
        // of up to a frequency
        spectrum =
            structure->computeLinearSpectrum(rawTestCase.getKPoint(), 20);
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
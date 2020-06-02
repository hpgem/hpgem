

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
        logger.assert_always(false, "Unknown structure");
        // Won't happen.
        return std::unique_ptr<BandStructure<DIM>>();
    }
}

const std::vector<std::string> meshes2D({"meshes/mesh_D2X10N1.hpgem",
                                         "meshes/mesh_D2X20N1.hpgem",
                                         "meshes/mesh_D2X40N1.hpgem"});

int main(int argc, char **argv) {
    Base::parse_options(argc, argv);
    initDGMaxLogging();
    DGMax::printArguments(argc, argv);

    const std::size_t numFrequencies = 10;

    std::unique_ptr<DGMax::RunnableEVTestCase<2>> testCase2;
    DGMax::EVTestCase<2> rawTestCase({0.5, 0.8}, 0, numFrequencies);
    if (method.getValue() == "DGMAX") {
        testCase2 = std::unique_ptr<DGMax::RunnableEVTestCase<2>>(
            new DGMax::DGMaxEVTestCase<2>(rawTestCase, meshes2D,
                                        0.0,  // No expecations
                                        1, 100, nullptr));
    } else if (method.getValue() == "DIVDGMAX") {
        typename DivDGMaxDiscretization<2>::Stab stab;
        stab.stab1 = 5;
        stab.stab2 = 0;
        stab.stab3 = 5;
        stab.setAllFluxeTypes(DivDGMaxDiscretization<2>::FluxType::BREZZI);

        testCase2 = std::unique_ptr<DGMax::RunnableEVTestCase<2>>(
            new DGMax::DivDGMaxEVTestCase<2>(rawTestCase, meshes2D,
                                           0.0,  // No expectations
                                           1, stab, nullptr));
    } else {
        DGMaxLogger(ERROR, "Unknown method %", method.getValue());
        return -1;
    }

    std::unique_ptr<BandStructure<2>> structure =
        createStructure<2>(rawTestCase.getStructureId());

    DGMax::EVConvergenceResult refinementResult = testCase2->runWithResults(false);

    // Compute Errors //
    ////////////////////

    std::vector<std::vector<double>> errors;

    // Probably overkill
    std::map<double, std::size_t> spectrum =
        structure->computeSpectrum(rawTestCase.getKPoint(), 20);
    // Create a linear spectrum for easy comparison.
    std::vector<double> linearSpectrum;
    for (auto const &entry : spectrum) {
        for (std::size_t i = 0; i < entry.second; ++i) {
            linearSpectrum.emplace_back(entry.first);
        }
    }
    // TODO: Add option for the analytical spectrum to compute at least N
    // frequencies
    logger.assert_always(linearSpectrum.size() >= numFrequencies,
                         "Not enough analytical frequencies");

    // Print results //
    ///////////////////

    refinementResult.filterResults(0.01, false);
    refinementResult.printFrequencyTable(linearSpectrum);
    refinementResult.printErrorTable(linearSpectrum);

    return 0;
}
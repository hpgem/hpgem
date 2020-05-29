

#include <iomanip>
#include "Base/CommandLineOptions.h"

#include "DGMaxLogger.h"
#include "DGMaxProgramUtils.h"
#include "Algorithms/DivDGMaxEigenValue.h"
#include "Algorithms/DGMaxEigenValue.h"
#include "Utils/HomogeneousBandStructure.h"

auto &method = Base::register_argument<std::string>(
    '\0', "method",
    "The method to be used, either 'DGMAX' or 'DIVDGMAX' (default)", false,
    "DIVDGMAX");
auto &order = Base::register_argument<std::size_t>(
    'p', "order", "The polynomial order to be used, defaults to p=1", false, 1);

auto &numEigenvalues = Base::register_argument<std::size_t>(
    'e', "eigenvalues", "The number of eigenvalues to compute", false, 10);

template <std::size_t DIM>
class Solver {
   public:
    Solver(std::size_t order) : order_(order){};
    virtual ~Solver() = default;
    virtual std::unique_ptr<BaseEigenvalueResult<DIM>> solve(
        const std::string &meshFileName, std::size_t structureIndex,
        const LinearAlgebra::SmallVector<DIM> &kpoint,
        std::size_t numEigenvalues) = 0;

   protected:
    // Order of the solver
    std::size_t order_;
    // K-vector
    // structure
    // dim?
};

template <std::size_t DIM>
class DivDGMaxSolver : public Solver<DIM> {
   public:
    DivDGMaxSolver(std::size_t order) : Solver<DIM>(order){};

    std::unique_ptr<BaseEigenvalueResult<DIM>> solve(
        const std::string &meshFileName, std::size_t structureIndex,
        const LinearAlgebra::SmallVector<DIM> &kpoint,
        std::size_t numEigenvalues) override;
};

template <std::size_t DIM>
class DGMaxSolver : public Solver<DIM> {
   public:
    DGMaxSolver(std::size_t order) : Solver<DIM>(order){};

    std::unique_ptr<BaseEigenvalueResult<DIM>> solve(
        const std::string &meshFileName, std::size_t structureIndex,
        const LinearAlgebra::SmallVector<DIM> &kpoint,
        std::size_t numEigenvalues) override;
};

/// Result of the convergence test
struct Result {
    /// The frequencies of the refinement levels
    std::vector<std::vector<double>> frequencies_;

    /// Frequencies computed in analytical sense.
    std::vector<double> theoreticalFrequencies_;

    std::size_t maxNumberOfFrequencies() {
        std::size_t max = 0;
        for (const auto &level : frequencies_) {
            max = std::max(max, level.size());
        }
        return max;
    }

    void printFrequencyTable(std::vector<double> theoretical);
    void printErrorTable(std::vector<double> theoretical);
};

void Result::printFrequencyTable(std::vector<double> theoretical) {
    std::size_t max = maxNumberOfFrequencies();

    // Header line
    std::cout << "|";
    for (std::size_t i = 0; i < max; ++i) {
        if (i < theoretical.size()) {
            std::cout << " " << std::setprecision(4) << std::setw(6)
                      << theoretical[i] << " |";
        } else {
            // When lacking an theoretical frequency, just print a placeholder
            std::cout << " f_ " << std::setw(4) << i << " |";
        }
    }
    std::cout << std::endl << "|";  // For the header bottom line
    for (std::size_t i = 0; i < max; ++i) {
        std::cout << "--------|";
    }
    std::cout << std::endl;
    // Frequencies
    for (auto &levelFrequency : frequencies_) {
        std::cout << "|";
        for (double &frequency : levelFrequency) {
            std::cout << " " << std::setprecision(4) << std::setw(6)
                      << frequency << " |";
        }
        std::cout << std::endl;
    }
    // Spacing newline
    std::cout << std::endl;
}

void Result::printErrorTable(std::vector<double> theoretical) {

    std::size_t max = maxNumberOfFrequencies();

    std::vector<std::vector<double>> errors;
    // Compute actual errors
    for (auto &levelFrequencies : frequencies_) {
        std::vector<double> temp;
        for (std::size_t i = 0;
             i < levelFrequencies.size() && i < theoretical.size(); ++i) {
            temp.emplace_back(std::abs(levelFrequencies[i] - theoretical[i]));
        }
        errors.emplace_back(temp);
    }

    // Printing
    // Header line
    std::cout << "|";
    for (std::size_t i = 0; i < max; ++i) {
        std::cout << " ";  // Separator from previous column
        if (i < theoretical.size()) {
            std::cout << std::setprecision(4) << std::setw(8) << theoretical[i];
        } else {
            // Place holder
            std::cout << "f_" << std::setw(8) << std::right << i;
        }
        std::cout << "        "  // Space of second column
                  << " |";
    }
    std::cout << std::endl << "|";  // For the header bottom line
    for (std::size_t i = 0; i < max; ++i) {
        std::cout << "----------|-------|";
    }
    std::cout << std::endl;
    // Frequencies
    for (std::size_t i = 0; i < errors.size(); ++i) {
        auto &meshErrors = errors[i];
        std::cout << "|";
        for (std::size_t j = 0; j < meshErrors.size(); ++j) {
            double &error = meshErrors[j];
            std::cout << " " << std::scientific << std::setprecision(2)
                      << std::setw(6) << error << " |";
            if (i == 0 || errors[i - 1].size() <= j) {
                // For first row there is no convergence speed
                // so print a placeholder instead.
                std::cout << "     - |";
            } else {
                double factor = errors[i - 1][j] / error;
                std::cout << " " << std::fixed << std::setprecision(2)
                          << std::setw(5) << factor << " |";
            }
        }
        std::cout << std::endl;
    }
    // Spacing newline
    std::cout << std::endl;
}

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

int main(int argc, char **argv) {
    Base::parse_options(argc, argv);
    initDGMaxLogging();
    DGMax::printArguments(argc, argv);

    std::string meshes[] = {"meshes/mesh_D2X10N1.hpgem",
                            "meshes/mesh_D2X20N1.hpgem",
                            "meshes/mesh_D2X40N1.hpgem"};

    const std::size_t structureIndex = 0;
    std::unique_ptr<Solver<2>> solver;
    if (method.getValue() == "DGMAX") {
        solver =
            std::unique_ptr<Solver<2>>(new DGMaxSolver<2>(order.getValue()));
    } else if (method.getValue() == "DIVDGMAX") {
        solver =
            std::unique_ptr<Solver<2>>(new DivDGMaxSolver<2>(order.getValue()));
    } else {
        DGMaxLogger(ERROR, "Unknown method %", method.getValue());
        return -1;
    }

    std::unique_ptr<BandStructure<2>> structure =
        createStructure<2>(structureIndex);
    const LinearAlgebra::SmallVector<2> kpoint(
        {0.5, 0.5});  // Should be multiplied by M_PI?

    Result refinementResult;
    std::size_t numFrequencies = numEigenvalues.getValue();

    for (std::string &meshFileName : meshes) {
        auto result =
            solver->solve(meshFileName, structureIndex, kpoint, numFrequencies);
        // TODO: Raw eigenvalues is better than the computed frequencies in case
        // something goes wrong

        const std::vector<double> &resultFreqs = result->frequencies(0);
        // Discard zero/negative eigenvalues
        std::size_t index = 0;
        while (index < resultFreqs.size() &&
               (std::isnan(resultFreqs[index]) ||
                std::abs(resultFreqs[index]) < 1e-3)) {
            index++;
        }

        // Copy the required number of frequencies
        std::vector<double> freqs;
        for (; index < numFrequencies && index < resultFreqs.size(); ++index) {
            freqs.emplace_back(resultFreqs[index]);
        }
        refinementResult.frequencies_.emplace_back(freqs);
    }

    // Compute Errors //
    ////////////////////

    std::vector<std::vector<double>> errors;

    // Probably overkill
    std::map<double, std::size_t> spectrum =
        structure->computeSpectrum(kpoint, 20);
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

    for (auto &meshFrequencies : refinementResult.frequencies_) {
        std::vector<double> tempErrors(meshFrequencies.size(), 0);
        for (std::size_t j = 0; j < meshFrequencies.size(); ++j) {
            tempErrors[j] = std::abs(meshFrequencies[j] - linearSpectrum[j]);
        }
        errors.emplace_back(tempErrors);
    }

    // Print results //
    ///////////////////

    refinementResult.printFrequencyTable(linearSpectrum);
    refinementResult.printErrorTable(linearSpectrum);
    return 0;
}

template <std::size_t DIM>
std::unique_ptr<BaseEigenvalueResult<DIM>> DivDGMaxSolver<DIM>::solve(
    const std::string &meshFileName, std::size_t structureIndex,
    const LinearAlgebra::SmallVector<DIM> &kpoint, std::size_t numEigenvalues) {
    Base::ConfigurationData configData(2, 1);
    auto mesh = DGMax::readMesh<DIM>(
        meshFileName, &configData, [&](const Geometry::PointPhysical<DIM> &p) {
            return jelmerStructure(p, structureIndex);
        });
    DGMaxLogger(INFO, "Loaded mesh % with % local elements.", meshFileName,
                mesh->getNumberOfElements());
    KSpacePath<DIM> path = KSpacePath<DIM>::singleStepPath(kpoint);
    EigenValueProblem<DIM> input(path, numEigenvalues);

    DivDGMaxEigenValue<DIM> solver(*mesh);
    typename DivDGMaxDiscretization<DIM>::Stab stab;  // Change this?
    stab.stab1 = 100;
    stab.stab2 = 0;
    stab.stab3 = 1;
    typename DivDGMaxEigenValue<DIM>::Result result =
        solver.solve(input, stab, this->order_);

    return std::unique_ptr<BaseEigenvalueResult<DIM>>(
        new typename DivDGMaxEigenValue<DIM>::Result(result));
}

template <std::size_t DIM>
std::unique_ptr<BaseEigenvalueResult<DIM>> DGMaxSolver<DIM>::solve(
    const std::string &meshFileName, std::size_t structureIndex,
    const LinearAlgebra::SmallVector<DIM> &kpoint, std::size_t numEigenvalues) {
    Base::ConfigurationData configData(1, 1);
    auto mesh = DGMax::readMesh<DIM>(
        meshFileName, &configData, [&](const Geometry::PointPhysical<DIM> &p) {
            return jelmerStructure(p, structureIndex);
        });
    DGMaxLogger(INFO, "Loaded mesh % with % local elements.", meshFileName,
                mesh->getNumberOfElements());
    KSpacePath<DIM> path = KSpacePath<DIM>::singleStepPath(kpoint);
    EigenValueProblem<DIM> input(path, numEigenvalues);

    // TODO Vary the order
    DGMaxEigenValue<DIM> solver(*mesh, this->order_);
    auto result = solver.solve(input, 100);
    return std::unique_ptr<BaseEigenvalueResult<DIM>>(
        new typename DGMaxEigenValue<DIM>::Result(result));
}
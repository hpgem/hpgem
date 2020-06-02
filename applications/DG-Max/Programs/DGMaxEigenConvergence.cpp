

#include <iomanip>
#include <utility>
#include "Base/CommandLineOptions.h"

#include "DGMaxLogger.h"
#include "DGMaxProgramUtils.h"
#include "Algorithms/DivDGMaxEigenValue.h"
#include "Algorithms/DGMaxEigenValue.h"
#include "Utils/HomogeneousBandStructure.h"
#include "Utils/Verification/EigenvalueResult.h"

auto &method = Base::register_argument<std::string>(
    '\0', "method",
    "The method to be used, either 'DGMAX' or 'DIVDGMAX' (default)", false,
    "DIVDGMAX");

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

template <std::size_t DIM>
struct DGMaxTestCase {
    DGMaxTestCase(LinearAlgebra::SmallVector<DIM> kpoint,
                  std::vector<std::string> meshes, size_t structureId,
                  size_t numberOfEigenvalues, size_t order,
                  std::shared_ptr<DGMax::EigenvalueResult> expected = nullptr)
        : kpoint_(std::move(kpoint)),
          meshes_(std::move(meshes)),
          structureId_(structureId),
          numberOfEigenvalues_(numberOfEigenvalues),
          order_(order),
          expected_(std::move(expected)) {}

    LinearAlgebra::SmallVector<DIM> kpoint_;
    std::vector<std::string> meshes_;
    std::size_t structureId_;
    std::size_t numberOfEigenvalues_;
    std::size_t order_;
    std::shared_ptr<DGMax::EigenvalueResult> expected_;
};

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

static DGMaxTestCase<2> simpleVacuum({0.5, 0.8},
                                     {"meshes/mesh_D2X10N1.hpgem",
                                      "meshes/mesh_D2X20N1.hpgem",
                                      "meshes/mesh_D2X40N1.hpgem"},
                                     0, 10, 1);

int main(int argc, char **argv) {
    Base::parse_options(argc, argv);
    initDGMaxLogging();
    DGMax::printArguments(argc, argv);
    ;

    DGMaxTestCase<2> testCase = simpleVacuum;

    std::unique_ptr<Solver<2>> solver;
    if (method.getValue() == "DGMAX") {
        solver =
            std::unique_ptr<Solver<2>>(new DGMaxSolver<2>(testCase.order_));
    } else if (method.getValue() == "DIVDGMAX") {
        solver =
            std::unique_ptr<Solver<2>>(new DivDGMaxSolver<2>(testCase.order_));
    } else {
        DGMaxLogger(ERROR, "Unknown method %", method.getValue());
        return -1;
    }

    std::unique_ptr<BandStructure<2>> structure =
        createStructure<2>(testCase.structureId_);

    DGMax::EigenvalueResult refinementResult;
    std::size_t numFrequencies = testCase.numberOfEigenvalues_;

    for (std::string &meshFileName : testCase.meshes_) {
        auto result = solver->solve(meshFileName, testCase.structureId_,
                                    testCase.kpoint_, numFrequencies);
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
        std::vector<double> frequencies;
        for (; index < numFrequencies && index < resultFreqs.size(); ++index) {
            frequencies.emplace_back(resultFreqs[index]);
        }
        refinementResult.addLevel(frequencies);
    }

    // Compute Errors //
    ////////////////////

    std::vector<std::vector<double>> errors;

    // Probably overkill
    std::map<double, std::size_t> spectrum =
        structure->computeSpectrum(testCase.kpoint_, 20);
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

    if (testCase.expected_) {
        // There is a case to compare
        if (!refinementResult.equals(*testCase.expected_, 1e-10)) {
            std::cout << "Expected results" << std::endl;
            testCase.expected_->printFrequencyTable({});
            refinementResult.printFrequencyTable({});
            logger.assert_always(false, "Differences with expected output");
        }
    } else {

        refinementResult.printFrequencyTable(linearSpectrum);
        refinementResult.printErrorTable(linearSpectrum);
    }

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
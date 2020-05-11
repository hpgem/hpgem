

#include <iomanip>
#include "Base/CommandLineOptions.h"

#include "DGMaxLogger.h"
#include "DGMaxProgramUtils.h"
#include "Algorithms/DivDGMaxEigenValue.h"
#include "Algorithms/DGMaxEigenValue.h"
#include "Utils/HomogeneousBandStructure.h"

class Solver
{
public:
    virtual ~Solver() = default;
    virtual std::unique_ptr<BaseEigenvalueResult<2>> solve(
            const std::string &meshFileName, std::size_t structureIndex,
            const LinearAlgebra::SmallVector<2>& kpoint) = 0;

private:
    std::size_t p;
    // K-vector
    // structure
    // dim?
};

class DivDGMaxSolver : public Solver
{
    std::unique_ptr<BaseEigenvalueResult<2>> solve(
            const std::string &meshFileName, std::size_t structureIndex,
            const LinearAlgebra::SmallVector<2>& kpoint) override;
};

template<std::size_t DIM>
std::unique_ptr<BandStructure<DIM>> createStructure(std::size_t structureIndex)
{
    // Using same structure indices as for the Jelmer structure
    if (structureIndex == 0)
    {
        std::array<LinearAlgebra::SmallVector<DIM>, DIM> reciprocalVectors;
        for (std::size_t i = 0; i < DIM; ++i)
        {
            for (std::size_t j = 0; j < DIM; ++j)
            {
                reciprocalVectors[i][j] = 0;
            }
            // Assuming unit cube as unit cell
            reciprocalVectors[i][i] = 2*M_PI;
        }
        return std::unique_ptr<BandStructure<DIM>>(new HomogeneousBandStructure<DIM>(reciprocalVectors));
    }
    else
    {
        logger.assert_always(false, "Unknown structure");
        // Won't happen.
        return std::unique_ptr<BandStructure<DIM>>();
    }
}

int main(int argc, char** argv)
{
    Base::parse_options(argc, argv);
    initDGMaxLogging();
    DGMax::printArguments(argc, argv);

    std::string meshes[] = {
            "meshes/mesh_D2X10N1.hpgem",
            "meshes/mesh_D2X20N1.hpgem",
            "meshes/mesh_D2X40N1.hpgem"
    };

    const std::size_t structureIndex = 0;
    std::unique_ptr<Solver> solver(new DivDGMaxSolver());
    std::unique_ptr<BandStructure<2>> structure = createStructure<2>(structureIndex);
    const LinearAlgebra::SmallVector<2> kpoint ({0.5, 0.5});// Should be multiplied by M_PI?

    std::vector<std::vector<double>> frequencies;
    std::size_t numFrequencies = 10;


    for (std::string& meshFileName : meshes)
    {
        auto result = solver->solve(meshFileName, structureIndex, kpoint);
        // TODO: Raw eigenvalues is better than the computed frequencies in case something goes wrong

        // Copy the required number of frequencies
        std::vector<double> freqs;
        const std::vector<double>& resultFreqs = result->frequencies(0);
        for (std::size_t i = 0; i < numFrequencies && i < resultFreqs.size(); ++i)
        {
            freqs.emplace_back(resultFreqs[i]);
        }
        frequencies.emplace_back(freqs);
    }

    // Compute Errors //
    ////////////////////

    std::vector<std::vector<double>> errors;

    // Probably overkill
    std::map<double, std::size_t> spectrum = structure->computeSpectrum(kpoint, 20);
    // Create a linear spectrum for easy comparison.
    std::vector<double> linearSpectrum;
    for (auto const &entry : spectrum)
    {
        for(std::size_t i = 0; i < entry.second; ++i)
        {
            linearSpectrum.emplace_back(entry.first);
        }
    }
    //TODO: Add option for the analytical spectrum to compute at least N frequencies
    logger.assert_always(linearSpectrum.size() >= numFrequencies, "Not enough analytical frequencies");

    for(auto & meshFrequencies : frequencies)
    {
        std::vector<double> tempErrors (meshFrequencies.size(), 0);
        for(std::size_t j = 0; j < meshFrequencies.size(); ++j)
        {
            tempErrors[j] = std::abs(meshFrequencies[j] - linearSpectrum[j]);
        }
        errors.emplace_back(tempErrors);
    }

    // Print results //
    ///////////////////

    // Header line
    std::cout << "|";
    for (std::size_t i = 0; i < numFrequencies; ++i)
    {
        std::cout << " "
            << std::setprecision(4)
            << std::setw(6)
            << linearSpectrum[i] << " |";
    }
    std::cout << std::endl
        << "|"; // For the header bottom line
    for (std::size_t i = 0; i < numFrequencies; ++i)
    {
        std::cout << "--------|";
    }
    std::cout << std::endl;
    // Frequencies
    for (auto &meshFrequencies : frequencies)
    {
        std::cout << "|";
        for (double &frequency : meshFrequencies)
        {
            std::cout << " "
                << std::setprecision(4)
                << std::setw(6) << frequency << " |";
        }
        std::cout << std::endl;
    }
    // Spacing newline
    std::cout << std::endl;

    // Errors
    // Header line
    std::cout << "|";
    for (std::size_t i = 0; i < numFrequencies; ++i)
    {
        std::cout << " "
                  << std::setprecision(4)
                  << std::setw(8)
                  << linearSpectrum[i]
                  << "        " // Space of second column
                  << " |";
    }
    std::cout << std::endl
              << "|"; // For the header bottom line
    for (std::size_t i = 0; i < numFrequencies; ++i)
    {
        std::cout << "----------|-------|";
    }
    std::cout << std::endl;
    // Frequencies
    for (std::size_t i = 0; i < errors.size(); ++i)
    {
        auto &meshErrors = errors[i];
        std::cout << "|";
        for (std::size_t j = 0; j < meshErrors.size(); ++j)
        {
            double &error = meshErrors[j];
            std::cout << " "
                      << std::scientific
                      << std::setprecision(2)
                      << std::setw(6) << error << " |";
            if (i == 0)
            {
                // For first row there is no convergence speed
                // so print a placeholder instead.
                std::cout << "     - |";
            }
            else
            {
                double factor = errors[i-1][j]/error;
                std::cout << " "
                    << std::fixed
                    << std::setprecision(2)
                    << std::setw(5)
                    << factor << " |";
            }
        }
        std::cout << std::endl;
    }
    // Spacing newline
    std::cout << std::endl;

}


std::unique_ptr<BaseEigenvalueResult<2>> DivDGMaxSolver::solve(
        const std::string &meshFileName, std::size_t structureIndex,
        const LinearAlgebra::SmallVector<2> &kpoint)
{
    Base::ConfigurationData configData (2, 1);
    auto mesh = DGMax::readMesh<2>(meshFileName, &configData,
            [&](const Geometry::PointPhysical<2>& p) {
        // TODO: Hardcoded structure
        return jelmerStructure(p, structureIndex);
    });
    DGMaxLogger(INFO, "Loaded mesh % with % local elements.", meshFileName, mesh->getNumberOfElements());
    KSpacePath<2> path = KSpacePath<2>::singleStepPath(kpoint);
    EigenValueProblem<2> input(path, 10); // TODO: Customize number of eigenvalues

    DivDGMaxEigenValue<2> solver (*mesh);
    typename DivDGMaxDiscretization<2>::Stab stab; // Change this?
    stab.stab1 = 100;
    stab.stab2 = 0;
    stab.stab3 = 1;
    // TODO: move p-value
    typename DivDGMaxEigenValue<2>::Result result = solver.solve(input, stab, 1);


    return std::unique_ptr<BaseEigenvalueResult<2>>(new DivDGMaxEigenValue<2>::Result(result));
}

#ifndef HPGEM_EIGENVALUERESULT_H
#define HPGEM_EIGENVALUERESULT_H

#include <vector>

namespace DGMax {

/// \brief Frequencies from eigenvalue computations on a set of meshes.
class EigenvalueResult {
   public:
    /// Print the raw frequencies
    ///
    /// \param theoretical Theoretical frequencies used as header.
    void printFrequencyTable(std::vector<double> theoretical);
    /// Print a table with errors and convergence rates
    ///
    /// \param theoretical The theoretical results, used for printing
    /// convergence order.
    void printErrorTable(std::vector<double> theoretical);

    /// Compare with another result to test whether the results are equal.
    /// \param other The other result to compare with
    /// \param tolerance The absolute tolerance in the frequency differences
    /// \return Whether this and other are the same up to tolerance.
    bool equals(const EigenvalueResult &other, double tolerance);

    /// Add the frequencies of a level.
    void addLevel(const std::vector<double> &frequencies) {
        frequencyLevels_.emplace_back(frequencies);
    }

    /// Maximum number of frequencies on any level.
    std::size_t maxNumberOfFrequencies() {
        std::size_t max = 0;
        for (const auto &level : frequencyLevels_) {
            max = std::max(max, level.size());
        }
        return max;
    }

   private:
    std::vector<std::vector<double>> frequencyLevels_;
};

}  // namespace DGMax

#endif  // HPGEM_EIGENVALUERESULT_H

#ifndef HPGEM_EVCONVERGENCERESULT_H
#define HPGEM_EVCONVERGENCERESULT_H

#include <vector>

#include "Logger.h"

namespace DGMax {

/// \brief Frequencies from eigenvalue computations on a set of meshes.
class EVConvergenceResult {
   public:
    /// Print the raw frequencies
    ///
    /// \param theoretical Theoretical frequencies used as header.
    void printFrequencyTable(std::vector<double> theoretical) const;
    /// Print a table with errors and convergence rates
    ///
    /// \param theoretical The theoretical results, used for printing
    /// convergence order.
    void printErrorTable(std::vector<double> theoretical) const;

    void filterResults(double minimum, bool removeNaN);

    std::size_t getNumberOfLevels() const { return frequencyLevels_.size(); }

    const std::vector<double> &getLevel(std::size_t level) const {
        logger.assert_debug(level < frequencyLevels_.size(), "Level too large");
        return frequencyLevels_[level];
    }

    /// Add the frequencies of a level.
    void addLevel(const std::vector<double> &frequencies) {
        frequencyLevels_.emplace_back(frequencies);
    }

    /// Maximum number of frequencies on any level.
    std::size_t maxNumberOfFrequencies() const {
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

#endif  // HPGEM_EVCONVERGENCERESULT_H

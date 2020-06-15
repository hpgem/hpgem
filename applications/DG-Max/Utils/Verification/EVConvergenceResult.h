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

#ifndef HPGEM_EVCONVERGENCERESULT_H
#define HPGEM_EVCONVERGENCERESULT_H

#include <vector>

#include "Logger.h"

namespace DGMax {

/// \brief Frequencies from eigenvalue computations on a set of meshes.
class EVConvergenceResult {
   public:
    EVConvergenceResult() = default;
    /// Constructor for hard coded results
    EVConvergenceResult(std::vector<std::vector<double>> rawResults)
        : frequencyLevels_(rawResults) {}
    EVConvergenceResult(
        const std::initializer_list<std::vector<double>> rawResults)
        : frequencyLevels_(rawResults.size())
    {
        std::size_t i = 0;
        for (const std::vector<double>& level : rawResults) {
            frequencyLevels_[i] = level;
            i++;
        }
    }

    /// Print the raw frequencies
    ///
    /// \param theoretical Theoretical frequencies used as header.
    void printFrequencyTable(std::vector<double> theoretical) const;
    /// Print a table with errors and convergence rates
    ///
    /// \param theoretical The theoretical results, used for printing
    /// convergence order.
    void printErrorTable(std::vector<double> theoretical) const;

    /// Print the results in a format that can be copied into code.
    /// \param precision The precision to print with
    void printResultCode(int precision) const;

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

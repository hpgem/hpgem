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

#include "EVConvergenceResult.h"

#include <iomanip>
#include <iostream>
#include <cmath>

namespace DGMax {

void EVConvergenceResult::printFrequencyTable(
    std::vector<double> theoretical) const {
    std::size_t max = maxNumberOfFrequencies();

    // Header line
    std::cout << "|";
    for (std::size_t i = 0; i < max; ++i) {
        if (i < theoretical.size()) {
            std::cout << " " << std::setprecision(4) << std::setw(6)
                      << theoretical[i] << " |";
        } else {
            // When lacking an theoretical frequency, just print a placeholder
            std::cout << " f_" << std::left << std::setw(4) << i << " |";
        }
    }
    std::cout << std::endl << "|";  // For the header bottom line
    for (std::size_t i = 0; i < max; ++i) {
        std::cout << "--------|";
    }
    std::cout << std::endl;
    // Frequencies
    for (auto &levelFrequency : frequencyLevels_) {
        std::cout << "|";
        for (const double &frequency : levelFrequency) {
            std::cout << " " << std::setprecision(4) << std::setw(6)
                      << frequency << " |";
        }
        std::cout << std::endl;
    }
    // Spacing newline
    std::cout << std::endl;
}

void EVConvergenceResult::printErrorTable(std::vector<double> theoretical) const {

    std::size_t max = maxNumberOfFrequencies();

    std::vector<std::vector<double>> errors;
    // Compute actual errors
    for (auto &levelFrequencies : frequencyLevels_) {
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
            std::cout << "f_" << std::setw(6) << std::left << i;
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

void EVConvergenceResult::printResultCode(int precision) const {
    std::cout << "{" << std::endl;
    for (std::size_t i = 0; i < frequencyLevels_.size(); ++i) {
        const std::vector<double> frequencies = frequencyLevels_[i];

        std::cout << "{";
        for (std::size_t j = 0; j < frequencies.size(); ++j) {
            std::cout << std::setprecision(precision) << frequencies[j];
            if (j != frequencies.size() - 1) {
                std::cout << ",";
            }
        }
        std::cout << "}";
        if (i != frequencyLevels_.size() - 1) {
            std::cout << ",";
        }
        std::cout << std::endl;
    }
    std::cout << "}" << std::endl;
}

void EVConvergenceResult::filterResults(double minimum, bool removeNaN) {
    for (std::vector<double>& frequencies : frequencyLevels_) {
        auto it = frequencies.begin();
        while (it != frequencies.end()) {
            if (*it < minimum || (removeNaN && std::isnan(*it))) {
                it = frequencies.erase(it);
            } else {
                it++;
            }
        }
    }
}

}  // namespace DGMax

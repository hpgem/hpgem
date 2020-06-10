#include "EVConvergenceResult.h"

#include <iomanip>
#include <iostream>
#include <cmath>

namespace DGMax {

void EVConvergenceResult::printFrequencyTable(std::vector<double> theoretical) const {
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

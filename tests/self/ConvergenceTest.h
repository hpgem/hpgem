/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2021, University of Twente
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
#ifndef HPGEM_CONVERGENCETEST_H
#define HPGEM_CONVERGENCETEST_H

#include "Logger.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace hpgem {

/**
 * Input of several meshes with expected errors for testing. This supports
 * multiple different error measures to be tested simultaneously.
 */
class ConvergenceTestSet {
   public:
    ConvergenceTestSet(std::vector<std::string> meshes,
                       double relativeAccuracy = 0.01)
        : meshes(std::move(meshes)),
          expectedErrors(),
          relativeAccuracy(relativeAccuracy){};

    // Backwards compatibility to test with a single measure
    ConvergenceTestSet(std::vector<std::string> meshes,
                       std::vector<double> expectedErrors,
                       double relativeAccuracy = 0.01)
        : meshes(std::move(meshes)),
          errorNames({"error"}),
          expectedErrors({std::move(expectedErrors)}),
          relativeAccuracy(relativeAccuracy){};

    /// Add a set of errors for testing
    void addExpectedErrors(std::string name, std::vector<double> errors) {
        errorNames.push_back(std::move(name));
        expectedErrors.push_back(std::move(errors));
    }

    std::vector<std::string> meshes;
    std::vector<std::string> errorNames;
    std::vector<std::vector<double>> expectedErrors;
    // Criterion for checking the errors:
    // |error - expectedError|  < relAccuracy * |expectedError|
    // Which allows slight numerical differences between platforms
    double relativeAccuracy;
};

/**
 * Run a test on multiple meshes. Primarily intended to do convergence tests,
 * and it will therefore print a simple convergence table with the ratios
 * between subsequent errors.
 *
 * @param testSet The meshes and errors for the test
 * @param ignoreFailures Whether to ignore failure (e.g. for recomputing the
 * errors)
 * @param solver Method that accepts the mesh, test input and mesh number and
 * produces the error.
 */
void runConvergenceTest(
    ConvergenceTestSet& testSet, bool ignoreFailures,
    std::function<std::vector<double>(std::string, std::size_t)> solver) {

    std::size_t numErrors = testSet.expectedErrors.size();
    std::vector<std::vector<double>> observedErrors(numErrors);

    for (std::size_t i = 0; i < testSet.meshes.size(); ++i) {
        std::vector<double> levelErrors = solver(testSet.meshes[i], i);
        for (std::size_t j = 0; j < numErrors; ++j) {
            double error = levelErrors[j];
            observedErrors[j].push_back(error);

            if (i < testSet.expectedErrors[j].size()) {
                double difference =
                    std::abs(error - testSet.expectedErrors[j][i]);
                logger.assert_always(
                    difference / testSet.expectedErrors[j][i] <
                            testSet.relativeAccuracy ||
                        ignoreFailures,
                    "Comparing to old results, expected %, got % "
                    "for % at the %-th level",
                    testSet.expectedErrors[j][i], testSet.errorNames[j], error,
                    i);
            } else if (!ignoreFailures) {
                // Require expected error information available when actually
                // testing
                logger.assert_always(false, "No error data for % at level %",
                                     testSet.errorNames[j], i);
            }
        }
    }

    // Print convergence table
    for (std::size_t j = 0; j < numErrors; ++j) {
        std::cout << testSet.errorNames[j] << std::endl;
        for (std::size_t i = 0; i < observedErrors[j].size(); ++i) {
            std::cout << std::setprecision(8) << std::setw(15)
                      << std::scientific << observedErrors[j][i];
            // Separator that allows copy pasting the result into the code
            // including the convergence rate. Note double space to match the
            // formatting standard.
            std::cout << ",  //";

            // Compute convergence rate
            if (i == 0) {
                std::cout << "------";
            } else {
                double rate = observedErrors[j][i - 1] / observedErrors[j][i];
                std::cout << std::setprecision(2) << std::setw(6) << std::fixed
                          << rate;
            }
            std::cout << std::endl;
        }
    }
}

void runConvergenceTest(
    ConvergenceTestSet& testSet, bool ignoreFailures,
    std::function<double(std::string, std::size_t)> solver) {
    runConvergenceTest(
        testSet, ignoreFailures,
        [&solver](std::string mesh, std::size_t level) -> std::vector<double> {
            return {solver(std::move(mesh), level)};
        });
}

}  // namespace hpgem

#endif  // HPGEM_CONVERGENCETEST_H

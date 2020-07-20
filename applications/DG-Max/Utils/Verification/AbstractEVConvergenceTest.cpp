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

#include "AbstractEVConvergenceTest.h"

#include "Logger.h"

namespace DGMax {

template <std::size_t DIM>
EVConvergenceResult AbstractEVConvergenceTest<DIM>::run(bool failOnDifference) {
    EVConvergenceResult convergenceResult;

    for (std::size_t level = 0; level < getNumberOfLevels(); ++level) {
        std::unique_ptr<AbstractEigenvalueResult<DIM>> result =
            runInternal(level);

        bool hasResult(result);
        logger.assert_always(hasResult, "Null result");

        convergenceResult.addLevel(result->frequencies(0));

        if (getExpected() != nullptr) {
            bool correct = compareWithExpected(level, *result);

            logger.assert_always(correct || !failOnDifference,
                                 "Results differ from expected");
        }
    }
    return convergenceResult;
}

template <std::size_t DIM>
bool AbstractEVConvergenceTest<DIM>::compareWithExpected(
    std::size_t level, const AbstractEigenvalueResult<DIM>& result) const {

    auto expected = getExpected();
    if (expected == nullptr || expected->getNumberOfLevels() <= level) {
        logger(WARN, "No expected results for level %", level);
        return true;
    }

    const std::vector<double>& expectedLevel = expected->getLevel(level);
    const std::vector<double>& resultLevel = result.frequencies(0);

    // Skip over the starting NaNs
    auto resultIter = resultLevel.begin();
    while (expectNaNsOrNegative() && resultIter != resultLevel.end() &&
           (std::isnan(*resultIter) || (*resultIter) < 0.0)) {
        ++resultIter;
    }
    auto expectedIter = expectedLevel.begin();
    while (expectNaNsOrNegative() && expectedIter != expectedLevel.end() &&
           (std::isnan(*expectedIter) || (*expectedIter) < 0.0)) {
        ++expectedIter;
    }
    // Check if the result is sorted as it should be
    if (!std::is_sorted(resultIter, resultLevel.end())) {
        logger(ERROR, "Non sorted result");
        return false;
    }
    logger.assert_always(std::is_sorted(expectedIter, expectedLevel.end()),
                         "Non sorted expected result");
    // Skip over any numerical zeros
    while (resultIter != resultLevel.end() && (*resultIter) >= 0.0 &&
           (*resultIter) < expectedNumericalZeroThreshold()) {
        ++resultIter;
    }
    while (expectedIter != expectedLevel.end() && (*expectedIter) >= 0.0 &&
           (*expectedIter) < expectedNumericalZeroThreshold()) {
        ++expectedIter;
    }

    // Check if there are enough values to compare
    logger.assert_always(
        expectedLevel.end() - expectedIter >= minimumNumberOfResults(level),
        "Not enough expected values for a comparison");
    if (resultLevel.end() - resultIter < minimumNumberOfResults(level)) {
        logger(ERROR,
               "Not enough result values for a comparison, minimum %, got %",
               minimumNumberOfResults(level), resultLevel.end() - resultIter);
        return false;
    }
    // Now the actual comparison

    bool same = true;
    while (resultIter != resultLevel.end() &&
           expectedIter != expectedLevel.end()) {
        double diff = std::abs(*resultIter - *expectedIter);
        if (std::abs(diff) >= getTolerance()) {
            logger(WARN, "Different eigenvalue detected: Expexted %, result %",
                   *expectedIter, *resultIter);
            same = false;
            break;
        }
        ++resultIter;
        ++expectedIter;
    }
    // Print complete output
    if (!same) {
        std::stringstream expectedString;
        expectedString << "Expected:";
        for (auto& expectedFreq : expectedLevel) {
            expectedString << " " << expectedFreq;
        }

        std::stringstream actualString;
        actualString << "Computed:";
        for (auto& resultFreq : resultLevel) {
            actualString << " " << resultFreq;
        }
        logger(ERROR, "Results\n%\n%\n", expectedString.str(),
               actualString.str());
    }
    return same;
}

template class AbstractEVConvergenceTest<2>;
template class AbstractEVConvergenceTest<3>;

}  // namespace DGMax

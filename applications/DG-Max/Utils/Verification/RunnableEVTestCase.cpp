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

#include "RunnableEVTestCase.h"

#include "Logger.h"

namespace DGMax {

template <std::size_t DIM>
EVConvergenceResult RunnableEVTestCase<DIM>::run(bool failOnDifference) {
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
bool RunnableEVTestCase<DIM>::compareWithExpected(
    std::size_t level, const AbstractEigenvalueResult<DIM>& result) const {

    auto expected = getExpected();
    if (expected == nullptr || expected->getNumberOfLevels() <= level) {
        return true;
    }
    const std::vector<double>& expectedLevel = expected->getLevel(level);
    const std::vector<double>& resultLevel = result.frequencies(0);

    if (expectedLevel.size() != resultLevel.size()) {
        logger(ERROR,
               "Different number of eigenvalues expected % got % at level %",
               expectedLevel.size(), resultLevel.size(), level);
        return false;
    }
    for (std::size_t i = 0; i < expectedLevel.size(); ++i) {
        double diff = std::abs(expectedLevel[i] - resultLevel[i]);
        if (diff >= getTolerance() ||
            (std::isnan(expectedLevel[i]) != std::isnan(resultLevel[i]))) {
            logger(ERROR,
                   "Different %-th eigenvalue at level %: Expected % got %", i,
                   level, expectedLevel[i], resultLevel[i]);
            return false;
        }
    }
    return true;
}

template class RunnableEVTestCase<2>;
template class RunnableEVTestCase<3>;

}  // namespace DGMax

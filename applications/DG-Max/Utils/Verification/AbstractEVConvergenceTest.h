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

#ifndef HPGEM_ABSTRACTEVCONVERGENCETEST_H
#define HPGEM_ABSTRACTEVCONVERGENCETEST_H

#include <memory>

#include "ProblemTypes/AbstractEigenvalueResult.h"
#include "EVTestPoint.h"

namespace DGMax {

/// Abstract convergence test for bandstructure eigenvalue solvers.
template <std::size_t DIM>
class AbstractEVConvergenceTest {
   public:
    virtual ~AbstractEVConvergenceTest() = default;

    /// Run the convergence test
    ///
    /// \param failOnDifference If there are differences from the expected
    /// result, create an assertion failure.
    /// \return The results
    EVConvergenceResult run(bool failOnDifference);

   protected:
    /// Number of mesh levels
    virtual std::size_t getNumberOfLevels() const = 0;
    /// Absolute tolerance for the frequencies in checking the results
    virtual double getTolerance() const = 0;
    /// Expected result (null if none)
    virtual const EVConvergenceResult* getExpected() const = 0;
    /// Run the actual algorithm on a single level.
    virtual std::unique_ptr<AbstractEigenvalueResult<DIM>> runInternal(
        std::size_t level) = 0;

   private:
    /// Compare the results with the expected results of a specific level.
    ///
    /// \param level The level in the expected results
    /// \param result The actual computed result
    /// \return  Whether the result is inconsistent with the expected data.
    bool compareWithExpected(std::size_t level,
                             const AbstractEigenvalueResult<DIM>& result) const;
};

}  // namespace DGMax

#endif  // HPGEM_ABSTRACTEVCONVERGENCETEST_H

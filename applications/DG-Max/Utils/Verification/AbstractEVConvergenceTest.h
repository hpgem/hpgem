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

#include "ProblemTypes/AbstractEigenvalueSolver.h"
#include "EVTestPoint.h"

namespace DGMax {

/// Abstract convergence test for bandstructure eigenvalue solvers.
template <std::size_t DIM>
class AbstractEVConvergenceTest {

   protected:
    // Driver that solves the eigenvalue problem at a single point and stores
    // the result.
    class Driver : public AbstractEigenvalueSolverDriver<DIM> {
       public:
        Driver(const LinearAlgebra::SmallVector<DIM>& kpoint,
               std::size_t numberOfEigenvalues)
            : point_(kpoint),
              numberOfEigenvalues_(numberOfEigenvalues),
              nextInvoked_(false),
              frequencyResult_(0){};

        bool stop() const final { return nextInvoked_; };

        void nextKPoint() final { nextInvoked_ = true; }

        LinearAlgebra::SmallVector<DIM> getCurrentKPoint() const final {
            return point_;
        }
        size_t getNumberOfKPoints() const override { return 1; }
        size_t getTargetNumberOfEigenvalues() const override {
            return numberOfEigenvalues_;
        }
        void handleResult(AbstractEigenvalueResult<DIM>& result) override {
            frequencyResult_ = result.getFrequencies();
        }

        std::vector<double> getResult() { return frequencyResult_; }

       private:
        LinearAlgebra::SmallVector<DIM> point_;
        std::size_t numberOfEigenvalues_;
        bool nextInvoked_;
        std::vector<double> frequencyResult_;
    };

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
    virtual void runInternal(Driver& driver, std::size_t level) = 0;

    /// Options for checking the result, by default perfect results

    /// Whether to expect NaNs and negative values at the lower end of the
    /// frequency spectrum.
    virtual bool expectNaNsOrNegative() const { return false; }
    /// The threshold for numerical zeros (use negative values for no expected
    /// zeros).
    virtual double expectedNumericalZeroThreshold() const { return -1.0; }

    /// The minimum number of actual results (not filtered by
    /// expectNaNsOrNegative nor by expectedNumericalZeroThreshold()) that
    /// should be included for a result to be consistent with the expected
    /// result.
    ///
    /// \param level The level
    /// \return The number of actual frequencies that should match.
    virtual std::size_t minimumNumberOfResults(std::size_t level) const {
        return getExpected()->getLevel(level).size();
    }

    virtual const LinearAlgebra::SmallVector<DIM>& getKPoint() const = 0;

    virtual std::size_t getTargetNumberOfEigenvalues() const = 0;

   private:
    /// Compare the results with the expected results of a specific level.
    ///
    /// \param level The level in the expected results
    /// \param result The actual computed frequencies
    /// \return  Whether the result is consistent with the expected data.
    bool compareWithExpected(std::size_t level,
                             const std::vector<double>& frequencies) const;
};

}  // namespace DGMax

#endif  // HPGEM_ABSTRACTEVCONVERGENCETEST_H

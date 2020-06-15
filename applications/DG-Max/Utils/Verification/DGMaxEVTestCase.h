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

#ifndef HPGEM_DGMAXEVTESTCASE_H
#define HPGEM_DGMAXEVTESTCASE_H

#include <memory>

#include "Utils/Verification/RunnableEVTestCase.h"
#include "Utils/Verification/EVTestCase.h"

namespace DGMax {

template <std::size_t DIM>
class DGMaxEVTestCase : public RunnableEVTestCase<DIM> {
   public:
    DGMaxEVTestCase(EVTestCase<DIM> testCase,
                    std::vector<std::string> meshFileNames, double tolerance,
                    std::size_t order, double stab,
                    EVConvergenceResult* expected)
        : testCase_(testCase),
          meshFileNames_(std::move(meshFileNames)),
          tolerance_(tolerance),
          order_(order),
          stab_(stab),
          expected_(expected){};

   protected:
    std::size_t getNumberOfLevels() const override {
        return meshFileNames_.size();
    }
    double getTolerance() const override { return tolerance_; }
    const EVConvergenceResult* getExpected() const override {
        return expected_;
    }
    /// Run the actual algorithm on a single level.
    std::unique_ptr<AbstractEigenvalueResult<DIM>> runInternal(
        std::size_t level) override;

   private:
    EVTestCase<DIM> testCase_;
    std::vector<std::string> meshFileNames_;
    double tolerance_;
    std::size_t order_;
    double stab_;
    // TODO: This is usually static data
    EVConvergenceResult* expected_;
};

}  // namespace DGMax

#endif  // HPGEM_DGMAXEVTESTCASE_H

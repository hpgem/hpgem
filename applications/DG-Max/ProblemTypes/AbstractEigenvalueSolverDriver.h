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
#ifndef HPGEM_EIGENVALUEDRIVER_H
#define HPGEM_EIGENVALUEDRIVER_H

#include "LinearAlgebra/SmallVector.h"
#include "AbstractEigenvalueResult.h"

/**
 * Driver for an AbstractEigenvalueSolver, that instructs the solver what
 * eigenvalue problems to solve and what to do with the results.
 *
 * @tparam DIM The dimension of the problem.
 */
template <std::size_t DIM>
class AbstractEigenvalueSolverDriver {
   public:
    virtual ~AbstractEigenvalueSolverDriver() = default;

    /**
     * @return Whether to stop solving
     */
    virtual bool stop() const = 0;
    /**
     * Advance to the next wave vector point for solving
     */
    virtual void nextKPoint() = 0;

    /**
     * @return The current k-point
     */
    virtual LinearAlgebra::SmallVector<DIM> getCurrentKPoint() const = 0;

    /**
     * The total number of wave vectors/k-points, if known a-priori.
     * @return The number of wave vector points, or 0 if unknown.
     */
    virtual std::size_t getNumberOfKPoints() const = 0;

    /**
     * The target number of eigenvalues to solve for at the current k-point.
     * Depending on the actual configuration of the solver the result may
     * include less eigenvalues (e.g. if a maximum iteration count was reached)
     * or slightly more (e.g. when a solver iteration provides more than
     * required).
     *
     * @return The number of eigenvalues.
     */
    virtual std::size_t getTargetNumberOfEigenvalues() const = 0;

    /**
     * Handle the result of a solved eigenvalue problem.
     * @param result The result handle, any references to results are only valid
     *   during the call to this method.
     */
    virtual void handleResult(AbstractEigenvalueResult<DIM>& result) = 0;
};

#endif  // HPGEM_EIGENVALUEDRIVER_H

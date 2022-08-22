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
#ifndef HPGEM_HARMONICERRORDRIVER_H
#define HPGEM_HARMONICERRORDRIVER_H

#include <ProblemTypes/AbstractHarmonicSolverDriver.h>
#include <Output/VTKSpecificTimeWriter.h>

template <std::size_t dim>
class HarmonicErrorDriver : public DGMax::AbstractHarmonicSolverDriver<dim> {
   private:
    using VecR = LinearAlgebra::SmallVector<dim>;
    using VecC = LinearAlgebra::SmallVectorC<dim>;

    struct Solution {
        VecC theory;
        VecC numerical;
    };

   public:
    HarmonicErrorDriver(ExactHarmonicProblem<dim>& problem)
        : problem_(&problem),
          nextCalled_(false),
          errorResult_(std::nan("")),
          plotter_(nullptr){};

    bool stop() const override { return nextCalled_; }
    void nextProblem() override { nextCalled_ = true; }

    const HarmonicProblem<dim>& currentProblem() const override {
        return *problem_;
    }
    void handleResult(DGMax::AbstractHarmonicResult<dim>& result) override {
        errorResult_ = result.computeL2Error(*problem_);
        if (plotter_ != nullptr) {
            result.writeVTK(*plotter_);
            // Plot the exact solution for comparison

            std::map<std::string, std::function<double(Solution&)>> scalars;
            scalars["Emag-sol"] = [](Solution& v) { return v.theory.l2Norm(); };
            scalars["Error"] = [](Solution& v) {
                return (v.theory - v.numerical).l2Norm();
            };
            std::map<std::string, std::function<VecR(Solution&)>> vectors;
            vectors["ESol-real"] = [](Solution& v) { return v.theory.real(); };
            vectors["ESol-imag"] = [](Solution& v) { return v.theory.imag(); };

            plotter_->template writeMultiple<Solution>(
                [this, &result](Base::Element* element,
                                const Geometry::PointReference<dim>& p,
                                std::size_t) {
                    Solution sol;
                    sol.theory = problem_->exactSolution(
                        element->template referenceToPhysical(p));
                    sol.numerical = result.computeField(element, p);
                    return sol;
                },
                scalars, vectors);
        }
    }

    double getError() { return errorResult_; }

    void setOutputPlotter(Output::VTKSpecificTimeWriter<dim>* plotter) {
        plotter_ = plotter;
    }

   private:
    ExactHarmonicProblem<dim>* problem_;
    bool nextCalled_;
    double errorResult_;
    Output::VTKSpecificTimeWriter<dim>* plotter_;
};

#endif  // HPGEM_HARMONICERRORDRIVER_H

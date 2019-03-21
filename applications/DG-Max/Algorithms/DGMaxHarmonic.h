/*
This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
http://hpgem.org/about-the-code/team

This code is distributed using BSD 3-Clause License. A copy of which can found below.


Copyright (c) 2018, Univesity of Twenete
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef ALGORITHMS_DGMAXHARMONIC_h
#define ALGORITHMS_DGMAXHARMONIC_h

#include "../BaseExtended.h"
#include "../ProblemTypes/HarmonicProblem.h"

#include "DGMaxDiscretization.h"

/// \brief Solver for a harmonic problem to find the fields.
template<std::size_t DIM>
class DGMaxHarmonic
{
public:
    explicit DGMaxHarmonic(hpGemUIExtentions<DIM>& base, std::size_t order);
    void solve(const HarmonicProblem<DIM>& harmonicProblem, double stab);

    std::map<typename DGMaxDiscretization<DIM>::NormType, double> computeError(
            const typename std::set<typename DGMaxDiscretization<DIM>::NormType>& norms,
            const typename DGMaxDiscretization<DIM>::InputFunction& exactSolution,
            const typename DGMaxDiscretization<DIM>::InputFunction& exactSolutionCurl
    ) const;

    std::map<typename DGMaxDiscretization<DIM>::NormType, double> computeError(
            const std::set<typename DGMaxDiscretization<DIM>::NormType>& norms,
            const ExactHarmonicProblem<DIM>& problem
    ) const;

    void writeTec(std::string fileName) const;

private:
    hpGemUIExtentions<DIM>& base_;
    DGMaxDiscretization<DIM> discretization;
};


#endif //ALGORITHMS_DGMAXHARMONIC_h

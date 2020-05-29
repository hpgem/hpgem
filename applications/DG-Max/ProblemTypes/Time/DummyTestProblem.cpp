/*
This file forms part of hpGEM. This package has been developed over a number of
years by various people at the University of Twente and a full list of
contributors can be found at http://hpgem.org/about-the-code/team

This code is distributed using BSD 3-Clause License. A copy of which can found
below.


Copyright (c) 2018, Univesity of Twenete
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "DummyTestProblem.h"

template <std::size_t DIM>
void DummyTestProblem<DIM>::initialConditionDerivative(
    const Geometry::PointPhysical<DIM> &point,
    LinearAlgebra::SmallVector<DIM> &result) const {
    for (std::size_t i = 0; i < DIM; i++) {
        result[i] = 0;
    }
}

template <std::size_t DIM>
void DummyTestProblem<DIM>::exactSolution(
    const Geometry::PointPhysical<DIM> &p, double t,
    LinearAlgebra::SmallVector<DIM> &ret) const {
    // TODO: Code literaly copied, no verification
    // ret[0]=sin(M_PI*2*p[1])*sin(M_PI*2*p[2]);
    // ret[1]=sin(M_PI*2*p[2])*sin(M_PI*2*p[0]);
    // ret[2]=sin(M_PI*2*p[0])*sin(M_PI*2*p[1]);
    // ret*=cos(sqrt(2)*2*M_PI*t);

    // for comparison with Freekjan's report.
    ret[0] = sin(M_PI * p[1]) * sin(M_PI * p[2]);
    ret[1] = sin(M_PI * p[2]) * sin(M_PI * p[0]);
    ret[2] = sin(M_PI * p[0]) * sin(M_PI * p[1]);
    ret *= cos(sqrt(2) * M_PI * t);
    // ret[0] = p[2];
    // ret[1] = p[0];
    // ret[2] = p[1];

    //     ret[0]=p[0]*(1-p[0]);
    //     ret[1]=0;
    // 	   ret[2]=0;
}

template <std::size_t DIM>
void DummyTestProblem<DIM>::exactSolutionCurl(
    const Geometry::PointPhysical<DIM> &p, double t,
    LinearAlgebra::SmallVector<DIM> &ret) const {
    // TODO: Code literaly copied, no verification
    // ret[0]=sin(M_PI*2*p[0])*(cos(M_PI*2*p[1])-cos(M_PI*2*p[2]));
    // ret[1]=sin(M_PI*2*p[1])*(cos(M_PI*2*p[2])-cos(M_PI*2*p[0]));
    // ret[2]=sin(M_PI*2*p[2])*(cos(M_PI*2*p[0])-cos(M_PI*2*p[1]));
    // ret*=cos(sqrt(2)*2*M_PI*t)*2*M_PI;

    // for comparison with Freekjan's report.
    ret[0] = sin(M_PI * p[0]) * (cos(M_PI * p[1]) - cos(M_PI * p[2]));
    ret[1] = sin(M_PI * p[1]) * (cos(M_PI * p[2]) - cos(M_PI * p[0]));
    ret[2] = sin(M_PI * p[2]) * (cos(M_PI * p[0]) - cos(M_PI * p[1]));
    ret *= cos(sqrt(2) * M_PI * t) * M_PI;

    // ret[0] = 1.0;
    // ret[1] = 1.0;
    // ret[2] = 1.0;

    //          ret[0]=0;ret[1]=0;ret[2]=0;
}

template <std::size_t DIM>
double DummyTestProblem<DIM>::timeScalingBoundary(double t) const {
    return 1.0;  // for comparison with Freekjan's report. The 1.0 is not
                 // physically meaningful. The source term is actually 0.0, but
                 // we want to avoid division by zero in timedependent code.
}

template <std::size_t DIM>
void DummyTestProblem<DIM>::sourceTermRef(
    const Geometry::PointPhysical<DIM> &point,
    LinearAlgebra::SmallVector<DIM> &result) const {
    exactSolution(point, 0, result);
    // 	ret*=-1;
    // ret*=M_PI*M_PI*8-1;
    result *= M_PI * M_PI * 2 - 1;
    // ret *= -1;
    // for comparison with the time-dependent code in Freekjan's report.
    // ret[0] = 0.0;
    // ret[1] = 0.0;
    // ret[2] = 0.0;
}

template <std::size_t DIM>
double DummyTestProblem<DIM>::timeScalingSource(double t) const {
    return 0;
}

template class DummyTestProblem<2>;
template class DummyTestProblem<3>;

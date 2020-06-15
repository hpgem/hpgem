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

#include "SampleHarmonicProblems.h"

template <std::size_t DIM>
SampleHarmonicProblems<DIM>::SampleHarmonicProblems(
    SampleHarmonicProblems<DIM>::Problem problem, double omega)
    : problem_(problem), omega_(omega) {
    // Note: the tolerances here are just a guess.
    switch (problem) {
        case CONSTANT: {
            logger.assert_debug(
                std::abs(omega_) > 1e-3,
                "Given omega is very close to the eigenvalue 0.");
        }
        case SARMANY2010: {
            logger.assert_debug(
                std::abs(omega_ * omega_ - 2 * M_PI * M_PI) > 1e-1,
                "Given omega^2 is very close to the eigenvalue 2 pi^2.");
            break;
        }
    }
}

template <std::size_t DIM>
double SampleHarmonicProblems<DIM>::omega() const {
    switch (problem_) {
        case CONSTANT:
            return omega_;
        case SARMANY2010:
            return omega_;
        default:
            logger.assert_debug(false, "Not implemented for this problem.");
            return -1;
    }
}

template <std::size_t DIM>
void SampleHarmonicProblems<DIM>::exactSolution(
    const Geometry::PointPhysical<DIM> &point,
    LinearAlgebra::SmallVector<DIM> &result) const {
    switch (problem_) {
        case CONSTANT: {
            result.set(1);
            break;
        }
        case SARMANY2010: {
            sarmanyx(point, result);
            break;
        }
        default:
            logger.assert_debug(false, "Not implemented for this problem.");
    }
}

template <std::size_t DIM>
void SampleHarmonicProblems<DIM>::exactSolutionCurl(
    const Geometry::PointPhysical<DIM> &point,
    LinearAlgebra::SmallVector<DIM> &result) const {
    switch (problem_) {
        case CONSTANT: {
            result.set(0);
            break;
        }
        case SARMANY2010: {
            double x = point[0], y = point[1], z = point[2];
            result[0] = sin(M_PI * x) * (cos(M_PI * y) - cos(M_PI * z));
            result[1] = sin(M_PI * y) * (cos(M_PI * z) - cos(M_PI * x));
            result[2] = sin(M_PI * z) * (cos(M_PI * x) - cos(M_PI * y));
            result *= M_PI;
            break;
        }
        default:
            logger.assert_debug(false, "Not implemented for this problem.");
    }
}

template <std::size_t DIM>
void SampleHarmonicProblems<DIM>::sourceTerm(
    const Geometry::PointPhysical<DIM> &point,
    LinearAlgebra::SmallVector<DIM> &result) const {
    switch (problem_) {
        case CONSTANT: {
            result.set(-omega_ * omega_);
            break;
        }
        case SARMANY2010: {
            sarmanyx(point, result);
            result *= (2 * M_PI * M_PI - omega_ * omega_);
            break;
        }
        default:
            logger.assert_debug(false, "Not implemented for this problem.");
    }
}

template <std::size_t DIM>
void SampleHarmonicProblems<DIM>::sarmanyx(
    const Geometry::PointPhysical<DIM> &point,
    LinearAlgebra::SmallVector<DIM> &result) const {
    if (DIM == 3) {
        double sx = sin(M_PI * point[0]), sy = sin(M_PI * point[1]),
               sz = sin(M_PI * point[2]);
        result[0] = sy * sz;
        result[1] = sz * sx;
        result[2] = sx * sy;
    } else {
        logger.assert_debug(DIM == 3, "Sarmany test case only works in 3D.");
    }
}

template class SampleHarmonicProblems<2>;
template class SampleHarmonicProblems<3>;
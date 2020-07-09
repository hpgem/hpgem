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

#include "SampleTestProblems.h"

using namespace hpgem;


template <std::size_t DIM>
SampleTestProblems<DIM>::SampleTestProblems(
    SampleTestProblems<DIM>::Problem problem)
    : problem_(problem) {}

template <std::size_t DIM>
void SampleTestProblems<DIM>::initialConditionDerivative(
    const Geometry::PointPhysical<DIM> &point,
    LinearAlgebra::SmallVector<DIM> &result) const {
    switch (problem_) {
        case CONSTANT: {
            result.set(0);
            break;
        }
        case LINEAR: {
            result.set(1);
            break;
        }
        case SINSIN: {
            sinx(point, result);
            // d sin(pi t)/dt at t = 0.
            result *= M_PI;
            break;
        }
        case SARMANY2013: {
            // Derivatives of all three cosines are 0.
            result.set(0);
            break;
        }
        default:
            logger.assert_debug(false, "Unknown problem");
    }
}

template <std::size_t DIM>
void SampleTestProblems<DIM>::sourceTermRef(
    const Geometry::PointPhysical<DIM> &point,
    LinearAlgebra::SmallVector<DIM> &result) const {
    switch (problem_) {
        case CONSTANT: {
            result.set(0);
            break;
        }
        case LINEAR: {
            result.set(0);
            break;
        }
        case SINSIN: {
            sinx(point, result);
            result *= -M_PI * M_PI;
            break;
        }
        case SARMANY2013: {
            sarmany2013x(point, result);
            break;
        }
        default:
            logger.assert_debug(false, "Unknown problem");
    }
}

template <std::size_t DIM>
void SampleTestProblems<DIM>::exactSolution(
    const Geometry::PointPhysical<DIM> &point, double t,
    LinearAlgebra::SmallVector<DIM> &result) const {
    switch (problem_) {
        case CONSTANT: {
            result.set(1);
            break;
        }
        case LINEAR: {
            result.set(t);
            break;
        }
        case SINSIN: {
            sinx(point, result);
            result *= sin(M_PI * t);
            break;
        }
        case SARMANY2013: {
            sarmany2013x(point, result);
            result *= (cos(t) + cos(t / 2) + cos(t / 3));
            break;
        }
        default:
            logger.assert_debug(false, "Unknown problem");
    }
}

template <std::size_t DIM>
void SampleTestProblems<DIM>::exactSolutionCurl(
    const Geometry::PointPhysical<DIM> &point, double t,
    LinearAlgebra::SmallVector<DIM> &result) const {
    switch (problem_) {
        case CONSTANT:
        case LINEAR:
        case SINSIN: {
            result.set(0);
            break;
        }
        case SARMANY2013: {
            if (DIM == 3) {
                double x = point[0], y = point[1], z = point[2];
                result[0] = sin(M_PI * x) * (cos(M_PI * y) - cos(M_PI * z));
                result[1] = sin(M_PI * y) * (cos(M_PI * z) - cos(M_PI * x));
                result[2] = sin(M_PI * z) * (cos(M_PI * x) - cos(M_PI * y));

                result *= M_PI * (cos(t) + cos(t / 2) + cos(t / 3));
            } else {
                logger.assert_debug(false,
                                    "Sarmany problem is restricted to 3D.");
            }
            break;
        }
        default:
            logger.assert_debug(false, "Unknown problem");
    }
}

template <std::size_t DIM>
double SampleTestProblems<DIM>::referenceTimeBoundary() const {
    switch (problem_) {
        case CONSTANT:
            return 0.0;
        case LINEAR:
            return 1.0;
        case SINSIN:
            return 0.5;
        case SARMANY2013:
            return 0.0;
        default:
            logger.assert_debug(false, "Unknown problem");
            return 0;
    }
}

template <std::size_t DIM>
double SampleTestProblems<DIM>::timeScalingBoundary(double t) const {
    switch (problem_) {
        case CONSTANT:
            return 1;
        case LINEAR:
            return t;
        case SINSIN:
            return sin(M_PI * t);  // TODO: For the moment we assume 0 boundary
                                   // conditions
        case SARMANY2013:
            // Division by three to ensure that at 0 it is 1.
            return (cos(t) + cos(t / 2) + cos(t / 3)) / 3.0;
        default:
            logger.assert_debug(false, "Unknown problem");
            return 0;
    }
}

template <std::size_t DIM>
double SampleTestProblems<DIM>::timeScalingSource(double t) const {
    switch (problem_) {
        case CONSTANT:
            return 0;
        case LINEAR:
            return t;
        case SINSIN:
            return sin(M_PI * t);
        case SARMANY2013: {
            double p22 = 2 * M_PI * M_PI;
            return (p22 - 1) * cos(t) + (p22 - 1 / 4.0) * cos(t / 2) +
                   (p22 - 1 / 9.0) * cos(t / 3);
        }
        default:
            logger.assert_debug(false, "Unknown problem");
            return 0;
    }
}

template <std::size_t DIM>
void SampleTestProblems<DIM>::sinx(
    const Geometry::PointPhysical<DIM> &point,
    LinearAlgebra::SmallVector<DIM> &result) const {
    for (int i = 0; i < DIM; ++i) {
        double v = sin(M_PI * point[i]);
        v *= v;
        result[i] = v;
    }
}

template <std::size_t DIM>
void SampleTestProblems<DIM>::sarmany2013x(
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

template class SampleTestProblems<2>;
template class SampleTestProblems<3>;

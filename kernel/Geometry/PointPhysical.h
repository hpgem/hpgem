/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2014, University of Twente
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

#ifndef HPGEM_KERNEL_POINTPHYSICAL_H
#define HPGEM_KERNEL_POINTPHYSICAL_H

#include "Point.h"
#include "PointPhysicalBase.h"
#include <complex>
namespace Geometry {
template <std::size_t DIM>
class PointPhysical;

template <std::size_t DIM>
class PointPhysical : public Point<DIM>, public PointPhysicalBase {

   public:
    PointPhysical() : Point<DIM>() {}

    PointPhysical(const PointPhysical& p) : Point<DIM>(p) {}

    explicit PointPhysical(const Point<DIM>& p) : Point<DIM>(p) {}

    PointPhysical(const LinearAlgebra::SmallVector<DIM>& coord)
        : Point<DIM>(coord) {}

    PointPhysical(std::initializer_list<double> data) : Point<DIM>(data) {}

    PointPhysical operator*(double right) const {
        return PointPhysical(Point<DIM>::coordinates_ * right);
    }

    PointPhysical operator/(double right) const {
        return PointPhysical(Point<DIM>::coordinates_ / right);
    }

    // please note that for type-safety this function cannot be removed in favor
    // of the Point::operator+
    PointPhysical operator+(const PointPhysical& right) const {
        return PointPhysical(Point<DIM>::coordinates_ + right.coordinates_);
    }

    PointPhysical operator-(const PointPhysical& right) const {
        return PointPhysical(Point<DIM>::coordinates_ - right.coordinates_);
    }

    PointPhysical operator-() const { return *this * -1.; }

    PointPhysical& operator=(const PointPhysical& right) {
        Point<DIM>::coordinates_ = right.coordinates_;
        return *this;
    }

    void axpy(const double& alpha, const PointPhysical& x) {
        Point<DIM>::coordinates_.axpy(alpha, x.coordinates_);
    }

    const double* data() const { return Point<DIM>::coordinates_.data(); }

    std::size_t size() const override { return Point<DIM>::size(); }
};

template <std::size_t DIM>
PointPhysical<DIM> operator*(double left, const PointPhysical<DIM>& right);
}  // namespace Geometry

#include "PointPhysical_Impl.h"

#endif // HPGEM_KERNEL_POINTPHYSICAL_H

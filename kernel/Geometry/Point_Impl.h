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
#include <iostream>
#include "Point.h"
namespace hpgem {
namespace Geometry {
template <std::size_t DIM>
Point<DIM>::Point() : coordinates_() {}

template <std::size_t DIM>
Point<DIM>::Point(double coords[]) : coordinates_(coords) {}

template <std::size_t DIM>
Point<DIM>::Point(std::initializer_list<double> data) : coordinates_(data) {}

template <std::size_t DIM>
Point<DIM>::Point(const Point<DIM>& other) : coordinates_(other.coordinates_) {}

template <std::size_t DIM>
Point<DIM>::Point(const LinearAlgebra::SmallVector<DIM>& coord)
    : coordinates_(coord) {}

template <std::size_t DIM>
bool Point<DIM>::operator==(const Point<DIM>& right) const {
    return coordinates_ == right.coordinates_;
}

template <std::size_t DIM>
bool Point<DIM>::operator!=(const Point<DIM>& right) const {
    return !(*this == right);
}

template <std::size_t DIM>
bool Point<DIM>::operator<(const Point<DIM>& right) const {
    return coordinates_ < right.coordinates_;
}

template <std::size_t DIM>
Point<DIM>& Point<DIM>::operator+=(const Point<DIM>& right) {
    coordinates_ += right.coordinates_;
    return *this;
}

template <std::size_t DIM>
Point<DIM>& Point<DIM>::operator-=(const Point<DIM>& right) {
    coordinates_ -= right.coordinates_;
    return *this;
}

template <std::size_t DIM>
Point<DIM>& Point<DIM>::operator*=(double right) {
    coordinates_.operator*=(right);
    return *this;
}

template <std::size_t DIM>
Point<DIM> Point<DIM>::operator*(double right) const {
    return Point(coordinates_ * right);
}

template <std::size_t DIM>
Point<DIM> Point<DIM>::operator*(double right) {
    return Point(coordinates_ * right);
}

template <std::size_t DIM>
Point<DIM> Point<DIM>::operator+(const Point<DIM>& right) const {
    return Point(coordinates_ + right.coordinates_);
}

template <std::size_t DIM>
Point<DIM> Point<DIM>::operator+(const Point<DIM>& right) {
    return Point(coordinates_ + right.coordinates_);
}

template <std::size_t DIM>
Point<DIM> Point<DIM>::operator-(const Point<DIM>& right) const {
    return Point(coordinates_ - right.coordinates_);
}

template <std::size_t DIM>
Point<DIM> Point<DIM>::operator-(const Point<DIM>& right) {
    return Point(coordinates_ - right.coordinates_);
}

template <std::size_t DIM>
std::size_t Point<DIM>::size() {
    return coordinates_.size();
}

template <std::size_t DIM>
std::size_t Point<DIM>::size() const {
    return coordinates_.size();
}

template<std::size_t DIM>
double Point<DIM>::l2Norm() const {
    return coordinates_.l2Norm();
}

template<std::size_t DIM>
double Point<DIM>::l2NormSquared() const {
    return coordinates_.l2NormSquared();
}

template <std::size_t DIM>
double Point<DIM>::getCoordinate(std::size_t n) const {
    logger.assert_debug(n < size(),
                        "In Point::getCoordinate, entry % is requested while "
                        "the dimension is %",
                        n, size());
    return coordinates_[n];
}

template <std::size_t DIM>
const typename LinearAlgebra::SmallVector<DIM>& Point<DIM>::getCoordinates()
    const {
    return coordinates_;
}

template <std::size_t DIM>
void Point<DIM>::setCoordinate(std::size_t n, const double& coord) {
    logger.assert_debug(n < size(),
                        "In Point::setCoordinate, trying to set entry % while "
                        "the dimension is %",
                        n, size());
    coordinates_[n] = coord;
}

template <std::size_t DIM>
void Point<DIM>::setCoordinates(const LinearAlgebra::SmallVector<DIM>& coord) {
    coordinates_ = coord;
}

template <std::size_t DIM>
double& Point<DIM>::operator[](std::size_t n) {
    logger.assert_debug(
        n < size(),
        "In Point::operator[], entry % is requested while the dimension is %",
        n, size());
    return coordinates_[n];
}

template <std::size_t DIM>
const double& Point<DIM>::operator[](std::size_t n) const {
    logger.assert_debug(n < size(),
                        "In Point::operator[] const, entry % is requested "
                        "while the dimension is %",
                        n, size());
    return coordinates_[n];
}

template <std::size_t DIM>
Point<DIM>& Point<DIM>::operator=(const Point<DIM>& rhs) {
    coordinates_ = rhs.coordinates_;
    return *this;
}

template <std::size_t DIM>
std::ostream& operator<<(std::ostream& os, const Point<DIM>& point) {
    os << "point={";
    for (std::size_t i = 0; i < point.size(); i++) {
        if (i < point.size() - 1)
            os << point[i] << ',';
        else
            os << point[i];
    }
    os << "}";
    return os;
}

template <std::size_t DIM>
Point<DIM> operator*(const double& left, const Point<DIM>& right) {
    return Point<DIM>(right * left);
}
}  // namespace Geometry
}  // namespace hpgem
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

#ifndef HPGEM_KERNEL_SMALLVECTOR_H
#define HPGEM_KERNEL_SMALLVECTOR_H

#include "Logger.h"
#include <array>
#include <cmath>
#include "MiddleSizeVector.h"
#include <algorithm>
#include <numeric>

namespace hpgem {
namespace LinearAlgebra {
class MiddleSizeVector;
/// \class SmallVector
/// \brief This is a vector of doubles
///
/// \details
/// This implements a vector of doubles and all the standard operators for it.
/// Note it is encapulating a std::array for its data storage.
template <std::size_t numberOfRows>
class SmallVector {

   public:
    SmallVector() : data_() {}

    /*SmallVector(std::array<double, numberOfRows> t)
        : data_(t)
    {
    }*/

    SmallVector(const SmallVector& other) : data_(other.data_) {}

    // this constructor is implicit because both vector types should allow for
    // the same mathematical operations, with the only potential difference
    // being in the implementation
    SmallVector(const MiddleSizeVector& other) : data_() {
        logger.assert_debug(
            other.size() == numberOfRows,
            "Cannot construct a vector of size % from a vector of size %",
            numberOfRows, other.size());
        for (std::size_t i = 0; i < numberOfRows; ++i) {
            logger.assert_debug(std::abs(other[i] - std::real(other[i])) < 1e-9,
                                "trying to construct a real vector from a "
                                "vector with nonzero complex component");
            data_[i] = std::real(other[i]);
        }
    }

    SmallVector(const std::vector<double>& vec){
         logger.assert_debug(
            vec.size() == numberOfRows,
            "Cannot construct a vector of size % from a std::vector of size %",
            numberOfRows, vec.size());
        std::copy_n(vec.begin(), numberOfRows, data_.begin());
    }

    SmallVector(SmallVector&& other) : data_(std::move(other.data_)) {}

    SmallVector(const double array[]) : data_() {
        std::copy(array, array + numberOfRows, data_.begin());
    }

    SmallVector(std::initializer_list<double> data) : data_() {
        logger.assert_debug(data.size() == numberOfRows,
                            "provided array has size %, but should have size %",
                            data.size(), numberOfRows);
        std::copy(data.begin(), data.end(), data_.begin());
    }

    SmallVector& operator=(const SmallVector& right) {
        std::copy(right.data_.begin(), right.data_.end(), data_.begin());
        return *this;
    }

    SmallVector& operator=(const std::array<double, numberOfRows> l) {
        std::copy(l.begin(), l.end(), data_.begin());
        return *this;
    }

    SmallVector operator+(const SmallVector& right) const {
        SmallVector result;
        std::transform(data_.begin(), data_.end(), right.data_.begin(),
                       result.data_.begin(), std::plus<double>());
        return result;
    }

    SmallVector operator-(const SmallVector& right) const {
        SmallVector result;
        std::transform(data_.begin(), data_.end(), right.data_.begin(),
                       result.data_.begin(), std::minus<double>());
        return result;
    }

    SmallVector operator*(const double& right) const {
        SmallVector result;
        std::transform(
            data_.begin(), data_.end(), result.data_.begin(),
            std::bind(std::multiplies<double>(), std::placeholders::_1, right));
        return result;
    }

    /// Computes inner product between two vectors.
    double operator*(const SmallVector& right) const {
        return std::inner_product(right.data_.begin(), right.data_.end(),
                                  data_.begin(), 0.);
    }

    SmallVector& operator/=(const double& right) {
        std::transform(
            data_.begin(), data_.end(), data_.begin(),
            std::bind(std::divides<double>(), std::placeholders::_1, right));
        return *this;
    }

    SmallVector operator/(const double& right) const {
        SmallVector result;
        std::transform(
            data_.begin(), data_.end(), result.data_.begin(),
            std::bind(std::divides<double>(), std::placeholders::_1, right));
        return result;
    }

    void axpy(double a, const SmallVector& x) {
        for (std::size_t i = 0; i < numberOfRows; ++i) {
            data_[i] += a * x[i];
        }
    }

    /// This function is dangerous to use, since it compares doubles without
    /// a tolerance interval to see if they are equal.
    bool operator==(const SmallVector& right) const {
        for (std::size_t i = 0; i < numberOfRows; ++i) {
            if (data_[i] != right[i]) {
                return false;
            }
        }
        return true;
    }

    /// This function is dangerous to use, since it compares doubles without
    /// a tolerance interval to see if they are equal.
    bool operator<(const SmallVector& right) const {
        for (std::size_t i = 0; i < numberOfRows; ++i) {
            if (data_[i] < right[i]) {
                return true;
            }
            if (data_[i] > right[i]) {
                return false;
            }
        }
        return false;
    }

    SmallVector& operator+=(const SmallVector& right) {
        std::transform(data_.begin(), data_.end(), right.data_.begin(),
                       data_.begin(), std::plus<double>());
        return *this;
    }

    SmallVector& operator-=(const SmallVector& right) {
        std::transform(data_.begin(), data_.end(), right.data_.begin(),
                       data_.begin(), std::minus<double>());
        return *this;
    }

    SmallVector& operator*=(const double& right) {
        std::transform(
            data_.begin(), data_.end(), data_.begin(),
            std::bind(std::multiplies<double>(), std::placeholders::_1, right));
        return *this;
    }

    double& operator[](std::size_t n) {
        logger.assert_debug(n < numberOfRows,
                            "Requested entry %, but there are only % entries",
                            n, numberOfRows);
        return data_[n];
    }

    const double& operator[](std::size_t n) const {
        logger.assert_debug(n < numberOfRows,
                            "Requested entry %, but there are only % entries",
                            n, numberOfRows);
        return data_[n];
    }

    double& operator()(std::size_t n) {
        logger.assert_debug(n < numberOfRows,
                            "Requested entry %, but there are only % entries",
                            n, numberOfRows);
        return data_[n];
    }

    const double& operator()(std::size_t n) const {
        logger.assert_debug(n < numberOfRows,
                            "Requested entry %, but there are only % entries",
                            n, numberOfRows);
        return data_[n];
    }

    std::size_t size() const { return numberOfRows; }
    const double* data() const { return data_.data(); }

    double* data() { return data_.data(); }

    SmallVector operator-() const { return *this * -1.; }

    void set(double value) {
        for (std::size_t i = 0; i < numberOfRows; ++i) {
            data_[i] = value;
        }
    }

    /// \brief Cross product for 2D and 3D vectors.
    ///
    /// Computes the tradition cross product for a right handed coordinate
    /// system.
    ///
    /// Note, for 2D vectors this will return the value in the first
    /// component of the vector with the second one zero.
    /// \param other The second argument of the cross product.
    /// \param result A vector to store the result in.
    void crossProduct(const SmallVector& other, SmallVector& result) const;

    SmallVector crossProduct(const SmallVector& other) const {
        SmallVector result;
        crossProduct(other, result);
        return result;
    }

    double l2NormSquared() const { return this->operator*(*this); }

    double l2Norm() const { return std::sqrt(l2NormSquared()); }

    /// \brief Append a number to the vector.
    ///
    /// \param value The value to append
    /// \return A new vector with the value appended to the values in this
    /// vector.
    SmallVector<numberOfRows + 1> append(double value) {
        SmallVector<numberOfRows + 1> result;
        for (std::size_t i = 0; i < numberOfRows; ++i) result[i] = data_[i];
        result[numberOfRows] = value;
        return result;
    }

   private:
    std::array<double, numberOfRows> data_;
};

template <std::size_t numberOfRows>
SmallVector<numberOfRows> operator*(const double& left,
                                    const SmallVector<numberOfRows>& right) {
    return right * left;
}

template <std::size_t numberOfRows>
std::ostream& operator<<(std::ostream& os, const SmallVector<numberOfRows>& A) {
    os << "(";
    for (std::size_t i = 0; i < numberOfRows; ++i) {
        os << A[i] << " ";
    }
    os << ")";
    return os;
}

template <>
inline void SmallVector<2>::crossProduct(
    const LinearAlgebra::SmallVector<2>& other,
    LinearAlgebra::SmallVector<2>& result) const {
    result[0] = data_[0] * other[1] - data_[1] * other[0];
    result[1] = 0;
}

template <>
inline void SmallVector<3>::crossProduct(
    const LinearAlgebra::SmallVector<3>& other,
    LinearAlgebra::SmallVector<3>& result) const {
    result[0] = data_[1] * other[2] - data_[2] * other[1];
    result[1] = data_[2] * other[0] - data_[0] * other[2];
    result[2] = data_[0] * other[1] - data_[1] * other[0];
}

template <std::size_t numberOfRows>
inline void SmallVector<numberOfRows>::crossProduct(
    const LinearAlgebra::SmallVector<numberOfRows>& other,
    LinearAlgebra::SmallVector<numberOfRows>& result) const {
    logger.assert_debug(
        false, "Cross product only defined for 2 and 3D vectors, not for % ",
        numberOfRows);
}

#ifdef HPGEM_USE_COMPLEX_PETSC
template <std::size_t nRows>
MiddleSizeVector::MiddleSizeVector(const SmallVector<nRows>& other)
    : data_(nRows) {
    logger(WARN,
           "Constructing middle size vector from small vector, consider using "
           "small vectors everywhere for fixed length vectors of size <= 4");
    for (std::size_t i = 0; i < nRows; ++i) {
        logger.assert_debug(std::abs(other[i] - std::real(other[i])) < 1e-9,
                            "attempting to construct a real vector from a "
                            "vector with a complex part");
        data_[i] = std::real(other[i]);
    }
}
#else
template <std::size_t numberOfRows>
MiddleSizeVector::MiddleSizeVector(const SmallVector<numberOfRows>& other)
    : data_(other.data(), other.data() + numberOfRows) {
    logger(
        WARN,
        "Constructing middle size vector from small vector, consider "
        "using small vectors everywhere for fixed length vectors of size <= 4. "
        "This warning might also be caused because of an operation between a "
        "small vector and a middle size vector.");
}
#endif

}  // namespace LinearAlgebra
}  // namespace hpgem
#endif  // HPGEM_KERNEL_SMALLVECTOR_H

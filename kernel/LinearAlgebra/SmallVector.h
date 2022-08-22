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

#ifndef HPGEM_KERNEL_GSmallVector_H
#define HPGEM_KERNEL_GSmallVector_H

// Note: there is an depenency between MiddleSizeVector and SmallVector to allow
// conversion.
#include "MiddleSizeVector.h"

#include "Logger.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <numeric>

namespace hpgem {
namespace LinearAlgebra {

namespace Detail {
// Several helper functions for GSmallVector.
//   - CrossProduct requires partial specialization
//
template <typename T>
void computeCrossProduct(const GSmallVector<2, T>&, const GSmallVector<2, T>&,
                         GSmallVector<2, T>&);
template <typename T>
void computeCrossProduct(const GSmallVector<3, T>&, const GSmallVector<3, T>&,
                         GSmallVector<3, T>&);
template <typename std::size_t nRows, typename T>
void computeCrossProduct(const GSmallVector<nRows, T>&,
                         const GSmallVector<nRows, T>&,
                         GSmallVector<nRows, T>&);

template <typename T1, typename T2>
struct ResultT {
    using Type = std::complex<double>;
};
template <>
struct ResultT<double, double> {
    using Type = double;
};

// Conversion operators to allow converting vectors especially needed for
// MiddleSizeVector conversions.
// Note: using an out parameter to allow overloading
template <typename T>
inline void convert(const T& in, T& out) {
    out = in;
}
template <typename T>
inline void convert(const T& in, std::complex<T>& out) {
    out = in;
}
template <typename T>
inline void convert(const std::complex<T>& in, T& out) {
    logger.assert_debug(std::abs(in.imag()) < 1e-9,
                        "Converting complex number with non zero complex part "
                        "into a real one.");
    out = in.real();
}

inline double conj(const double& v) { return v; }
inline std::complex<double> conj(const std::complex<double>& v) {
    return std::conj(v);
}

inline std::size_t hash(double v) noexcept { return std::hash<double>()(v); }
inline std::size_t hash(std::complex<double> v) noexcept {
    // Same as
    // https://stackoverflow.com/questions/1646807/quick-and-simple-hash-code-combinations
    std::size_t hash = 17;
    hash = hash * 31 + std::hash<double>()(v.real());
    hash = hash * 31 + std::hash<double>()(v.imag());

    return hash;
}

inline bool less(const double& v1, const double& v2) noexcept {
    return v1 < v2;
}
inline bool less(const std::complex<double>& v1,
                 const std::complex<double>& v2) noexcept {
    if (v1.real() != v2.real()) {
        return v1.real() < v2.real();
    } else {
        return v1.imag() < v2.imag();
    }
}

}  // namespace Detail

// Note specialization to SmallVector/SmallVectorC in Types.h

/// \class GSmallVector
/// \brief This is a fixed size vector of generic type (supported: double,
/// std::complex<double>)
///
/// \tparam numberOfRows The number of rows or columns in this vector
/// \tparam EntryT The value type
template <std::size_t numberOfRows, typename EntryT>
class GSmallVector {

   public:
    // This will actually zero the memory, unlike =default;
    GSmallVector() : data_(){};

    /*GSmallVector(std::array<double, numberOfRows> t)
        : data_(t)
    {
    }*/

    GSmallVector(const GSmallVector& other) : data_(other.data_) {}

    template <typename EntryTOther>
    GSmallVector(const GSmallVector<numberOfRows, EntryTOther>& other)
        : data_() {
        for (std::size_t i = 0; i < numberOfRows; ++i) {
            Detail::convert(other[i], data_[i]);
        }
    }

    // this constructor is implicit because both vector types should allow for
    // the same mathematical operations, with the only potential difference
    // being in the implementation
    GSmallVector(const MiddleSizeVector& other) : data_() {
        logger.assert_debug(
            other.size() == numberOfRows,
            "Cannot construct a vector of size % from a vector of size %",
            numberOfRows, other.size());
        for (std::size_t i = 0; i < numberOfRows; ++i) {
            Detail::convert(other[i], data_[i]);
        }
    }

    GSmallVector(const std::vector<EntryT>& vec) {
        logger.assert_debug(
            vec.size() == numberOfRows,
            "Cannot construct a vector of size % from a std::vector of size %",
            numberOfRows, vec.size());
        std::copy_n(vec.begin(), numberOfRows, data_.begin());
    }

    GSmallVector(GSmallVector&& other) : data_(std::move(other.data_)) {}

    GSmallVector(const EntryT array[]) : data_() {
        std::copy(array, array + numberOfRows, data_.begin());
    }

    GSmallVector(std::initializer_list<EntryT> data) : data_() {
        logger.assert_debug(data.size() == numberOfRows,
                            "provided array has size %, but should have size %",
                            data.size(), numberOfRows);
        std::copy(data.begin(), data.end(), data_.begin());
    }

    static GSmallVector constant(EntryT entry) {
        GSmallVector result;
        result.set(entry);
        return result;
    }

    GSmallVector& operator=(const GSmallVector& right) {
        std::copy(right.data_.begin(), right.data_.end(), data_.begin());
        return *this;
    }

    GSmallVector& operator=(const std::array<EntryT, numberOfRows> l) {
        std::copy(l.begin(), l.end(), data_.begin());
        return *this;
    }

    GSmallVector operator+(const GSmallVector& right) const {
        GSmallVector result;
        std::transform(data_.begin(), data_.end(), right.data_.begin(),
                       result.data_.begin(), std::plus<EntryT>());
        return result;
    }

    GSmallVector operator-(const GSmallVector& right) const {
        GSmallVector result;
        std::transform(data_.begin(), data_.end(), right.data_.begin(),
                       result.data_.begin(), std::minus<EntryT>());
        return result;
    }

    template <typename EntryT2>
    GSmallVector<numberOfRows, typename Detail::ResultT<EntryT, EntryT2>::Type>
        operator*(const EntryT2& right) const {
        GSmallVector<numberOfRows,
                     typename Detail::ResultT<EntryT, EntryT2>::Type>
            result;
        for (std::size_t i = 0; i < numberOfRows; ++i) {
            result[i] = right * data_[i];
        }
        return result;
    }

    /// Computes inner product between two vectors.
    template <typename EntryT2>
    auto operator*(const GSmallVector<numberOfRows, EntryT2>& right) const {
        typename Detail::ResultT<EntryT, EntryT2>::Type result = 0.0;
        for (std::size_t i = 0; i < numberOfRows; ++i) {
            result += data_[i] * Detail::conj(right[i]);
        }
        return result;
    }

    GSmallVector& operator/=(const EntryT& right) {
        std::transform(
            data_.begin(), data_.end(), data_.begin(),
            std::bind(std::divides<EntryT>(), std::placeholders::_1, right));
        return *this;
    }

    GSmallVector operator/(const EntryT& right) const {
        GSmallVector result;
        std::transform(
            data_.begin(), data_.end(), result.data_.begin(),
            std::bind(std::divides<EntryT>(), std::placeholders::_1, right));
        return result;
    }

    void axpy(EntryT a, const GSmallVector& x) {
        for (std::size_t i = 0; i < numberOfRows; ++i) {
            data_[i] += a * x[i];
        }
    }

    /// This function is dangerous to use, since it compares doubles without
    /// a tolerance interval to see if they are equal.
    bool operator==(const GSmallVector& right) const {
        for (std::size_t i = 0; i < numberOfRows; ++i) {
            if (data_[i] != right[i]) {
                return false;
            }
        }
        return true;
    }

    /// Lexicographic ordering of vectors, note that it does not use any
    /// tolerances.
    bool operator<(const GSmallVector& right) const {
        for (std::size_t i = 0; i < numberOfRows; ++i) {
            if (data_[i] != right.data_[i]) {
                return Detail::less(data_[i], right.data_[i]);
            }
        }
        // Equals
        return false;
    }

    GSmallVector& operator+=(const GSmallVector& right) {
        std::transform(data_.begin(), data_.end(), right.data_.begin(),
                       data_.begin(), std::plus<EntryT>());
        return *this;
    }

    GSmallVector& operator-=(const GSmallVector& right) {
        std::transform(data_.begin(), data_.end(), right.data_.begin(),
                       data_.begin(), std::minus<EntryT>());
        return *this;
    }

    GSmallVector& operator*=(const EntryT& right) {
        std::transform(
            data_.begin(), data_.end(), data_.begin(),
            std::bind(std::multiplies<EntryT>(), std::placeholders::_1, right));
        return *this;
    }

    EntryT& operator[](std::size_t n) {
        logger.assert_debug(n < numberOfRows,
                            "Requested entry %, but there are only % entries",
                            n, numberOfRows);
        return data_[n];
    }

    const EntryT& operator[](std::size_t n) const {
        logger.assert_debug(n < numberOfRows,
                            "Requested entry %, but there are only % entries",
                            n, numberOfRows);
        return data_[n];
    }

    EntryT& operator()(std::size_t n) {
        logger.assert_debug(n < numberOfRows,
                            "Requested entry %, but there are only % entries",
                            n, numberOfRows);
        return data_[n];
    }

    const EntryT& operator()(std::size_t n) const {
        logger.assert_debug(n < numberOfRows,
                            "Requested entry %, but there are only % entries",
                            n, numberOfRows);
        return data_[n];
    }

    std::size_t size() const { return numberOfRows; }
    const EntryT* data() const { return data_.data(); }

    EntryT* data() { return data_.data(); }

    GSmallVector operator-() const { return *this * -1.; }

    void set(EntryT value) {
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
    void crossProduct(const GSmallVector<numberOfRows, EntryT>& other,
                      GSmallVector<numberOfRows, EntryT>& result) const {
        Detail::computeCrossProduct(*this, other, result);
    }

    GSmallVector crossProduct(
        const GSmallVector<numberOfRows, EntryT>& other) const {
        GSmallVector result;
        crossProduct(other, result);
        return result;
    }

    // For implementation of the cross product without having to specialize
    // for N=2 and N=3
    template <typename T>
    friend void Detail::computeCrossProduct(const GSmallVector<2, T>&,
                                            const GSmallVector<2, T>&,
                                            GSmallVector<2, T>&);
    template <typename T>
    friend void Detail::computeCrossProduct(const GSmallVector<3, T>&,
                                            const GSmallVector<3, T>&,
                                            GSmallVector<3, T>&);
    template <std::size_t n, typename T>
    friend void Detail::computeCrossProduct(const GSmallVector<n, T>&,
                                            const GSmallVector<n, T>&,
                                            GSmallVector<n, T>&);

    double l2NormSquared() const {
        double result = 0.0;
        for (std::size_t i = 0; i < numberOfRows; ++i) {
            result += std::norm(data_[i]);
        }
        return result;
    }

    double l2Norm() const { return std::sqrt(l2NormSquared()); }

    /// \brief Append a number to the vector.
    ///
    /// \param value The value to append
    /// \return A new vector with the value appended to the values in this
    /// vector.
    GSmallVector<numberOfRows + 1, EntryT> append(EntryT value) {
        GSmallVector<numberOfRows + 1, EntryT> result;
        for (std::size_t i = 0; i < numberOfRows; ++i) result[i] = data_[i];
        result[numberOfRows] = value;
        return result;
    }

    std::size_t hash() const noexcept {
        // Similar to
        // https://stackoverflow.com/questions/1646807/quick-and-simple-hash-code-combinations

        std::size_t hash = 17;
        for (std::size_t i = 0; i < numberOfRows; ++i) {
            // Basic hash combining function
            std::size_t v = Detail::hash(data_[i]);
            hash = hash * 31 + v;
        }
        return hash;
    }

    SmallVector<numberOfRows> real() const {
        SmallVector<numberOfRows> res;
        for (std::size_t i = 0; i < numberOfRows; ++i) {
            res[i] = std::real(data_[i]);
        }
        return res;
    }

    SmallVector<numberOfRows> imag() const {
        SmallVector<numberOfRows> res;
        for (std::size_t i = 0; i < numberOfRows; ++i) {
            res[i] = std::imag(data_[i]);
        }
        return res;
    }

    GSmallVector conj() const {
        GSmallVector res;
        for (std::size_t i = 0; i < numberOfRows; ++i) {
            res.data_[i] = Detail::conj(data_[i]);
        }
        return res;
    }

   private:
    std::array<EntryT, numberOfRows> data_;
};

// Multiplication of number by vector
template <std::size_t numberOfRows, typename EntryT>
GSmallVector<numberOfRows, typename Detail::ResultT<EntryT, double>::Type>
    operator*(const double& left,
              const GSmallVector<numberOfRows, EntryT>& right) {
    return right * left;
}

template <std::size_t numberOfRows, typename EntryT>
GSmallVector<numberOfRows, std::complex<double>> operator*(
    const std::complex<double>& left,
    const GSmallVector<numberOfRows, EntryT>& right) {
    return right * left;
}

/// 2D cross product for the double cross product of the form
/// A x (B x C) with left = A, right = B x C (possibly B = nabla)
/// Needed due to storing the cross product as 'x' entry
template <typename EntryT>
GSmallVector<2, EntryT> leftDoubledCrossProduct(
    const GSmallVector<2, EntryT>& left, const GSmallVector<2, EntryT>& right) {
    GSmallVector<2, EntryT> result;
    result[0] = left[1] * right[0];
    result[1] = -left[0] * right[0];
    return result;
}
template <std::size_t numberOfRows, typename EntryT>
GSmallVector<numberOfRows, EntryT> leftDoubledCrossProduct(
    const GSmallVector<numberOfRows, EntryT>& left,
    const GSmallVector<numberOfRows, EntryT>& right) {
    return left.crossProduct(right);
}

template <std::size_t numberOfRows, typename EntryT>
std::ostream& operator<<(std::ostream& os,
                         const GSmallVector<numberOfRows, EntryT>& A) {
    os << "(";
    for (std::size_t i = 0; i < numberOfRows; ++i) {
        os << A[i] << " ";
    }
    os << ")";
    return os;
}

namespace Detail {

template <typename T>
void computeCrossProduct(const GSmallVector<2, T>& v1,
                         const GSmallVector<2, T>& v2,
                         GSmallVector<2, T>& result) {
    result[0] = v1.data_[0] * v2.data_[1] - v1.data_[1] * v2.data_[0];
    result[1] = 0.0;
}
template <typename T>
void computeCrossProduct(const GSmallVector<3, T>& v1,
                         const GSmallVector<3, T>& v2,
                         GSmallVector<3, T>& result) {
    result[0] = v1.data_[1] * v2.data_[2] - v1.data_[2] * v2.data_[1];
    result[1] = v1.data_[2] * v2.data_[0] - v1.data_[0] * v2.data_[2];
    result[2] = v1.data_[0] * v2.data_[1] - v1.data_[1] * v2.data_[0];
}
template <typename std::size_t nRows, typename T>
void computeCrossProduct(const GSmallVector<nRows, T>&,
                         const GSmallVector<nRows, T>&,
                         GSmallVector<nRows, T>&) {
    logger.assert_debug(
        false, "Cross product only defined for 2 and 3D vectors, not for % ",
        nRows);
}

}  // namespace Detail

template <std::size_t nRows, typename EntryT>
MiddleSizeVector::MiddleSizeVector(const GSmallVector<nRows, EntryT>& other)
    : data_(nRows) {
    logger(WARN,
           "Constructing middle size vector from small vector, consider using "
           "small vectors everywhere for fixed length vectors of size <= 4");
    for (std::size_t i = 0; i < nRows; ++i) {
        Detail::convert(other[i], data_[i]);
    }
}

// Standard instances
extern template class GSmallVector<0, double>;
extern template class GSmallVector<1, double>;
extern template class GSmallVector<2, double>;
extern template class GSmallVector<3, double>;
extern template class GSmallVector<4, double>;

extern template class GSmallVector<0, std::complex<double>>;
extern template class GSmallVector<1, std::complex<double>>;
extern template class GSmallVector<2, std::complex<double>>;
extern template class GSmallVector<3, std::complex<double>>;
extern template class GSmallVector<4, std::complex<double>>;

}  // namespace LinearAlgebra
}  // namespace hpgem

namespace std {

template <std::size_t numberOfRows>
struct hash<hpgem::LinearAlgebra::GSmallVector<numberOfRows, double>> {
    std::size_t operator()(
        const hpgem::LinearAlgebra::GSmallVector<numberOfRows, double>& v)
        const noexcept {
        return v.hash();
    }
};

}  // namespace std

#endif  // HPGEM_KERNEL_GSmallVector_H

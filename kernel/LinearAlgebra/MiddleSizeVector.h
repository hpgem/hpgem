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

#ifndef MIDDLESIZEVECTOR_H_
#define MIDDLESIZEVECTOR_H_

#include "Logger.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <complex>

namespace LinearAlgebra {
template <std::size_t numberOfRows>
class SmallVector;

/// \class MiddleSizeVector
/// \brief This is a vector of doubles
///
/// \details
/// This implements a vector of doubles and all the standard operators for it.
class MiddleSizeVector {

   public:
    // we will need an extra macro hpgem_use_complex_numbers if complex numbers
    // get different use cases
#ifdef HPGEM_USE_COMPLEX_PETSC
    using type = std::complex<double>;
#else
    using type = double;
#endif

    MiddleSizeVector();

    explicit MiddleSizeVector(std::size_t m);

    MiddleSizeVector(std::initializer_list<type> t);

    // this constructor is implicit because both vector types should allow for
    // the same mathematical operations, with the only potential difference
    // being in the implementation
    MiddleSizeVector(const MiddleSizeVector& other);

    MiddleSizeVector(MiddleSizeVector&& other);

    // implemented with SmallVector for dependency reasons
    template <std::size_t numberOfRows>
    MiddleSizeVector(const SmallVector<numberOfRows>& other);

    MiddleSizeVector(const type array[], std::size_t size);

    void resize(std::size_t size);

    MiddleSizeVector& operator=(const MiddleSizeVector& right);

    MiddleSizeVector& operator=(const std::initializer_list<type> l);

    MiddleSizeVector operator+(const MiddleSizeVector& right) const;

    MiddleSizeVector operator-(const MiddleSizeVector& right) const;

    MiddleSizeVector operator*(const type& right) const;

    /// Computes inner product between two vectors.
    type operator*(const MiddleSizeVector& right) const;

    MiddleSizeVector& operator/=(const type& right);

    MiddleSizeVector operator/(const type& right) const;

    void axpy(type a, const MiddleSizeVector& x);

    /// This function is dangerous to use, since it compares doubles without
    /// a tolerance interval to see if they are equal.
    bool operator==(const MiddleSizeVector& right) const;

    /// This function is dangerous to use, since it compares doubles without
    /// a tolerance interval to see if they are equal.
    bool operator<(const MiddleSizeVector& right) const;

    MiddleSizeVector& operator+=(const MiddleSizeVector& right);

    MiddleSizeVector& operator-=(const MiddleSizeVector& right);

    MiddleSizeVector& operator*=(const type& right);

    type& operator[](std::size_t n);

    const type& operator[](std::size_t n) const {
        logger.assert_debug(n < data_.size(),
                            "Requested entry %, but there are only % entries",
                            n, data_.size());
        return data_[n];
    }

    type& operator()(std::size_t n) {
        logger.assert_debug(n < data_.size(),
                            "Requested entry %, but there are only % entries",
                            n, data_.size());
        return data_[n];
    }

    const type& operator()(std::size_t n) const {
        logger.assert_debug(n < data_.size(),
                            "Requested entry %, but there are only % entries",
                            n, data_.size());
        return data_[n];
    }

    std::size_t size() const { return data_.size(); }

    const type* data() const { return data_.data(); }

    type* data() { return data_.data(); }

   private:
    std::vector<type> data_;
};

MiddleSizeVector operator*(const MiddleSizeVector::type& left,
                           const MiddleSizeVector& right);

MiddleSizeVector operator-(const MiddleSizeVector& right);

std::ostream& operator<<(std::ostream& os, const MiddleSizeVector& A);

}  // namespace LinearAlgebra
#endif

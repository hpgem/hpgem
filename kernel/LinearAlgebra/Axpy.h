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

#ifndef HPGEM_KERNEL_AXPY_H
#define HPGEM_KERNEL_AXPY_H
#include "MiddleSizeMatrix.h"
#include "MiddleSizeVector.h"

#include <type_traits>

namespace LinearAlgebra {
///\brief Implementation for the AXPY operation, i.e. Y = alpha * X + Y.
/// If you need a variant that's not working on a scalar or class that
/// implemented axpy itself, please implement it in this file.

//!\brief axpy operation for scalar values, for example a double or an int.
// for typical use cases S is most likely to be double or std::complex<double>
template <typename T, typename S>
inline typename std::enable_if<std::is_arithmetic<T>::value, void>::type axpy(
    S alpha, const T& x, T& y) {
    y += alpha * x;
}

//! \brief axpy operation for classes that have implemented axpy themselves.
/*! Please note that the syntax here is different than it is normally for axpy.
 * If a class has not implemented axpy, the compiler will automatically look
 * for another option of this templated function.
 * The user is responsible for checking if x and y have the same size.
 */
template <typename T, typename S>
inline typename std::enable_if<std::is_class<T>::value, void>::type axpy(
    S a, const T& x, T& y) {
    y.axpy(a, x);
}

template <typename S>
inline void axpy(S a, const std::complex<double>& x, std::complex<double>& y) {
    y += a * x;
}

}  // namespace LinearAlgebra

#endif  // HPGEM_KERNEL_AXPY_H

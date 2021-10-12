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

#include "SmallMatrix.h"

namespace hpgem {

namespace LinearAlgebra {

template <std::size_t numberOfRows, std::size_t numberOfColumns>
SmallVector<numberOfRows> SmallMatrix<numberOfRows, numberOfColumns>::operator*(
    SmallVector<numberOfColumns>& right) {

    typename SmallVector<numberOfColumns>::EigenType& rightData = right;
    typename SmallVector<numberOfRows>::EigenType res = data2_ * rightData;
    return res;
}

template <std::size_t numberOfRows, std::size_t numberOfColumns>
SmallVector<numberOfRows> SmallMatrix<numberOfRows, numberOfColumns>::operator*(
    SmallVector<numberOfColumns>& right) const {

    typename SmallVector<numberOfColumns>::EigenType& rightData = right;
    typename SmallVector<numberOfRows>::EigenType res = data2_ * rightData;
    return res;
}

template <std::size_t numberOfRows, std::size_t numberOfColumns>
template <std::size_t K>
SmallMatrix<numberOfRows, K>
    SmallMatrix<numberOfRows, numberOfColumns>::operator*(
        const SmallMatrix<numberOfColumns, K>& other) {

    const typename SmallMatrix<numberOfColumns, K>::EigenType& otherMat = other;
    typename SmallMatrix<numberOfRows, K>::EigenType res;
    res = data2_ * otherMat;
    return res;
}

template <std::size_t numberOfRows, std::size_t numberOfColumns>
template <std::size_t K>
SmallMatrix<numberOfRows, K>
    SmallMatrix<numberOfRows, numberOfColumns>::operator*(
        const SmallMatrix<numberOfColumns, K>& other) const {

    const typename SmallMatrix<numberOfColumns, K>::EigenType& otherMat = other;
    typename SmallMatrix<numberOfRows, K>::EigenType res;
    res = data2_ * otherMat;
    return res;
}

template <std::size_t numberOfRows, std::size_t numberOfColumns>
SmallMatrix<numberOfRows, numberOfColumns>&
    SmallMatrix<numberOfRows, numberOfColumns>::operator*=(
        const SmallMatrix<numberOfColumns, numberOfColumns>& other) {
    data2_ *= other.data2_;
    return (*this);
}

template <std::size_t numberOfRows, std::size_t numberOfColumns>
SmallVector<numberOfRows> SmallMatrix<
    numberOfRows, numberOfColumns>::computeWedgeStuffVector() const {
    // copied from MiddleSizeMatrix to prevent constructing a temporary
    // MiddleSizeMatrix
    logger.assert_debug(
        numberOfColumns == numberOfRows - 1,
        "Matrix has wrong dimensions to construct the wedge stuff vector");
    SmallVector<numberOfRows> result;

    switch (numberOfRows) {
        case 2:
            result[0] = -(*this)(1, 0);
            result[1] = +(*this)(0, 0);
            break;
        case 3:
            result[0] =
                (*this)(1, 0) * (*this)(2, 1) - (*this)(2, 0) * (*this)(1, 1);
            result[1] =
                (*this)(0, 1) * (*this)(2, 0) -
                (*this)(0, 0) * (*this)(2, 1);  // includes minus sign already!
            result[2] =
                (*this)(0, 0) * (*this)(1, 1) - (*this)(1, 0) * (*this)(0, 1);
            break;
        case 4:
            result[0] = (*this)(1, 0) * (-(*this)(2, 1) * (*this)(3, 2) +
                                         (*this)(3, 1) * (*this)(2, 2)) +
                        (*this)(2, 0) * ((*this)(1, 1) * (*this)(3, 2) -
                                         (*this)(3, 1) * (*this)(1, 2)) +
                        (*this)(3, 0) * (-(*this)(1, 1) * (*this)(2, 2) +
                                         (*this)(2, 1) * (*this)(1, 2));

            result[1] = (*this)(0, 0) * ((*this)(2, 1) * (*this)(3, 2) -
                                         (*this)(3, 1) * (*this)(2, 2)) +
                        (*this)(2, 0) * (-(*this)(0, 1) * (*this)(3, 2) +
                                         (*this)(3, 1) * (*this)(0, 2)) +
                        (*this)(3, 0) * ((*this)(0, 1) * (*this)(2, 2) -
                                         (*this)(2, 1) * (*this)(0, 2));
            result[2] = (*this)(0, 0) * (-(*this)(1, 1) * (*this)(3, 2) +
                                         (*this)(3, 1) * (*this)(1, 2)) +
                        (*this)(1, 0) * ((*this)(0, 1) * (*this)(3, 2) -
                                         (*this)(3, 1) * (*this)(0, 2)) +
                        (*this)(3, 0) * (-(*this)(0, 1) * (*this)(1, 2) +
                                         (*this)(1, 1) * (*this)(0, 2));
            result[3] = (*this)(0, 0) * ((*this)(1, 1) * (*this)(2, 2) -
                                         (*this)(2, 1) * (*this)(1, 2)) +
                        (*this)(1, 0) * (-(*this)(0, 1) * (*this)(2, 2) +
                                         (*this)(2, 1) * (*this)(0, 2)) +
                        (*this)(2, 0) * ((*this)(0, 1) * (*this)(1, 2) -
                                         (*this)(1, 1) * (*this)(0, 2));
            break;
        default:
            logger(ERROR, "Wedge product not implemented for this dimension");
    }  // end switch

    return (result);
}

// class template specialization for this one function is a waste of code
// duplication just let the compiler figure out which case it needs
template <std::size_t numberOfRows, std::size_t numberOfColumns>
double SmallMatrix<numberOfRows, numberOfColumns>::determinant() const {
    return data2_.determinant();
}

template <std::size_t numberOfRows, std::size_t numberOfColumns>
SmallMatrix<numberOfRows, numberOfColumns>
    SmallMatrix<numberOfRows, numberOfColumns>::inverse() const {
    EigenType res = data2_.inverse();
    return res;
}

template <std::size_t numberOfRows, std::size_t numberOfColumns>
template <std::size_t numberOfRightHandSideColumns>
void SmallMatrix<numberOfRows, numberOfColumns>::solve(
    SmallMatrix<numberOfRows, numberOfRightHandSideColumns>& B) const {
    logger.assert_debug(numberOfRows == numberOfColumns,
                        "can only solve for square matrixes");

    const typename SmallMatrix<
        numberOfRows, numberOfRightHandSideColumns>::EigenType& rawVec = B;
    typename SmallMatrix<numberOfRows, numberOfRightHandSideColumns>::EigenType
        res;
    res = data2_.fullPivLu().solve(rawVec);
    // To prevent aliasing
    B = res;
}

template <std::size_t numberOfRows, std::size_t numberOfColumns>
void SmallMatrix<numberOfRows, numberOfColumns>::solve(
    SmallVector<numberOfRows>& b) const {
    logger.assert_debug(numberOfRows == numberOfColumns,
                        "can only solve for square matrixes");

    const typename SmallVector<numberOfRows>::EigenType& rawVec = b;
    typename SmallVector<numberOfRows>::EigenType res;
    // To prevent aliasing
    res = data2_.fullPivLu().solve(rawVec);
    b = res;
}

template <std::size_t numberOfRows, std::size_t numberOfColumns>
SmallVector<numberOfColumns> operator*(
    SmallVector<numberOfRows>& vec,
    SmallMatrix<numberOfRows, numberOfColumns>& mat) {
    const typename SmallVector<numberOfRows>::EigenType& rawVec = vec;
    const typename SmallMatrix<numberOfRows, numberOfColumns>::EigenType
        rawMat = mat;
    typename SmallVector<numberOfColumns>::EigenType res =
        rawVec.adjoint() * rawMat;
    return res;
}
}  // namespace LinearAlgebra
}  // namespace hpgem

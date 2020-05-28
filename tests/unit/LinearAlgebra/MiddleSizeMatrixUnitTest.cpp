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
#include <cstdlib>
#include <cmath>
#include "LinearAlgebra/MiddleSizeMatrix.h"
#include "LinearAlgebra/MiddleSizeVector.h"
#include "Logger.h"

using LinearAlgebra::MiddleSizeMatrix;
using LinearAlgebra::MiddleSizeVector;

int main(int argc, char** argv) {
    // constructors
    LinearAlgebra::MiddleSizeVector vec0{0., 1.}, vec1{2., 3.}, vec2{4., 5.},
        vec3{6., 7.};
    MiddleSizeMatrix A0, A22(2, 2), A23(2, 3), A32(3, 2), destroy(3, 3, 1),
        count0({vec0, vec1}), count1({vec2, vec3}), copy(count0), bla({count0}),
        merge({count0, count1});
    logger.assert_always(destroy.getNumberOfRows() == 3, "Rows in a matrix");
    logger.assert_always(destroy.getNumberOfColumns() == 3,
                         "Columns in a matrix");
    logger.assert_always(destroy.size() == 9, "Size of a matrix");
    for (std::size_t i = 0; i < destroy.size(); ++i) {
        logger.assert_always(std::abs(destroy[i] - 1.) < 1e-12,
                             "Entry of a matrix");
    }
    for (std::size_t i = 0; i < destroy.getNumberOfRows(); ++i) {
        for (std::size_t j = 0; j < destroy.getNumberOfColumns(); ++j) {
            logger.assert_always(std::abs(destroy(i, j) - 1.) < 1e-12,
                                 "Entry of a matrix");
        }
    }
    MiddleSizeMatrix moved(std::move(destroy));
    logger.assert_always(moved.getNumberOfRows() == 3, "Rows in a matrix");
    logger.assert_always(moved.getNumberOfColumns() == 3,
                         "Columns in a matrix");
    logger.assert_always(moved.size() == 9, "Size of a matrix");
    for (std::size_t i = 0; i < moved.size(); ++i) {
        logger.assert_always(std::abs(moved[i] - 1.) < 1e-12,
                             "Entry of a matrix");
    }
    for (std::size_t i = 0; i < moved.getNumberOfRows(); ++i) {
        for (std::size_t j = 0; j < moved.getNumberOfColumns(); ++j) {
            logger.assert_always(std::abs(moved(i, j) - 1.) < 1e-12,
                                 "Entry of a matrix");
        }
    }
    logger.assert_always(A0.getNumberOfRows() == 0, "Rows in a matrix");
    logger.assert_always(A0.getNumberOfColumns() == 0, "Columns in a matrix");
    logger.assert_always(A0.size() == 0, "Size of a matrix");
    logger.assert_always(A22.getNumberOfRows() == 2, "Rows in a matrix");
    logger.assert_always(A22.getNumberOfColumns() == 2, "Columns in a matrix");
    logger.assert_always(A22.size() == 4, "Size of a matrix");
    logger.assert_always(A23.getNumberOfRows() == 2, "Rows in a matrix");
    logger.assert_always(A23.getNumberOfColumns() == 3, "Columns in a matrix");
    logger.assert_always(A23.size() == 6, "Size of a matrix");
    logger.assert_always(A32.getNumberOfRows() == 3, "Rows in a matrix");
    logger.assert_always(A32.getNumberOfColumns() == 2, "Columns in a matrix");
    logger.assert_always(A32.size() == 6, "Size of a matrix");
    logger.assert_always(count0.getNumberOfRows() == 2, "Rows in a matrix");
    logger.assert_always(count0.getNumberOfColumns() == 2,
                         "Columns in a matrix");
    logger.assert_always(count0.size() == 4, "Size of a matrix");
    for (std::size_t i = 0; i < count0.size(); ++i) {
        logger.assert_always(std::abs(count0[i] - double(i)) < 1e-12,
                             "Entry of a matrix");
    }
    logger.assert_always(std::abs(count0(0, 0) - 0.) < 1e-12,
                         "Entry of a matrix");
    logger.assert_always(std::abs(count0(1, 0) - 1.) < 1e-12,
                         "Entry of a matrix");
    logger.assert_always(std::abs(count0(0, 1) - 2.) < 1e-12,
                         "Entry of a matrix");
    logger.assert_always(std::abs(count0(1, 1) - 3.) < 1e-12,
                         "Entry of a matrix");
    logger.assert_always(count1.getNumberOfRows() == 2, "Rows in a matrix");
    logger.assert_always(count1.getNumberOfColumns() == 2,
                         "Columns in a matrix");
    logger.assert_always(count1.size() == 4, "Size of a matrix");
    for (std::size_t i = 0; i < count1.size(); ++i) {
        logger.assert_always(std::abs(count1[i] - 4. - double(i)) < 1e-12,
                             "Entry of a matrix");
    }
    logger.assert_always(copy.getNumberOfRows() == 2, "Rows in a matrix");
    logger.assert_always(copy.getNumberOfColumns() == 2, "Columns in a matrix");
    logger.assert_always(copy.size() == 4, "Size of a matrix");
    for (std::size_t i = 0; i < copy.size(); ++i) {
        logger.assert_always(std::abs(copy[i] - double(i)) < 1e-12,
                             "Entry of a matrix");
    }
    logger.assert_always(bla.getNumberOfRows() == 2, "Rows in a matrix");
    logger.assert_always(bla.getNumberOfColumns() == 2, "Columns in a matrix");
    logger.assert_always(bla.size() == 4, "Size of a matrix");
    for (std::size_t i = 0; i < bla.size(); ++i) {
        logger.assert_always(std::abs(bla[i] - double(i)) < 1e-12,
                             "Entry of a matrix");
    }
    logger.assert_always(merge.getNumberOfRows() == 2, "Rows in a matrix");
    logger.assert_always(merge.getNumberOfColumns() == 4,
                         "Columns in a matrix");
    logger.assert_always(merge.size() == 8, "Size of a matrix");
    for (std::size_t i = 0; i < merge.size(); ++i) {
        logger.assert_always(std::abs(merge[i] - double(i)) < 1e-12,
                             "Entry of a matrix");
    }
    A23(0, 0) = 0.0;
    A23(1, 0) = 0.1;
    A23(0, 1) = 0.2;
    A23(1, 1) = 0.3;
    A23(0, 2) = 0.4;
    A23(1, 2) = 0.5;
    A32[0] = 0.0;
    A32[1] = 0.1;
    A32[2] = 0.2;
    A32[3] = 0.3;
    A32[4] = 0.4;
    A32[5] = 0.5;

    // out-of-place operators
    logger.assert_always(std::abs((count0 * count1)(0, 0) - 10.) < 1e-12,
                         "multiply");
    logger.assert_always(std::abs((count0 * count1)(1, 0) - 19.) < 1e-12,
                         "multiply");
    logger.assert_always(std::abs((count0 * count1)(0, 1) - 14.) < 1e-12,
                         "multiply");
    logger.assert_always(std::abs((count0 * count1)(1, 1) - 27.) < 1e-12,
                         "multiply");
    logger.assert_always(std::abs((count1 * count0)(0, 0) - 6.) < 1e-12,
                         "multiply");
    logger.assert_always(std::abs((count1 * count0)(1, 0) - 7.) < 1e-12,
                         "multiply");
    logger.assert_always(std::abs((count1 * count0)(0, 1) - 26.) < 1e-12,
                         "multiply");
    logger.assert_always(std::abs((count1 * count0)(1, 1) - 31.) < 1e-12,
                         "multiply");
    logger.assert_always(std::abs((count0 * A23)(0, 0) - .2) < 1e-12,
                         "multiply");
    logger.assert_always(std::abs((count0 * A23)(1, 0) - .3) < 1e-12,
                         "multiply");
    logger.assert_always(std::abs((count0 * A23)(0, 1) - .6) < 1e-12,
                         "multiply");
    logger.assert_always(std::abs((count0 * A23)(1, 1) - 1.1) < 1e-12,
                         "multiply");
    logger.assert_always(std::abs((count0 * A23)(0, 2) - 1.) < 1e-12,
                         "multiply");
    logger.assert_always(std::abs((count0 * A23)(1, 2) - 1.9) < 1e-12,
                         "multiply");
    logger.assert_always(std::abs((A32 * count0)(0, 0) - .3) < 1e-12,
                         "multiply");
    logger.assert_always(std::abs((A32 * count0)(1, 0) - .4) < 1e-12,
                         "multiply");
    logger.assert_always(std::abs((A32 * count0)(2, 0) - .5) < 1e-12,
                         "multiply");
    logger.assert_always(std::abs((A32 * count0)(0, 1) - .9) < 1e-12,
                         "multiply");
    logger.assert_always(std::abs((A32 * count0)(1, 1) - 1.4) < 1e-12,
                         "multiply");
    logger.assert_always(std::abs((A32 * count0)(2, 1) - 1.9) < 1e-12,
                         "multiply");
    logger.assert_always(std::abs((A32 * A23)(0, 0) - .03) < 1e-12, "multiply");
    logger.assert_always(std::abs((A32 * A23)(1, 0) - .04) < 1e-12, "multiply");
    logger.assert_always(std::abs((A32 * A23)(2, 0) - .05) < 1e-12, "multiply");
    logger.assert_always(std::abs((A32 * A23)(0, 1) - .09) < 1e-12, "multiply");
    logger.assert_always(std::abs((A32 * A23)(1, 1) - .14) < 1e-12, "multiply");
    logger.assert_always(std::abs((A32 * A23)(2, 1) - .19) < 1e-12, "multiply");
    logger.assert_always(std::abs((A32 * A23)(0, 2) - .15) < 1e-12, "multiply");
    logger.assert_always(std::abs((A32 * A23)(1, 2) - .24) < 1e-12, "multiply");
    logger.assert_always(std::abs((A32 * A23)(2, 2) - .33) < 1e-12, "multiply");
    logger.assert_always(std::abs((A23 * A32)(0, 0) - .1) < 1e-12, "multiply");
    logger.assert_always(std::abs((A23 * A32)(1, 0) - .13) < 1e-12, "multiply");
    logger.assert_always(std::abs((A23 * A32)(0, 1) - .28) < 1e-12, "multiply");
    logger.assert_always(std::abs((A23 * A32)(1, 1) - .40) < 1e-12, "multiply");
    logger.assert_always(std::abs((vec1 * count0) * vec1 - 45.) < 1e-12,
                         "multiply");
    logger.assert_always(std::abs(vec1 * (count0 * vec1) - 45.) < 1e-12,
                         "multiply");
    MiddleSizeVector size3 = {3., 4., 5.};
    logger.assert_always(std::abs(size3 * (A32 * vec0) - 5.) < 1e-12,
                         "multiply");
    logger.assert_always(std::abs((size3 * A32) * vec0 - 5.) < 1e-12,
                         "multiply");
    logger.assert_always(std::abs((2 * count0 * 2)(0, 0) - 0.) < 1e-12,
                         "multiply");
    logger.assert_always(std::abs((2 * count0 * 2)(1, 0) - 4.) < 1e-12,
                         "multiply");
    logger.assert_always(std::abs((2 * count0 * 2)(0, 1) - 8.) < 1e-12,
                         "multiply");
    logger.assert_always(std::abs((2 * count0 * 2)(1, 1) - 12.) < 1e-12,
                         "multiply");
    logger.assert_always(std::abs((count0 / 2.)(0, 0) - 0.) < 1e-12, "divide");
    logger.assert_always(std::abs((count0 / 2.)(1, 0) - .5) < 1e-12, "divide");
    logger.assert_always(std::abs((count0 / 2.)(0, 1) - 1.) < 1e-12, "divide");
    logger.assert_always(std::abs((count0 / 2.)(1, 1) - 1.5) < 1e-12, "divide");
    logger.assert_always(std::abs((count0 + count1)(0, 0) - 4.) < 1e-12, "add");
    logger.assert_always(std::abs((count0 + count1)(1, 0) - 6.) < 1e-12, "add");
    logger.assert_always(std::abs((count0 + count1)(0, 1) - 8.) < 1e-12, "add");
    logger.assert_always(std::abs((count0 + count1)(1, 1) - 10.) < 1e-12,
                         "add");
    logger.assert_always(std::abs((count1 + count0)(0, 0) - 4.) < 1e-12, "add");
    logger.assert_always(std::abs((count1 + count0)(1, 0) - 6.) < 1e-12, "add");
    logger.assert_always(std::abs((count1 + count0)(0, 1) - 8.) < 1e-12, "add");
    logger.assert_always(std::abs((count1 + count0)(1, 1) - 10.) < 1e-12,
                         "add");
    logger.assert_always(std::abs((count0 - count1)(0, 0) + 4.) < 1e-12,
                         "subtract");
    logger.assert_always(std::abs((count0 - count1)(1, 0) + 4.) < 1e-12,
                         "subtract");
    logger.assert_always(std::abs((count0 - count1)(0, 1) + 4.) < 1e-12,
                         "subtract");
    logger.assert_always(std::abs((count0 - count1)(1, 1) + 4.) < 1e-12,
                         "subtract");
    logger.assert_always(std::abs((count1 - count0)(0, 0) - 4.) < 1e-12,
                         "subtract");
    logger.assert_always(std::abs((count1 - count0)(1, 0) - 4.) < 1e-12,
                         "subtract");
    logger.assert_always(std::abs((count1 - count0)(0, 1) - 4.) < 1e-12,
                         "subtract");
    logger.assert_always(std::abs((count1 - count0)(1, 1) - 4.) < 1e-12,
                         "subtract");

    // assignent operators
    destroy = 4;
    logger.assert_always(destroy.getNumberOfRows() == 1, "Rows in a matrix");
    logger.assert_always(destroy.getNumberOfColumns() == 1,
                         "Columns in a matrix");
    logger.assert_always(destroy.size() == 1, "Size of a matrix");
    logger.assert_always(std::abs(destroy[0] - 4.) < 1e-12,
                         "Entry of a matrix");
    MiddleSizeMatrix extra = copy;
    logger.assert_always(extra.getNumberOfRows() == 2, "Rows in a matrix");
    logger.assert_always(extra.getNumberOfColumns() == 2,
                         "Columns in a matrix");
    logger.assert_always(extra.size() == 4, "Size of a matrix");
    for (std::size_t i = 0; i < extra.size(); ++i) {
        logger.assert_always(std::abs(extra[i] - double(i)) < 1e-12,
                             "Entry of a matrix");
    }
    copy = count1;
    logger.assert_always(copy.getNumberOfRows() == 2, "Rows in a matrix");
    logger.assert_always(copy.getNumberOfColumns() == 2, "Columns in a matrix");
    logger.assert_always(copy.size() == 4, "Size of a matrix");
    for (std::size_t i = 0; i < copy.size(); ++i) {
        logger.assert_always(std::abs(copy[i] - double(i) - 4.) < 1e-12,
                             "Entry of a matrix");
    }
    A0 = std::move(destroy);
    logger.assert_always(A0.getNumberOfRows() == 1, "Rows in a matrix");
    logger.assert_always(A0.getNumberOfColumns() == 1, "Columns in a matrix");
    logger.assert_always(A0.size() == 1, "Size of a matrix");
    logger.assert_always(std::abs(A0[0] - 4.) < 1e-12, "Entry of a matrix");
    count0 *= count1;
    logger.assert_always(count0.getNumberOfRows() == 2, "Rows in a matrix");
    logger.assert_always(count0.getNumberOfColumns() == 2,
                         "Columns in a matrix");
    logger.assert_always(count0.size() == 4, "Size of a matrix");
    logger.assert_always(std::abs((count0)(0, 0) - 10.) < 1e-12, "multiply");
    logger.assert_always(std::abs((count0)(1, 0) - 19.) < 1e-12, "multiply");
    logger.assert_always(std::abs((count0)(0, 1) - 14.) < 1e-12, "multiply");
    logger.assert_always(std::abs((count0)(1, 1) - 27.) < 1e-12, "multiply");
    count1 *= A23;
    logger.assert_always(count1.getNumberOfRows() == 2, "Rows in a matrix");
    logger.assert_always(count1.getNumberOfColumns() == 3,
                         "Columns in a matrix");
    logger.assert_always(count1.size() == 6, "Size of a matrix");
    logger.assert_always(std::abs((count1)(0, 0) - .6) < 1e-12, "multiply");
    logger.assert_always(std::abs((count1)(1, 0) - .7) < 1e-12, "multiply");
    logger.assert_always(std::abs((count1)(0, 1) - 2.6) < 1e-12, "multiply");
    logger.assert_always(std::abs((count1)(1, 1) - 3.1) < 1e-12, "multiply");
    logger.assert_always(std::abs((count1)(0, 2) - 4.6) < 1e-12, "multiply");
    logger.assert_always(std::abs((count1)(1, 2) - 5.5) < 1e-12, "multiply");
    count0 *= 4;
    logger.assert_always(count0.getNumberOfRows() == 2, "Rows in a matrix");
    logger.assert_always(count0.getNumberOfColumns() == 2,
                         "Columns in a matrix");
    logger.assert_always(count0.size() == 4, "Size of a matrix");
    logger.assert_always(std::abs((count0)(0, 0) - 40.) < 1e-12, "multiply");
    logger.assert_always(std::abs((count0)(1, 0) - 76.) < 1e-12, "multiply");
    logger.assert_always(std::abs((count0)(0, 1) - 56.) < 1e-12, "multiply");
    logger.assert_always(std::abs((count0)(1, 1) - 108.) < 1e-12, "multiply");
    count0 /= 2;
    logger.assert_always(count0.getNumberOfRows() == 2, "Rows in a matrix");
    logger.assert_always(count0.getNumberOfColumns() == 2,
                         "Columns in a matrix");
    logger.assert_always(count0.size() == 4, "Size of a matrix");
    logger.assert_always(std::abs((count0)(0, 0) - 20.) < 1e-12, "divide");
    logger.assert_always(std::abs((count0)(1, 0) - 38.) < 1e-12, "divide");
    logger.assert_always(std::abs((count0)(0, 1) - 28.) < 1e-12, "divide");
    logger.assert_always(std::abs((count0)(1, 1) - 54.) < 1e-12, "divide");
    copy += extra;
    logger.assert_always(copy.getNumberOfRows() == 2, "Rows in a matrix");
    logger.assert_always(copy.getNumberOfColumns() == 2, "Columns in a matrix");
    logger.assert_always(copy.size() == 4, "Size of a matrix");
    logger.assert_always(std::abs((copy)(0, 0) - 4.) < 1e-12, "add");
    logger.assert_always(std::abs((copy)(1, 0) - 6.) < 1e-12, "add");
    logger.assert_always(std::abs((copy)(0, 1) - 8.) < 1e-12, "add");
    logger.assert_always(std::abs((copy)(1, 1) - 10.) < 1e-12, "add");
    extra -= count0;
    logger.assert_always(extra.getNumberOfRows() == 2, "Rows in a matrix");
    logger.assert_always(extra.getNumberOfColumns() == 2,
                         "Columns in a matrix");
    logger.assert_always(extra.size() == 4, "Size of a matrix");
    logger.assert_always(std::abs((extra)(0, 0) + 20.) < 1e-12, "subtract");
    logger.assert_always(std::abs((extra)(1, 0) + 37.) < 1e-12, "subtract");
    logger.assert_always(std::abs((extra)(0, 1) + 26.) < 1e-12, "subtract");
    logger.assert_always(std::abs((extra)(1, 1) + 51.) < 1e-12, "subtract");
    extra.axpy(3., copy);
    logger.assert_always(extra.getNumberOfRows() == 2, "Rows in a matrix");
    logger.assert_always(extra.getNumberOfColumns() == 2,
                         "Columns in a matrix");
    logger.assert_always(extra.size() == 4, "Size of a matrix");
    logger.assert_always(std::abs((extra)(0, 0) + 8.) < 1e-12, "ax+y");
    logger.assert_always(std::abs((extra)(1, 0) + 19.) < 1e-12, "ax+y");
    logger.assert_always(std::abs((extra)(0, 1) + 2.) < 1e-12, "ax+y");
    logger.assert_always(std::abs((extra)(1, 1) + 21.) < 1e-12, "ax+y");
    A0.resize(3, 7);
    logger.assert_always(A0.getNumberOfRows() == 3, "Rows in a matrix");
    logger.assert_always(A0.getNumberOfColumns() == 7, "Columns in a matrix");
    logger.assert_always(A0.size() == 21, "Size of a matrix");
    logger.assert_always(std::abs(A0[0] - 4.) < 1e-12, "Entry of a matrix");

    // wedge stuff
    logger.assert_always(
        std::abs((MiddleSizeMatrix(vec0).computeWedgeStuffVector()) *
                     (MiddleSizeMatrix(vec0).computeWedgeStuffVector()) -
                 vec0 * vec0) < 1e-12,
        "norm of wedge stuff vector");
    logger.assert_always(
        std::abs((MiddleSizeMatrix(vec0).computeWedgeStuffVector()) * vec0) <
            1e-12,
        "direction of wedge stuff vector");
    logger.assert_always(
        std::abs((MiddleSizeMatrix(vec1).computeWedgeStuffVector()) *
                     (MiddleSizeMatrix(vec1).computeWedgeStuffVector()) -
                 vec1 * vec1) < 1e-12,
        "norm of wedge stuff vector");
    logger.assert_always(
        std::abs((MiddleSizeMatrix(vec1).computeWedgeStuffVector()) * vec1) <
            1e-12,
        "direction of wedge stuff vector");
    logger.assert_always(
        std::abs((MiddleSizeMatrix(vec2).computeWedgeStuffVector()) *
                     (MiddleSizeMatrix(vec2).computeWedgeStuffVector()) -
                 vec2 * vec2) < 1e-12,
        "norm of wedge stuff vector");
    logger.assert_always(
        std::abs((MiddleSizeMatrix(vec2).computeWedgeStuffVector()) * vec2) <
            1e-12,
        "direction of wedge stuff vector");
    logger.assert_always(
        std::abs((MiddleSizeMatrix(vec3).computeWedgeStuffVector()) *
                     (MiddleSizeMatrix(vec3).computeWedgeStuffVector()) -
                 vec3 * vec3) < 1e-12,
        "norm of wedge stuff vector");
    logger.assert_always(
        std::abs((MiddleSizeMatrix(vec3).computeWedgeStuffVector()) * vec3) <
            1e-12,
        "direction of wedge stuff vector");
    MiddleSizeVector vec3D0{0., 1., 2.}, vec3D1{3., 4., 5.},
        vec3D2{0., -1., 2.};
    ///\todo test that the norm of the 3D wedge stuff vector equals the area of
    ///the triangle formed by nodes {0, 0, 0}, v1 and v2
    logger.assert_always(
        std::abs(
            (MiddleSizeMatrix({vec3D0, vec3D1}).computeWedgeStuffVector()) *
            vec3D0) < 1e-12,
        "direction of wedge stuff vector");
    logger.assert_always(
        std::abs(
            (MiddleSizeMatrix({vec3D0, vec3D1}).computeWedgeStuffVector()) *
            vec3D1) < 1e-12,
        "direction of wedge stuff vector");
    logger.assert_always(
        std::abs(
            (MiddleSizeMatrix({vec3D1, vec3D2}).computeWedgeStuffVector()) *
            vec3D1) < 1e-12,
        "direction of wedge stuff vector");
    logger.assert_always(
        std::abs(
            (MiddleSizeMatrix({vec3D1, vec3D2}).computeWedgeStuffVector()) *
            vec3D2) < 1e-12,
        "direction of wedge stuff vector");
    logger.assert_always(
        std::abs(
            (MiddleSizeMatrix({vec3D0, vec3D2}).computeWedgeStuffVector()) *
            vec3D0) < 1e-12,
        "direction of wedge stuff vector");
    logger.assert_always(
        std::abs(
            (MiddleSizeMatrix({vec3D0, vec3D2}).computeWedgeStuffVector()) *
            vec3D2) < 1e-12,
        "direction of wedge stuff vector");
    MiddleSizeVector vec4D0{0., 1., 2., 3.}, vec4D1{4., 5., 6., 7.},
        vec4D2{0., -1., 2., -3.}, vec4D3{0., -1., -2., 3.};
    ///\todo test that the norm of the 4D wedge stuff vector equals the area of
    ///the tetrahedron formed by nodes {0, 0, 0}, v1 and v2 and v3
    logger.assert_always(std::abs((MiddleSizeMatrix({vec4D0, vec4D1, vec4D2})
                                       .computeWedgeStuffVector()) *
                                  vec4D0) < 1e-12,
                         "direction of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix({vec4D0, vec4D1, vec4D2})
                                       .computeWedgeStuffVector()) *
                                  vec4D1) < 1e-12,
                         "direction of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix({vec4D0, vec4D1, vec4D2})
                                       .computeWedgeStuffVector()) *
                                  vec4D2) < 1e-12,
                         "direction of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix({vec4D0, vec4D1, vec4D3})
                                       .computeWedgeStuffVector()) *
                                  vec4D0) < 1e-12,
                         "direction of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix({vec4D0, vec4D1, vec4D3})
                                       .computeWedgeStuffVector()) *
                                  vec4D1) < 1e-12,
                         "direction of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix({vec4D0, vec4D1, vec4D3})
                                       .computeWedgeStuffVector()) *
                                  vec4D3) < 1e-12,
                         "direction of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix({vec4D0, vec4D2, vec4D3})
                                       .computeWedgeStuffVector()) *
                                  vec4D0) < 1e-12,
                         "direction of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix({vec4D0, vec4D2, vec4D3})
                                       .computeWedgeStuffVector()) *
                                  vec4D2) < 1e-12,
                         "direction of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix({vec4D0, vec4D2, vec4D3})
                                       .computeWedgeStuffVector()) *
                                  vec4D3) < 1e-12,
                         "direction of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix({vec4D1, vec4D2, vec4D3})
                                       .computeWedgeStuffVector()) *
                                  vec4D1) < 1e-12,
                         "direction of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix({vec4D1, vec4D2, vec4D3})
                                       .computeWedgeStuffVector()) *
                                  vec4D2) < 1e-12,
                         "direction of wedge stuff vector");
    logger.assert_always(std::abs((MiddleSizeMatrix({vec4D1, vec4D2, vec4D3})
                                       .computeWedgeStuffVector()) *
                                  vec4D3) < 1e-12,
                         "direction of wedge stuff vector");

    copy.concatenate(extra);
    logger.assert_always(copy.getNumberOfRows() == 4, "Rows in a matrix");
    logger.assert_always(copy.getNumberOfColumns() == 2, "Columns in a matrix");
    logger.assert_always(copy.size() == 8, "Size of a matrix");
    logger.assert_always(std::abs((copy)(0, 0) - 4.) < 1e-12, "concatenate");
    logger.assert_always(std::abs((copy)(1, 0) - 6.) < 1e-12, "concatenate");
    logger.assert_always(std::abs((copy)(0, 1) - 8.) < 1e-12, "concatenate");
    logger.assert_always(std::abs((copy)(1, 1) - 10.) < 1e-12, "concatenate");
    logger.assert_always(std::abs((copy)(2, 0) + 8.) < 1e-12, "concatenate");
    logger.assert_always(std::abs((copy)(3, 0) + 19.) < 1e-12, "concatenate");
    logger.assert_always(std::abs((copy)(2, 1) + 2.) < 1e-12, "concatenate");
    logger.assert_always(std::abs((copy)(3, 1) + 21.) < 1e-12, "concatenate");

    logger.assert_always(
        MiddleSizeMatrix({vec3D0, vec3D1}).getColumn(0).size() == 3,
        "getColumn");
    logger.assert_always(
        std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getColumn(0)) -
                  vec3D0)[0]) < 1e-12,
        "getColumn");
    logger.assert_always(
        std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getColumn(0)) -
                  vec3D0)[1]) < 1e-12,
        "getColumn");
    logger.assert_always(
        std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getColumn(0)) -
                  vec3D0)[2]) < 1e-12,
        "getColumn");
    logger.assert_always(
        MiddleSizeMatrix({vec3D0, vec3D1}).getColumn(1).size() == 3,
        "getColumn");
    logger.assert_always(
        std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getColumn(1)) -
                  vec3D1)[0]) < 1e-12,
        "getColumn");
    logger.assert_always(
        std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getColumn(1)) -
                  vec3D1)[1]) < 1e-12,
        "getColumn");
    logger.assert_always(
        std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getColumn(1)) -
                  vec3D1)[2]) < 1e-12,
        "getColumn");
    logger.assert_always(
        MiddleSizeMatrix({vec3D0, vec3D1}).getRow(0).size() == 2, "getRow");
    logger.assert_always(
        std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getRow(0)))[0] - 0.) <
            1e-12,
        "getRow");
    logger.assert_always(
        std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getRow(0)))[1] - 3.) <
            1e-12,
        "getRow");
    logger.assert_always(
        MiddleSizeMatrix({vec3D0, vec3D1}).getRow(1).size() == 2, "getRow");
    logger.assert_always(
        std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getRow(1)))[0] - 1.) <
            1e-12,
        "getRow");
    logger.assert_always(
        std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getRow(1)))[1] - 4.) <
            1e-12,
        "getRow");
    logger.assert_always(
        MiddleSizeMatrix({vec3D0, vec3D1}).getRow(2).size() == 2, "getRow");
    logger.assert_always(
        std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getRow(2)))[0] - 2.) <
            1e-12,
        "getRow");
    logger.assert_always(
        std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getRow(2)))[1] - 5.) <
            1e-12,
        "getRow");

    ///\todo figure out a way to test a LU factorisation

    MiddleSizeVector duplicate = vec2;
    count0.solve(vec2);
    logger.assert_always(
        std::abs((vec2 - count0.inverse() * duplicate)[0]) < 1e-12,
        "inverse and solve");
    logger.assert_always(
        std::abs((vec2 - count0.inverse() * duplicate)[1]) < 1e-12,
        "inverse and solve");
    logger.assert_always(
        std::abs((count0 - count0.inverse().inverse())[0]) < 1e-12, "inverse");
    logger.assert_always(
        std::abs((count0 - count0.inverse().inverse())[1]) < 1e-12, "inverse");
    logger.assert_always(
        std::abs((count0 - count0.inverse().inverse())[2]) < 1e-12, "inverse");
    logger.assert_always(
        std::abs((count0 - count0.inverse().inverse())[3]) < 1e-12, "inverse");
    logger.assert_always(
        std::abs((count0 - count0.transpose().transpose())[0]) < 1e-12,
        "transpose");
    logger.assert_always(
        std::abs((count0 - count0.transpose().transpose())[1]) < 1e-12,
        "transpose");
    logger.assert_always(
        std::abs((count0 - count0.transpose().transpose())[2]) < 1e-12,
        "transpose");
    logger.assert_always(
        std::abs((count0 - count0.transpose().transpose())[3]) < 1e-12,
        "transpose");
    logger.assert_always(std::abs(A23(0, 0) - A23.transpose()(0, 0)) < 1e-12,
                         "transpose");
    logger.assert_always(std::abs(A23(1, 0) - A23.transpose()(0, 1)) < 1e-12,
                         "transpose");
    logger.assert_always(std::abs(A23(0, 1) - A23.transpose()(1, 0)) < 1e-12,
                         "transpose");
    logger.assert_always(std::abs(A23(1, 1) - A23.transpose()(1, 1)) < 1e-12,
                         "transpose");
    logger.assert_always(std::abs(A23(0, 2) - A23.transpose()(2, 0)) < 1e-12,
                         "transpose");
    logger.assert_always(std::abs(A23(1, 2) - A23.transpose()(2, 1)) < 1e-12,
                         "transpose");

    auto data = A23.data();
    logger.assert_always(std::abs(data[0] - A23[0]) < 1e-12, "data");
    logger.assert_always(std::abs(data[1] - A23[1]) < 1e-12, "data");
    logger.assert_always(std::abs(data[2] - A23[2]) < 1e-12, "data");
    logger.assert_always(std::abs(data[3] - A23[3]) < 1e-12, "data");
    logger.assert_always(std::abs(data[4] - A23[4]) < 1e-12, "data");
    logger.assert_always(std::abs(data[5] - A23[5]) < 1e-12, "data");
    data[2] = 17.3;
    logger.assert_always(std::abs(A23[2] - 17.3) < 1e-12, "data");

    std::cout << A32 << std::endl;

    return 0;
}

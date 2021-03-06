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

#include "../catch.hpp"

using namespace hpgem;
using LinearAlgebra::MiddleSizeMatrix;
using LinearAlgebra::MiddleSizeVector;

TEST_CASE("MiddleSizeMatrixUnitTest", "[MiddleSizeMatrixUnitTest]") {
    // constructors
    LinearAlgebra::MiddleSizeVector vec0{0., 1.}, vec1{2., 3.}, vec2{4., 5.},
        vec3{6., 7.};
    MiddleSizeMatrix A0, A22(2, 2), A23(2, 3), A32(3, 2), destroy(3, 3, 1),
        count0({vec0, vec1}), count1({vec2, vec3}), copy(count0), bla({count0}),
        merge({count0, count1});
    INFO("Rows in a matrix");
    CHECK(destroy.getNumberOfRows() == 3);
    INFO("Columns in a matrix");
    CHECK(destroy.getNumberOfColumns() == 3);
    INFO("Size of a matrix");
    CHECK(destroy.size() == 9);
    for (std::size_t i = 0; i < destroy.size(); ++i) {
        INFO("Entry of a matrix");
        CHECK(std::abs(destroy[i] - 1.) < 1e-12);
    }
    for (std::size_t i = 0; i < destroy.getNumberOfRows(); ++i) {
        for (std::size_t j = 0; j < destroy.getNumberOfColumns(); ++j) {
            INFO("Entry of a matrix");
            CHECK(std::abs(destroy(i, j) - 1.) < 1e-12);
        }
    }
    MiddleSizeMatrix moved(std::move(destroy));
    INFO("Rows in a matrix");
    CHECK(moved.getNumberOfRows() == 3);
    INFO("Columns in a matrix");
    CHECK(moved.getNumberOfColumns() == 3);
    INFO("Size of a matrix");
    CHECK(moved.size() == 9);
    for (std::size_t i = 0; i < moved.size(); ++i) {
        INFO("Entry of a matrix");
        CHECK(std::abs(moved[i] - 1.) < 1e-12);
    }
    for (std::size_t i = 0; i < moved.getNumberOfRows(); ++i) {
        for (std::size_t j = 0; j < moved.getNumberOfColumns(); ++j) {
            INFO("Entry of a matrix");
            CHECK(std::abs(moved(i, j) - 1.) < 1e-12);
        }
    }
    INFO("Rows in a matrix");
    CHECK(A0.getNumberOfRows() == 0);
    INFO("Columns in a matrix");
    CHECK(A0.getNumberOfColumns() == 0);
    INFO("Size of a matrix");
    CHECK(A0.size() == 0);
    INFO("Rows in a matrix");
    CHECK(A22.getNumberOfRows() == 2);
    INFO("Columns in a matrix");
    CHECK(A22.getNumberOfColumns() == 2);
    INFO("Size of a matrix");
    CHECK(A22.size() == 4);
    INFO("Rows in a matrix");
    CHECK(A23.getNumberOfRows() == 2);
    INFO("Columns in a matrix");
    CHECK(A23.getNumberOfColumns() == 3);
    INFO("Size of a matrix");
    CHECK(A23.size() == 6);
    INFO("Rows in a matrix");
    CHECK(A32.getNumberOfRows() == 3);
    INFO("Columns in a matrix");
    CHECK(A32.getNumberOfColumns() == 2);
    INFO("Size of a matrix");
    CHECK(A32.size() == 6);
    INFO("Rows in a matrix");
    CHECK(count0.getNumberOfRows() == 2);
    INFO("Columns in a matrix");
    CHECK(count0.getNumberOfColumns() == 2);
    INFO("Size of a matrix");
    CHECK(count0.size() == 4);
    for (std::size_t i = 0; i < count0.size(); ++i) {
        INFO("Entry of a matrix");
        CHECK(std::abs(count0[i] - double(i)) < 1e-12);
    }
    INFO("Entry of a matrix");
    CHECK(std::abs(count0(0, 0) - 0.) < 1e-12);
    INFO("Entry of a matrix");
    CHECK(std::abs(count0(1, 0) - 1.) < 1e-12);
    INFO("Entry of a matrix");
    CHECK(std::abs(count0(0, 1) - 2.) < 1e-12);
    INFO("Entry of a matrix");
    CHECK(std::abs(count0(1, 1) - 3.) < 1e-12);
    INFO("Rows in a matrix");
    CHECK(count1.getNumberOfRows() == 2);
    INFO("Columns in a matrix");
    CHECK(count1.getNumberOfColumns() == 2);
    INFO("Size of a matrix");
    CHECK(count1.size() == 4);
    for (std::size_t i = 0; i < count1.size(); ++i) {
        INFO("Entry of a matrix");
        CHECK(std::abs(count1[i] - 4. - double(i)) < 1e-12);
    }
    INFO("Rows in a matrix");
    CHECK(copy.getNumberOfRows() == 2);
    INFO("Columns in a matrix");
    CHECK(copy.getNumberOfColumns() == 2);
    INFO("Size of a matrix");
    CHECK(copy.size() == 4);
    for (std::size_t i = 0; i < copy.size(); ++i) {
        INFO("Entry of a matrix");
        CHECK(std::abs(copy[i] - double(i)) < 1e-12);
    }
    INFO("Rows in a matrix");
    CHECK(bla.getNumberOfRows() == 2);
    INFO("Columns in a matrix");
    CHECK(bla.getNumberOfColumns() == 2);
    INFO("Size of a matrix");
    CHECK(bla.size() == 4);
    for (std::size_t i = 0; i < bla.size(); ++i) {
        INFO("Entry of a matrix");
        CHECK(std::abs(bla[i] - double(i)) < 1e-12);
    }
    INFO("Rows in a matrix");
    CHECK(merge.getNumberOfRows() == 2);
    INFO("Columns in a matrix");
    CHECK(merge.getNumberOfColumns() == 4);
    INFO("Size of a matrix");
    CHECK(merge.size() == 8);
    for (std::size_t i = 0; i < merge.size(); ++i) {
        INFO("Entry of a matrix");
        CHECK(std::abs(merge[i] - double(i)) < 1e-12);
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
    INFO("multiply");
    CHECK(std::abs((count0 * count1)(0, 0) - 10.) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((count0 * count1)(1, 0) - 19.) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((count0 * count1)(0, 1) - 14.) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((count0 * count1)(1, 1) - 27.) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((count1 * count0)(0, 0) - 6.) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((count1 * count0)(1, 0) - 7.) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((count1 * count0)(0, 1) - 26.) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((count1 * count0)(1, 1) - 31.) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((count0 * A23)(0, 0) - .2) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((count0 * A23)(1, 0) - .3) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((count0 * A23)(0, 1) - .6) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((count0 * A23)(1, 1) - 1.1) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((count0 * A23)(0, 2) - 1.) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((count0 * A23)(1, 2) - 1.9) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((A32 * count0)(0, 0) - .3) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((A32 * count0)(1, 0) - .4) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((A32 * count0)(2, 0) - .5) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((A32 * count0)(0, 1) - .9) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((A32 * count0)(1, 1) - 1.4) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((A32 * count0)(2, 1) - 1.9) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((A32 * A23)(0, 0) - .03) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((A32 * A23)(1, 0) - .04) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((A32 * A23)(2, 0) - .05) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((A32 * A23)(0, 1) - .09) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((A32 * A23)(1, 1) - .14) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((A32 * A23)(2, 1) - .19) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((A32 * A23)(0, 2) - .15) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((A32 * A23)(1, 2) - .24) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((A32 * A23)(2, 2) - .33) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((A23 * A32)(0, 0) - .1) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((A23 * A32)(1, 0) - .13) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((A23 * A32)(0, 1) - .28) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((A23 * A32)(1, 1) - .40) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((vec1 * count0) * vec1 - 45.) < 1e-12);
    INFO("multiply");
    CHECK(std::abs(vec1 * (count0 * vec1) - 45.) < 1e-12);
    MiddleSizeVector size3 = {3., 4., 5.};
    INFO("multiply");
    CHECK(std::abs(size3 * (A32 * vec0) - 5.) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((size3 * A32) * vec0 - 5.) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((2 * count0 * 2)(0, 0) - 0.) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((2 * count0 * 2)(1, 0) - 4.) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((2 * count0 * 2)(0, 1) - 8.) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((2 * count0 * 2)(1, 1) - 12.) < 1e-12);
    INFO("divide");
    CHECK(std::abs((count0 / 2.)(0, 0) - 0.) < 1e-12);
    INFO("divide");
    CHECK(std::abs((count0 / 2.)(1, 0) - .5) < 1e-12);
    INFO("divide");
    CHECK(std::abs((count0 / 2.)(0, 1) - 1.) < 1e-12);
    INFO("divide");
    CHECK(std::abs((count0 / 2.)(1, 1) - 1.5) < 1e-12);
    INFO("add");
    CHECK(std::abs((count0 + count1)(0, 0) - 4.) < 1e-12);
    INFO("add");
    CHECK(std::abs((count0 + count1)(1, 0) - 6.) < 1e-12);
    INFO("add");
    CHECK(std::abs((count0 + count1)(0, 1) - 8.) < 1e-12);
    INFO("add");
    CHECK(std::abs((count0 + count1)(1, 1) - 10.) < 1e-12);
    INFO("add");
    CHECK(std::abs((count1 + count0)(0, 0) - 4.) < 1e-12);
    INFO("add");
    CHECK(std::abs((count1 + count0)(1, 0) - 6.) < 1e-12);
    INFO("add");
    CHECK(std::abs((count1 + count0)(0, 1) - 8.) < 1e-12);
    INFO("add");
    CHECK(std::abs((count1 + count0)(1, 1) - 10.) < 1e-12);
    INFO("subtract");
    CHECK(std::abs((count0 - count1)(0, 0) + 4.) < 1e-12);
    INFO("subtract");
    CHECK(std::abs((count0 - count1)(1, 0) + 4.) < 1e-12);
    INFO("subtract");
    CHECK(std::abs((count0 - count1)(0, 1) + 4.) < 1e-12);
    INFO("subtract");
    CHECK(std::abs((count0 - count1)(1, 1) + 4.) < 1e-12);
    INFO("subtract");
    CHECK(std::abs((count1 - count0)(0, 0) - 4.) < 1e-12);
    INFO("subtract");
    CHECK(std::abs((count1 - count0)(1, 0) - 4.) < 1e-12);
    INFO("subtract");
    CHECK(std::abs((count1 - count0)(0, 1) - 4.) < 1e-12);
    INFO("subtract");
    CHECK(std::abs((count1 - count0)(1, 1) - 4.) < 1e-12);

    // assignent operators
    destroy = 4;
    INFO("Rows in a matrix");
    CHECK(destroy.getNumberOfRows() == 1);
    INFO("Columns in a matrix");
    CHECK(destroy.getNumberOfColumns() == 1);
    INFO("Size of a matrix");
    CHECK(destroy.size() == 1);
    INFO("Entry of a matrix");
    CHECK(std::abs(destroy[0] - 4.) < 1e-12);
    MiddleSizeMatrix extra = copy;
    INFO("Rows in a matrix");
    CHECK(extra.getNumberOfRows() == 2);
    INFO("Columns in a matrix");
    CHECK(extra.getNumberOfColumns() == 2);
    INFO("Size of a matrix");
    CHECK(extra.size() == 4);
    for (std::size_t i = 0; i < extra.size(); ++i) {
        INFO("Entry of a matrix");
        CHECK(std::abs(extra[i] - double(i)) < 1e-12);
    }
    copy = count1;
    INFO("Rows in a matrix");
    CHECK(copy.getNumberOfRows() == 2);
    INFO("Columns in a matrix");
    CHECK(copy.getNumberOfColumns() == 2);
    INFO("Size of a matrix");
    CHECK(copy.size() == 4);
    for (std::size_t i = 0; i < copy.size(); ++i) {
        INFO("Entry of a matrix");
        CHECK(std::abs(copy[i] - double(i) - 4.) < 1e-12);
    }
    A0 = std::move(destroy);
    INFO("Rows in a matrix");
    CHECK(A0.getNumberOfRows() == 1);
    INFO("Columns in a matrix");
    CHECK(A0.getNumberOfColumns() == 1);
    INFO("Size of a matrix");
    CHECK(A0.size() == 1);
    INFO("Entry of a matrix");
    CHECK(std::abs(A0[0] - 4.) < 1e-12);
    count0 *= count1;
    INFO("Rows in a matrix");
    CHECK(count0.getNumberOfRows() == 2);
    INFO("Columns in a matrix");
    CHECK(count0.getNumberOfColumns() == 2);
    INFO("Size of a matrix");
    CHECK(count0.size() == 4);
    INFO("multiply");
    CHECK(std::abs((count0)(0, 0) - 10.) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((count0)(1, 0) - 19.) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((count0)(0, 1) - 14.) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((count0)(1, 1) - 27.) < 1e-12);
    count1 *= A23;
    INFO("Rows in a matrix");
    CHECK(count1.getNumberOfRows() == 2);
    INFO("Columns in a matrix");
    CHECK(count1.getNumberOfColumns() == 3);
    INFO("Size of a matrix");
    CHECK(count1.size() == 6);
    INFO("multiply");
    CHECK(std::abs((count1)(0, 0) - .6) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((count1)(1, 0) - .7) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((count1)(0, 1) - 2.6) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((count1)(1, 1) - 3.1) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((count1)(0, 2) - 4.6) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((count1)(1, 2) - 5.5) < 1e-12);
    count0 *= 4;
    INFO("Rows in a matrix");
    CHECK(count0.getNumberOfRows() == 2);
    INFO("Columns in a matrix");
    CHECK(count0.getNumberOfColumns() == 2);
    INFO("Size of a matrix");
    CHECK(count0.size() == 4);
    INFO("multiply");
    CHECK(std::abs((count0)(0, 0) - 40.) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((count0)(1, 0) - 76.) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((count0)(0, 1) - 56.) < 1e-12);
    INFO("multiply");
    CHECK(std::abs((count0)(1, 1) - 108.) < 1e-12);
    count0 /= 2;
    INFO("Rows in a matrix");
    CHECK(count0.getNumberOfRows() == 2);
    INFO("Columns in a matrix");
    CHECK(count0.getNumberOfColumns() == 2);
    INFO("Size of a matrix");
    CHECK(count0.size() == 4);
    INFO("divide");
    CHECK(std::abs((count0)(0, 0) - 20.) < 1e-12);
    INFO("divide");
    CHECK(std::abs((count0)(1, 0) - 38.) < 1e-12);
    INFO("divide");
    CHECK(std::abs((count0)(0, 1) - 28.) < 1e-12);
    INFO("divide");
    CHECK(std::abs((count0)(1, 1) - 54.) < 1e-12);
    copy += extra;
    INFO("Rows in a matrix");
    CHECK(copy.getNumberOfRows() == 2);
    INFO("Columns in a matrix");
    CHECK(copy.getNumberOfColumns() == 2);
    INFO("Size of a matrix");
    CHECK(copy.size() == 4);
    INFO("add");
    CHECK(std::abs((copy)(0, 0) - 4.) < 1e-12);
    INFO("add");
    CHECK(std::abs((copy)(1, 0) - 6.) < 1e-12);
    INFO("add");
    CHECK(std::abs((copy)(0, 1) - 8.) < 1e-12);
    INFO("add");
    CHECK(std::abs((copy)(1, 1) - 10.) < 1e-12);
    extra -= count0;
    INFO("Rows in a matrix");
    CHECK(extra.getNumberOfRows() == 2);
    INFO("Columns in a matrix");
    CHECK(extra.getNumberOfColumns() == 2);
    INFO("Size of a matrix");
    CHECK(extra.size() == 4);
    INFO("subtract");
    CHECK(std::abs((extra)(0, 0) + 20.) < 1e-12);
    INFO("subtract");
    CHECK(std::abs((extra)(1, 0) + 37.) < 1e-12);
    INFO("subtract");
    CHECK(std::abs((extra)(0, 1) + 26.) < 1e-12);
    INFO("subtract");
    CHECK(std::abs((extra)(1, 1) + 51.) < 1e-12);
    extra.axpy(3., copy);
    INFO("Rows in a matrix");
    CHECK(extra.getNumberOfRows() == 2);
    INFO("Columns in a matrix");
    CHECK(extra.getNumberOfColumns() == 2);
    INFO("Size of a matrix");
    CHECK(extra.size() == 4);
    INFO("ax+y");
    CHECK(std::abs((extra)(0, 0) + 8.) < 1e-12);
    INFO("ax+y");
    CHECK(std::abs((extra)(1, 0) + 19.) < 1e-12);
    INFO("ax+y");
    CHECK(std::abs((extra)(0, 1) + 2.) < 1e-12);
    INFO("ax+y");
    CHECK(std::abs((extra)(1, 1) + 21.) < 1e-12);
    A0.resize(3, 7);
    INFO("Rows in a matrix");
    CHECK(A0.getNumberOfRows() == 3);
    INFO("Columns in a matrix");
    CHECK(A0.getNumberOfColumns() == 7);
    INFO("Size of a matrix");
    CHECK(A0.size() == 21);
    INFO("Entry of a matrix");
    CHECK(std::abs(A0[0] - 4.) < 1e-12);

    // wedge stuff
    INFO("norm of wedge stuff vector");
    CHECK(std::abs((MiddleSizeMatrix(vec0).computeWedgeStuffVector()) *
                       (MiddleSizeMatrix(vec0).computeWedgeStuffVector()) -
                   vec0 * vec0) < 1e-12);
    INFO("direction of wedge stuff vector");
    CHECK(std::abs((MiddleSizeMatrix(vec0).computeWedgeStuffVector()) * vec0) <
          1e-12);
    INFO("norm of wedge stuff vector");
    CHECK(std::abs((MiddleSizeMatrix(vec1).computeWedgeStuffVector()) *
                       (MiddleSizeMatrix(vec1).computeWedgeStuffVector()) -
                   vec1 * vec1) < 1e-12);
    INFO("direction of wedge stuff vector");
    CHECK(std::abs((MiddleSizeMatrix(vec1).computeWedgeStuffVector()) * vec1) <
          1e-12);
    INFO("norm of wedge stuff vector");
    CHECK(std::abs((MiddleSizeMatrix(vec2).computeWedgeStuffVector()) *
                       (MiddleSizeMatrix(vec2).computeWedgeStuffVector()) -
                   vec2 * vec2) < 1e-12);
    INFO("direction of wedge stuff vector");
    CHECK(std::abs((MiddleSizeMatrix(vec2).computeWedgeStuffVector()) * vec2) <
          1e-12);
    INFO("norm of wedge stuff vector");
    CHECK(std::abs((MiddleSizeMatrix(vec3).computeWedgeStuffVector()) *
                       (MiddleSizeMatrix(vec3).computeWedgeStuffVector()) -
                   vec3 * vec3) < 1e-12);
    INFO("direction of wedge stuff vector");
    CHECK(std::abs((MiddleSizeMatrix(vec3).computeWedgeStuffVector()) * vec3) <
          1e-12);
    MiddleSizeVector vec3D0{0., 1., 2.}, vec3D1{3., 4., 5.},
        vec3D2{0., -1., 2.};
    ///\todo test that the norm of the 3D wedge stuff vector equals the area of
    /// the triangle formed by nodes {0, 0, 0}, v1 and v2
    INFO("direction of wedge stuff vector");
    CHECK(std::abs(
              (MiddleSizeMatrix({vec3D0, vec3D1}).computeWedgeStuffVector()) *
              vec3D0) < 1e-12);
    INFO("direction of wedge stuff vector");
    CHECK(std::abs(
              (MiddleSizeMatrix({vec3D0, vec3D1}).computeWedgeStuffVector()) *
              vec3D1) < 1e-12);
    INFO("direction of wedge stuff vector");
    CHECK(std::abs(
              (MiddleSizeMatrix({vec3D1, vec3D2}).computeWedgeStuffVector()) *
              vec3D1) < 1e-12);
    INFO("direction of wedge stuff vector");
    CHECK(std::abs(
              (MiddleSizeMatrix({vec3D1, vec3D2}).computeWedgeStuffVector()) *
              vec3D2) < 1e-12);
    INFO("direction of wedge stuff vector");
    CHECK(std::abs(
              (MiddleSizeMatrix({vec3D0, vec3D2}).computeWedgeStuffVector()) *
              vec3D0) < 1e-12);
    INFO("direction of wedge stuff vector");
    CHECK(std::abs(
              (MiddleSizeMatrix({vec3D0, vec3D2}).computeWedgeStuffVector()) *
              vec3D2) < 1e-12);
    MiddleSizeVector vec4D0{0., 1., 2., 3.}, vec4D1{4., 5., 6., 7.},
        vec4D2{0., -1., 2., -3.}, vec4D3{0., -1., -2., 3.};
    ///\todo test that the norm of the 4D wedge stuff vector equals the area of
    /// the tetrahedron formed by nodes {0, 0, 0}, v1 and v2 and v3
    INFO("direction of wedge stuff vector");
    CHECK(std::abs((MiddleSizeMatrix({vec4D0, vec4D1, vec4D2})
                        .computeWedgeStuffVector()) *
                   vec4D0) < 1e-12);
    INFO("direction of wedge stuff vector");
    CHECK(std::abs((MiddleSizeMatrix({vec4D0, vec4D1, vec4D2})
                        .computeWedgeStuffVector()) *
                   vec4D1) < 1e-12);
    INFO("direction of wedge stuff vector");
    CHECK(std::abs((MiddleSizeMatrix({vec4D0, vec4D1, vec4D2})
                        .computeWedgeStuffVector()) *
                   vec4D2) < 1e-12);
    INFO("direction of wedge stuff vector");
    CHECK(std::abs((MiddleSizeMatrix({vec4D0, vec4D1, vec4D3})
                        .computeWedgeStuffVector()) *
                   vec4D0) < 1e-12);
    INFO("direction of wedge stuff vector");
    CHECK(std::abs((MiddleSizeMatrix({vec4D0, vec4D1, vec4D3})
                        .computeWedgeStuffVector()) *
                   vec4D1) < 1e-12);
    INFO("direction of wedge stuff vector");
    CHECK(std::abs((MiddleSizeMatrix({vec4D0, vec4D1, vec4D3})
                        .computeWedgeStuffVector()) *
                   vec4D3) < 1e-12);
    INFO("direction of wedge stuff vector");
    CHECK(std::abs((MiddleSizeMatrix({vec4D0, vec4D2, vec4D3})
                        .computeWedgeStuffVector()) *
                   vec4D0) < 1e-12);
    INFO("direction of wedge stuff vector");
    CHECK(std::abs((MiddleSizeMatrix({vec4D0, vec4D2, vec4D3})
                        .computeWedgeStuffVector()) *
                   vec4D2) < 1e-12);
    INFO("direction of wedge stuff vector");
    CHECK(std::abs((MiddleSizeMatrix({vec4D0, vec4D2, vec4D3})
                        .computeWedgeStuffVector()) *
                   vec4D3) < 1e-12);
    INFO("direction of wedge stuff vector");
    CHECK(std::abs((MiddleSizeMatrix({vec4D1, vec4D2, vec4D3})
                        .computeWedgeStuffVector()) *
                   vec4D1) < 1e-12);
    INFO("direction of wedge stuff vector");
    CHECK(std::abs((MiddleSizeMatrix({vec4D1, vec4D2, vec4D3})
                        .computeWedgeStuffVector()) *
                   vec4D2) < 1e-12);
    INFO("direction of wedge stuff vector");
    CHECK(std::abs((MiddleSizeMatrix({vec4D1, vec4D2, vec4D3})
                        .computeWedgeStuffVector()) *
                   vec4D3) < 1e-12);

    copy.concatenate(extra);
    INFO("Rows in a matrix");
    CHECK(copy.getNumberOfRows() == 4);
    INFO("Columns in a matrix");
    CHECK(copy.getNumberOfColumns() == 2);
    INFO("Size of a matrix");
    CHECK(copy.size() == 8);
    INFO("concatenate");
    CHECK(std::abs((copy)(0, 0) - 4.) < 1e-12);
    INFO("concatenate");
    CHECK(std::abs((copy)(1, 0) - 6.) < 1e-12);
    INFO("concatenate");
    CHECK(std::abs((copy)(0, 1) - 8.) < 1e-12);
    INFO("concatenate");
    CHECK(std::abs((copy)(1, 1) - 10.) < 1e-12);
    INFO("concatenate");
    CHECK(std::abs((copy)(2, 0) + 8.) < 1e-12);
    INFO("concatenate");
    CHECK(std::abs((copy)(3, 0) + 19.) < 1e-12);
    INFO("concatenate");
    CHECK(std::abs((copy)(2, 1) + 2.) < 1e-12);
    INFO("concatenate");
    CHECK(std::abs((copy)(3, 1) + 21.) < 1e-12);

    INFO("getColumn");
    CHECK(MiddleSizeMatrix({vec3D0, vec3D1}).getColumn(0).size() == 3);
    INFO("getColumn");
    CHECK(std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getColumn(0)) -
                    vec3D0)[0]) < 1e-12);
    INFO("getColumn");
    CHECK(std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getColumn(0)) -
                    vec3D0)[1]) < 1e-12);
    INFO("getColumn");
    CHECK(std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getColumn(0)) -
                    vec3D0)[2]) < 1e-12);
    INFO("getColumn");
    CHECK(MiddleSizeMatrix({vec3D0, vec3D1}).getColumn(1).size() == 3);
    INFO("getColumn");
    CHECK(std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getColumn(1)) -
                    vec3D1)[0]) < 1e-12);
    INFO("getColumn");
    CHECK(std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getColumn(1)) -
                    vec3D1)[1]) < 1e-12);
    INFO("getColumn");
    CHECK(std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getColumn(1)) -
                    vec3D1)[2]) < 1e-12);
    INFO("getRow");
    CHECK(MiddleSizeMatrix({vec3D0, vec3D1}).getRow(0).size() == 2);
    INFO("getRow");
    CHECK(std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getRow(0)))[0] - 0.) <
          1e-12);
    INFO("getRow");
    CHECK(std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getRow(0)))[1] - 3.) <
          1e-12);
    INFO("getRow");
    CHECK(MiddleSizeMatrix({vec3D0, vec3D1}).getRow(1).size() == 2);
    INFO("getRow");
    CHECK(std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getRow(1)))[0] - 1.) <
          1e-12);
    INFO("getRow");
    CHECK(std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getRow(1)))[1] - 4.) <
          1e-12);
    INFO("getRow");
    CHECK(MiddleSizeMatrix({vec3D0, vec3D1}).getRow(2).size() == 2);
    INFO("getRow");
    CHECK(std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getRow(2)))[0] - 2.) <
          1e-12);
    INFO("getRow");
    CHECK(std::abs(((MiddleSizeMatrix({vec3D0, vec3D1}).getRow(2)))[1] - 5.) <
          1e-12);

    ///\todo figure out a way to test a LU factorisation

    MiddleSizeVector duplicate = vec2;
    count0.solve(vec2);
    INFO("inverse and solve");
    CHECK(std::abs((vec2 - count0.inverse() * duplicate)[0]) < 1e-12);
    INFO("inverse and solve");
    CHECK(std::abs((vec2 - count0.inverse() * duplicate)[1]) < 1e-12);
    INFO("inverse");
    CHECK(std::abs((count0 - count0.inverse().inverse())[0]) < 1e-12);
    INFO("inverse");
    CHECK(std::abs((count0 - count0.inverse().inverse())[1]) < 1e-12);
    INFO("inverse");
    CHECK(std::abs((count0 - count0.inverse().inverse())[2]) < 1e-12);
    INFO("inverse");
    CHECK(std::abs((count0 - count0.inverse().inverse())[3]) < 1e-12);
    INFO("transpose");
    CHECK(std::abs((count0 - count0.transpose().transpose())[0]) < 1e-12);
    INFO("transpose");
    CHECK(std::abs((count0 - count0.transpose().transpose())[1]) < 1e-12);
    INFO("transpose");
    CHECK(std::abs((count0 - count0.transpose().transpose())[2]) < 1e-12);
    INFO("transpose");
    CHECK(std::abs((count0 - count0.transpose().transpose())[3]) < 1e-12);
    INFO("transpose");
    CHECK(std::abs(A23(0, 0) - A23.transpose()(0, 0)) < 1e-12);
    INFO("transpose");
    CHECK(std::abs(A23(1, 0) - A23.transpose()(0, 1)) < 1e-12);
    INFO("transpose");
    CHECK(std::abs(A23(0, 1) - A23.transpose()(1, 0)) < 1e-12);
    INFO("transpose");
    CHECK(std::abs(A23(1, 1) - A23.transpose()(1, 1)) < 1e-12);
    INFO("transpose");
    CHECK(std::abs(A23(0, 2) - A23.transpose()(2, 0)) < 1e-12);
    INFO("transpose");
    CHECK(std::abs(A23(1, 2) - A23.transpose()(2, 1)) < 1e-12);

    auto data = A23.data();
    INFO("data");
    CHECK(std::abs(data[0] - A23[0]) < 1e-12);
    INFO("data");
    CHECK(std::abs(data[1] - A23[1]) < 1e-12);
    INFO("data");
    CHECK(std::abs(data[2] - A23[2]) < 1e-12);
    INFO("data");
    CHECK(std::abs(data[3] - A23[3]) < 1e-12);
    INFO("data");
    CHECK(std::abs(data[4] - A23[4]) < 1e-12);
    INFO("data");
    CHECK(std::abs(data[5] - A23[5]) < 1e-12);
    data[2] = 17.3;
    INFO("data");
    CHECK(std::abs(A23[2] - 17.3) < 1e-12);

    std::cout << A32 << std::endl;
}

TEST_CASE("MiddleSizeMatrix solve lower", "[MiddleSizeMatrixUnitTest]") {
    LinearAlgebra::MiddleSizeMatrix L(2, 2);
    // Define a test matrix
    // [1 0; 1 2]
    L(0, 0) = 1.0;
    L(1, 0) = 1.0;
    L(1, 1) = 2.0;

    LinearAlgebra::MiddleSizeVector vec({3.0, 4.0});

    SECTION("Solve Lx = b") {
        L.solveLowerTriangular(vec, LinearAlgebra::Side::OP_LEFT,
                               LinearAlgebra::Transpose::NOT);
        CHECK(vec(0) == 3.0);
        CHECK(vec(1) == 0.5);
    }

    SECTION("Solve xL = b") {
        L.solveLowerTriangular(vec, hpgem::LinearAlgebra::Side::OP_RIGHT,
                               hpgem::LinearAlgebra::Transpose::NOT);
        CHECK(vec(0) == 1.0);
        CHECK(vec(1) == 2.0);
    }

    SECTION("Solve L^T x = b") {
        L.solveLowerTriangular(vec, hpgem::LinearAlgebra::Side::OP_LEFT,
                               hpgem::LinearAlgebra::Transpose::TRANSPOSE);
        CHECK(vec(0) == 1.0);
        CHECK(vec(1) == 2.0);
    }

    SECTION("Solve x L^T = b") {
        L.solveLowerTriangular(vec, hpgem::LinearAlgebra::Side::OP_RIGHT,
                               hpgem::LinearAlgebra::Transpose::TRANSPOSE);
        CHECK(vec(0) == 3.0);
        CHECK(vec(1) == 0.5);
    }
    // Reset for next
}

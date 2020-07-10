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

// naming convention: <Digit><ClassName>_UnitTest.cpp where <Digit> is a number
// that will make sure the unit tests are ordered such that the first failing
// unit test indicate the culprit class and other 'unit' tests may assume
// correct execution of all prior unit tests
#include "Geometry/Jacobian.h"
#include "Logger.h"

#include <cmath>

#define CATCH_CONFIG_MAIN
#include "../catch.hpp"

using namespace hpgem;
using Geometry::Jacobian;

TEST_CASE("210Jacobian_UnitTest", "[210Jacobian_UnitTest]") {

    // square jacobians and determinants

    Jacobian<1, 1> dim1;

    dim1(0, 0) = 3.14;
    INFO("1D determinant");
    CHECK((std::abs(dim1.determinant() - 3.14) < 1e-12));

    dim1(0, 0) = 0;
    INFO("1D determinant");
    CHECK((std::abs(dim1.determinant()) < 1e-12));

    dim1(0, 0) = -2.81;
    INFO("1D determinant");
    CHECK((std::abs(dim1.determinant() + 2.81) < 1e-12));

    Jacobian<2, 2> dim2;

    dim2(0, 0) = 1.;
    dim2(1, 0) = 2.;
    dim2(0, 1) = 3.;
    dim2(1, 1) = 4.;
    INFO("2D determinant");
    CHECK((std::abs(dim2.determinant() + 2.) < 1e-12));

    dim2(0, 0) = -1.;
    dim2(1, 0) = 2.;
    dim2(0, 1) = 3.;
    dim2(1, 1) = 4.;
    INFO("2D determinant");
    CHECK((std::abs(dim2.determinant() + 10.) < 1e-12));

    Jacobian<2, 2> copy(dim2);

    dim2(0, 0) = 1.;
    dim2(1, 0) = -2.;
    dim2(0, 1) = 3.;
    dim2(1, 1) = 4.;
    INFO("2D determinant");
    CHECK((std::abs(dim2.determinant() - 10.) < 1e-12));

    dim2(0, 0) = 1.;
    dim2(1, 0) = 2.;
    dim2(0, 1) = -3.;
    dim2(1, 1) = 4.;
    INFO("2D determinant");
    CHECK((std::abs(dim2.determinant() - 10.) < 1e-12));

    dim2(0, 0) = 1.;
    dim2(1, 0) = 2.;
    dim2(0, 1) = 3.;
    dim2(1, 1) = -4.;
    INFO("2D determinant");
    CHECK((std::abs(dim2.determinant() + 10.) < 1e-12));

    dim2(0, 0) = -1.;
    dim2(1, 0) = -2.;
    dim2(0, 1) = 3.;
    dim2(1, 1) = 4.;
    INFO("2D determinant");
    CHECK((std::abs(dim2.determinant() - 2.) < 1e-12));

    dim2(0, 0) = -1.;
    dim2(1, 0) = 2.;
    dim2(0, 1) = -3.;
    dim2(1, 1) = 4.;
    INFO("2D determinant");
    CHECK((std::abs(dim2.determinant() - 2.) < 1e-12));

    dim2(0, 0) = -1.;
    dim2(1, 0) = 2.;
    dim2(0, 1) = 3.;
    dim2(1, 1) = -4.;
    INFO("2D determinant");
    CHECK((std::abs(dim2.determinant() + 2.) < 1e-12));

    dim2(0, 0) = 1.;
    dim2(1, 0) = 2.;
    dim2(0, 1) = 2.;
    dim2(1, 1) = 4.;
    INFO("2D determinant");
    CHECK((std::abs(dim2.determinant()) < 1e-12));
    INFO("copy constructor");
    CHECK((std::abs(copy.determinant() + 10.) < 1e-12));

    Jacobian<2, 2> product;

    product = dim2.multiplyJacobiansInto(copy);
    INFO("multiply JacobiansInto - square matrixes");
    CHECK((product.getNCols() == 2 && product.getNRows() == 2 &&
           std::abs(product.determinant()) < 1e-12));

    Jacobian<3, 3> dim3;
    dim3(0, 0) = 1.;
    dim3(1, 0) = 2.;
    dim3(2, 0) = 3.;
    dim3(0, 1) = 4.;
    dim3(1, 1) = 5.;
    dim3(2, 1) = 6.;
    dim3(0, 2) = 7.;
    dim3(1, 2) = 8.;
    dim3(2, 2) = 9.;
    INFO("3D determinant");
    CHECK((std::abs(dim3.determinant()) < 1e-12));

    dim3(0, 0) = -1.;
    dim3(1, 0) = 2.;
    dim3(2, 0) = 3.;
    dim3(0, 1) = 4.;
    dim3(1, 1) = 5.;
    dim3(2, 1) = 6.;
    dim3(0, 2) = 7.;
    dim3(1, 2) = 8.;
    dim3(2, 2) = 9.;
    INFO("3D determinant");
    CHECK((std::abs(dim3.determinant() - 6.) < 1e-12));

    dim3(0, 0) = 1.;
    dim3(1, 0) = -2.;
    dim3(2, 0) = 3.;
    dim3(0, 1) = 4.;
    dim3(1, 1) = 5.;
    dim3(2, 1) = 6.;
    dim3(0, 2) = 7.;
    dim3(1, 2) = 8.;
    dim3(2, 2) = 9.;
    INFO("3D determinant");
    CHECK((std::abs(dim3.determinant() + 24.) < 1e-12));

    dim3(0, 0) = 1.;
    dim3(1, 0) = 2.;
    dim3(2, 0) = -3.;
    dim3(0, 1) = 4.;
    dim3(1, 1) = 5.;
    dim3(2, 1) = 6.;
    dim3(0, 2) = 7.;
    dim3(1, 2) = 8.;
    dim3(2, 2) = 9.;
    INFO("3D determinant");
    CHECK((std::abs(dim3.determinant() - 18.) < 1e-12));

    dim3(0, 0) = 1.;
    dim3(1, 0) = 2.;
    dim3(2, 0) = 3.;
    dim3(0, 1) = -4.;
    dim3(1, 1) = 5.;
    dim3(2, 1) = 6.;
    dim3(0, 2) = 7.;
    dim3(1, 2) = 8.;
    dim3(2, 2) = 9.;
    INFO("3D determinant");
    CHECK((std::abs(dim3.determinant() + 48.) < 1e-12));

    dim3(0, 0) = 1.;
    dim3(1, 0) = 2.;
    dim3(2, 0) = 3.;
    dim3(0, 1) = 4.;
    dim3(1, 1) = -5.;
    dim3(2, 1) = 6.;
    dim3(0, 2) = 7.;
    dim3(1, 2) = 8.;
    dim3(2, 2) = 9.;
    INFO("3D determinant");
    CHECK((std::abs(dim3.determinant() - 120.) < 1e-12));

    dim3(0, 0) = 1.;
    dim3(1, 0) = 2.;
    dim3(2, 0) = 3.;
    dim3(0, 1) = 4.;
    dim3(1, 1) = 5.;
    dim3(2, 1) = -6.;
    dim3(0, 2) = 7.;
    dim3(1, 2) = 8.;
    dim3(2, 2) = 9.;
    INFO("3D determinant");
    CHECK((std::abs(dim3.determinant() + 72.) < 1e-12));

    dim3(0, 0) = 1.;
    dim3(1, 0) = 2.;
    dim3(2, 0) = 3.;
    dim3(0, 1) = 4.;
    dim3(1, 1) = 5.;
    dim3(2, 1) = 6.;
    dim3(0, 2) = -7.;
    dim3(1, 2) = 8.;
    dim3(2, 2) = 9.;
    INFO("3D determinant");
    CHECK((std::abs(dim3.determinant() - 42.) < 1e-12));

    dim3(0, 0) = 1.;
    dim3(1, 0) = 2.;
    dim3(2, 0) = 3.;
    dim3(0, 1) = 4.;
    dim3(1, 1) = 5.;
    dim3(2, 1) = 6.;
    dim3(0, 2) = 7.;
    dim3(1, 2) = -8.;
    dim3(2, 2) = 9.;
    INFO("3D determinant");
    CHECK((std::abs(dim3.determinant() + 96.) < 1e-12));

    dim3(0, 0) = 1.;
    dim3(1, 0) = 2.;
    dim3(2, 0) = 3.;
    dim3(0, 1) = 4.;
    dim3(1, 1) = 5.;
    dim3(2, 1) = 6.;
    dim3(0, 2) = 7.;
    dim3(1, 2) = 8.;
    dim3(2, 2) = -9.;
    INFO("3D determinant");
    CHECK((std::abs(dim3.determinant() - 54.) < 1e-12));

    Jacobian<1, 2> rectangle21;
    Jacobian<2, 1> rectangle12;
    Jacobian<2, 3> rectangle32;
    Jacobian<3, 2> rectangle23;

    rectangle21(0, 0) = 1.;
    rectangle21(1, 0) = 2.;
    rectangle12(0, 0) = 3.;
    rectangle12(0, 1) = 4.;

    rectangle32(0, 0) = 1.;
    rectangle32(0, 1) = 2.;
    rectangle32(1, 0) = 3.;
    rectangle32(1, 1) = 4.;
    rectangle32(2, 0) = 5.;
    rectangle32(2, 1) = 6.;
    rectangle23(0, 0) = 7.;
    rectangle23(0, 1) = 8.;
    rectangle23(0, 2) = 9.;
    rectangle23(1, 0) = 10.;
    rectangle23(1, 1) = 11.;
    rectangle23(1, 2) = 12.;

    auto product2 = rectangle21.multiplyJacobiansInto(dim1);
    INFO("multiply JacobiansInto - square & rectangular matrixes");
    CHECK((product2.getNumberOfColumns() == 1 &&
           product2.getNumberOfRows() == 2));

    auto product4 = rectangle12.multiplyJacobiansInto(dim2);
    INFO("multiply JacobiansInto - square & rectangular matrixes");
    CHECK((product4.getNumberOfColumns() == 2 &&
           product4.getNumberOfRows() == 1));

    auto product6 = rectangle23.multiplyJacobiansInto(dim3);
    INFO("multiply JacobiansInto - square & rectangular matrixes");
    CHECK((product6.getNumberOfColumns() == 3 &&
           product6.getNumberOfRows() == 2));

    auto product8 = rectangle32.multiplyJacobiansInto(dim2);
    INFO("multiply JacobiansInto - square & rectangular matrixes");
    CHECK((product8.getNumberOfColumns() == 2 &&
           product8.getNumberOfRows() == 3));
}

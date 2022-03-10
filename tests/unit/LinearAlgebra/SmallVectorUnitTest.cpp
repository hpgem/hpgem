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
#include "LinearAlgebra/SmallVector.h"
#include "LinearAlgebra/SmallMatrix.h"
#include "Logger.h"

#include "../catch.hpp"

using namespace hpgem;
using LinearAlgebra::SmallVector;
using LinearAlgebra::SmallVectorC;
using namespace std::complex_literals;

TEST_CASE("Empty vector is zero", "[SmallVectorUnitTest]") {
    // Test that the default constructor actually zeros the entries. When the
    // entries are not zero-ed, they will contain whatever data is left by the
    // previous user of that piece of memory.
    {
        // To make it more likely that we use polluted memory (i.e. memory that
        // is not accidentally containing zeros), we first create the same
        // objects and fill them with non-zero values.
        LinearAlgebra::SmallVector<1> vecR;
        LinearAlgebra::SmallVectorC<1> vecC;
        vecR[0] = 1.0;
        vecC[0] = 1.0;
    }
    {
        // Test that the default initialization will zero the entries.
        LinearAlgebra::SmallVector<1> vecR;
        LinearAlgebra::SmallVectorC<1> vecC;
        for (std::size_t i = 0; i < 1; ++i) {
            REQUIRE(vecR[i] == 0.0);
            REQUIRE(vecC[i] == 0.0);
        }
    }
}

TEST_CASE("Cross product 3D", "[SmallVectorUnitTest]") {
    SmallVector<3> x({1, 0, 0}), y({0, 1, 0}), z({0, 0, 1});

    // Checking the orientation.
    INFO("Cross xy = z");
    CHECK((x.crossProduct(y) - z).l2Norm() < 1e-12);
    INFO("Cross yz = x");
    CHECK((y.crossProduct(z) - x).l2Norm() < 1e-12);
    INFO("Cross zx = y");
    CHECK((z.crossProduct(x) - y).l2Norm() < 1e-12);

    // Checking sign inversion
    INFO("Cross yx = -z");
    CHECK((y.crossProduct(x) + z).l2Norm() < 1e-12);
    INFO("Cross zy = -x");
    CHECK((z.crossProduct(y) + x).l2Norm() < 1e-12);
    INFO("Cross xz = -y");
    CHECK((x.crossProduct(z) + y).l2Norm() < 1e-12);

    // Checking parallel vectors
    SmallVector<3> w({1, 3, -1});
    INFO("Cross: w x w = 0");
    CHECK(w.crossProduct(w).l2Norm() < 1e-12);
    INFO("Cross: w x -2w = 0");
    CHECK(w.crossProduct(-2 * w).l2Norm() < 1e-12);

    // Checking random example from internet
    SmallVector<3> a({3, -3, 1}), b({4, 9, 2}), ab({-15, -2, 39});
    INFO("Cross example 1");
    CHECK((a.crossProduct(b) - ab).l2Norm() < 1e-12);

    // Checking equivalence with the wedge stuff factor from SmallMatrix
    LinearAlgebra::SmallMatrix<3, 2> stuff({x, y});
    INFO("cross equals matrix wedge stuff factor");
    CHECK((x.crossProduct(y) - stuff.computeWedgeStuffVector()).l2Norm() <
          1e-12);

    stuff = LinearAlgebra::SmallMatrix<3, 2>({a, b});
    INFO("cross equals matrix wedge stuff factor");
    CHECK((a.crossProduct(b) - stuff.computeWedgeStuffVector()).l2Norm() <
          1e-12);
}

// Run 3D case first, as this uses that
TEST_CASE("Cross product 2D", "[SmallVectorUnitTest]") {
    SmallVector<2> x({1, 0}), y({0, 1});

    // Checking basic unit vectors crosses
    INFO("Cross xy = 1");
    CHECK((x.crossProduct(y) - x).l2Norm() < 1e-12);
    INFO("Cross yx = -1");
    CHECK((y.crossProduct(x) + x).l2Norm() < 1e-12);
    INFO("Cross xx = 0");
    CHECK(x.crossProduct(x).l2Norm() < 1e-12);

    // Check that it matches with 3D cross products
    SmallVector<2> a({2, 4}), b({-3, 5});
    SmallVector<2> abCross2 = a.crossProduct(b);
    SmallVector<3> abCross3 = a.append(0).crossProduct(b.append(0));
    INFO("Cross in 2D and 3D match");
    CHECK(std::abs(abCross2[0] - abCross3[2]) < 1e-12);
}

TEST_CASE("SmallVectorUnitTest", "[SmallVectorUnitTest]") {
    //(operator[], operator() and size() don't have a stand-alone test, but are
    // used throughout the test in multiple assertions

    // various constructors
    double data[] = {1., 2., 3.};
    SmallVector<0> x0;
    SmallVector<1> x1, y1(x1);
    SmallVector<3> fromArray(data);
    std::vector<double> vec{1.0, 3, 5};
    SmallVector<3> fromVector(vec);
    SmallVector<4> destroy, convenient({5., 6., 8., 9.});
    data[0] = 4;
    INFO("Constructor creates a vector of the wrong size");
    CHECK(x0.size() == 0);
    INFO("Constructor creates a vector of the wrong size");
    CHECK(fromVector.size() == 3);
    INFO("Constructor creates a vector of the wrong size");
    CHECK(x1.size() == 1);
    INFO("Constructor creates a vector of the wrong size");
    CHECK(y1.size() == 1);
    INFO("Copied array does not pass equality test");
    CHECK(y1 == x1);
    INFO("Copy constructor does not copy values");
    CHECK(std::abs(y1[0] - x1[0]) < 1e-12);
    INFO("Constructor creates a vector of the wrong size");
    CHECK(destroy.size() == 4);
    INFO("Constructor creates a vector of the wrong size");
    CHECK(convenient.size() == 4);
    INFO("Constructor creates a vector of the wrong size");
    CHECK(fromArray.size() == 3);
    destroy[2] = 4;
    SmallVector<4> moved(std::move(destroy));
    INFO("Constructor creates a vector of the wrong size");
    CHECK(moved.size() == 4);
    INFO("Constructor from array does not copy!");
    CHECK(std::abs(moved[2] - 4.) < 1e-12);
    INFO("Constructor from array does not copy!");
    CHECK(std::abs(fromArray(0) - 1.) < 1e-12);
    INFO("Constructor from array does not copy!");
    CHECK(std::abs(fromArray(1) - 2.) < 1e-12);
    INFO("Constructor from array does not copy!");
    CHECK(std::abs(fromArray[2] - 3.) < 1e-12);
    INFO("Initializer list constructor does not copy!");
    CHECK(std::abs(convenient[0] - 5.) < 1e-12);
    INFO("Initializer list constructor does not copy!");
    CHECK(std::abs(convenient[1] - 6.) < 1e-12);
    INFO("Initializer list constructor does not copy!");
    CHECK(std::abs(convenient[2] - 8.) < 1e-12);
    INFO("Initializer list constructor does not copy!");
    CHECK(std::abs(convenient[3] - 9.) < 1e-12);
    SmallVector<4> assigned = moved;
    INFO("Constructor creates a vector of the wrong size");
    CHECK(assigned.size() == 4);
    INFO("Constructor from array does not copy!");
    CHECK(std::abs(assigned[2] - 4.) < 1e-12);
    destroy = convenient;
    INFO("Assignment operator creates a vector of the wrong size");
    CHECK(destroy.size() == 4);
    INFO("Assignment operator from array does not copy!");
    CHECK(std::abs(convenient(0) - 5.) < 1e-12);
    INFO("Assignment operator from array does not copy!");
    CHECK(std::abs(convenient[1] - 6.) < 1e-12);
    INFO("Assignment operator from array does not copy!");
    CHECK(std::abs(convenient(2) - 8.) < 1e-12);
    INFO("Assignment operator from array does not copy!");
    CHECK(std::abs(convenient(3) - 9.) < 1e-12);

    SmallVector<2> p2 = {{0.8, 1.8}};
    SmallVector<2> pc2;
    SmallVector<2> pv2 = {{0.8, 0.8}};
    SmallVector<2> pw2 = {{1.8, 1.8}};
    SmallVector<2> px2 = {{1.8, 0.8}};
    SmallVector<2> py2 = {{0.6, 0.7}};
    const SmallVector<2> pr2 = pc2 = p2;
    INFO("Constructor creates a vector of the wrong size");
    CHECK(pr2.size() == 2);
    INFO("Constructor creates a vector of the wrong size");
    CHECK(pc2.size() == 2);
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("assignment operator");
        CHECK((std::abs(pc2[i] - 0.8 - i) < 1e-12));
        INFO("assignment operator");
        CHECK((std::abs(pr2[i] - 0.8 - i) < 1e-12));
    }
    INFO("equality operator");
    CHECK((pr2 == pc2 && pc2 == pr2 && pc2 == p2 &&
           !(pr2 == pv2 || pv2 == pr2 || p2 == pv2)));
    INFO("equality operator");
    CHECK(!(pr2 == pw2 || pw2 == pr2 || p2 == pw2 || pr2 == px2 || px2 == pr2 ||
            p2 == px2 || pr2 == py2 || py2 == pr2 || p2 == py2));
    pc2 += p2;
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("increment operator");
        CHECK((std::abs(pc2[i] - 1.6 - 2 * i) < 1e-12));
    }
    pc2 -= pv2;
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("decrement operator");
        CHECK((std::abs(pc2[i] - 0.8 - 2 * i) < 1e-12));
    }
    pc2 *= 4.;
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("multiply operator");
        CHECK((std::abs(pc2[i] - 3.2 - 8 * i) < 1e-12));
    }
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("multiplication");
        CHECK((std::abs((pc2 * -0.25)[i] + 0.8 + 2 * i) < 1e-12));
    }
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("addition");
        CHECK((std::abs((pc2 + pv2)[i] - 4. - 8 * i) < 1e-12));
    }
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("subtraction");
        CHECK((std::abs((pv2 - pc2)[i] + 2.4 + 8 * i) < 1e-12));
    }
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("unary -");
        CHECK((std::abs((-pc2)[i] + 3.2 + 8 * i) < 1e-12));
    }
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("left multiplication");
        CHECK((std::abs((0.25 * pc2)[i] - 0.8 - 2 * i) < 1e-12));
    }
    INFO("in-product");
    CHECK((std::abs((pv2 * pc2) - 1.44 * 8) < 1e-12));
    pc2 /= 4.;
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("divide operator");
        CHECK((std::abs(pc2[i] - 0.8 - 2 * i) < 1e-12));
    }
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("division");
        CHECK((std::abs((pc2 / -0.25)[i] + 3.2 + 8 * i) < 1e-12));
    }
    pc2.axpy(-3., p2);
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("combined multiply and addition (y = ax + y)");
        CHECK((std::abs(pc2[i] + 1.6 + i) < 1e-12));
    }

    x0.data();

    SmallVector<4> magicCopy(convenient.data());
    INFO("Assignment operator from array does not copy!");
    CHECK(std::abs(magicCopy(0) - 5.) < 1e-12);
    INFO("Assignment operator from array does not copy!");
    CHECK(std::abs(magicCopy(1) - 6.) < 1e-12);
    INFO("Assignment operator from array does not copy!");
    CHECK(std::abs(magicCopy(2) - 8.) < 1e-12);
    INFO("Assignment operator from array does not copy!");
    CHECK(std::abs(magicCopy(3) - 9.) < 1e-12);
    convenient.data()[3] = 2.;
    INFO("Assignment operator from array does not copy!");
    CHECK(std::abs(magicCopy(3) - 9.) < 1e-12);
    INFO("Assignment operator from array does not copy!");
    CHECK(std::abs(convenient(3) - 2.) < 1e-12);

    std::cout << pc2 << convenient << std::endl;
}

TEST_CASE("L2Norm", "[SmallVectorUnitTest]") {
    // Following tests use integers and should therefore be exact.
    SmallVector<2> v1 = {2.0, 3.0};
    CHECK(v1.l2NormSquared() == 13);
    SmallVectorC<2> vc1 = {2.0 + 3.0i, 4.0 + 5.0i};
    CHECK(vc1.l2NormSquared() == 4 + 9 + 16 + 25);
}

TEST_CASE("Complex Inner product", "[SmallVectorUnitTest]") {
    SmallVectorC<2> v1 = {1.0 + 2.0i, 1.0};
    SmallVectorC<2> v2 = {1.0i, 3.0 + 2.0i};
    // Expected value: 5-3i
    // (1+2i) * conj(1i) = 2-1i
    // 1 * conj(3+2i)    = 3-2i
    CHECK(v1 * v2 == 5.0 - 3.0i);
    // Flipping the order conjugates the output
    CHECK(v2 * v1 == 5.0 + 3.0i);

    // Consistency with l2 norm squared
    INFO("l2 norm <-> inner product equivalence");
    CHECK(v1 * v1 == 1.0 + 4.0 + 1.0);
    CHECK(v2 * v2 == 1.0 + 4.0 + 9.0);
}

TEST_CASE("Conversion complex <-> real", "[SmallVectorUnitTest]") {
    SmallVector<1> vreal = {1.0};
    SmallVectorC<1> vcreal = {1.0 + 0i};

    INFO("Conversion real -> complex");
    CHECK(SmallVectorC<1>(vreal) == vcreal);

    INFO("Conversion complex -> real");
    CHECK(SmallVector<1>(vcreal) == vreal);
}

TEST_CASE("Vector real imag", "[SmallVectorUnitTest]") {
    SmallVectorC<2> cvector = {1.0 + 2.0i, 3.0 + 4.0i};
    SmallVector<2> realPart = {1.0, 3.0};
    SmallVector<2> imagPart = {2.0, 4.0};
    SmallVector<2> zeroVec = {0.0, 0.0};

    {
        INFO("Check complex vector")
        CHECK(cvector.real() == realPart);
        CHECK(cvector.imag() == imagPart);
    }
    {
        INFO("Check real vector");
        CHECK(realPart.real() == realPart);
        CHECK(realPart.imag() == zeroVec);
    }
}
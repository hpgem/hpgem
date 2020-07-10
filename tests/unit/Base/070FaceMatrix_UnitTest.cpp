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

#include "Logger.h"

// Test the class FaceMatrix.
#include <iostream>

#include "Base/FaceMatrix.h"

#define CATCH_CONFIG_MAIN
#include "../catch.hpp"

using namespace hpgem;
TEST_CASE("070FaceMatrix_UnitTest", "[070FaceMatrix_UnitTest]") {
    // std::cout << "Test if FaceMatrix works.\n";

    const Base::Side sL = Base::Side::LEFT;
    const Base::Side sR = Base::Side::RIGHT;

    std::size_t nDOFLeft = 3;
    std::size_t nDOFRight = 5;

    Base::FaceMatrix F(nDOFLeft, nDOFRight);
    for (std::size_t i = 0; i < nDOFLeft + nDOFRight; i++) {
        for (std::size_t j = 0; j < nDOFLeft + nDOFRight; j++) {
            F(i, j) = i * 10 + j;
        }
    }

    Base::FaceMatrix F2(F);
    F2.axpy(100, F);

    Base::FaceMatrix F3;
    F3.resize(nDOFLeft, nDOFRight);
    F3 = F;
    F3 *= 100;
    F3 += F;

    Base::FaceMatrix F4;
    F4.resize(nDOFLeft, nDOFRight);
    F4.setEntireMatrix(F3.getEntireMatrix());

    Base::FaceMatrix F5;
    F5.resize(F4.getNumberOfDegreesOfFreedom(sL),
              F4.getNumberOfDegreesOfFreedom(sR));
    F5.setElementMatrix(F4.getElementMatrix(sL, sL), sL, sL);
    F5.setElementMatrix(F4.getElementMatrix(sL, sR), sL, sR);
    F5.setElementMatrix(F4.getElementMatrix(sR, sL), sR, sL);
    F5.setElementMatrix(F4.getElementMatrix(sR, sR), sR, sR);

    Base::Side iS = sL;
    Base::Side jS = sR;
    std::size_t iVB = 0;
    std::size_t jVB = 0;
    LinearAlgebra::MiddleSizeMatrix::type z = 0;
    for (std::size_t i = 0; i < nDOFLeft + nDOFRight; i++) {
        if (i < nDOFLeft) {
            iS = Base::Side::LEFT;
            iVB = i;
        } else {
            iS = Base::Side::RIGHT;
            iVB = i - nDOFLeft;
        }

        for (std::size_t j = 0; j < nDOFLeft + nDOFRight; j++) {
            if (j < nDOFLeft) {
                jS = Base::Side::LEFT;
                jVB = j;
            } else {
                jS = Base::Side::RIGHT;
                jVB = j - nDOFLeft;
            }

            z = i * 10 + j;
            z = z * 100. + z;

            logger.assert_always(F2.getElementMatrix(iS, jS)(iVB, jVB) == z,
                                 "Face matrix is wrong");
            logger.assert_always(F3(i, j) == z, "Face matrix is wrong");
            logger.assert_always(F4.getEntireMatrix()(i, j) == z,
                                 "Face matrix is wrong");
            logger.assert_always(F5(i, j) == z, "Face matrix is wrong");
        }
    }

    // std::cout << "Matrix F2:\n" << F2.getEntireMatrix() << "\n";
}

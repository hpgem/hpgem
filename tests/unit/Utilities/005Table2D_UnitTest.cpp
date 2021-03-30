/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2020, University of Twente
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

#include "Utilities/Table2D.h"
#include "../catch.hpp"

using namespace hpgem;

TEST_CASE("Resizing", "[Table2D]") {
    Utilities::Table2D<int> dummy;
    REQUIRE(dummy.getNumberOfRows() == 0);
    REQUIRE(dummy.getNumberOfColumns() == 0);
    REQUIRE(dummy.getSize() == 0);

    dummy.resize(2, 3);
    INFO("Resizing row count");
    CHECK(dummy.getNumberOfRows() == 2);
    INFO("Resizing column count");
    CHECK(dummy.getNumberOfColumns() == 3);
    INFO("Resizing total size")
    CHECK(dummy.getSize() == 2 * 3);
}

TEST_CASE("Entry values", "[Table2D]") {
    const std::size_t ROWS = 2;
    const std::size_t COLS = 3;
    const int START_VALUE = 42;

    Utilities::Table2D<int> table(ROWS, COLS, START_VALUE);
    REQUIRE(table.getNumberOfRows() == ROWS);
    REQUIRE(table.getNumberOfColumns() == COLS);
    for (std::size_t i = 0; i < table.getSize(); ++i) {
        REQUIRE(table[i] == START_VALUE);
    }

    SECTION("Filling") {
        table.fill(31);
        for (std::size_t i = 0; i < table.getSize(); ++i) {
            REQUIRE(table[i] == 31);
        }
    }

    SECTION("Filling 2") {
        table = 20;
        for (std::size_t i = 0; i < table.getSize(); ++i) {
            REQUIRE(table[i] == 20);
        }
    }

    SECTION("Setting single value") {
        table(1, 1) = 2;
        for (std::size_t i = 0; i < ROWS; ++i) {
            for (std::size_t j = 0; j < COLS; ++j) {
                if (i == 1 && j == 1) {
                    REQUIRE(table(i, j) == 2);
                } else {
                    REQUIRE(table(i, j) == START_VALUE);
                }
            }
        }
    }
}

TEST_CASE("Bool table", "[Table2D]") {
    // Had problems with assigning boolean values in the past due to how vector
    // stores booleans as bitfields.
    Utilities::Table2D<bool> boolTable(2, 2, true);
    boolTable(0, 1) = true;
    REQUIRE(boolTable(0, 1) == true);
}
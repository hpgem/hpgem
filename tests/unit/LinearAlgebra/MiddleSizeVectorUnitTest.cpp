/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, University of Twente
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <iostream>
#include <cstdlib>
#include <cmath>
#include "LinearAlgebra/MiddleSizeVector.h"
#include "Logger.h"

using LinearAlgebra::MiddleSizeVector;

int main(int argc, char* argv[])
{
    //(operator[], operator() and size() don't have a stand-alone test, but are used throughout the test in multiple assertions

    //various constructors
    double data[] = {1., 2., 3.};
    MiddleSizeVector x0, x1(1), y1(x1), destroy(4), fromArray(data,3), convenient({5.,6.,7.,8.,9.});
    data[0] = 4;
    logger.assert_always(x0.size() == 0, "Constructor creates a vector of the wrong size");
    logger.assert_always(x1.size() == 1, "Constructor creates a vector of the wrong size");
    logger.assert_always(y1.size() == 1, "Constructor creates a vector of the wrong size");
    logger.assert_always(y1 == x1, "Copied array does not pass equality test");
    logger.assert_always(std::abs(y1[0] - x1[0]) < 1e-12, "Copy constructor does not copy values");
    logger.assert_always(destroy.size() == 4, "Constructor creates a vector of the wrong size");
    logger.assert_always(convenient.size() == 5, "Constructor creates a vector of the wrong size");
    logger.assert_always(fromArray.size() == 3, "Constructor creates a vector of the wrong size");
    destroy[2] = 4;
    MiddleSizeVector moved(std::move(destroy));
    logger.assert_always(moved.size() == 4, "Constructor creates a vector of the wrong size");
    logger.assert_always(std::abs(moved[2] - 4.) < 1e-12, "Constructor from array does not copy!");
    logger.assert_always(std::abs(fromArray(0) - 1.) < 1e-12, "Constructor from array does not copy!");
    logger.assert_always(std::abs(fromArray(1) - 2.) < 1e-12, "Constructor from array does not copy!");
    logger.assert_always(std::abs(fromArray[2] - 3.) < 1e-12, "Constructor from array does not copy!");
    logger.assert_always(std::abs(convenient[0] - 5.) < 1e-12, "Initializer list constructor does not copy!");
    logger.assert_always(std::abs(convenient[1] - 6.) < 1e-12, "Initializer list constructor does not copy!");
    logger.assert_always(std::abs(convenient[2] - 7.) < 1e-12, "Initializer list constructor does not copy!");
    logger.assert_always(std::abs(convenient[3] - 8.) < 1e-12, "Initializer list constructor does not copy!");
    logger.assert_always(std::abs(convenient[4] - 9.) < 1e-12, "Initializer list constructor does not copy!");
    MiddleSizeVector assigned = moved;
    logger.assert_always(assigned.size() == 4, "Constructor creates a vector of the wrong size");
    logger.assert_always(std::abs(assigned[2] - 4.) < 1e-12, "Constructor from array does not copy!");
    destroy = convenient;
    logger.assert_always(destroy.size() == 5, "Assignment operator creates a vector of the wrong size");
    logger.assert_always(std::abs(convenient(0) - 5.) < 1e-12, "Assignment operator from array does not copy!");
    logger.assert_always(std::abs(convenient(1) - 6.) < 1e-12, "Assignment operator from array does not copy!");
    logger.assert_always(std::abs(convenient(2) - 7.) < 1e-12, "Assignment operator from array does not copy!");
    logger.assert_always(std::abs(convenient(3) - 8.) < 1e-12, "Assignment operator from array does not copy!");
    logger.assert_always(std::abs(convenient[4] - 9.) < 1e-12, "Assignment operator from array does not copy!");

    //resize
    x1.resize(5);
    moved.resize(3);
    logger.assert_always(std::abs(moved[2] - 4.) < 1e-12, "Resizing destroys values");
    moved(2) = 3;
    logger.assert_always(moved.size() == 3, "Resizing does not resize properly");
    logger.assert_always(x1.size() == 5, "Resizing does not resize properly");
    logger.assert_always(y1.size() == 1, "Constructor creates a vector of the wrong size");
    logger.assert_always(assigned.size() == 4, "Constructor creates a vector of the wrong size");
    logger.assert_always(std::abs(assigned[2] - 4.) < 1e-12, "Resizing destroys values");

    MiddleSizeVector p2 = {0.8, 1.8};
    MiddleSizeVector pc2(7);
    MiddleSizeVector pv2 = {0.8, 0.8};
    MiddleSizeVector pw2 = {1.8, 1.8};
    MiddleSizeVector px2 = {1.8, 0.8};
    MiddleSizeVector py2 = {0.6, 0.7};
    const MiddleSizeVector pr2 = pc2 = p2;
    logger.assert_always(pr2.size() == 2, "Constructor creates a vector of the wrong size");
    logger.assert_always(pc2.size() == 2, "Constructor creates a vector of the wrong size");
    for (std::size_t i = 0; i < 2; ++i)
    {
        logger.assert_always((std::abs(pc2[i] - 0.8 - i) < 1e-12), "assignment operator");
        logger.assert_always((std::abs(pr2[i] - 0.8 - i) < 1e-12), "assignment operator");
    }
    logger.assert_always((pr2 == pc2 && pc2 == pr2 && pc2 == p2 && !(pr2 == pv2 || pv2 == pr2 || p2 == pv2 || p2 == convenient || convenient == p2)), "equality operator");
    logger.assert_always(!(pr2 == pw2 || pw2 == pr2 || p2 == pw2 || pr2 == px2 || px2 == pr2 || p2 == px2 || pr2 == py2 || py2 == pr2 || p2 == py2), "equality operator");
    pc2 += p2;
    for (std::size_t i = 0; i < 2; ++i)
    {
        logger.assert_always((std::abs(pc2[i] - 1.6 - 2 * i) < 1e-12), "increment operator");
    }
    pc2 -= pv2;
    for (std::size_t i = 0; i < 2; ++i)
    {
        logger.assert_always((std::abs(pc2[i] - 0.8 - 2 * i) < 1e-12), "decrement operator");
    }
    pc2 *= 4.;
    for (std::size_t i = 0; i < 2; ++i)
    {
        logger.assert_always((std::abs(pc2[i] - 3.2 - 8 * i) < 1e-12), "multiply operator");
    }
    for (std::size_t i = 0; i < 2; ++i)
    {
        logger.assert_always((std::abs((pc2 * -0.25)[i] + 0.8 + 2 * i) < 1e-12), "multiplication");
    }
    for (std::size_t i = 0; i < 2; ++i)
    {
        logger.assert_always((std::abs((pc2 + pv2)[i] - 4. - 8 * i) < 1e-12), "addition");
    }
    for (std::size_t i = 0; i < 2; ++i)
    {
        logger.assert_always((std::abs((pv2 - pc2)[i] + 2.4 + 8 * i) < 1e-12), "subtraction");
    }
    for (std::size_t i = 0; i < 2; ++i)
    {
        logger.assert_always((std::abs((-pc2)[i] + 3.2 + 8 * i) < 1e-12), "unary -");
    }
    for (std::size_t i = 0; i < 2; ++i)
    {
        logger.assert_always((std::abs((0.25 * pc2)[i] - 0.8 - 2 * i) < 1e-12), "left multiplication");
    }
    logger.assert_always((std::abs((pv2 * pc2) - 1.44 * 8) < 1e-12), "in-product");
    pc2 /= 4.;
    for (std::size_t i = 0; i < 2; ++i)
    {
        logger.assert_always((std::abs(pc2[i] - 0.8 - 2 * i) < 1e-12), "divide operator");
    }
    for (std::size_t i = 0; i < 2; ++i)
    {
        logger.assert_always((std::abs((pc2 / -0.25)[i] + 3.2 + 8 * i) < 1e-12), "division");
    }
    pc2.axpy(-3., p2);
    for (std::size_t i = 0; i < 2; ++i)
    {
        logger.assert_always((std::abs(pc2[i] + 1.6 + i) < 1e-12), "combined multiply and addition (y = ax + y)");
    }

    x0.data();

    MiddleSizeVector magicCopy(convenient.data(), convenient.size());
    logger.assert_always(std::abs(magicCopy(0) - 5.) < 1e-12, "Assignment operator from array does not copy!");
    logger.assert_always(std::abs(magicCopy(1) - 6.) < 1e-12, "Assignment operator from array does not copy!");
    logger.assert_always(std::abs(magicCopy(2) - 7.) < 1e-12, "Assignment operator from array does not copy!");
    logger.assert_always(std::abs(magicCopy(3) - 8.) < 1e-12, "Assignment operator from array does not copy!");
    logger.assert_always(std::abs(magicCopy[4] - 9.) < 1e-12, "Assignment operator from array does not copy!");
    convenient.data()[3] = 2.;
    logger.assert_always(std::abs(convenient(3) - 2.) < 1e-12, "Assignment operator from array does not copy!");


    std::cout << pc2 << convenient << std::endl;
}


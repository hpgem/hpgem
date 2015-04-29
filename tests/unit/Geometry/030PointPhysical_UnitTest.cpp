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

//naming convention: <Digit><ClassName>_UnitTest.cpp where <Digit> is a number that will make sure
//the unit tests are ordered such that the first failing unit test indicate the culprit class and
//other 'unit' tests may assume correct execution of all prior unit tests
#include "Geometry/PointPhysical.h"
#include <iostream>
#include "Logger.h"
#include <cmath>
using Geometry::PointPhysical;

int main()
{
    
    double *coord0 = nullptr;
    double coord1[] = {1.1};
    double coord2[] = {1.2, 2.2};
    double coord3[] = {1.3, 2.3, 3.3};
    double coord4[] = {1.4, 2.4, 3.4, 4.4};
    
    LinearAlgebra::NumericalVector vec0(coord0, 0), vec1(coord1, 1), vec2(coord2, 2), vec3(coord3, 3), vec4(coord4, 4);
    
    Geometry::Point pb0(0), pb1(1), pb2(2), pb3(3), pb4(4);
    
    //testing constructors up to DIM=4
    PointPhysical p0(0), p1(1), p2(2), p3(3), p4(4), pp0(pb0), pp1(pb1), pp2(pb2), pp3(pb3), pp4(pb4), pv0(vec0), pv1(vec1), pv2(vec2), pv3(vec3), pv4(vec4); //OMG no int[] constructor
            
    logger.assert_always((p1[0] == 0.), "1D default constructor");
    for (std::size_t i = 0; i < 2; ++i)
    {
        logger.assert_always((p2[i] == 0.), "2D default constructor");
    }
    for (std::size_t i = 0; i < 3; ++i)
    {
        logger.assert_always((p3[i] == 0.), "3D default constructor");
    }
    for (std::size_t i = 0; i < 4; ++i)
    {
        logger.assert_always((p4[i] == 0.), "4D default constructor");
    }
    
    logger.assert_always((pp1[0] == 0.), "1D copy constructor");
    for (std::size_t i = 0; i < 2; ++i)
    {
        logger.assert_always((pp2[i] == 0.), "2D copy constructor");
    }
    for (std::size_t i = 0; i < 3; ++i)
    {
        logger.assert_always((pp3[i] == 0.), "3D copy constructor");
    }
    for (std::size_t i = 0; i < 4; ++i)
    {
        logger.assert_always((pp4[i] == 0.), "4D copy constructor");
    }
    
    logger.assert_always((std::abs(pv1[0] - 1.1) < 1e-12), "1D from NumericalVector constructor");
    for (std::size_t i = 0; i < 2; ++i)
    {
        logger.assert_always((std::abs(pv2[i] - 1.2 - i) < 1e-12), "2D from NumericalVector constructor");
    }
    for (std::size_t i = 0; i < 3; ++i)
    {
        logger.assert_always((std::abs(pv3[i] - 1.3 - i) < 1e-12), "3D from NumericalVector constructor");
    }
    for (std::size_t i = 0; i < 4; ++i)
    {
        logger.assert_always((std::abs(pv4[i] - 1.4 - i) < 1e-12), "4D from NumericalVector constructor");
    }
    
    //testing operators
    
    const PointPhysical pr0 = pp0 = pv0;
    const PointPhysical pr1 = pp1 = pv1;
    const PointPhysical pr2 = pp2 = pv2;
    const PointPhysical pr3 = pp3 = pv3;
    const PointPhysical pr4 = pp4 = pv4;
    
    pv0 * 6.;
    logger.assert_always((std::abs((pv1 * 5.)[0] - 5.5) < 1e-12), "1D multiplication");
    for (std::size_t i = 0; i < 2; ++i)
    {
        logger.assert_always((std::abs((pv2 * 4.)[i] - 4.8 - 4 * i) < 1e-12), "2D multiplication");
    }
    for (std::size_t i = 0; i < 3; ++i)
    {
        logger.assert_always((std::abs((pv3 * 3.)[i] - 3.9 - 3 * i) < 1e-12), "3D multiplication");
    }
    for (std::size_t i = 0; i < 4; ++i)
    {
        logger.assert_always((std::abs((pv4 * 2.)[i] - 2.8 - 2 * i) < 1e-12), "4D multiplication");
    }
    
    logger.assert_always(((pr0 * 0.) == p0), "0D multiplication");
    logger.assert_always(((pr1 * 0.) == p1), "1D multiplication");
    logger.assert_always(((pr2 * 0.) == p2), "2D multiplication");
    logger.assert_always(((pr3 * 0.) == p3), "3D multiplication");
    logger.assert_always(((pr4 * 0.) == p4), "4D multiplication");
    
    pp0 + pv0;
    logger.assert_always((std::abs((pp1 + pv1)[0] - 2.2) < 1e-12), "1D addition");
    for (std::size_t i = 0; i < 2; ++i)
    {
        logger.assert_always((std::abs((pp2 + pv2)[i] - 2.4 - 2 * i) < 1e-12), "2D addition");
    }
    for (std::size_t i = 0; i < 3; ++i)
    {
        logger.assert_always((std::abs((pp3 + pv3)[i] - 2.6 - 2 * i) < 1e-12), "3D addition");
    }
    for (std::size_t i = 0; i < 4; ++i)
    {
        logger.assert_always((std::abs((pp4 + pv4)[i] - 2.8 - 2 * i) < 1e-12), "4D addition");
    }
    
    pr0 + pv0;
    logger.assert_always((std::abs((pr1 + pv1)[0] - 2.2) < 1e-12), "1D addition");
    for (std::size_t i = 0; i < 2; ++i)
    {
        logger.assert_always((std::abs((pr2 + pv2)[i] - 2.4 - 2 * i) < 1e-12), "2D addition");
    }
    for (std::size_t i = 0; i < 3; ++i)
    {
        logger.assert_always((std::abs((pr3 + pv3)[i] - 2.6 - 2 * i) < 1e-12), "3D addition");
    }
    for (std::size_t i = 0; i < 4; ++i)
    {
        logger.assert_always((std::abs((pr4 + pv4)[i] - 2.8 - 2 * i) < 1e-12), "4D addition");
    }
    
    pp0 - pv0;
    logger.assert_always((std::abs((pp1 - pv1)[0]) < 1e-12), "1D subtraction");
    for (std::size_t i = 0; i < 2; ++i)
    {
        logger.assert_always((std::abs((pp2 - pv2)[i]) < 1e-12), "2D subtraction");
    }
    for (std::size_t i = 0; i < 3; ++i)
    {
        logger.assert_always((std::abs((pp3 - pv3)[i]) < 1e-12), "3D subtraction");
    }
    for (std::size_t i = 0; i < 4; ++i)
    {
        logger.assert_always((std::abs((pp4 - pv4)[i]) < 1e-12), "4D subtraction");
    }
    
    pr0 - pv0;
    logger.assert_always((std::abs((pr1 - pv1)[0]) < 1e-12), "1D subtraction");
    for (std::size_t i = 0; i < 2; ++i)
    {
        logger.assert_always((std::abs((pr2 - pv2)[i]) < 1e-12), "2D subtraction");
    }
    for (std::size_t i = 0; i < 3; ++i)
    {
        logger.assert_always((std::abs((pr3 - pv3)[i]) < 1e-12), "3D subtraction");
    }
    for (std::size_t i = 0; i < 4; ++i)
    {
        logger.assert_always((std::abs((pr4 - pv4)[i]) < 1e-12), "4D subtraction");
    }
    
    //and point already works so done
    
}


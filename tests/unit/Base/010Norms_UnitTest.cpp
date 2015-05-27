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
#include "Base/L2Norm.h"

#include <cmath>
#include <iostream>
#include "LinearAlgebra/MiddleSizeVector.h"
#include "Geometry/PointPhysical.h"
#include "Logger.h"

int main()
{
    
    double *test0(nullptr), test1[1], test2[2], test3[3];
    
    LinearAlgebra::MiddleSizeVector vec0D(test0, 0);
    Geometry::PointPhysical<0> point0D(vec0D);
    
    logger.assert_always(Base::L2Norm(vec0D) == 0, "0D case");
    logger.assert_always(Base::L2Norm(point0D) == 0, "0D case");
    
    test1[0] = 1;
    
    LinearAlgebra::MiddleSizeVector vec1D(test1, 1);
    Geometry::PointPhysical<1> point1D(vec1D);
    
    logger.assert_always(std::abs(Base::L2Norm(vec1D) - 1) < 1e-12, "1D case, positive");
    logger.assert_always(std::abs(Base::L2Norm(point1D) - 1) < 1e-12, "1D case, positive");
    
    vec1D[0] = -1;
    point1D[0] = -1;
    
    logger.assert_always(std::abs(Base::L2Norm(vec1D) - 1) < 1e-12, "1D case, negative");
    logger.assert_always(std::abs(Base::L2Norm(point1D) - 1) < 1e-12, "1D case, negative");
    
    vec1D[0] = 4.38573895783677438;
    point1D[0] = 4.38573895783677438;
    
    logger.assert_always(std::abs(Base::L2Norm(vec1D) - 4.38573895783677438) < 1e-12, "non-unit data");
    logger.assert_always(std::abs(Base::L2Norm(point1D) - 4.38573895783677438) < 1e-12, "non-unit data");
    
    test2[0] = 1;
    test2[1] = 1;
    
    LinearAlgebra::MiddleSizeVector vec2D(test2, 2);
    LinearAlgebra::SmallVector<2> smallVec2D(test2);
    Geometry::PointPhysical<2> point2D(vec2D);

    logger.assert_always(std::abs(Base::L2Norm(vec2D) - std::sqrt(2.)) < 1e-12, "2D case, positive");
    logger.assert_always(std::abs(Base::L2Norm(smallVec2D) - std::sqrt(2.)) < 1e-12, "2D case, positive");
    logger.assert_always(std::abs(Base::L2Norm(point2D) - std::sqrt(2.)) < 1e-12, "2D case, positive");

    vec2D[0] = -1;
    smallVec2D[0] = -1;
    point2D[0] = -1;

    logger.assert_always(std::abs(Base::L2Norm(vec2D) - std::sqrt(2.)) < 1e-12, "2D case, mix");
    logger.assert_always(std::abs(Base::L2Norm(smallVec2D) - std::sqrt(2.)) < 1e-12, "2D case, mix");
    logger.assert_always(std::abs(Base::L2Norm(point2D) - std::sqrt(2.)) < 1e-12, "2D case, mix");

    vec2D[1] = -1;
    smallVec2D[1] = -1;
    point2D[1] = -1;

    logger.assert_always(std::abs(Base::L2Norm(vec2D) - std::sqrt(2.)) < 1e-12, "2D case, negative");
    logger.assert_always(std::abs(Base::L2Norm(smallVec2D) - std::sqrt(2.)) < 1e-12, "2D case, negative");
    logger.assert_always(std::abs(Base::L2Norm(point2D) - std::sqrt(2.)) < 1e-12, "2D case, negative");
    
    test3[0] = 1;
    test3[1] = 1;
    test3[2] = 2;
    
    LinearAlgebra::MiddleSizeVector vec3D(test3, 3);
    Geometry::PointPhysical<3> point3D(vec3D);
    
    logger.assert_always(std::abs(Base::L2Norm(vec3D) - std::sqrt(6.)) < 1e-12, "3D case, positive");
    logger.assert_always(std::abs(Base::L2Norm(point3D) - std::sqrt(6.)) < 1e-12, "3D case, positive");
    
    vec3D[0] = -1;
    point3D[0] = -1;
    
    logger.assert_always(std::abs(Base::L2Norm(vec3D) - std::sqrt(6.)) < 1e-12, "3D case, mix");
    logger.assert_always(std::abs(Base::L2Norm(point3D) - std::sqrt(6.)) < 1e-12, "3D case, mix");
    
    vec3D[1] = -1;
    point3D[1] = -1;
    
    logger.assert_always(std::abs(Base::L2Norm(vec3D) - std::sqrt(6.)) < 1e-12, "3D case, mix");
    logger.assert_always(std::abs(Base::L2Norm(point3D) - std::sqrt(6.)) < 1e-12, "3D case, mix");
    
    vec3D[2] = -2;
    point3D[2] = -2;
    
    logger.assert_always(std::abs(Base::L2Norm(vec3D) - std::sqrt(6.)) < 1e-12, "3D case, negative");
    logger.assert_always(std::abs(Base::L2Norm(point3D) - std::sqrt(6.)) < 1e-12, "3D case, negative");
    
    return 0;
}


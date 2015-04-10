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
#include "Base/AssembleBasisFunctionSet.h"

#include "Base/BasisFunctionSet.h"
#include "Geometry/PointReference.h"
#include "LinearAlgebra/NumericalVector.h"
#include "Base/BaseBasisFunction.h"
#include "Logger.h"

int main()
{
    
    Base::BasisFunctionSet all1DbasisFunctions(5);
    Base::AssembleBasisFunctionSet_1D_Ord5_A0(all1DbasisFunctions);
    Geometry::PointReference point1D(1);
    for (std::size_t i = 0; i < all1DbasisFunctions.size(); ++i)
    {
        const Base::BaseBasisFunction* test = all1DbasisFunctions[i];
        for (point1D[0] = -1.5; point1D[0] < 1.51; point1D[0] += 0.1)
        {
            logger.assert_always((test->eval(point1D) == all1DbasisFunctions.eval(i, point1D)), "eval");
            logger.assert_always((test->evalDeriv0(point1D) == all1DbasisFunctions.evalDeriv(i, 0, point1D)), "derivative");
        }
    }
    
    Base::BasisFunctionSet all2DbasisFunctions(5);
    Base::AssembleBasisFunctionSet_2D_Ord5_A0(all2DbasisFunctions);
    Geometry::PointReference point2D(2);
    for (std::size_t i = 0; i < all2DbasisFunctions.size(); ++i)
    {
        const Base::BaseBasisFunction* test = all2DbasisFunctions[i];
        for (point2D[0] = -1.5; point2D[0] < 1.51; point2D[0] += 0.1)
        {
            for (point2D[1] = -1.5; point2D[1] < 1.51; point2D[1] += 0.1)
            {
                logger.assert_always((test->eval(point2D) == all2DbasisFunctions.eval(i, point2D)), "eval");
                logger.assert_always((test->evalDeriv0(point2D) == all2DbasisFunctions.evalDeriv(i, 0, point2D)), "derivative");
                logger.assert_always((test->evalDeriv1(point2D) == all2DbasisFunctions.evalDeriv(i, 1, point2D)), "derivative");
            }
        }
    }
    
    Base::BasisFunctionSet all3DbasisFunctions(5);
    Base::AssembleBasisFunctionSet_3D_Ord5_A0(all3DbasisFunctions);
    Geometry::PointReference point3D(3);
    for (std::size_t i = 0; i < all3DbasisFunctions.size(); ++i)
    {
        const Base::BaseBasisFunction* test = all3DbasisFunctions[i];
        for (point3D[0] = -1.5; point3D[0] < 1.51; point3D[0] += 0.1)
        {
            for (point3D[1] = -1.5; point3D[1] < 1.51; point3D[1] += 0.1)
            {
                for (point3D[2] = -1.5; point3D[2] < 1.51; point3D[2] += 0.1)
                {
                    logger.assert_always((test->eval(point3D) == all3DbasisFunctions.eval(i, point3D)), "eval");
                    logger.assert_always((test->evalDeriv0(point3D) == all3DbasisFunctions.evalDeriv(i, 0, point3D)), "derivative");
                    logger.assert_always((test->evalDeriv1(point3D) == all3DbasisFunctions.evalDeriv(i, 1, point3D)), "derivative");
                    logger.assert_always((test->evalDeriv2(point3D) == all3DbasisFunctions.evalDeriv(i, 2, point3D)), "derivative");
                }
            }
        }
    }
}


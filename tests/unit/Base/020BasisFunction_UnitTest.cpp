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

//I can only check that basisFunctions return finite numbers for some of the points that you can pass to them
//actually validating the correct implementation of the promised basisfunction requires computing said basisfunction
//however, the 'unit' tests for the quadrature rules and most of the self tests rely quite heavily on the correct functionality
//of the basisfunctions so some assurance can be gained from the entire test suite.
//Testing derivatives is much easier since it relies on using the same numerical approximation for all basis-function over and over again
//(but of course this is not as accurate as the actual derivative should be) -FB

#include "Base/AssembleBasisFunctionSet.h"

#include "Base/BasisFunctionSet.h"
#include "Geometry/PointReference.h"
#include "LinearAlgebra/MiddleSizeVector.h"
#include "Base/BaseBasisFunction.h"
#include "Logger.h"

#include <cmath>
int main()
{
    
    // 1D
    
    Base::BasisFunctionSet all1DbasisFunctions(5); //WARNING: this breaks the ordering of the unit tests, but it is basically the only way to collect all basisfunctions in an indexable way
    Base::AssembleBasisFunctionSet_1D_Ord5_A0(all1DbasisFunctions);
    Geometry::PointReference<1> point1D;
    LinearAlgebra::SmallVector<1> ret;
    for (std::size_t i = 0; i < all1DbasisFunctions.size(); ++i)
    {
        const Base::BaseBasisFunction* test = all1DbasisFunctions[i];
        for (point1D[0] = -1.5; point1D[0] < 1.51; point1D[0] += 0.1)
        {
            point1D[0] += -1.e-8;
            double x0 = test->eval((point1D));
            point1D[0] += 2.e-8;
            double x1 = test->eval((point1D));
            
            point1D[0] += -1e-8;
            ret = test->evalDeriv((point1D));
            logger.assert_always(std::abs(ret[0] - 5.e7 * (x1 - x0)) < 1e-5, "derivative");
            logger.assert_always(std::abs(test->evalDeriv0((point1D)) - 5.e7 * (x1 - x0)) < 1e-5, "derivative");
        }
    }
    
    // 2D
    
    Base::BasisFunctionSet all2DbasisFunctions(5); //WARNING: this breaks the ordering of the unit tests, but it is basically the only way to collect all basisfunctions in an indexable way
    Base::AssembleBasisFunctionSet_2D_Ord5_A1(all2DbasisFunctions);
    Geometry::PointReference<2> point2D;
    LinearAlgebra::SmallVector<2> ret2;
    for (std::size_t i = 0; i < all2DbasisFunctions.size(); ++i)
    {
        const Base::BaseBasisFunction* test = all2DbasisFunctions[i];
        for (point2D[0] = -1.5; point2D[0] < 1.51; point2D[0] += 0.3)
        {
            for (point2D[1] = -1.5; point2D[1] < 1.51; point2D[1] += 0.3)
            {
                point2D[0] += -1.e-8;
                double x0 = test->eval((point2D));
                point2D[0] += 2.e-8;
                double x1 = test->eval((point2D));
                
                point2D[0] += -1e-8;
                ret2 = test->evalDeriv((point2D));
                logger.assert_always(std::abs(ret2[0] - 5.e7 * (x1 - x0)) < 1e-5, "derivative");
                logger.assert_always(std::abs(test->evalDeriv0((point2D)) - 5.e7 * (x1 - x0)) < 1e-5, "derivative");
                
                point2D[1] += -1.e-8;
                x0 = test->eval((point2D));
                point2D[1] += 2.e-8;
                x1 = test->eval((point2D));
                
                point2D[1] += -1e-8;
                logger.assert_always(std::abs(ret2[1] - 5.e7 * (x1 - x0)) < 1e-5, "derivative");
                logger.assert_always(std::abs(test->evalDeriv1((point2D)) - 5.e7 * (x1 - x0)) < 1e-5, "derivative");
                
            }
        }
    }
    
    //3D
    
    Base::BasisFunctionSet all3DbasisFunctions(5); //WARNING: this breaks the ordering of the unit tests, but it is basically the only way to collect all basisfunctions in an indexable way
    Base::AssembleBasisFunctionSet_3D_Ord5_A1(all3DbasisFunctions);
    Geometry::PointReference<3> point3D;
    LinearAlgebra::SmallVector<3> ret3;
    for (std::size_t i = 0; i < all3DbasisFunctions.size(); ++i)
    {
        const Base::BaseBasisFunction* test = all3DbasisFunctions[i];
        for (point3D[0] = -1.5; point3D[0] < 1.51; point3D[0] += 0.6)
        {
            for (point3D[1] = -1.5; point3D[1] < 1.51; point3D[1] += 0.7)
            {
                for (point3D[2] = -1.5; point3D[2] < 1.51; point3D[2] += 0.8)
                {
                    point3D[0] += -1.e-8;
                    double x0 = test->eval((point3D));
                    point3D[0] += 2.e-8;
                    double x1 = test->eval((point3D));
                    
                    point3D[0] += -1e-8;
                    ret3 = test->evalDeriv((point3D));
                    logger.assert_always(std::abs(ret3[0] - 5.e7 * (x1 - x0)) < 1e-5, "derivative");
                    logger.assert_always(std::abs(test->evalDeriv0((point3D)) - 5.e7 * (x1 - x0)) < 1e-5, "derivative");
                    
                    point3D[1] += -1.e-8;
                    x0 = test->eval((point3D));
                    point3D[1] += 2.e-8;
                    x1 = test->eval((point3D));
                    
                    point3D[1] += -1e-8;
                    logger.assert_always(std::abs(ret3[1] - 5.e7 * (x1 - x0)) < 1e-5, "derivative");
                    logger.assert_always(std::abs(test->evalDeriv1((point3D)) - 5.e7 * (x1 - x0)) < 1e-5, "derivative");
                    
                    point3D[2] += -1.e-8;
                    x0 = test->eval((point3D));
                    point3D[2] += 2.e-8;
                    x1 = test->eval((point3D));
                    
                    point3D[2] += -1e-8;
                    logger.assert_always(std::abs(ret3[2] - 5.e7 * (x1 - x0)) < 1e-5, "derivative");
                    logger.assert_always(std::abs(test->evalDeriv2((point3D)) - 5.e7 * (x1 - x0)) < 1e-5, "derivative");
                    
                }
            }
        }
    }
    
    return 0;
}

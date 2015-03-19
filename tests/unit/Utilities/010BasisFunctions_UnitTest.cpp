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
//see also Base/???BasisFunction_UnitTest.cpp - running the same series of checks on different basisfunctions
#include "Utilities/BasisFunctions1DH1ConformingLine.h"
#include "Utilities/BasisFunctions2DH1ConformingSquare.h"
#include "Utilities/BasisFunctions2DH1ConformingTriangle.h"
#include "Utilities/BasisFunctions3DH1ConformingCube.h"
#include "Utilities/BasisFunctions3DH1ConformingTetrahedron.h"
#include "Utilities/BasisFunctions3DH1ConformingPrism.h"
#include "Logger.h"

#include "Base/BasisFunctionSet.h"
#include "Base/L2Norm.h"
#include "Geometry/PointReference.h"
#include "LinearAlgebra/NumericalVector.h"
#include "Geometry/ReferenceGeometry.h"
#include "Geometry/ReferenceLine.h"
#include <cmath>
using Base::L2Norm;

int main()
{
    
    // 1D
    
    Base::BasisFunctionSet *all1DbasisFunctions = Utilities::createDGBasisFunctionSet1DH1Line(5);
    Geometry::PointReference point1D(1);
    LinearAlgebra::NumericalVector ret(1);
    for (std::size_t i = 0; i < all1DbasisFunctions->size(); ++i)
    {
        const Base::BaseBasisFunction* test = (*all1DbasisFunctions)[i];
        for (point1D[0] = -1.5; point1D[0] < 1.51; point1D[0] += 0.1)
        {
            test->eval(point1D, ret);
            logger.assert_always(((test->eval(point1D) - ret[0]) < 1e-12), "eval");
            
            point1D[0] += -1.e-8;
            double x0 = test->eval(point1D);
            point1D[0] += 2.e-8;
            double x1 = test->eval(point1D);
            
            point1D[0] += -1e-8;
            ret = test->evalDeriv(point1D);
            logger.assert_always((std::abs(ret[0] - 5.e7 * (x1 - x0)) < 1e-5), "derivative");
            logger.assert_always((std::abs(test->evalDeriv0(point1D) - 5.e7 * (x1 - x0)) < 1e-5), "derivative");
        }
    }
    
    delete all1DbasisFunctions;
    
    // 2D
    
    Base::BasisFunctionSet *all2DbasisFunctions = Utilities::createDGBasisFunctionSet2DH1Square(5);
    Geometry::PointReference point2D(2);
    for (std::size_t i = 0; i < all2DbasisFunctions->size(); ++i)
    {
        const Base::BaseBasisFunction* test = (*all2DbasisFunctions)[i];
        for (point2D[0] = -1.5; point2D[0] < 1.51; point2D[0] += 0.1)
        {
            for (point2D[1] = -1.5; point2D[1] < 1.51; point2D[1] += 0.1)
            {
                test->eval(point2D, ret);
                logger.assert_always(((test->eval(point2D) - ret[0]) < 1e-12), "eval");
                
                point2D[0] += -1.e-8;
                double x0 = test->eval(point2D);
                point2D[0] += 2.e-8;
                double x1 = test->eval(point2D);
                
                point2D[0] += -1e-8;
                ret.resize(2);
                ret = test->evalDeriv(point2D);
                logger.assert_always((std::abs(ret[0] - 5.e7 * (x1 - x0)) < 1e-5), "derivative");
                logger.assert_always((std::abs(test->evalDeriv0(point2D) - 5.e7 * (x1 - x0)) < 1e-5), "derivative");
                
                point2D[1] += -1.e-8;
                x0 = test->eval(point2D);
                point2D[1] += 2.e-8;
                x1 = test->eval(point2D);
                
                point2D[1] += -1e-8;
                logger.assert_always((std::abs(ret[1] - 5.e7 * (x1 - x0)) < 1e-5), "derivative");
                logger.assert_always((std::abs(test->evalDeriv1(point2D) - 5.e7 * (x1 - x0)) < 1e-5), "derivative");
                
                ret.resize(1);
            }
        }
    }
    
    delete all2DbasisFunctions;
    
    all2DbasisFunctions = Utilities::createDGBasisFunctionSet2DH1Triangle(5);
    for (std::size_t i = 0; i < all2DbasisFunctions->size(); ++i)
    {
        const Base::BaseBasisFunction* test = (*all2DbasisFunctions)[i];
        for (point2D[0] = -1.5; point2D[0] < 1.51; point2D[0] += 0.1)
        {
            for (point2D[1] = -1.5; point2D[1] < 1.51; point2D[1] += 0.1)
            {
                test->eval(point2D, ret);
                logger.assert_always(((test->eval(point2D) - ret[0]) < 1e-12), "eval");
                
                point2D[0] += -1.e-8;
                double x0 = test->eval(point2D);
                point2D[0] += 2.e-8;
                double x1 = test->eval(point2D);
                
                point2D[0] += -1e-8;
                ret.resize(2);
                ret = test->evalDeriv(point2D); //exact to within absolute OR relative tolerace of 1e-5
                logger.assert_always((std::abs(ret[0] - 5.e7 * (x1 - x0)) < 1e-5 || (L2Norm(ret) > 1 && std::abs(ret[0] - 5.e7 * (x1 - x0)) < 1e-5 * L2Norm(ret))), "derivative");
                logger.assert_always((std::abs(test->evalDeriv0(point2D) - 5.e7 * (x1 - x0)) < 1e-5 || (L2Norm(ret) > 1 && std::abs(ret[0] - 5.e7 * (x1 - x0)) < 1e-5 * L2Norm(ret))), "derivative");
                
                point2D[1] += -1.e-8;
                x0 = test->eval(point2D);
                point2D[1] += 2.e-8;
                x1 = test->eval(point2D);
                
                point2D[1] += -1e-8;
                logger.assert_always((std::abs(ret[1] - 5.e7 * (x1 - x0)) < 1e-5 || (L2Norm(ret) > 1 && std::abs(ret[1] - 5.e7 * (x1 - x0)) < 1e-5 * L2Norm(ret))), "derivative");
                logger.assert_always((std::abs(test->evalDeriv1(point2D) - 5.e7 * (x1 - x0)) < 1e-5 || (L2Norm(ret) > 1 && std::abs(ret[1] - 5.e7 * (x1 - x0)) < 1e-5 * L2Norm(ret))), "derivative");
                
                ret.resize(1);
            }
        }
    }
    
    delete all2DbasisFunctions;
    
    //3D
    
    Base::BasisFunctionSet *all3DbasisFunctions = Utilities::createDGBasisFunctionSet3DH1Cube(5);
    Geometry::PointReference point3D(3);
    for (std::size_t i = 0; i < all3DbasisFunctions->size(); ++i)
    {
        const Base::BaseBasisFunction* test = (*all3DbasisFunctions)[i];
        for (point3D[0] = -1.5; point3D[0] < 1.51; point3D[0] += 0.25)
        {
            for (point3D[1] = -1.5; point3D[1] < 1.51; point3D[1] += 0.25)
            {
                for (point3D[2] = -1.5; point3D[2] < 1.51; point3D[2] += 0.25)
                {
                    test->eval(point3D, ret);
                    logger.assert_always(((test->eval(point3D) - ret[0]) < 1e-12), "eval");
                    
                    point3D[0] += -1.e-8;
                    double x0 = test->eval(point3D);
                    point3D[0] += 2.e-8;
                    double x1 = test->eval(point3D);
                    
                    point3D[0] += -1e-8;
                    ret.resize(3);
                    ret = test->evalDeriv(point3D);
                    logger.assert_always((std::abs(ret[0] - 5.e7 * (x1 - x0)) < 1e-5 || (L2Norm(ret) > 1 && std::abs(ret[0] - 5.e7 * (x1 - x0)) < 1e-5 * L2Norm(ret))), "derivative");
                    logger.assert_always((std::abs(test->evalDeriv0(point3D) - 5.e7 * (x1 - x0)) < 1e-5 || (L2Norm(ret) > 1 && std::abs(ret[0] - 5.e7 * (x1 - x0)) < 1e-5 * L2Norm(ret))), "derivative");
                    
                    point3D[1] += -1.e-8;
                    x0 = test->eval(point3D);
                    point3D[1] += 2.e-8;
                    x1 = test->eval(point3D);
                    
                    point3D[1] += -1e-8;
                    logger.assert_always((std::abs(ret[1] - 5.e7 * (x1 - x0)) < 1e-5 || (L2Norm(ret) > 1 && std::abs(ret[1] - 5.e7 * (x1 - x0)) < 1e-5 * L2Norm(ret))), "derivative");
                    logger.assert_always((std::abs(test->evalDeriv1(point3D) - 5.e7 * (x1 - x0)) < 1e-5 || (L2Norm(ret) > 1 && std::abs(ret[1] - 5.e7 * (x1 - x0)) < 1e-5 * L2Norm(ret))), "derivative");
                    
                    point3D[2] += -1.e-8;
                    x0 = test->eval(point3D);
                    point3D[2] += 2.e-8;
                    x1 = test->eval(point3D);
                    
                    point3D[2] += -1e-8;
                    logger.assert_always((std::abs(ret[2] - 5.e7 * (x1 - x0)) < 1e-5 || (L2Norm(ret) > 1 && std::abs(ret[2] - 5.e7 * (x1 - x0)) < 1e-5 * L2Norm(ret))), "derivative");
                    logger.assert_always((std::abs(test->evalDeriv2(point3D) - 5.e7 * (x1 - x0)) < 1e-5 || (L2Norm(ret) > 1 && std::abs(ret[2] - 5.e7 * (x1 - x0)) < 1e-5 * L2Norm(ret))), "derivative");
                    
                    ret.resize(1);
                }
            }
        }
    }
    
    delete all3DbasisFunctions;
    
    all3DbasisFunctions = Utilities::createDGBasisFunctionSet3DH1Tetrahedron(5);
    for (std::size_t i = 0; i < all3DbasisFunctions->size(); ++i)
    {
        const Base::BaseBasisFunction* test = (*all3DbasisFunctions)[i];
        for (point3D[0] = -1.5; point3D[0] < 1.51; point3D[0] += 0.25)
        {
            for (point3D[1] = -1.5; point3D[1] < 1.51; point3D[1] += 0.25)
            {
                for (point3D[2] = -1.5; point3D[2] < 1.51; point3D[2] += 0.25)
                {
                    test->eval(point3D, ret);
                    logger.assert_always(((test->eval(point3D) - ret[0]) < 1e-12), "eval");
                    
                    point3D[0] += -1.e-8;
                    double x0 = test->eval(point3D);
                    point3D[0] += 2.e-8;
                    double x1 = test->eval(point3D);
                    
                    point3D[0] += -1e-8;
                    ret.resize(3);
                    ret = test->evalDeriv(point3D);
                    logger.assert_always((std::abs(ret[0] - 5.e7 * (x1 - x0)) < 1e-5 || (L2Norm(ret) > 1 && std::abs(ret[0] - 5.e7 * (x1 - x0)) < 1e-5 * L2Norm(ret))), "derivative");
                    logger.assert_always((std::abs(test->evalDeriv0(point3D) - 5.e7 * (x1 - x0)) < 1e-5 || (L2Norm(ret) > 1 && std::abs(ret[0] - 5.e7 * (x1 - x0)) < 1e-5 * L2Norm(ret))), "derivative");
                    
                    point3D[1] += -1.e-8;
                    x0 = test->eval(point3D);
                    point3D[1] += 2.e-8;
                    x1 = test->eval(point3D);
                    
                    point3D[1] += -1e-8;
                    logger.assert_always((std::abs(ret[1] - 5.e7 * (x1 - x0)) < 1e-5 || (L2Norm(ret) > 1 && std::abs(ret[1] - 5.e7 * (x1 - x0)) < 1e-5 * L2Norm(ret))), "derivative");
                    logger.assert_always((std::abs(test->evalDeriv1(point3D) - 5.e7 * (x1 - x0)) < 1e-5 || (L2Norm(ret) > 1 && std::abs(ret[1] - 5.e7 * (x1 - x0)) < 1e-5 * L2Norm(ret))), "derivative");
                    
                    point3D[2] += -1.e-8;
                    x0 = test->eval(point3D);
                    point3D[2] += 2.e-8;
                    x1 = test->eval(point3D);
                    
                    point3D[2] += -1e-8;
                    logger.assert_always((std::abs(ret[2] - 5.e7 * (x1 - x0)) < 1e-5 || (L2Norm(ret) > 1 && std::abs(ret[2] - 5.e7 * (x1 - x0)) < 1e-5 * L2Norm(ret))), "derivative");
                    logger.assert_always((std::abs(test->evalDeriv2(point3D) - 5.e7 * (x1 - x0)) < 1e-5 || (L2Norm(ret) > 1 && std::abs(ret[2] - 5.e7 * (x1 - x0)) < 1e-5 * L2Norm(ret))), "derivative");
                    
                    ret.resize(1);
                }
            }
        }
    }
    
    delete all3DbasisFunctions;
    
    all3DbasisFunctions = Utilities::createDGBasisFunctionSet3DH1ConformingPrism(5);
    for (std::size_t i = 0; i < all3DbasisFunctions->size(); ++i)
    {
        const Base::BaseBasisFunction* test = (*all3DbasisFunctions)[i];
        for (point3D[0] = -1.5; point3D[0] < 1.51; point3D[0] += 0.25)
        {
            for (point3D[1] = -1.5; point3D[1] < 1.51; point3D[1] += 0.25)
            {
                for (point3D[2] = -1.5; point3D[2] < 1.51; point3D[2] += 0.25)
                {
                    test->eval(point3D, ret);
                    logger.assert_always(((test->eval(point3D) - ret[0]) < 1e-12), "eval");
                    
                    point3D[0] += -1.e-8;
                    double x0 = test->eval(point3D);
                    point3D[0] += 2.e-8;
                    double x1 = test->eval(point3D);
                    
                    point3D[0] += -1e-8;
                    ret.resize(3);
                    ret = test->evalDeriv(point3D);
                    logger.assert_always((std::abs(ret[0] - 5.e7 * (x1 - x0)) < 1e-5 || (L2Norm(ret) > 1 && std::abs(ret[0] - 5.e7 * (x1 - x0)) < 1e-5 * L2Norm(ret))), "derivative");
                    logger.assert_always((std::abs(test->evalDeriv0(point3D) - 5.e7 * (x1 - x0)) < 1e-5 || (L2Norm(ret) > 1 && std::abs(ret[0] - 5.e7 * (x1 - x0)) < 1e-5 * L2Norm(ret))), "derivative");
                    
                    point3D[1] += -1.e-8;
                    x0 = test->eval(point3D);
                    point3D[1] += 2.e-8;
                    x1 = test->eval(point3D);
                    
                    point3D[1] += -1e-8;
                    logger.assert_always((std::abs(ret[1] - 5.e7 * (x1 - x0)) < 1e-5 || (L2Norm(ret) > 1 && std::abs(ret[1] - 5.e7 * (x1 - x0)) < 1e-5 * L2Norm(ret))), "derivative");
                    logger.assert_always((std::abs(test->evalDeriv1(point3D) - 5.e7 * (x1 - x0)) < 1e-5 || (L2Norm(ret) > 1 && std::abs(ret[1] - 5.e7 * (x1 - x0)) < 1e-5 * L2Norm(ret))), "derivative");
                    
                    point3D[2] += -1.e-8;
                    x0 = test->eval(point3D);
                    point3D[2] += 2.e-8;
                    x1 = test->eval(point3D);
                    
                    point3D[2] += -1e-8;
                    logger.assert_always((std::abs(ret[2] - 5.e7 * (x1 - x0)) < 1e-5 || (L2Norm(ret) > 1 && std::abs(ret[2] - 5.e7 * (x1 - x0)) < 1e-5 * L2Norm(ret))), "derivative");
                    logger.assert_always((std::abs(test->evalDeriv2(point3D) - 5.e7 * (x1 - x0)) < 1e-5 || (L2Norm(ret) > 1 && std::abs(ret[2] - 5.e7 * (x1 - x0)) < 1e-5 * L2Norm(ret))), "derivative");
                    
                    ret.resize(1);
                }
            }
        }
    }
    
    delete all3DbasisFunctions;
    
    return 0;
}


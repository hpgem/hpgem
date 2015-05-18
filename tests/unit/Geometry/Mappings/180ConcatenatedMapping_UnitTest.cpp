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
#include "Geometry/Mappings/ConcatenatedMapping.h"
#include "Geometry/Mappings/MappingToRefPointToPoint.h"
#include "Geometry/Mappings/MappingToRefPointToLine.h"
#include "Geometry/Mappings/MappingToRefLineToLine.h"
#include "Geometry/Mappings/MappingToRefLineToSquare.h"
#include "Geometry/Mappings/MappingToRefSquareToSquare.h"
#include "Geometry/Mappings/MappingToRefSquareToCube.h"
#include "Geometry/Mappings/MappingToRefCubeToCube.h"
#include "Geometry/Mappings/MappingToRefCubeToHypercube.h"
#include "Logger.h"

#include "Geometry/ReferencePoint.h"
#include "Geometry/ReferenceLine.h"
#include "Geometry/ReferenceSquare.h"
#include "Geometry/ReferenceCube.h"
#include "Geometry/ReferenceHypercube.h"
#include "Geometry/PointReference.h"
#include "Geometry/Jacobian.h"
#include "LinearAlgebra/MiddleSizeVector.h"
#include <cmath>

int main()
{
    
    Geometry::Point point0D(0), point1D(1), point2D(2), point3D(3), point4D(4), compare0D(0), compare1D(1), compare2D(2), compare3D(3), compare4D(4), orig0D(0), orig1D(1), orig2D(2), orig3D(3), orig4D(4);
    
    //0->0->0
    
    Geometry::ReferenceGeometry* source = &Geometry::ReferencePoint::Instance();
    Geometry::ReferenceGeometry* target = &Geometry::ReferencePoint::Instance();
    
    Geometry::Jacobian jac(0, 0);
    
    std::vector<std::size_t> nodesAfterTransformation(1);
    
    Geometry::ConcatenatedMapping* test = new Geometry::ConcatenatedMapping(Geometry::MappingToRefPointToPoint::Instance(), Geometry::MappingToRefPointToPoint::Instance());
    nodesAfterTransformation[0] = 0;
    
    point0D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig0D));
    logger.assert_always((source->isInternalPoint(*Geometry::PointReferenceFactory::instance()->makePoint(orig0D)) == target->isInternalPoint(*Geometry::PointReferenceFactory::instance()->makePoint(point0D))), "transform");
    
    jac = test->calcJacobian(*Geometry::PointReferenceFactory::instance()->makePoint(orig0D));
    
    for (std::size_t i = 0; i < source->getNumberOfNodes(); ++i)
    {
        orig0D = source->getNode(i);
        compare0D = target->getNode(nodesAfterTransformation[i]);
        point0D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig0D));
    }
    
    logger.assert_always((test->getTargetDimension() == 0), "getTargetDimension");
    
    //0->0->1
    
    target = &Geometry::ReferenceLine::Instance();
    jac.resize(0, 1);
    delete test;
    test = new Geometry::ConcatenatedMapping(Geometry::MappingToRefPointToPoint::Instance(), Geometry::MappingToRefPointToLine0::Instance());
    
    point1D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig0D));
    logger.assert_always((source->isInternalPoint(*Geometry::PointReferenceFactory::instance()->makePoint(orig0D)) == target->isInternalPoint(*Geometry::PointReferenceFactory::instance()->makePoint(point1D))), "transform");
    
    jac = test->calcJacobian(*Geometry::PointReferenceFactory::instance()->makePoint(orig0D));
    
    for (std::size_t i = 0; i < source->getNumberOfNodes(); ++i)
    {
        orig0D = source->getNode(i);
        compare1D = target->getNode(nodesAfterTransformation[i]);
        point1D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig0D));
    }
    
    logger.assert_always((test->getTargetDimension() == 1), "getTargetDimension");
    
    //0->1->1
    
    delete test;
    test = new Geometry::ConcatenatedMapping(Geometry::MappingToRefPointToLine1::Instance(), Geometry::MappingToRefLineToLine1::Instance());
    
    point1D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig0D));
    logger.assert_always((source->isInternalPoint(*Geometry::PointReferenceFactory::instance()->makePoint(orig0D)) == target->isInternalPoint(*Geometry::PointReferenceFactory::instance()->makePoint(point1D))), "transform");
    
    jac = test->calcJacobian(*Geometry::PointReferenceFactory::instance()->makePoint(orig0D));
    
    for (std::size_t i = 0; i < source->getNumberOfNodes(); ++i)
    {
        orig0D = source->getNode(i);
        compare1D = target->getNode(nodesAfterTransformation[i]);
        point1D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig0D));
    }
    
    logger.assert_always((test->getTargetDimension() == 1), "getTargetDimension");
    
    //0->1->2
    
    target = &Geometry::ReferenceSquare::Instance();
    jac.resize(0, 2);
    delete test;
    test = new Geometry::ConcatenatedMapping(Geometry::MappingToRefPointToLine1::Instance(), Geometry::MappingToRefLineToSquare3::Instance());
    nodesAfterTransformation[0] = 3;
    
    point2D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig0D));
    logger.assert_always((source->isInternalPoint(*Geometry::PointReferenceFactory::instance()->makePoint(orig0D)) == target->isInternalPoint(*Geometry::PointReferenceFactory::instance()->makePoint(point2D))), "transform");
    
    jac = test->calcJacobian(*Geometry::PointReferenceFactory::instance()->makePoint(orig0D));
    
    for (std::size_t i = 0; i < source->getNumberOfNodes(); ++i)
    {
        orig0D = source->getNode(i);
        compare2D = target->getNode(nodesAfterTransformation[i]);
        point2D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig0D));
    }
    
    logger.assert_always((test->getTargetDimension() == 2), "getTargetDimension");
    
    //1->1->1
    
    source = &Geometry::ReferenceLine::Instance();
    target = &Geometry::ReferenceLine::Instance();
    jac.resize(1, 1);
    nodesAfterTransformation.resize(2);
    delete test;
    test = new Geometry::ConcatenatedMapping(Geometry::MappingToRefLineToLine1::Instance(), Geometry::MappingToRefLineToLine0::Instance());
    nodesAfterTransformation[0] = 1;
    nodesAfterTransformation[1] = 0;
    
    for (orig1D[0] = -2.8189; orig1D[0] < 3.141; orig1D[0] += 0.1)
    {
        point1D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig1D));
        logger.assert_always((source->isInternalPoint(*Geometry::PointReferenceFactory::instance()->makePoint(orig1D)) == target->isInternalPoint(*Geometry::PointReferenceFactory::instance()->makePoint(point1D))), "transform");
        
        orig1D[0] += -1.e-8;
        compare1D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig1D));
        orig1D[0] += 2.e-8;
        point1D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig1D));
        
        orig1D[0] += -1e-8;
        jac = test->calcJacobian(*Geometry::PointReferenceFactory::instance()->makePoint(orig1D));
        logger.assert_always((std::abs(jac[0] - 5.e7 * (point1D[0] - compare1D[0])) < 1e-5), "jacobian"); //estimate is a bit rough, but should work for most mappings
    }
    
    for (std::size_t i = 0; i < source->getNumberOfNodes(); ++i)
    {
        orig1D = source->getNode(i);
        compare1D = target->getNode(nodesAfterTransformation[i]);
        point1D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig1D));
        logger.assert_always((std::abs(point1D[0] - compare1D[0]) < 1e-12), "transform");
    }
    
    logger.assert_always((test->getTargetDimension() == 1), "getTargetDimension");
    
    //1->1->2
    
    target = &Geometry::ReferenceSquare::Instance();
    jac.resize(2, 1);
    delete test;
    test = new Geometry::ConcatenatedMapping(Geometry::MappingToRefLineToLine1::Instance(), Geometry::MappingToRefLineToSquare3::Instance());
    nodesAfterTransformation[0] = 3;
    nodesAfterTransformation[1] = 2;
    
    for (orig1D[0] = -2.8189; orig1D[0] < 3.141; orig1D[0] += 0.1)
    {
        point2D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig1D));
        logger.assert_always((source->isInternalPoint(*Geometry::PointReferenceFactory::instance()->makePoint(orig1D)) == target->isInternalPoint(*Geometry::PointReferenceFactory::instance()->makePoint(point2D))), "transform");
        
        orig1D[0] += -1.e-8;
        compare2D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig1D));
        orig1D[0] += 2.e-8;
        point2D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig1D));
        
        orig1D[0] += -1e-8;
        jac = test->calcJacobian(*Geometry::PointReferenceFactory::instance()->makePoint(orig1D));
        logger.assert_always((std::abs(jac[0] - 5.e7 * (point2D[0] - compare2D[0])) < 1e-5), "jacobian"); //estimate is a bit rough, but should work for most mappings
        logger.assert_always((std::abs(jac[1] - 5.e7 * (point2D[1] - compare2D[1])) < 1e-5), "jacobian"); //implementations are very strongly recommended to be more accurate
    }
    
    for (std::size_t i = 0; i < source->getNumberOfNodes(); ++i)
    {
        orig1D = source->getNode(i);
        compare2D = target->getNode(nodesAfterTransformation[i]);
        point2D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig1D));
        logger.assert_always((std::abs(point2D[0] - compare2D[0]) < 1e-12), "transform");
        logger.assert_always((std::abs(point2D[1] - compare2D[1]) < 1e-12), "transform");
    }
    
    logger.assert_always((test->getTargetDimension() == 2), "getTargetDimension");
    
    //1->2->2
    
    delete test;
    test = new Geometry::ConcatenatedMapping(Geometry::MappingToRefLineToSquare3::Instance(), Geometry::MappingToRefSquareToSquare7::Instance());
    nodesAfterTransformation[0] = 1;
    nodesAfterTransformation[1] = 3;
    
    for (orig1D[0] = -2.8189; orig1D[0] < 3.141; orig1D[0] += 0.1)
    {
        point2D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig1D));
        logger.assert_always((source->isInternalPoint(*Geometry::PointReferenceFactory::instance()->makePoint(orig1D)) == target->isInternalPoint(*Geometry::PointReferenceFactory::instance()->makePoint(point2D))), "transform");
        
        orig1D[0] += -1.e-8;
        compare2D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig1D));
        orig1D[0] += 2.e-8;
        point2D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig1D));
        
        orig1D[0] += -1e-8;
        jac = test->calcJacobian(*Geometry::PointReferenceFactory::instance()->makePoint(orig1D));
        //std::cout<<jac<<std::endl<<5.e7*(point2D-compare2D)<<std::endl;
        logger.assert_always((std::abs(jac[0] - 5.e7 * (point2D[0] - compare2D[0])) < 1e-5), "jacobian"); //estimate is a bit rough, but should work for most mappings
        logger.assert_always((std::abs(jac[1] - 5.e7 * (point2D[1] - compare2D[1])) < 1e-5), "jacobian"); //implementations are very strongly recommended to be more accurate
    }
    
    for (std::size_t i = 0; i < source->getNumberOfNodes(); ++i)
    {
        orig1D = source->getNode(i);
        compare2D = target->getNode(nodesAfterTransformation[i]);
        point2D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig1D));
        logger.assert_always((std::abs(point2D[0] - compare2D[0]) < 1e-12), "transform");
        logger.assert_always((std::abs(point2D[1] - compare2D[1]) < 1e-12), "transform");
    }
    
    logger.assert_always((test->getTargetDimension() == 2), "getTargetDimension");
    
    //1->2->3
    
    target = &Geometry::ReferenceCube::Instance();
    jac.resize(3, 1);
    delete test;
    test = new Geometry::ConcatenatedMapping(Geometry::MappingToRefLineToSquare1::Instance(), Geometry::MappingToRefSquareToCube2::Instance());
    nodesAfterTransformation[0] = 0;
    nodesAfterTransformation[1] = 4;
    
    for (orig1D[0] = -2.8189; orig1D[0] < 3.141; orig1D[0] += 0.1)
    {
        point3D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig1D));
        logger.assert_always((source->isInternalPoint(*Geometry::PointReferenceFactory::instance()->makePoint(orig1D)) == target->isInternalPoint(*Geometry::PointReferenceFactory::instance()->makePoint(point3D))), "transform");
        
        orig1D[0] += -1.e-8;
        compare3D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig1D));
        orig1D[0] += 2.e-8;
        point3D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig1D));
        
        orig1D[0] += -1e-8;
        jac = test->calcJacobian(*Geometry::PointReferenceFactory::instance()->makePoint(orig1D));
        logger.assert_always((std::abs(jac[0] - 5.e7 * (point3D[0] - compare3D[0])) < 1e-5), "jacobian"); //estimate is a bit rough, but should work for most mappings
        logger.assert_always((std::abs(jac[1] - 5.e7 * (point3D[1] - compare3D[1])) < 1e-5), "jacobian"); //implementations are very strongly recommended to be more accurate
    }
    
    for (std::size_t i = 0; i < source->getNumberOfNodes(); ++i)
    {
        orig1D = source->getNode(i);
        compare3D = target->getNode(nodesAfterTransformation[i]);
        point3D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig1D));
        logger.assert_always((std::abs(point3D[0] - compare3D[0]) < 1e-12), "transform");
        logger.assert_always((std::abs(point3D[1] - compare3D[1]) < 1e-12), "transform");
    }
    
    logger.assert_always((test->getTargetDimension() == 3), "getTargetDimension");
    
    //2->2->2
    
    source = &Geometry::ReferenceSquare::Instance();
    target = &Geometry::ReferenceSquare::Instance();
    nodesAfterTransformation.resize(4);
    jac.resize(2, 2);
    delete test;
    test = new Geometry::ConcatenatedMapping(Geometry::MappingToRefSquareToSquare3::Instance(), Geometry::MappingToRefSquareToSquare2::Instance());
    nodesAfterTransformation[0] = 1;
    nodesAfterTransformation[1] = 3;
    nodesAfterTransformation[2] = 0;
    nodesAfterTransformation[3] = 2;
    
    for (orig2D[0] = -1.51; orig2D[0] < 1.51; orig2D[0] += 0.2)
    {
        for (orig2D[1] = -1.51; orig2D[1] < 1.51; orig2D[1] += 0.2)
        {
            point2D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D));
            logger.assert_always((source->isInternalPoint(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D)) == target->isInternalPoint(*Geometry::PointReferenceFactory::instance()->makePoint(point2D))), "transform");
            
            orig2D[0] += -1.e-8;
            compare2D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D));
            orig2D[0] += 2.e-8;
            point2D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D));
            
            orig2D[0] += -1e-8;
            jac = test->calcJacobian(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D));
            logger.assert_always((std::abs(jac[0] - 5.e7 * (point2D[0] - compare2D[0])) < 1e-5), "jacobian"); //estimate is a bit rough, but should work for most mappings
            logger.assert_always((std::abs(jac[1] - 5.e7 * (point2D[1] - compare2D[1])) < 1e-5), "jacobian"); //implementations are very strongly recommended to be more accurate
                    
            orig2D[1] += -1.e-8;
            compare2D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D));
            orig2D[1] += 2.e-8;
            point2D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D));
            
            orig2D[1] += -1e-8;
            jac = test->calcJacobian(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D));
            logger.assert_always((std::abs(jac[2] - 5.e7 * (point2D[0] - compare2D[0])) < 1e-5), "jacobian");
            logger.assert_always((std::abs(jac[3] - 5.e7 * (point2D[1] - compare2D[1])) < 1e-5), "jacobian");
        }
    }
    
    for (std::size_t i = 0; i < source->getNumberOfNodes(); ++i)
    {
        orig2D = source->getNode(i);
        compare2D = target->getNode(nodesAfterTransformation[i]);
        point2D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D));
        logger.assert_always((std::abs(point2D[0] - compare2D[0]) < 1e-12), "transform");
        logger.assert_always((std::abs(point2D[1] - compare2D[1]) < 1e-12), "transform");
    }
    
    logger.assert_always((test->getTargetDimension() == 2), "getTargetDimension");
    
    //2->2->3
    
    target = &Geometry::ReferenceCube::Instance();
    jac.resize(3, 2);
    delete test;
    test = new Geometry::ConcatenatedMapping(Geometry::MappingToRefSquareToSquare3::Instance(), Geometry::MappingToRefSquareToCube2::Instance());
    nodesAfterTransformation[0] = 4;
    nodesAfterTransformation[1] = 0;
    nodesAfterTransformation[2] = 6;
    nodesAfterTransformation[3] = 2;
    
    for (orig2D[0] = -1.51; orig2D[0] < 1.51; orig2D[0] += 0.2)
    {
        for (orig2D[1] = -1.51; orig2D[1] < 1.51; orig2D[1] += 0.2)
        {
            point3D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D));
            logger.assert_always((source->isInternalPoint(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D)) == target->isInternalPoint(*Geometry::PointReferenceFactory::instance()->makePoint(point3D))), "transform");
            
            orig2D[0] += -1.e-8;
            compare3D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D));
            orig2D[0] += 2.e-8;
            point3D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D));
            
            orig2D[0] += -1e-8;
            jac = test->calcJacobian(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D));
            logger.assert_always((std::abs(jac[0] - 5.e7 * (point3D[0] - compare3D[0])) < 1e-5), "jacobian"); //estimate is a bit rough, but should work for most mappings
            logger.assert_always((std::abs(jac[1] - 5.e7 * (point3D[1] - compare3D[1])) < 1e-5), "jacobian"); //implementations are very strongly recommended to be more accurate
            logger.assert_always((std::abs(jac[2] - 5.e7 * (point3D[2] - compare3D[2])) < 1e-5), "jacobian");
            
            orig2D[1] += -1.e-8;
            compare3D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D));
            orig2D[1] += 2.e-8;
            point3D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D));
            
            orig2D[1] += -1e-8;
            jac = test->calcJacobian(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D));
            logger.assert_always((std::abs(jac[3] - 5.e7 * (point3D[0] - compare3D[0])) < 1e-5), "jacobian");
            logger.assert_always((std::abs(jac[4] - 5.e7 * (point3D[1] - compare3D[1])) < 1e-5), "jacobian");
            logger.assert_always((std::abs(jac[5] - 5.e7 * (point3D[2] - compare3D[2])) < 1e-5), "jacobian");
        }
    }
    
    for (std::size_t i = 0; i < source->getNumberOfNodes(); ++i)
    {
        orig2D = source->getNode(i);
        compare3D = target->getNode(nodesAfterTransformation[i]);
        point3D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D));
        logger.assert_always((std::abs(point3D[0] - compare3D[0]) < 1e-12), "transform");
        logger.assert_always((std::abs(point3D[1] - compare3D[1]) < 1e-12), "transform");
        logger.assert_always((std::abs(point3D[2] - compare3D[2]) < 1e-12), "transform");
    }
    
    logger.assert_always((test->getTargetDimension() == 3), "getTargetDimension");
    
    //2->3->3
    
    delete test;
    test = new Geometry::ConcatenatedMapping(Geometry::MappingToRefSquareToCube3::Instance(), Geometry::MappingToRefCubeToCube2::Instance());
    nodesAfterTransformation[0] = 4;
    nodesAfterTransformation[1] = 0;
    nodesAfterTransformation[2] = 6;
    nodesAfterTransformation[3] = 2;
    
    for (orig2D[0] = -1.51; orig2D[0] < 1.51; orig2D[0] += 0.2)
    {
        for (orig2D[1] = -1.51; orig2D[1] < 1.51; orig2D[1] += 0.2)
        {
            point3D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D));
            logger.assert_always((source->isInternalPoint(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D)) == target->isInternalPoint(*Geometry::PointReferenceFactory::instance()->makePoint(point3D))), "transform");
            
            orig2D[0] += -1.e-8;
            compare3D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D));
            orig2D[0] += 2.e-8;
            point3D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D));
            
            orig2D[0] += -1e-8;
            jac = test->calcJacobian(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D));
            logger.assert_always((std::abs(jac[0] - 5.e7 * (point3D[0] - compare3D[0])) < 1e-5), "jacobian"); //estimate is a bit rough, but should work for most mappings
            logger.assert_always((std::abs(jac[1] - 5.e7 * (point3D[1] - compare3D[1])) < 1e-5), "jacobian"); //implementations are very strongly recommended to be more accurate
            logger.assert_always((std::abs(jac[2] - 5.e7 * (point3D[2] - compare3D[2])) < 1e-5), "jacobian");
            
            orig2D[1] += -1.e-8;
            compare3D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D));
            orig2D[1] += 2.e-8;
            point3D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D));
            
            orig2D[1] += -1e-8;
            jac = test->calcJacobian(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D));
            logger.assert_always((std::abs(jac[3] - 5.e7 * (point3D[0] - compare3D[0])) < 1e-5), "jacobian");
            logger.assert_always((std::abs(jac[4] - 5.e7 * (point3D[1] - compare3D[1])) < 1e-5), "jacobian");
            logger.assert_always((std::abs(jac[5] - 5.e7 * (point3D[2] - compare3D[2])) < 1e-5), "jacobian");
        }
    }
    
    /*for(std::size_t i=0;i<source->getNumberOfNodes();++i){///\TODO figure out how the cube->cube maps actually map their nodes
     source->getNode(i,orig2D);
     target->getNode(nodesAfterTransformation[i],compare3D);
     test->transform(orig2D,point3D);
     logger.assert_always(("transform",std::abs(point3D[0]-compare3D[0])<1e-12));
     logger.assert_always(("transform",std::abs(point3D[1]-compare3D[1])<1e-12));
     logger.assert_always(("transform",std::abs(point3D[2]-compare3D[2])<1e-12));
     }*/

    logger.assert_always((test->getTargetDimension() == 3), "getTargetDimension");
    
    //2->3->4
    jac.resize(4, 2);
    
    target = &Geometry::ReferenceHypercube::Instance();
    delete test;
    test = new Geometry::ConcatenatedMapping(Geometry::MappingToRefSquareToCube3::Instance(), Geometry::MappingToRefCubeToHypercube4::Instance());
    nodesAfterTransformation[0] = 3;
    nodesAfterTransformation[1] = 7;
    nodesAfterTransformation[2] = 11;
    nodesAfterTransformation[3] = 15;
    
    for (orig2D[0] = -1.51; orig2D[0] < 1.51; orig2D[0] += 0.2)
    {
        for (orig2D[1] = -1.51; orig2D[1] < 1.51; orig2D[1] += 0.2)
        {
            point4D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D));
            logger.assert_always((source->isInternalPoint(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D)) == target->isInternalPoint(*Geometry::PointReferenceFactory::instance()->makePoint(point4D))), "transform");
            
            orig2D[0] += -1.e-8;
            compare4D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D));
            orig2D[0] += 2.e-8;
            point4D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D));
            
            orig2D[0] += -1e-8;
            jac = test->calcJacobian(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D));
            logger.assert_always((std::abs(jac[0] - 5.e7 * (point4D[0] - compare4D[0])) < 1e-5), "jacobian"); //estimate is a bit rough, but should work for most mappings
            logger.assert_always((std::abs(jac[1] - 5.e7 * (point4D[1] - compare4D[1])) < 1e-5), "jacobian"); //implementations are very strongly recommended to be more accurate
            logger.assert_always((std::abs(jac[2] - 5.e7 * (point4D[2] - compare4D[2])) < 1e-5), "jacobian");
            logger.assert_always((std::abs(jac[3] - 5.e7 * (point4D[3] - compare4D[3])) < 1e-5), "jacobian");
            
            orig2D[1] += -1.e-8;
            compare4D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D));
            orig2D[1] += 2.e-8;
            point4D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D));
            
            orig2D[1] += -1e-8;
            jac = test->calcJacobian(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D));
            logger.assert_always((std::abs(jac[4] - 5.e7 * (point4D[0] - compare4D[0])) < 1e-5), "jacobian");
            logger.assert_always((std::abs(jac[5] - 5.e7 * (point4D[1] - compare4D[1])) < 1e-5), "jacobian");
            logger.assert_always((std::abs(jac[6] - 5.e7 * (point4D[2] - compare4D[2])) < 1e-5), "jacobian");
            logger.assert_always((std::abs(jac[7] - 5.e7 * (point4D[3] - compare4D[3])) < 1e-5), "jacobian");
        }
    }
    
    for (std::size_t i = 0; i < source->getNumberOfNodes(); ++i)
    {
        orig2D = source->getNode(i);
        compare4D = target->getNode(nodesAfterTransformation[i]);
        point4D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig2D));
        logger.assert_always((std::abs(point4D[0] - compare4D[0]) < 1e-12), "transform");
        logger.assert_always((std::abs(point4D[1] - compare4D[1]) < 1e-12), "transform");
        logger.assert_always((std::abs(point4D[2] - compare4D[2]) < 1e-12), "transform");
        logger.assert_always((std::abs(point4D[3] - compare4D[3]) < 1e-12), "transform");
    }
    
    logger.assert_always((test->getTargetDimension() == 4), "getTargetDimension");
    
    //3->3->3
    
    source = &Geometry::ReferenceCube::Instance();
    target = &Geometry::ReferenceCube::Instance();
    nodesAfterTransformation.resize(8);
    jac.resize(3, 3);
    delete test;
    test = new Geometry::ConcatenatedMapping(Geometry::MappingToRefCubeToCube3::Instance(), Geometry::MappingToRefCubeToCube7::Instance());
    nodesAfterTransformation[0] = 4;
    nodesAfterTransformation[1] = 0;
    nodesAfterTransformation[2] = 6;
    nodesAfterTransformation[3] = 2;
    
    for (orig3D[0] = -1.51; orig3D[0] < 1.51; orig3D[0] += 0.6)
    {
        for (orig3D[1] = -1.51; orig3D[1] < 1.51; orig3D[1] += 0.7)
        {
            for (orig3D[2] = -1.51; orig3D[2] < 1.51; orig3D[2] += 0.8)
            {
                point3D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig3D));
                logger.assert_always((source->isInternalPoint(*Geometry::PointReferenceFactory::instance()->makePoint(orig3D)) == target->isInternalPoint(*Geometry::PointReferenceFactory::instance()->makePoint(point3D))), "transform");
                
                orig3D[0] += -1.e-8;
                compare3D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig3D));
                orig3D[0] += 2.e-8;
                point3D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig3D));
                
                orig3D[0] += -1e-8;
                jac = test->calcJacobian(*Geometry::PointReferenceFactory::instance()->makePoint(orig3D));
                logger.assert_always((std::abs(jac[0] - 5.e7 * (point3D[0] - compare3D[0])) < 1e-5), "jacobian"); //estimate is a bit rough, but should work for most mappings
                logger.assert_always((std::abs(jac[1] - 5.e7 * (point3D[1] - compare3D[1])) < 1e-5), "jacobian"); //implementations are very strongly recommended to be more accurate
                logger.assert_always((std::abs(jac[2] - 5.e7 * (point3D[2] - compare3D[2])) < 1e-5), "jacobian");
                
                orig3D[1] += -1.e-8;
                compare3D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig3D));
                orig3D[1] += 2.e-8;
                point3D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig3D));
                
                orig3D[1] += -1e-8;
                jac = test->calcJacobian(*Geometry::PointReferenceFactory::instance()->makePoint(orig3D));
                logger.assert_always((std::abs(jac[3] - 5.e7 * (point3D[0] - compare3D[0])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[4] - 5.e7 * (point3D[1] - compare3D[1])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[5] - 5.e7 * (point3D[2] - compare3D[2])) < 1e-5), "jacobian");
                
                orig3D[2] += -1.e-8;
                compare3D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig3D));
                orig3D[2] += 2.e-8;
                point3D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig3D));
                
                orig3D[2] += -1e-8;
                jac = test->calcJacobian(*Geometry::PointReferenceFactory::instance()->makePoint(orig3D));
                logger.assert_always((std::abs(jac[6] - 5.e7 * (point3D[0] - compare3D[0])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[7] - 5.e7 * (point3D[1] - compare3D[1])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[8] - 5.e7 * (point3D[2] - compare3D[2])) < 1e-5), "jacobian");
            }
        }
    }
    
    /*for(std::size_t i=0;i<source->getNumberOfNodes();++i){///\TODO figure out how the cube->cube maps actually map their nodes
     source->getNode(i,orig3D);
     target->getNode(nodesAfterTransformation[i],compare3D);
     test->transform(orig3D,point3D);
     logger.assert_always(("transform",std::abs(point3D[0]-compare3D[0])<1e-12));
     logger.assert_always(("transform",std::abs(point3D[1]-compare3D[1])<1e-12));
     logger.assert_always(("transform",std::abs(point3D[2]-compare3D[2])<1e-12));
     }*/

    logger.assert_always((test->getTargetDimension() == 3), "getTargetDimension");
    
    //3->3->4
    
    target = &Geometry::ReferenceHypercube::Instance();
    jac.resize(4, 3);
    delete test;
    test = new Geometry::ConcatenatedMapping(Geometry::MappingToRefCubeToCube3::Instance(), Geometry::MappingToRefCubeToHypercube7::Instance());
    nodesAfterTransformation[0] = 4;
    nodesAfterTransformation[1] = 0;
    nodesAfterTransformation[2] = 6;
    nodesAfterTransformation[3] = 2;
    
    for (orig3D[0] = -1.51; orig3D[0] < 1.51; orig3D[0] += 0.6)
    {
        for (orig3D[1] = -1.51; orig3D[1] < 1.51; orig3D[1] += 0.7)
        {
            for (orig3D[2] = -1.51; orig3D[2] < 1.51; orig3D[2] += 0.8)
            {
                point4D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig3D));
                logger.assert_always((source->isInternalPoint(*Geometry::PointReferenceFactory::instance()->makePoint(orig3D)) == target->isInternalPoint(*Geometry::PointReferenceFactory::instance()->makePoint(point4D))), "transform");
                
                orig3D[0] += -1.e-8;
                compare4D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig3D));
                orig3D[0] += 2.e-8;
                point4D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig3D));
                
                orig3D[0] += -1e-8;
                jac = test->calcJacobian(*Geometry::PointReferenceFactory::instance()->makePoint(orig3D));
                logger.assert_always((std::abs(jac[0] - 5.e7 * (point4D[0] - compare4D[0])) < 1e-5), "jacobian"); //estimate is a bit rough, but should work for most mappings
                logger.assert_always((std::abs(jac[1] - 5.e7 * (point4D[1] - compare4D[1])) < 1e-5), "jacobian"); //implementations are very strongly recommended to be more accurate
                logger.assert_always((std::abs(jac[2] - 5.e7 * (point4D[2] - compare4D[2])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[3] - 5.e7 * (point4D[3] - compare4D[3])) < 1e-5), "jacobian");
                
                orig3D[1] += -1.e-8;
                compare4D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig3D));
                orig3D[1] += 2.e-8;
                point4D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig3D));
                
                orig3D[1] += -1e-8;
                jac = test->calcJacobian(*Geometry::PointReferenceFactory::instance()->makePoint(orig3D));
                logger.assert_always((std::abs(jac[4] - 5.e7 * (point4D[0] - compare4D[0])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[5] - 5.e7 * (point4D[1] - compare4D[1])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[6] - 5.e7 * (point4D[2] - compare4D[2])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[7] - 5.e7 * (point4D[3] - compare4D[3])) < 1e-5), "jacobian");
                
                orig3D[2] += -1.e-8;
                compare4D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig3D));
                orig3D[2] += 2.e-8;
                point4D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig3D));
                
                orig3D[2] += -1e-8;
                jac = test->calcJacobian(*Geometry::PointReferenceFactory::instance()->makePoint(orig3D));
                logger.assert_always((std::abs(jac[8] - 5.e7 * (point4D[0] - compare4D[0])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[9] - 5.e7 * (point4D[1] - compare4D[1])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[10] - 5.e7 * (point4D[2] - compare4D[2])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[11] - 5.e7 * (point4D[3] - compare4D[3])) < 1e-5), "jacobian");
            }
        }
    }
    
    /*for(std::size_t i=0;i<source->getNumberOfNodes();++i){///\TODO figure out how the cube->cube maps actually map their nodes
     source->getNode(i,orig3D);
     target->getNode(nodesAfterTransformation[i],compare4D);
     test->transform(orig3D,point4D);
     logger.assert_always(("transform",std::abs(point4D[0]-compare4D[0])<1e-12));
     logger.assert_always(("transform",std::abs(point4D[1]-compare4D[1])<1e-12));
     logger.assert_always(("transform",std::abs(point4D[2]-compare4D[2])<1e-12));
     logger.assert_always(("transform",std::abs(point4D[3]-compare4D[3])<1e-12));
     }*/

    logger.assert_always((test->getTargetDimension() == 4), "getTargetDimension");
    
    //0(->1->)2(->3->)4 (chaining)
    
    source = &Geometry::ReferencePoint::Instance();
    jac.resize(4, 0);
    nodesAfterTransformation.resize(1);
    delete test;
    
    Geometry::ConcatenatedMapping map1(Geometry::MappingToRefPointToLine1::Instance(), Geometry::MappingToRefLineToSquare1::Instance());
    Geometry::ConcatenatedMapping map2(Geometry::MappingToRefSquareToCube1::Instance(), Geometry::MappingToRefCubeToHypercube1::Instance());
    
    test = new Geometry::ConcatenatedMapping(map1, map2);
    nodesAfterTransformation[0] = 8;
    
    point4D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig0D));
    logger.assert_always((source->isInternalPoint(*Geometry::PointReferenceFactory::instance()->makePoint(orig0D)) == target->isInternalPoint(*Geometry::PointReferenceFactory::instance()->makePoint(point4D))), "transform");
    
    jac = test->calcJacobian(*Geometry::PointReferenceFactory::instance()->makePoint(orig0D));
    
    for (std::size_t i = 0; i < source->getNumberOfNodes(); ++i)
    { ///\TODO figure out how the cube->cube maps actually map their nodes
        orig0D = source->getNode(i);
        compare4D = target->getNode(nodesAfterTransformation[i]);
        point4D = test->transform(*Geometry::PointReferenceFactory::instance()->makePoint(orig0D));
        logger.assert_always((std::abs(point4D[0] - compare4D[0]) < 1e-12), "transform");
        logger.assert_always((std::abs(point4D[1] - compare4D[1]) < 1e-12), "transform");
        logger.assert_always((std::abs(point4D[2] - compare4D[2]) < 1e-12), "transform");
        logger.assert_always((std::abs(point4D[3] - compare4D[3]) < 1e-12), "transform");
    }
    
    logger.assert_always((test->getTargetDimension() == 4), "getTargetDimension");
    return 0;
}


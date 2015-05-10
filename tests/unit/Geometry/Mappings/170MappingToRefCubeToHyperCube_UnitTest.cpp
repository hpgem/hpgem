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
#include "Geometry/Mappings/MappingToRefCubeToHypercube.h"
#include "Logger.h"

#include "Geometry/ReferenceCube.h"
#include "Geometry/ReferenceHypercube.h"
#include "Geometry/PointReference.h"
#include "Geometry/Jacobian.h"
#include "LinearAlgebra/NumericalVector.h"
#include <cmath>
int main()
{
    
    Geometry::PointReference refPoint(3), point(4), compare(4);
    
    Geometry::ReferenceHypercube& eGeom = Geometry::ReferenceHypercube::Instance();
    Geometry::ReferenceCube& fGeom = Geometry::ReferenceCube::Instance();
    
    Geometry::Jacobian jac(4, 3);
    
    std::vector<std::size_t> nodesAfterTransformation(8);
    
    const Geometry::MappingReferenceToReference* test = &Geometry::MappingToRefCubeToHypercube0::Instance();
    nodesAfterTransformation[0] = 0;
    nodesAfterTransformation[1] = 1;
    nodesAfterTransformation[2] = 2;
    nodesAfterTransformation[3] = 3;
    nodesAfterTransformation[4] = 4;
    nodesAfterTransformation[5] = 5;
    nodesAfterTransformation[6] = 6;
    nodesAfterTransformation[7] = 7;
    
    for (refPoint[0] = -2.8189; refPoint[0] < 3.141; refPoint[0] += 0.1)
    {
        for (refPoint[1] = -2.8189; refPoint[1] < 3.141; refPoint[1] += 0.1)
        {
            for (refPoint[2] = -2.8189; refPoint[2] < 3.141; refPoint[2] += 0.1)
            {
                point = test->transform(refPoint);
                logger.assert_always((fGeom.isInternalPoint(refPoint) == eGeom.isInternalPoint(point)), "transform");
                
                refPoint[0] += -1.e-8;
                compare = test->transform(refPoint);
                refPoint[0] += 2.e-8;
                point = test->transform(refPoint);
                
                refPoint[0] += -1e-8;
                jac = test->calcJacobian(refPoint);
                logger.assert_always((std::abs(jac[0] - 5.e7 * (point[0] - compare[0])) < 1e-5), "jacobian"); //estimate is a bit rough, but should work for most mappings
                logger.assert_always((std::abs(jac[1] - 5.e7 * (point[1] - compare[1])) < 1e-5), "jacobian"); //implementations are very strongly recommended to be more accurate
                logger.assert_always((std::abs(jac[2] - 5.e7 * (point[2] - compare[2])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[3] - 5.e7 * (point[3] - compare[3])) < 1e-5), "jacobian");
                
                refPoint[1] += -1.e-8;
                compare = test->transform(refPoint);
                refPoint[1] += 2.e-8;
                point = test->transform(refPoint);
                
                refPoint[1] += -1e-8;
                jac = test->calcJacobian(refPoint);
                logger.assert_always((std::abs(jac[4] - 5.e7 * (point[0] - compare[0])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[5] - 5.e7 * (point[1] - compare[1])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[6] - 5.e7 * (point[2] - compare[2])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[7] - 5.e7 * (point[3] - compare[3])) < 1e-5), "jacobian");
                
                refPoint[2] += -1.e-8;
                compare = test->transform(refPoint);
                refPoint[2] += 2.e-8;
                point = test->transform(refPoint);
                
                refPoint[2] += -1e-8;
                jac = test->calcJacobian(refPoint);
                logger.assert_always((std::abs(jac[8] - 5.e7 * (point[0] - compare[0])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[9] - 5.e7 * (point[1] - compare[1])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[10] - 5.e7 * (point[2] - compare[2])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[11] - 5.e7 * (point[3] - compare[3])) < 1e-5), "jacobian");
            }
        }
    }
    
    for (std::size_t i = 0; i < fGeom.getNumberOfNodes(); ++i)
    {
        refPoint = fGeom.getNode(i);
        compare = eGeom.getNode(nodesAfterTransformation[i]);
        point = test->transform(refPoint);
        logger.assert_always((std::abs(point[0] - compare[0]) < 1e-12), "transform");
        logger.assert_always((std::abs(point[1] - compare[1]) < 1e-12), "transform");
        logger.assert_always((std::abs(point[2] - compare[2]) < 1e-12), "transform");
        logger.assert_always((std::abs(point[3] - compare[3]) < 1e-12), "transform");
    }
    
    logger.assert_always((test->getTargetDimension() == 4), "getTargetDimension");
    
    test = &Geometry::MappingToRefCubeToHypercube1::Instance();
    nodesAfterTransformation[0] = 0;
    nodesAfterTransformation[1] = 1;
    nodesAfterTransformation[2] = 2;
    nodesAfterTransformation[3] = 3;
    nodesAfterTransformation[4] = 8;
    nodesAfterTransformation[5] = 9;
    nodesAfterTransformation[6] = 10;
    nodesAfterTransformation[7] = 11;
    
    for (refPoint[0] = -2.8189; refPoint[0] < 3.141; refPoint[0] += 0.1)
    {
        for (refPoint[1] = -2.8189; refPoint[1] < 3.141; refPoint[1] += 0.1)
        {
            for (refPoint[2] = -2.8189; refPoint[2] < 3.141; refPoint[2] += 0.1)
            {
                point = test->transform(refPoint);
                logger.assert_always((fGeom.isInternalPoint(refPoint) == eGeom.isInternalPoint(point)), "transform");
                
                refPoint[0] += -1.e-8;
                compare = test->transform(refPoint);
                refPoint[0] += 2.e-8;
                point = test->transform(refPoint);
                
                refPoint[0] += -1e-8;
                jac = test->calcJacobian(refPoint);
                logger.assert_always((std::abs(jac[0] - 5.e7 * (point[0] - compare[0])) < 1e-5), "jacobian"); //estimate is a bit rough, but should work for most mappings
                logger.assert_always((std::abs(jac[1] - 5.e7 * (point[1] - compare[1])) < 1e-5), "jacobian"); //implementations are very strongly recommended to be more accurate
                logger.assert_always((std::abs(jac[2] - 5.e7 * (point[2] - compare[2])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[3] - 5.e7 * (point[3] - compare[3])) < 1e-5), "jacobian");
                
                refPoint[1] += -1.e-8;
                compare = test->transform(refPoint);
                refPoint[1] += 2.e-8;
                point = test->transform(refPoint);
                
                refPoint[1] += -1e-8;
                jac = test->calcJacobian(refPoint);
                logger.assert_always((std::abs(jac[4] - 5.e7 * (point[0] - compare[0])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[5] - 5.e7 * (point[1] - compare[1])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[6] - 5.e7 * (point[2] - compare[2])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[7] - 5.e7 * (point[3] - compare[3])) < 1e-5), "jacobian");
                
                refPoint[2] += -1.e-8;
                compare = test->transform(refPoint);
                refPoint[2] += 2.e-8;
                point = test->transform(refPoint);
                
                refPoint[2] += -1e-8;
                jac = test->calcJacobian(refPoint);
                logger.assert_always((std::abs(jac[8] - 5.e7 * (point[0] - compare[0])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[9] - 5.e7 * (point[1] - compare[1])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[10] - 5.e7 * (point[2] - compare[2])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[11] - 5.e7 * (point[3] - compare[3])) < 1e-5), "jacobian");
            }
        }
    }
    
    for (std::size_t i = 0; i < fGeom.getNumberOfNodes(); ++i)
    {
        refPoint = fGeom.getNode(i);
        compare = eGeom.getNode(nodesAfterTransformation[i]);
        point = test->transform(refPoint);
        logger.assert_always((std::abs(point[0] - compare[0]) < 1e-12), "transform");
        logger.assert_always((std::abs(point[1] - compare[1]) < 1e-12), "transform");
        logger.assert_always((std::abs(point[2] - compare[2]) < 1e-12), "transform");
        logger.assert_always((std::abs(point[3] - compare[3]) < 1e-12), "transform");
    }
    
    logger.assert_always((test->getTargetDimension() == 4), "getTargetDimension");
    
    test = &Geometry::MappingToRefCubeToHypercube2::Instance();
    nodesAfterTransformation[0] = 0;
    nodesAfterTransformation[1] = 1;
    nodesAfterTransformation[2] = 4;
    nodesAfterTransformation[3] = 5;
    nodesAfterTransformation[4] = 8;
    nodesAfterTransformation[5] = 9;
    nodesAfterTransformation[6] = 12;
    nodesAfterTransformation[7] = 13;
    
    for (refPoint[0] = -2.8189; refPoint[0] < 3.141; refPoint[0] += 0.1)
    {
        for (refPoint[1] = -2.8189; refPoint[1] < 3.141; refPoint[1] += 0.1)
        {
            for (refPoint[2] = -2.8189; refPoint[2] < 3.141; refPoint[2] += 0.1)
            {
                point = test->transform(refPoint);
                logger.assert_always((fGeom.isInternalPoint(refPoint) == eGeom.isInternalPoint(point)), "transform");
                
                refPoint[0] += -1.e-8;
                compare = test->transform(refPoint);
                refPoint[0] += 2.e-8;
                point = test->transform(refPoint);
                
                refPoint[0] += -1e-8;
                jac = test->calcJacobian(refPoint);
                logger.assert_always((std::abs(jac[0] - 5.e7 * (point[0] - compare[0])) < 1e-5), "jacobian"); //estimate is a bit rough, but should work for most mappings
                logger.assert_always((std::abs(jac[1] - 5.e7 * (point[1] - compare[1])) < 1e-5), "jacobian"); //implementations are very strongly recommended to be more accurate
                logger.assert_always((std::abs(jac[2] - 5.e7 * (point[2] - compare[2])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[3] - 5.e7 * (point[3] - compare[3])) < 1e-5), "jacobian");
                
                refPoint[1] += -1.e-8;
                compare = test->transform(refPoint);
                refPoint[1] += 2.e-8;
                point = test->transform(refPoint);
                
                refPoint[1] += -1e-8;
                jac = test->calcJacobian(refPoint);
                logger.assert_always((std::abs(jac[4] - 5.e7 * (point[0] - compare[0])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[5] - 5.e7 * (point[1] - compare[1])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[6] - 5.e7 * (point[2] - compare[2])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[7] - 5.e7 * (point[3] - compare[3])) < 1e-5), "jacobian");
                
                refPoint[2] += -1.e-8;
                compare = test->transform(refPoint);
                refPoint[2] += 2.e-8;
                point = test->transform(refPoint);
                
                refPoint[2] += -1e-8;
                jac = test->calcJacobian(refPoint);
                logger.assert_always((std::abs(jac[8] - 5.e7 * (point[0] - compare[0])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[9] - 5.e7 * (point[1] - compare[1])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[10] - 5.e7 * (point[2] - compare[2])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[11] - 5.e7 * (point[3] - compare[3])) < 1e-5), "jacobian");
            }
        }
    }
    
    for (std::size_t i = 0; i < fGeom.getNumberOfNodes(); ++i)
    {
        refPoint = fGeom.getNode(i);
        compare = eGeom.getNode(nodesAfterTransformation[i]);
        point = test->transform(refPoint);
        logger.assert_always((std::abs(point[0] - compare[0]) < 1e-12), "transform");
        logger.assert_always((std::abs(point[1] - compare[1]) < 1e-12), "transform");
        logger.assert_always((std::abs(point[2] - compare[2]) < 1e-12), "transform");
        logger.assert_always((std::abs(point[3] - compare[3]) < 1e-12), "transform");
    }
    
    logger.assert_always((test->getTargetDimension() == 4), "getTargetDimension");
    
    test = &Geometry::MappingToRefCubeToHypercube3::Instance();
    nodesAfterTransformation[0] = 0;
    nodesAfterTransformation[1] = 2;
    nodesAfterTransformation[2] = 4;
    nodesAfterTransformation[3] = 6;
    nodesAfterTransformation[4] = 8;
    nodesAfterTransformation[5] = 10;
    nodesAfterTransformation[6] = 12;
    nodesAfterTransformation[7] = 14;
    
    for (refPoint[0] = -2.8189; refPoint[0] < 3.141; refPoint[0] += 0.1)
    {
        for (refPoint[1] = -2.8189; refPoint[1] < 3.141; refPoint[1] += 0.1)
        {
            for (refPoint[2] = -2.8189; refPoint[2] < 3.141; refPoint[2] += 0.1)
            {
                point = test->transform(refPoint);
                logger.assert_always((fGeom.isInternalPoint(refPoint) == eGeom.isInternalPoint(point)), "transform");
                
                refPoint[0] += -1.e-8;
                compare = test->transform(refPoint);
                refPoint[0] += 2.e-8;
                point = test->transform(refPoint);
                
                refPoint[0] += -1e-8;
                jac = test->calcJacobian(refPoint);
                refPoint[0] += -1e-8;
                jac = test->calcJacobian(refPoint);
                logger.assert_always((std::abs(jac[0] - 5.e7 * (point[0] - compare[0])) < 1e-5), "jacobian"); //estimate is a bit rough, but should work for most mappings
                logger.assert_always((std::abs(jac[1] - 5.e7 * (point[1] - compare[1])) < 1e-5), "jacobian"); //implementations are very strongly recommended to be more accurate
                logger.assert_always((std::abs(jac[2] - 5.e7 * (point[2] - compare[2])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[3] - 5.e7 * (point[3] - compare[3])) < 1e-5), "jacobian");
                
                refPoint[1] += -1.e-8;
                compare = test->transform(refPoint);
                refPoint[1] += 2.e-8;
                point = test->transform(refPoint);
                
                refPoint[1] += -1e-8;
                jac = test->calcJacobian(refPoint);
                logger.assert_always((std::abs(jac[4] - 5.e7 * (point[0] - compare[0])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[5] - 5.e7 * (point[1] - compare[1])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[6] - 5.e7 * (point[2] - compare[2])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[7] - 5.e7 * (point[3] - compare[3])) < 1e-5), "jacobian");
                
                refPoint[2] += -1.e-8;
                compare = test->transform(refPoint);
                refPoint[2] += 2.e-8;
                point = test->transform(refPoint);
                
                refPoint[2] += -1e-8;
                jac = test->calcJacobian(refPoint);
                logger.assert_always((std::abs(jac[8] - 5.e7 * (point[0] - compare[0])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[9] - 5.e7 * (point[1] - compare[1])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[10] - 5.e7 * (point[2] - compare[2])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[11] - 5.e7 * (point[3] - compare[3])) < 1e-5), "jacobian");
                
            }
        }
    }
    
    for (std::size_t i = 0; i < fGeom.getNumberOfNodes(); ++i)
    {
        refPoint = fGeom.getNode(i);
        compare = eGeom.getNode(nodesAfterTransformation[i]);
        point = test->transform(refPoint);
        logger.assert_always((std::abs(point[0] - compare[0]) < 1e-12), "transform");
        logger.assert_always((std::abs(point[1] - compare[1]) < 1e-12), "transform");
        logger.assert_always((std::abs(point[2] - compare[2]) < 1e-12), "transform");
        logger.assert_always((std::abs(point[3] - compare[3]) < 1e-12), "transform");
    }
    
    logger.assert_always((test->getTargetDimension() == 4), "getTargetDimension");
    
    test = &Geometry::MappingToRefCubeToHypercube4::Instance();
    nodesAfterTransformation[0] = 1;
    nodesAfterTransformation[1] = 3;
    nodesAfterTransformation[2] = 5;
    nodesAfterTransformation[3] = 7;
    nodesAfterTransformation[4] = 9;
    nodesAfterTransformation[5] = 11;
    nodesAfterTransformation[6] = 13;
    nodesAfterTransformation[7] = 15;
    
    for (refPoint[0] = -2.8189; refPoint[0] < 3.141; refPoint[0] += 0.1)
    {
        for (refPoint[1] = -2.8189; refPoint[1] < 3.141; refPoint[1] += 0.1)
        {
            for (refPoint[2] = -2.8189; refPoint[2] < 3.141; refPoint[2] += 0.1)
            {
                point = test->transform(refPoint);
                logger.assert_always((fGeom.isInternalPoint(refPoint) == eGeom.isInternalPoint(point)), "transform");
                
                refPoint[0] += -1.e-8;
                compare = test->transform(refPoint);
                refPoint[0] += 2.e-8;
                point = test->transform(refPoint);
                
                refPoint[0] += -1e-8;
                jac = test->calcJacobian(refPoint);
                logger.assert_always((std::abs(jac[0] - 5.e7 * (point[0] - compare[0])) < 1e-5), "jacobian"); //estimate is a bit rough, but should work for most mappings
                logger.assert_always((std::abs(jac[1] - 5.e7 * (point[1] - compare[1])) < 1e-5), "jacobian"); //implementations are very strongly recommended to be more accurate
                logger.assert_always((std::abs(jac[2] - 5.e7 * (point[2] - compare[2])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[3] - 5.e7 * (point[3] - compare[3])) < 1e-5), "jacobian");
                
                refPoint[1] += -1.e-8;
                compare = test->transform(refPoint);
                refPoint[1] += 2.e-8;
                point = test->transform(refPoint);
                
                refPoint[1] += -1e-8;
                jac = test->calcJacobian(refPoint);
                logger.assert_always((std::abs(jac[4] - 5.e7 * (point[0] - compare[0])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[5] - 5.e7 * (point[1] - compare[1])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[6] - 5.e7 * (point[2] - compare[2])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[7] - 5.e7 * (point[3] - compare[3])) < 1e-5), "jacobian");
                
                refPoint[2] += -1.e-8;
                compare = test->transform(refPoint);
                refPoint[2] += 2.e-8;
                point = test->transform(refPoint);
                
                refPoint[2] += -1e-8;
                jac = test->calcJacobian(refPoint);
                logger.assert_always((std::abs(jac[8] - 5.e7 * (point[0] - compare[0])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[9] - 5.e7 * (point[1] - compare[1])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[10] - 5.e7 * (point[2] - compare[2])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[11] - 5.e7 * (point[3] - compare[3])) < 1e-5), "jacobian");
            }
        }
    }
    
    for (std::size_t i = 0; i < fGeom.getNumberOfNodes(); ++i)
    {
        refPoint = fGeom.getNode(i);
        compare = eGeom.getNode(nodesAfterTransformation[i]);
        point = test->transform(refPoint);
        logger.assert_always((std::abs(point[0] - compare[0]) < 1e-12), "transform");
        logger.assert_always((std::abs(point[1] - compare[1]) < 1e-12), "transform");
        logger.assert_always((std::abs(point[2] - compare[2]) < 1e-12), "transform");
        logger.assert_always((std::abs(point[3] - compare[3]) < 1e-12), "transform");
    }
    
    logger.assert_always((test->getTargetDimension() == 4), "getTargetDimension");
    
    test = &Geometry::MappingToRefCubeToHypercube5::Instance();
    nodesAfterTransformation[0] = 2;
    nodesAfterTransformation[1] = 3;
    nodesAfterTransformation[2] = 6;
    nodesAfterTransformation[3] = 7;
    nodesAfterTransformation[4] = 10;
    nodesAfterTransformation[5] = 11;
    nodesAfterTransformation[6] = 14;
    nodesAfterTransformation[7] = 15;
    
    for (refPoint[0] = -2.8189; refPoint[0] < 3.141; refPoint[0] += 0.1)
    {
        for (refPoint[1] = -2.8189; refPoint[1] < 3.141; refPoint[1] += 0.1)
        {
            for (refPoint[2] = -2.8189; refPoint[2] < 3.141; refPoint[2] += 0.1)
            {
                point = test->transform(refPoint);
                logger.assert_always((fGeom.isInternalPoint(refPoint) == eGeom.isInternalPoint(point)), "transform");
                
                refPoint[0] += -1.e-8;
                compare = test->transform(refPoint);
                refPoint[0] += 2.e-8;
                point = test->transform(refPoint);
                
                refPoint[0] += -1e-8;
                jac = test->calcJacobian(refPoint);
                logger.assert_always((std::abs(jac[0] - 5.e7 * (point[0] - compare[0])) < 1e-5), "jacobian"); //estimate is a bit rough, but should work for most mappings
                logger.assert_always((std::abs(jac[1] - 5.e7 * (point[1] - compare[1])) < 1e-5), "jacobian"); //implementations are very strongly recommended to be more accurate
                logger.assert_always((std::abs(jac[2] - 5.e7 * (point[2] - compare[2])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[3] - 5.e7 * (point[3] - compare[3])) < 1e-5), "jacobian");
                
                refPoint[1] += -1.e-8;
                compare = test->transform(refPoint);
                refPoint[1] += 2.e-8;
                point = test->transform(refPoint);
                
                refPoint[1] += -1e-8;
                jac = test->calcJacobian(refPoint);
                logger.assert_always((std::abs(jac[4] - 5.e7 * (point[0] - compare[0])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[5] - 5.e7 * (point[1] - compare[1])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[6] - 5.e7 * (point[2] - compare[2])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[7] - 5.e7 * (point[3] - compare[3])) < 1e-5), "jacobian");
                
                refPoint[2] += -1.e-8;
                compare = test->transform(refPoint);
                refPoint[2] += 2.e-8;
                point = test->transform(refPoint);
                
                refPoint[2] += -1e-8;
                jac = test->calcJacobian(refPoint);
                logger.assert_always((std::abs(jac[8] - 5.e7 * (point[0] - compare[0])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[9] - 5.e7 * (point[1] - compare[1])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[10] - 5.e7 * (point[2] - compare[2])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[11] - 5.e7 * (point[3] - compare[3])) < 1e-5), "jacobian");
            }
        }
    }
    
    for (std::size_t i = 0; i < fGeom.getNumberOfNodes(); ++i)
    {
        refPoint = fGeom.getNode(i);
        compare = eGeom.getNode(nodesAfterTransformation[i]);
        point = test->transform(refPoint);
        logger.assert_always((std::abs(point[0] - compare[0]) < 1e-12), "transform");
        logger.assert_always((std::abs(point[1] - compare[1]) < 1e-12), "transform");
        logger.assert_always((std::abs(point[2] - compare[2]) < 1e-12), "transform");
        logger.assert_always((std::abs(point[3] - compare[3]) < 1e-12), "transform");
    }
    
    logger.assert_always((test->getTargetDimension() == 4), "getTargetDimension");
    
    test = &Geometry::MappingToRefCubeToHypercube6::Instance();
    nodesAfterTransformation[0] = 4;
    nodesAfterTransformation[1] = 5;
    nodesAfterTransformation[2] = 6;
    nodesAfterTransformation[3] = 7;
    nodesAfterTransformation[4] = 12;
    nodesAfterTransformation[5] = 13;
    nodesAfterTransformation[6] = 14;
    nodesAfterTransformation[7] = 15;
    
    for (refPoint[0] = -2.8189; refPoint[0] < 3.141; refPoint[0] += 0.1)
    {
        for (refPoint[1] = -2.8189; refPoint[1] < 3.141; refPoint[1] += 0.1)
        {
            for (refPoint[2] = -2.8189; refPoint[2] < 3.141; refPoint[2] += 0.1)
            {
                point = test->transform(refPoint);
                logger.assert_always((fGeom.isInternalPoint(refPoint) == eGeom.isInternalPoint(point)), "transform");
                
                refPoint[0] += -1.e-8;
                compare = test->transform(refPoint);
                refPoint[0] += 2.e-8;
                point = test->transform(refPoint);
                
                refPoint[0] += -1e-8;
                jac = test->calcJacobian(refPoint);
                logger.assert_always((std::abs(jac[0] - 5.e7 * (point[0] - compare[0])) < 1e-5), "jacobian"); //estimate is a bit rough, but should work for most mappings
                logger.assert_always((std::abs(jac[1] - 5.e7 * (point[1] - compare[1])) < 1e-5), "jacobian"); //implementations are very strongly recommended to be more accurate
                logger.assert_always((std::abs(jac[2] - 5.e7 * (point[2] - compare[2])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[3] - 5.e7 * (point[3] - compare[3])) < 1e-5), "jacobian");
                
                refPoint[1] += -1.e-8;
                compare = test->transform(refPoint);
                refPoint[1] += 2.e-8;
                point = test->transform(refPoint);
                
                refPoint[1] += -1e-8;
                jac = test->calcJacobian(refPoint);
                logger.assert_always((std::abs(jac[4] - 5.e7 * (point[0] - compare[0])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[5] - 5.e7 * (point[1] - compare[1])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[6] - 5.e7 * (point[2] - compare[2])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[7] - 5.e7 * (point[3] - compare[3])) < 1e-5), "jacobian");
                
                refPoint[2] += -1.e-8;
                compare = test->transform(refPoint);
                refPoint[2] += 2.e-8;
                point = test->transform(refPoint);
                
                refPoint[2] += -1e-8;
                jac = test->calcJacobian(refPoint);
                logger.assert_always((std::abs(jac[8] - 5.e7 * (point[0] - compare[0])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[9] - 5.e7 * (point[1] - compare[1])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[10] - 5.e7 * (point[2] - compare[2])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[11] - 5.e7 * (point[3] - compare[3])) < 1e-5), "jacobian");
            }
        }
    }
    
    for (std::size_t i = 0; i < fGeom.getNumberOfNodes(); ++i)
    {
        refPoint = fGeom.getNode(i);
        compare = eGeom.getNode(nodesAfterTransformation[i]);
        point = test->transform(refPoint);
        logger.assert_always((std::abs(point[0] - compare[0]) < 1e-12), "transform");
        logger.assert_always((std::abs(point[1] - compare[1]) < 1e-12), "transform");
        logger.assert_always((std::abs(point[2] - compare[2]) < 1e-12), "transform");
        logger.assert_always((std::abs(point[3] - compare[3]) < 1e-12), "transform");
    }
    
    logger.assert_always((test->getTargetDimension() == 4), "getTargetDimension");
    
    test = &Geometry::MappingToRefCubeToHypercube7::Instance();
    nodesAfterTransformation[0] = 8;
    nodesAfterTransformation[1] = 9;
    nodesAfterTransformation[2] = 10;
    nodesAfterTransformation[3] = 11;
    nodesAfterTransformation[4] = 12;
    nodesAfterTransformation[5] = 13;
    nodesAfterTransformation[6] = 14;
    nodesAfterTransformation[7] = 15;
    
    for (refPoint[0] = -2.8189; refPoint[0] < 3.141; refPoint[0] += 0.1)
    {
        for (refPoint[1] = -2.8189; refPoint[1] < 3.141; refPoint[1] += 0.1)
        {
            for (refPoint[2] = -2.8189; refPoint[2] < 3.141; refPoint[2] += 0.1)
            {
                point = test->transform(refPoint);
                logger.assert_always((fGeom.isInternalPoint(refPoint) == eGeom.isInternalPoint(point)), "transform");
                
                refPoint[0] += -1.e-8;
                compare = test->transform(refPoint);
                refPoint[0] += 2.e-8;
                point = test->transform(refPoint);
                
                refPoint[0] += -1e-8;
                jac = test->calcJacobian(refPoint);
                logger.assert_always((std::abs(jac[0] - 5.e7 * (point[0] - compare[0])) < 1e-5), "jacobian"); //estimate is a bit rough, but should work for most mappings
                logger.assert_always((std::abs(jac[1] - 5.e7 * (point[1] - compare[1])) < 1e-5), "jacobian"); //implementations are very strongly recommended to be more accurate
                logger.assert_always((std::abs(jac[2] - 5.e7 * (point[2] - compare[2])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[3] - 5.e7 * (point[3] - compare[3])) < 1e-5), "jacobian");
                
                refPoint[1] += -1.e-8;
                compare = test->transform(refPoint);
                refPoint[1] += 2.e-8;
                point = test->transform(refPoint);
                
                refPoint[1] += -1e-8;
                jac = test->calcJacobian(refPoint);
                logger.assert_always((std::abs(jac[4] - 5.e7 * (point[0] - compare[0])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[5] - 5.e7 * (point[1] - compare[1])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[6] - 5.e7 * (point[2] - compare[2])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[7] - 5.e7 * (point[3] - compare[3])) < 1e-5), "jacobian");
                
                refPoint[2] += -1.e-8;
                compare = test->transform(refPoint);
                refPoint[2] += 2.e-8;
                point = test->transform(refPoint);
                
                refPoint[2] += -1e-8;
                jac = test->calcJacobian(refPoint);
                logger.assert_always((std::abs(jac[8] - 5.e7 * (point[0] - compare[0])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[9] - 5.e7 * (point[1] - compare[1])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[10] - 5.e7 * (point[2] - compare[2])) < 1e-5), "jacobian");
                logger.assert_always((std::abs(jac[11] - 5.e7 * (point[3] - compare[3])) < 1e-5), "jacobian");
            }
        }
    }
    
    for (std::size_t i = 0; i < fGeom.getNumberOfNodes(); ++i)
    {
        refPoint = fGeom.getNode(i);
        compare = eGeom.getNode(nodesAfterTransformation[i]);
        point = test->transform(refPoint);
        logger.assert_always((std::abs(point[0] - compare[0]) < 1e-12), "transform");
        logger.assert_always((std::abs(point[1] - compare[1]) < 1e-12), "transform");
        logger.assert_always((std::abs(point[2] - compare[2]) < 1e-12), "transform");
        logger.assert_always((std::abs(point[3] - compare[3]) < 1e-12), "transform");
    }
    
    logger.assert_always((test->getTargetDimension() == 4), "getTargetDimension");
    
    return 0;
}


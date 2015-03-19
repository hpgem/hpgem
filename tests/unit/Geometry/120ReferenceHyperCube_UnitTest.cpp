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
#include "Geometry/ReferenceHypercube.h"
#include "Geometry/ReferenceCube.h"
#include "Geometry/ReferenceSquare.h"
#include "Geometry/ReferenceLine.h"
#include "Geometry/ReferencePoint.h"
#include <iostream>
#include "Logger.h"

#include "Geometry/PointReference.h"
#include "Geometry/Mappings/MappingToRefCubeToHypercube.h"
#include "Integration/QuadratureRules/GaussQuadratureRule.h"
#include <cmath>
using Geometry::ReferenceHypercube;

int main()
{
    ReferenceHypercube& test = ReferenceHypercube::Instance();
    
    Geometry::PointReference pTest(4);
    
    //testing basic functionality
    
    for (pTest[0] = -3.141; pTest[0] < -1.; pTest[0] += 0.1)
    {
        for (pTest[1] = -3.1416; pTest[1] < 3.1416; pTest[1] += 0.1)
        {
            for (pTest[2] = -3.1416; pTest[2] < 3.1416; pTest[2] += 0.1)
            {
                for (pTest[3] = -3.1416; pTest[3] < 3.1416; pTest[3] += 0.1)
                {
                    logger.assert_always((!test.isInternalPoint(pTest)), "isInternalPoint");
                }
            }
        }
    }
    for (; pTest[0] < 1; pTest[0] += 0.1)
    {
        for (pTest[1] = -3.1417; pTest[1] < -1.; pTest[1] += 0.1)
        {
            for (pTest[2] = -3.1416; pTest[2] < 3.1416; pTest[2] += 0.1)
            {
                for (pTest[3] = -3.1416; pTest[3] < 3.1416; pTest[3] += 0.1)
                {
                    logger.assert_always((!test.isInternalPoint(pTest)), "isInternalPoint");
                }
            }
        }
        for (; pTest[1] < 1.; pTest[1] += 0.1)
        {
            for (pTest[2] = -3.1416; pTest[2] < -1.; pTest[2] += 0.1)
            {
                for (pTest[3] = -3.1416; pTest[3] < 3.1416; pTest[3] += 0.1)
                {
                    logger.assert_always((!test.isInternalPoint(pTest)), "isInternalPoint");
                }
            }
            for (; pTest[2] < 1.; pTest[2] += 0.1)
            {
                for (pTest[3] = -3.1416; pTest[3] < -1.; pTest[3] += 0.1)
                {
                    logger.assert_always((!test.isInternalPoint(pTest)), "isInternalPoint");
                }
                for (; pTest[3] < 1.; pTest[3] += 0.1)
                {
                    logger.assert_always((test.isInternalPoint(pTest)), "isInternalPoint");
                }
                for (; pTest[3] < 3.1416; pTest[3] += 0.1)
                {
                    logger.assert_always((!test.isInternalPoint(pTest)), "isInternalPoint");
                }
            }
            for (; pTest[2] < 3.141; pTest[2] += 0.1)
            {
                for (pTest[3] = -3.1416; pTest[3] < 3.1416; pTest[3] += 0.1)
                {
                    logger.assert_always((!test.isInternalPoint(pTest)), "isInternalPoint");
                }
            }
        }
        for (; pTest[1] < 3.141; pTest[1] += 0.1)
        {
            for (pTest[2] = -3.1416; pTest[2] < 3.1416; pTest[2] += 0.1)
            {
                for (pTest[3] = -3.1416; pTest[3] < 3.1416; pTest[3] += 0.1)
                {
                    logger.assert_always((!test.isInternalPoint(pTest)), "isInternalPoint");
                }
            }
        }
    }
    for (; pTest[0] < 3.141; pTest[0] += 0.1)
    {
        for (pTest[1] = -3.1416; pTest[1] < 3.1416; pTest[1] += 0.1)
        {
            for (pTest[2] = -3.1416; pTest[2] < 3.1416; pTest[2] += 0.1)
            {
                for (pTest[3] = -3.1416; pTest[3] < 3.1416; pTest[3] += 0.1)
                {
                    logger.assert_always((!test.isInternalPoint(pTest)), "isInternalPoint");
                }
            }
        }
    }
    
    pTest = test.getCenter();
    logger.assert_always((test.isInternalPoint(pTest) && std::abs(pTest[0]) < 1e-12 && std::abs(pTest[1]) < 1e-12) && std::abs(pTest[2]) < 1e-12 && std::abs(pTest[3]) < 1e-12, "getCenter");
    pTest = test.getNode(0);
    logger.assert_always((std::abs(pTest[0] + 1) < 1e-12 && std::abs(pTest[1] + 1) < 1e-12 && std::abs(pTest[2] + 1) < 1e-12 && std::abs(pTest[3] + 1) < 1e-12), "getNode 0");
    pTest = test.getNode(1);
    logger.assert_always((std::abs(pTest[0] - 1) < 1e-12 && std::abs(pTest[1] + 1) < 1e-12 && std::abs(pTest[2] + 1) < 1e-12 && std::abs(pTest[3] + 1) < 1e-12), "getNode 1");
    pTest = test.getNode(2);
    logger.assert_always((std::abs(pTest[0] + 1) < 1e-12 && std::abs(pTest[1] - 1) < 1e-12 && std::abs(pTest[2] + 1) < 1e-12 && std::abs(pTest[3] + 1) < 1e-12), "getNode 2");
    pTest = test.getNode(3);
    logger.assert_always((std::abs(pTest[0] - 1) < 1e-12 && std::abs(pTest[1] - 1) < 1e-12 && std::abs(pTest[2] + 1) < 1e-12 && std::abs(pTest[3] + 1) < 1e-12), "getNode 3");
    pTest = test.getNode(4);
    logger.assert_always((std::abs(pTest[0] + 1) < 1e-12 && std::abs(pTest[1] + 1) < 1e-12 && std::abs(pTest[2] - 1) < 1e-12 && std::abs(pTest[3] + 1) < 1e-12), "getNode 4");
    pTest = test.getNode(5);
    logger.assert_always((std::abs(pTest[0] - 1) < 1e-12 && std::abs(pTest[1] + 1) < 1e-12 && std::abs(pTest[2] - 1) < 1e-12 && std::abs(pTest[3] + 1) < 1e-12), "getNode 5");
    pTest = test.getNode(6);
    logger.assert_always((std::abs(pTest[0] + 1) < 1e-12 && std::abs(pTest[1] - 1) < 1e-12 && std::abs(pTest[2] - 1) < 1e-12 && std::abs(pTest[3] + 1) < 1e-12), "getNode 6");
    pTest = test.getNode(7);
    logger.assert_always((std::abs(pTest[0] - 1) < 1e-12 && std::abs(pTest[1] - 1) < 1e-12 && std::abs(pTest[2] - 1) < 1e-12 && std::abs(pTest[3] + 1) < 1e-12), "getNode 7");
    pTest = test.getNode(8);
    logger.assert_always((std::abs(pTest[0] + 1) < 1e-12 && std::abs(pTest[1] + 1) < 1e-12 && std::abs(pTest[2] + 1) < 1e-12 && std::abs(pTest[3] - 1) < 1e-12), "getNode 8");
    pTest = test.getNode(9);
    logger.assert_always((std::abs(pTest[0] - 1) < 1e-12 && std::abs(pTest[1] + 1) < 1e-12 && std::abs(pTest[2] + 1) < 1e-12 && std::abs(pTest[3] - 1) < 1e-12), "getNode 9");
    pTest = test.getNode(10);
    logger.assert_always((std::abs(pTest[0] + 1) < 1e-12 && std::abs(pTest[1] - 1) < 1e-12 && std::abs(pTest[2] + 1) < 1e-12 && std::abs(pTest[3] - 1) < 1e-12), "getNode 10");
    pTest = test.getNode(11);
    logger.assert_always((std::abs(pTest[0] - 1) < 1e-12 && std::abs(pTest[1] - 1) < 1e-12 && std::abs(pTest[2] + 1) < 1e-12 && std::abs(pTest[3] - 1) < 1e-12), "getNode 11");
    pTest = test.getNode(12);
    logger.assert_always((std::abs(pTest[0] + 1) < 1e-12 && std::abs(pTest[1] + 1) < 1e-12 && std::abs(pTest[2] - 1) < 1e-12 && std::abs(pTest[3] - 1) < 1e-12), "getNode 12");
    pTest = test.getNode(13);
    logger.assert_always((std::abs(pTest[0] - 1) < 1e-12 && std::abs(pTest[1] + 1) < 1e-12 && std::abs(pTest[2] - 1) < 1e-12 && std::abs(pTest[3] - 1) < 1e-12), "getNode 13");
    pTest = test.getNode(14);
    logger.assert_always((std::abs(pTest[0] + 1) < 1e-12 && std::abs(pTest[1] - 1) < 1e-12 && std::abs(pTest[2] - 1) < 1e-12 && std::abs(pTest[3] - 1) < 1e-12), "getNode 14");
    pTest = test.getNode(15);
    logger.assert_always((std::abs(pTest[0] - 1) < 1e-12 && std::abs(pTest[1] - 1) < 1e-12 && std::abs(pTest[2] - 1) < 1e-12 && std::abs(pTest[3] - 1) < 1e-12), "getNode 15");
    std::cout << test.getName();
    
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 0) == 0), "getLocalNodeIndex 0"); //the nodes of the face must always be specified IN THIS SPECIFIC ORDER
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 1) == 1), "getLocalNodeIndex 0"); //im not sure if I like this myself, but this should at least verify
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 2) == 2), "getLocalNodeIndex 0"); //that all face nodes are specified, none are specified twice
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 3) == 3), "getLocalNodeIndex 0"); //and only face nodes are specified and the ordering of the nodes is consistent
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 4) == 4), "getLocalNodeIndex 0"); //across function calls
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 5) == 5), "getLocalNodeIndex 0");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 6) == 6), "getLocalNodeIndex 0");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 7) == 7), "getLocalNodeIndex 0");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 0) == 0), "getLocalNodeIndex 1");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 1) == 1), "getLocalNodeIndex 1");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 2) == 2), "getLocalNodeIndex 1");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 3) == 3), "getLocalNodeIndex 1");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 4) == 8), "getLocalNodeIndex 1");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 5) == 9), "getLocalNodeIndex 1");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 6) == 10), "getLocalNodeIndex 1");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 7) == 11), "getLocalNodeIndex 1");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 0) == 0), "getLocalNodeIndex 2");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 1) == 1), "getLocalNodeIndex 2");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 2) == 4), "getLocalNodeIndex 2");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 3) == 5), "getLocalNodeIndex 2");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 4) == 8), "getLocalNodeIndex 2");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 5) == 9), "getLocalNodeIndex 2");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 6) == 12), "getLocalNodeIndex 2");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 7) == 13), "getLocalNodeIndex 2");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 0) == 0), "getLocalNodeIndex 3");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 1) == 2), "getLocalNodeIndex 3");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 2) == 4), "getLocalNodeIndex 3");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 3) == 6), "getLocalNodeIndex 3");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 4) == 8), "getLocalNodeIndex 3");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 5) == 10), "getLocalNodeIndex 3");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 6) == 12), "getLocalNodeIndex 3");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 7) == 14), "getLocalNodeIndex 3");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 0) == 1), "getLocalNodeIndex 4");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 1) == 3), "getLocalNodeIndex 4");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 2) == 5), "getLocalNodeIndex 4");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 3) == 7), "getLocalNodeIndex 4");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 4) == 9), "getLocalNodeIndex 4");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 5) == 11), "getLocalNodeIndex 4");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 6) == 13), "getLocalNodeIndex 4");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 7) == 15), "getLocalNodeIndex 4");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 0) == 2), "getLocalNodeIndex 5");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 1) == 3), "getLocalNodeIndex 5");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 2) == 6), "getLocalNodeIndex 5");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 3) == 7), "getLocalNodeIndex 5");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 4) == 10), "getLocalNodeIndex 5");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 5) == 11), "getLocalNodeIndex 5");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 6) == 14), "getLocalNodeIndex 5");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 7) == 15), "getLocalNodeIndex 5");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(6, 0) == 4), "getLocalNodeIndex 6");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(6, 1) == 5), "getLocalNodeIndex 6");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(6, 2) == 6), "getLocalNodeIndex 6");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(6, 3) == 7), "getLocalNodeIndex 6");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(6, 4) == 12), "getLocalNodeIndex 6");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(6, 5) == 13), "getLocalNodeIndex 6");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(6, 6) == 14), "getLocalNodeIndex 6");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(6, 7) == 15), "getLocalNodeIndex 6");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(7, 0) == 8), "getLocalNodeIndex 7");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(7, 1) == 9), "getLocalNodeIndex 7");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(7, 2) == 10), "getLocalNodeIndex 7");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(7, 3) == 11), "getLocalNodeIndex 7");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(7, 4) == 12), "getLocalNodeIndex 7");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(7, 5) == 13), "getLocalNodeIndex 7");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(7, 6) == 14), "getLocalNodeIndex 7");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(7, 7) == 15), "getLocalNodeIndex 7");
    
    std::cout << test;
    
    //testing mappings and quadrature rules
    
    std::vector<std::size_t> faceIndices(8);
    //there is no 5D element so codim0mappings are not needed
    
    logger.assert_always((test.getNrOfCodim1Entities() == 8 && test.getNrOfCodim2Entities() == 24) && test.getNrOfCodim3Entities() == 32, "higher codimensional entities");
    logger.assert_always((test.getCodim1ReferenceGeometry(0) == &Geometry::ReferenceCube::Instance() && test.getCodim1ReferenceGeometry(1) == &Geometry::ReferenceCube::Instance() && test.getCodim1ReferenceGeometry(2) == &Geometry::ReferenceCube::Instance() && test.getCodim1ReferenceGeometry(3) == &Geometry::ReferenceCube::Instance() && test.getCodim1ReferenceGeometry(4) == &Geometry::ReferenceCube::Instance() && test.getCodim1ReferenceGeometry(5) == &Geometry::ReferenceCube::Instance() && test.getCodim1ReferenceGeometry(6) == &Geometry::ReferenceCube::Instance() && test.getCodim1ReferenceGeometry(7) == &Geometry::ReferenceCube::Instance()), "getCodim1ReferenceGeometry");
    logger.assert_always((test.getCodim2ReferenceGeometry(0) == &Geometry::ReferenceSquare::Instance() && test.getCodim2ReferenceGeometry(1) == &Geometry::ReferenceSquare::Instance() && test.getCodim2ReferenceGeometry(2) == &Geometry::ReferenceSquare::Instance() && test.getCodim2ReferenceGeometry(3) == &Geometry::ReferenceSquare::Instance() && test.getCodim2ReferenceGeometry(4) == &Geometry::ReferenceSquare::Instance() && test.getCodim2ReferenceGeometry(5) == &Geometry::ReferenceSquare::Instance() && test.getCodim2ReferenceGeometry(6) == &Geometry::ReferenceSquare::Instance() && test.getCodim2ReferenceGeometry(7) == &Geometry::ReferenceSquare::Instance() && test.getCodim2ReferenceGeometry(8) == &Geometry::ReferenceSquare::Instance() && test.getCodim2ReferenceGeometry(9) == &Geometry::ReferenceSquare::Instance() && test.getCodim2ReferenceGeometry(10) == &Geometry::ReferenceSquare::Instance() && test.getCodim2ReferenceGeometry(11) == &Geometry::ReferenceSquare::Instance() && test.getCodim2ReferenceGeometry(12) == &Geometry::ReferenceSquare::Instance() && test.getCodim2ReferenceGeometry(13) == &Geometry::ReferenceSquare::Instance() && test.getCodim2ReferenceGeometry(14) == &Geometry::ReferenceSquare::Instance() && test.getCodim2ReferenceGeometry(15) == &Geometry::ReferenceSquare::Instance() && test.getCodim2ReferenceGeometry(16) == &Geometry::ReferenceSquare::Instance() && test.getCodim2ReferenceGeometry(17) == &Geometry::ReferenceSquare::Instance() && test.getCodim2ReferenceGeometry(18) == &Geometry::ReferenceSquare::Instance() && test.getCodim2ReferenceGeometry(19) == &Geometry::ReferenceSquare::Instance() && test.getCodim2ReferenceGeometry(20) == &Geometry::ReferenceSquare::Instance() && test.getCodim2ReferenceGeometry(21) == &Geometry::ReferenceSquare::Instance() && test.getCodim2ReferenceGeometry(22) == &Geometry::ReferenceSquare::Instance() && test.getCodim2ReferenceGeometry(23) == &Geometry::ReferenceSquare::Instance()), "getCodim2ReferenceGeometry");
    logger.assert_always((test.getCodim1MappingPtr(0) == &Geometry::MappingToRefCubeToHypercube0::Instance()), "getCodim1MappingPtr");
    logger.assert_always((test.getCodim1MappingPtr(1) == &Geometry::MappingToRefCubeToHypercube1::Instance()), "getCodim1MappingPtr");
    logger.assert_always((test.getCodim1MappingPtr(2) == &Geometry::MappingToRefCubeToHypercube2::Instance()), "getCodim1MappingPtr");
    logger.assert_always((test.getCodim1MappingPtr(3) == &Geometry::MappingToRefCubeToHypercube3::Instance()), "getCodim1MappingPtr");
    logger.assert_always((test.getCodim1MappingPtr(4) == &Geometry::MappingToRefCubeToHypercube4::Instance()), "getCodim1MappingPtr");
    logger.assert_always((test.getCodim1MappingPtr(5) == &Geometry::MappingToRefCubeToHypercube5::Instance()), "getCodim1MappingPtr");
    logger.assert_always((test.getCodim1MappingPtr(6) == &Geometry::MappingToRefCubeToHypercube6::Instance()), "getCodim1MappingPtr");
    logger.assert_always((test.getCodim1MappingPtr(7) == &Geometry::MappingToRefCubeToHypercube7::Instance()), "getCodim1MappingPtr");
    faceIndices = test.getCodim1EntityLocalIndices(0);
    logger.assert_always((faceIndices[0] == test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 0)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[1] == test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 1)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[2] == test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 2)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[3] == test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 3)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[4] == test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 4)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[5] == test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 5)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[6] == test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 6)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[7] == test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 7)), "getCodim1EntityLocalIndices");
    faceIndices = test.getCodim1EntityLocalIndices(1);
    logger.assert_always((faceIndices[0] == test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 0)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[1] == test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 1)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[2] == test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 2)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[3] == test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 3)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[4] == test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 4)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[5] == test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 5)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[6] == test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 6)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[7] == test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 7)), "getCodim1EntityLocalIndices");
    faceIndices = test.getCodim1EntityLocalIndices(2);
    logger.assert_always((faceIndices[0] == test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 0)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[1] == test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 1)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[2] == test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 2)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[3] == test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 3)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[4] == test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 4)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[5] == test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 5)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[6] == test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 6)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[7] == test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 7)), "getCodim1EntityLocalIndices");
    faceIndices = test.getCodim1EntityLocalIndices(3);
    logger.assert_always((faceIndices[0] == test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 0)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[1] == test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 1)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[2] == test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 2)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[3] == test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 3)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[4] == test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 4)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[5] == test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 5)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[6] == test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 6)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[7] == test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 7)), "getCodim1EntityLocalIndices");
    faceIndices = test.getCodim1EntityLocalIndices(4);
    logger.assert_always((faceIndices[0] == test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 0)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[1] == test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 1)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[2] == test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 2)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[3] == test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 3)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[4] == test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 4)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[5] == test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 5)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[6] == test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 6)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[7] == test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 7)), "getCodim1EntityLocalIndices");
    faceIndices = test.getCodim1EntityLocalIndices(5);
    logger.assert_always((faceIndices[0] == test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 0)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[1] == test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 1)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[2] == test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 2)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[3] == test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 3)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[4] == test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 4)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[5] == test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 5)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[6] == test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 6)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[7] == test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 7)), "getCodim1EntityLocalIndices");
    faceIndices = test.getCodim1EntityLocalIndices(6);
    logger.assert_always((faceIndices[0] == test.getLocalNodeIndexFromFaceAndIndexOnFace(6, 0)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[1] == test.getLocalNodeIndexFromFaceAndIndexOnFace(6, 1)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[2] == test.getLocalNodeIndexFromFaceAndIndexOnFace(6, 2)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[3] == test.getLocalNodeIndexFromFaceAndIndexOnFace(6, 3)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[4] == test.getLocalNodeIndexFromFaceAndIndexOnFace(6, 4)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[5] == test.getLocalNodeIndexFromFaceAndIndexOnFace(6, 5)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[6] == test.getLocalNodeIndexFromFaceAndIndexOnFace(6, 6)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[7] == test.getLocalNodeIndexFromFaceAndIndexOnFace(6, 7)), "getCodim1EntityLocalIndices");
    faceIndices = test.getCodim1EntityLocalIndices(7);
    logger.assert_always((faceIndices[0] == test.getLocalNodeIndexFromFaceAndIndexOnFace(7, 0)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[1] == test.getLocalNodeIndexFromFaceAndIndexOnFace(7, 1)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[2] == test.getLocalNodeIndexFromFaceAndIndexOnFace(7, 2)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[3] == test.getLocalNodeIndexFromFaceAndIndexOnFace(7, 3)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[4] == test.getLocalNodeIndexFromFaceAndIndexOnFace(7, 4)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[5] == test.getLocalNodeIndexFromFaceAndIndexOnFace(7, 5)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[6] == test.getLocalNodeIndexFromFaceAndIndexOnFace(7, 6)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[7] == test.getLocalNodeIndexFromFaceAndIndexOnFace(7, 7)), "getCodim1EntityLocalIndices");
    
    //other codimensions are not implemented
    
    logger.assert_always((test.getGaussQuadratureRule(3)->order() >= 3), "quadrature rules");
    //assert(("quadrature rules",test.getGaussQuadratureRule(5)->order()>=5));///\TODO implement more quadrature rules
    //assert(("quadrature rules",test.getGaussQuadratureRule(7)->order()>=7));
    //assert(("quadrature rules",test.getGaussQuadratureRule(9)->order()>=9));
    //assert(("quadrature rules",test.getGaussQuadratureRule(11)->order()>=11));
    
    //testing functionality of abstract parent classes
    
    logger.assert_always((test.getNumberOfNodes() == 16), "number of nodes");
    logger.assert_always((test.getGeometryType() == Geometry::ReferenceGeometryType::HYPERCUBE), "type of geometry");
    
    ///\TODO if it is decided that getBasisFunctionValue and getBasisFucntionDerivative remain here, test them
    
    ///\TODO testing that the refinement maps behave exactly like the forwarded calls of this class
}


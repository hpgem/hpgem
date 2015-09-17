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
#include "Geometry/ReferenceTriangularPrism.h"
#include "Geometry/ReferenceTriangle.h"
#include "Geometry/ReferenceSquare.h"
#include "Geometry/ReferenceLine.h"
#include "Geometry/ReferencePoint.h"
#include <iostream>
#include "Logger.h"

#include "Geometry/PointReference.h"
#include "Geometry/Mappings/MappingToRefFaceToTriangularPrism.h"
#include "Integration/QuadratureRules/GaussQuadratureRule.h"
#include <cmath>
using Geometry::ReferenceTriangularPrism;

int main()
{
    ReferenceTriangularPrism& test = ReferenceTriangularPrism::Instance();
    
    Geometry::Point<3> pTest;
    
    //testing basic functionality
    
    for (pTest[0] = -1.51; pTest[0] < 0.; pTest[0] += 0.2)
    {
        for (pTest[1] = -1.51; pTest[1] < 1.51; pTest[1] += 0.2)
        {
            for (pTest[2] = -1.51; pTest[2] < 1.51; pTest[2] += 0.2)
            {
                logger.assert_always((!test.isInternalPoint(*Geometry::PointReferenceFactory<3>::instance()->makePoint(pTest))), "isInternalPoint");
            }
        }
    }
    for (; pTest[0] < 1; pTest[0] += 0.2)
    {
        for (pTest[1] = -1.51; pTest[1] < 0.; pTest[1] += 0.2)
        {
            for (pTest[2] = -1.51; pTest[2] < 1.51; pTest[2] += 0.2)
            {
                logger.assert_always((!test.isInternalPoint(*Geometry::PointReferenceFactory<3>::instance()->makePoint(pTest))), "isInternalPoint");
            }
        }
        for (; pTest[1] < 1. - pTest[0]; pTest[1] += 0.2)
        {
            for (pTest[2] = -1.51; pTest[2] < -1.; pTest[2] += 0.2)
            {
                logger.assert_always((!test.isInternalPoint(*Geometry::PointReferenceFactory<3>::instance()->makePoint(pTest))), "isInternalPoint");
            }
            for (; pTest[2] < 1.; pTest[2] += 0.2)
            {
                logger.assert_always((test.isInternalPoint(*Geometry::PointReferenceFactory<3>::instance()->makePoint(pTest))), "isInternalPoint");
            }
            for (; pTest[2] < 1.51; pTest[2] += 0.2)
            {
                logger.assert_always((!test.isInternalPoint(*Geometry::PointReferenceFactory<3>::instance()->makePoint(pTest))), "isInternalPoint");
            }
        }
        for (; pTest[1] < 1.51; pTest[1] += 0.2)
        {
            for (pTest[2] = -1.51; pTest[2] < 1.51; pTest[2] += 0.2)
            {
                logger.assert_always((!test.isInternalPoint(*Geometry::PointReferenceFactory<3>::instance()->makePoint(pTest))), "isInternalPoint");
            }
        }
    }
    for (; pTest[0] < 1.51; pTest[0] += 0.2)
    {
        for (pTest[1] = -1.51; pTest[1] < 1.51; pTest[1] += 0.2)
        {
            for (pTest[2] = -1.51; pTest[2] < 1.51; pTest[2] += 0.2)
            {
                logger.assert_always((!test.isInternalPoint(*Geometry::PointReferenceFactory<3>::instance()->makePoint(pTest))), "isInternalPoint");
            }
        }
    }
    
    pTest = test.getCenter();
    logger.assert_always((test.isInternalPoint(*Geometry::PointReferenceFactory<3>::instance()->makePoint(pTest)) && std::abs(pTest[0] - 1. / 3.) < 1e-12 && std::abs(pTest[1] - 1. / 3.) < 1e-12) && std::abs(pTest[2]) < 1e-12, "getCenter");
    pTest = test.getReferenceNodeCoordinate(0);
    logger.assert_always((std::abs(pTest[0]) < 1e-12 && std::abs(pTest[1]) < 1e-12 && std::abs(pTest[2] + 1) < 1e-12), "getNode 0");
    pTest = test.getReferenceNodeCoordinate(1);
    logger.assert_always((std::abs(pTest[0] - 1) < 1e-12 && std::abs(pTest[1]) < 1e-12 && std::abs(pTest[2] + 1) < 1e-12), "getNode 1");
    pTest = test.getReferenceNodeCoordinate(2);
    logger.assert_always((std::abs(pTest[0]) < 1e-12 && std::abs(pTest[1] - 1) < 1e-12 && std::abs(pTest[2] + 1) < 1e-12), "getNode 2");
    pTest = test.getReferenceNodeCoordinate(3);
    logger.assert_always((std::abs(pTest[0]) < 1e-12 && std::abs(pTest[1]) < 1e-12 && std::abs(pTest[2] - 1) < 1e-12), "getNode 3");
    pTest = test.getReferenceNodeCoordinate(4);
    logger.assert_always((std::abs(pTest[0] - 1) < 1e-12 && std::abs(pTest[1]) < 1e-12 && std::abs(pTest[2] - 1) < 1e-12), "getNode 4");
    pTest = test.getReferenceNodeCoordinate(5);
    logger.assert_always((std::abs(pTest[0]) < 1e-12 && std::abs(pTest[1] - 1) < 1e-12 && std::abs(pTest[2] - 1) < 1e-12), "getNode 5");
    std::cout << test.getName();
    
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 0) == 0), "getLocalNodeIndex 0"); //the nodes of the face must always be specified IN THIS SPECIFIC ORDER
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 1) == 2), "getLocalNodeIndex 0"); //im not sure if I like this myself, but this should at least verify
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 2) == 1), "getLocalNodeIndex 0"); //that all face nodes are specified, none are specified twice
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 0) == 3), "getLocalNodeIndex 1"); //and only face nodes are specified and the ordering of the nodes is consistent
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 1) == 4), "getLocalNodeIndex 1"); //across function calls
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 2) == 5), "getLocalNodeIndex 1");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 0) == 2), "getLocalNodeIndex 2");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 1) == 0), "getLocalNodeIndex 2");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 2) == 5), "getLocalNodeIndex 2");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 3) == 3), "getLocalNodeIndex 2");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 0) == 0), "getLocalNodeIndex 3");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 1) == 1), "getLocalNodeIndex 3");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 2) == 3), "getLocalNodeIndex 3");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 3) == 4), "getLocalNodeIndex 3");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 0) == 1), "getLocalNodeIndex 4");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 1) == 2), "getLocalNodeIndex 4");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 2) == 4), "getLocalNodeIndex 4");
    logger.assert_always((test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 3) == 5), "getLocalNodeIndex 4");
    
    std::cout << test;
    
    //testing mappings and quadrature rules
    
    std::vector<std::size_t> faceIndices(3);
    //codim0maps dont exist so they dont need to be found properly
    
    logger.assert_always((test.getNrOfCodim1Entities() == 5 && test.getNrOfCodim2Entities() == 9) && test.getNrOfCodim3Entities() == 6, "higher codimensional entities");
    logger.assert_always((test.getCodim1ReferenceGeometry(0) == &Geometry::ReferenceTriangle::Instance() && test.getCodim1ReferenceGeometry(1) == &Geometry::ReferenceTriangle::Instance() && test.getCodim1ReferenceGeometry(2) == &Geometry::ReferenceSquare::Instance() && test.getCodim1ReferenceGeometry(3) == &Geometry::ReferenceSquare::Instance() && test.getCodim1ReferenceGeometry(4) == &Geometry::ReferenceSquare::Instance()), "getCodim1ReferenceGeometry");
    logger.assert_always((test.getCodim2ReferenceGeometry(0) == &Geometry::ReferenceLine::Instance() && test.getCodim2ReferenceGeometry(1) == &Geometry::ReferenceLine::Instance() && test.getCodim2ReferenceGeometry(2) == &Geometry::ReferenceLine::Instance() && test.getCodim2ReferenceGeometry(3) == &Geometry::ReferenceLine::Instance() && test.getCodim2ReferenceGeometry(4) == &Geometry::ReferenceLine::Instance() && test.getCodim2ReferenceGeometry(5) == &Geometry::ReferenceLine::Instance() && test.getCodim2ReferenceGeometry(6) == &Geometry::ReferenceLine::Instance() && test.getCodim2ReferenceGeometry(7) == &Geometry::ReferenceLine::Instance() && test.getCodim2ReferenceGeometry(8) == &Geometry::ReferenceLine::Instance()), "getCodim2ReferenceGeometry");
    logger.assert_always((test.getCodim1MappingPtr(0) == &Geometry::MappingToRefFaceToTriangularPrism0::Instance()), "getCodim1MappingPtr");
    logger.assert_always((test.getCodim1MappingPtr(1) == &Geometry::MappingToRefFaceToTriangularPrism1::Instance()), "getCodim1MappingPtr");
    logger.assert_always((test.getCodim1MappingPtr(2) == &Geometry::MappingToRefFaceToTriangularPrism2::Instance()), "getCodim1MappingPtr");
    logger.assert_always((test.getCodim1MappingPtr(3) == &Geometry::MappingToRefFaceToTriangularPrism3::Instance()), "getCodim1MappingPtr");
    logger.assert_always((test.getCodim1MappingPtr(4) == &Geometry::MappingToRefFaceToTriangularPrism4::Instance()), "getCodim1MappingPtr");
    faceIndices = test.getCodim1EntityLocalIndices(0);
    logger.assert_always((faceIndices[0] == test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 0)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[1] == test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 1)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[2] == test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 2)), "getCodim1EntityLocalIndices");
    faceIndices = test.getCodim1EntityLocalIndices(1);
    logger.assert_always((faceIndices[0] == test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 0)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[1] == test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 1)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[2] == test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 2)), "getCodim1EntityLocalIndices");
    faceIndices.resize(4);
    faceIndices = test.getCodim1EntityLocalIndices(2);
    logger.assert_always((faceIndices[0] == test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 0)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[1] == test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 1)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[2] == test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 2)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[3] == test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 3)), "getCodim1EntityLocalIndices");
    faceIndices = test.getCodim1EntityLocalIndices(3);
    logger.assert_always((faceIndices[0] == test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 0)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[1] == test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 1)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[2] == test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 2)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[3] == test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 3)), "getCodim1EntityLocalIndices");
    faceIndices = test.getCodim1EntityLocalIndices(4);
    logger.assert_always((faceIndices[0] == test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 0)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[1] == test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 1)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[2] == test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 2)), "getCodim1EntityLocalIndices");
    logger.assert_always((faceIndices[3] == test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 3)), "getCodim1EntityLocalIndices");
    faceIndices.resize(2);
    faceIndices = test.getCodim2EntityLocalIndices(0);
    logger.assert_always((faceIndices[0] == 0), "getCodim2EntityLocalIndices");
    logger.assert_always((faceIndices[1] == 1), "getCodim2EntityLocalIndices");
    faceIndices = test.getCodim2EntityLocalIndices(1);
    logger.assert_always((faceIndices[0] == 0), "getCodim2EntityLocalIndices");
    logger.assert_always((faceIndices[1] == 2), "getCodim2EntityLocalIndices");
    faceIndices = test.getCodim2EntityLocalIndices(2);
    logger.assert_always((faceIndices[0] == 1), "getCodim2EntityLocalIndices");
    logger.assert_always((faceIndices[1] == 2), "getCodim2EntityLocalIndices");
    faceIndices = test.getCodim2EntityLocalIndices(3);
    logger.assert_always((faceIndices[0] == 3), "getCodim2EntityLocalIndices");
    logger.assert_always((faceIndices[1] == 4), "getCodim2EntityLocalIndices");
    faceIndices = test.getCodim2EntityLocalIndices(4);
    logger.assert_always((faceIndices[0] == 3), "getCodim2EntityLocalIndices");
    logger.assert_always((faceIndices[1] == 5), "getCodim2EntityLocalIndices");
    faceIndices = test.getCodim2EntityLocalIndices(5);
    logger.assert_always((faceIndices[0] == 4), "getCodim2EntityLocalIndices");
    logger.assert_always((faceIndices[1] == 5), "getCodim2EntityLocalIndices");
    faceIndices = test.getCodim2EntityLocalIndices(6);
    logger.assert_always((faceIndices[0] == 0), "getCodim2EntityLocalIndices");
    logger.assert_always((faceIndices[1] == 3), "getCodim2EntityLocalIndices");
    faceIndices = test.getCodim2EntityLocalIndices(7);
    logger.assert_always((faceIndices[0] == 1), "getCodim2EntityLocalIndices");
    logger.assert_always((faceIndices[1] == 4), "getCodim2EntityLocalIndices");
    faceIndices = test.getCodim2EntityLocalIndices(8);
    logger.assert_always((faceIndices[0] == 2), "getCodim2EntityLocalIndices");
    logger.assert_always((faceIndices[1] == 5), "getCodim2EntityLocalIndices");
    faceIndices.resize(1);
    faceIndices = test.getCodim3EntityLocalIndices(0);
    logger.assert_always((faceIndices[0] == 0), "getCodim3EntityLocalIndices");
    faceIndices = test.getCodim3EntityLocalIndices(1);
    logger.assert_always((faceIndices[0] == 1), "getCodim3EntityLocalIndices");
    faceIndices = test.getCodim3EntityLocalIndices(2);
    logger.assert_always((faceIndices[0] == 2), "getCodim3EntityLocalIndices");
    faceIndices = test.getCodim3EntityLocalIndices(3);
    logger.assert_always((faceIndices[0] == 3), "getCodim3EntityLocalIndices");
    faceIndices = test.getCodim3EntityLocalIndices(4);
    logger.assert_always((faceIndices[0] == 4), "getCodim3EntityLocalIndices");
    faceIndices = test.getCodim3EntityLocalIndices(5);
    logger.assert_always((faceIndices[0] == 5), "getCodim3EntityLocalIndices");
    
    logger.assert_always((test.getGaussQuadratureRule(3)->order() >= 3), "quadrature rules");
    logger.assert_always((test.getGaussQuadratureRule(5)->order() >= 5), "quadrature rules");
    logger.assert_always((test.getGaussQuadratureRule(7)->order() >= 7), "quadrature rules");
    ///\todo add extra quadrature rules
    //assert(("quadrature rules",test.getGaussQuadratureRule(9)->order()>=9));
    //assert(("quadrature rules",test.getGaussQuadratureRule(11)->order()>=11));
    
    //testing functionality of abstract parent classes
    
    logger.assert_always((test.getNumberOfNodes() == 6), "number of nodes");
    logger.assert_always((test.getGeometryType() == Geometry::ReferenceGeometryType::TRIANGULARPRISM), "type of geometry");
    
    ///\todo testing that the refinement maps behave exactly like the forwarded calls of this class
}


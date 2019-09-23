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

#include "Base/MeshManipulator.h"
#include "Geometry/PointPhysical.h"
#include "Base/Element.h"
#include "Geometry/PhysicalGeometry.h"
#include "Base/CommandLineOptions.h"
#include <vector>
#include "CMakeDefinitions.h"
#include <sstream>
#include "Output/TecplotDiscontinuousSolutionWriter.h"

using namespace std;

class Dummy2
{
public:
    Dummy2(){}
    void operator()(const Base::Element* el, const Geometry::PointReference<2>& p, ostream& os)
    {
    }
};

class Dummy3
{
public:
    Dummy3(){}
    void operator()(const Base::Element* el, const Geometry::PointReference<3>& p, ostream& os)
    {
    }
};

template<std::size_t DIM>
void testMesh(Base::MeshManipulator<DIM>* test)
{
    std::unordered_set<std::size_t> elementIDs, faceIDs, edgeIDs, nodeIDs;
    std::cout << test->getElementsList(Base::IteratorType::GLOBAL).size() << std::endl;
    std::cout << test->getElementsList(Base::IteratorType::LOCAL).size() << std::endl;
    for (Base::Element* element : test->getElementsList())
    {
        logger.assert_always((elementIDs.find(element->getID()) == elementIDs.end()), "duplicate element ID");
        elementIDs.insert(element->getID());
        logger.assert_always((element->getNumberOfFaces() == element->getReferenceGeometry()->getNumberOfCodim1Entities()), "confusion about the number of faces");
        if (test->dimension() == 2)
        {
            logger.assert_always((element->getNumberOfEdges() == 0), "confusion about the number of edges");
        }
        else
        {
            logger.assert_always((element->getNumberOfEdges() == element->getReferenceGeometry()->getNumberOfCodim2Entities()), "confusion about the number of edges");
        }
        logger.assert_always((element->getNumberOfNodes() == element->getReferenceGeometry()->getNumberOfNodes()), "confusion about the number of vertices");
        for (std::size_t i = 0; i < element->getNumberOfFaces(); ++i)
        {
            logger.assert_always((element->getFace(i) != nullptr), "missing Face no. % in element %", i, element->getID());
        }
        for (std::size_t i = 0; i < element->getNumberOfEdges(); ++i)
        {
            logger.assert_always((element->getEdge(i) != nullptr), "missing Face");
        }
        for (std::size_t i = 0; i < element->getNumberOfNodes(); ++i)
        {
            logger.assert_always((element->getNode(i) != nullptr), "missing Face");
        }
    }
    for (Base::Face* face : test->getFacesList())
    {
        logger.assert_always((faceIDs.find(face->getID()) == faceIDs.end()), "duplicate face ID");
        faceIDs.insert(face->getID());
        logger.assert_always((face->getPtrElementLeft()->getFace(face->localFaceNumberLeft()) == face), "element<->face matching");
        if (face->isInternal())
        {
            std::vector<std::size_t> leftNodes(face->getReferenceGeometry()->getNumberOfNodes()), rightNodes(leftNodes);
            leftNodes = face->getPtrElementLeft()->getPhysicalGeometry()->getGlobalFaceNodeIndices(face->localFaceNumberLeft());
            rightNodes = face->getPtrElementRight()->getPhysicalGeometry()->getGlobalFaceNodeIndices(face->localFaceNumberRight());
            for (std::size_t i = 0; i < leftNodes.size(); ++i)
            {
                bool found = false;
                for (std::size_t j = 0; j < rightNodes.size(); ++j)
                {
                    found |= leftNodes[i] == rightNodes[j];
                }
                logger.assert_always((found), "face positioning1");
            }
            logger.assert_always((leftNodes.size() == rightNodes.size()), "face positioning2");
            logger.assert_always((face->getPtrElementRight()->getFace(face->localFaceNumberRight()) == face), "element<->face matching");
        }
    }
    for (Base::Edge* edge : test->getEdgesList())
    {
        logger.assert_always((edgeIDs.find(edge->getID()) == edgeIDs.end()), "duplicate edge ID");
        edgeIDs.insert(edge->getID());
        logger.assert_always((edge->getElement(0)->getEdge(edge->getEdgeNumber(0)) == edge), "element<->edge matching");
        std::vector<std::size_t> firstNodes(edge->getElement(0)->getReferenceGeometry()->getCodim2ReferenceGeometry(edge->getEdgeNumber(0))->getNumberOfNodes()), otherNodes(firstNodes);
        firstNodes = edge->getElement(0)->getReferenceGeometry()->getCodim2EntityLocalIndices(edge->getEdgeNumber(0));
        for (std::size_t i = 0; i < firstNodes.size(); ++i)
        {
            firstNodes[i] = edge->getElement(0)->getPhysicalGeometry()->getNodeIndex(firstNodes[i]);
        }
        for (std::size_t i = 1; i < edge->getNumberOfElements(); ++i)
        {
            logger.assert_always((edge->getElement(i)->getEdge(edge->getEdgeNumber(i)) == edge), "element<->edge matching");
            otherNodes = edge->getElement(i)->getReferenceGeometry()->getCodim2EntityLocalIndices(edge->getEdgeNumber(i));
            for (std::size_t j = 0; j < otherNodes.size(); ++j)
            {
                otherNodes[j] = edge->getElement(i)->getPhysicalGeometry()->getNodeIndex(otherNodes[j]);
            }
            for (std::size_t k = 0; k < firstNodes.size(); ++k)
            {
                bool found = false;
                for (std::size_t j = 0; j < otherNodes.size(); ++j)
                {
                    found |= firstNodes[k] == otherNodes[j];
                }
                logger.assert_always((found), "edge positioning");
            }
            logger.assert_always((firstNodes.size() == otherNodes.size()), "edge positioning");
        }
    }
    for (Base::Node* node : test->getNodesList())
    {
        logger.assert_always((nodeIDs.find(node->getID()) == nodeIDs.end()), "duplicate node ID");
        nodeIDs.insert(node->getID());
        logger.assert_always((node->getElement(0)->getNode(node->getNodeNumber(0)) == node), "element<->node matching");
        Geometry::PointPhysical<DIM> pFirst, pOther;
        pFirst = node->getElement(0)->getPhysicalGeometry()->getLocalNodeCoordinates(node->getNodeNumber(0));
        for (std::size_t i = 1; i < node->getNumberOfElements(); ++i)
        {
            logger.assert_always((node->getElement(i)->getNode(node->getNodeNumber(i)) == node), "element<->node matching");
            pOther = node->getElement(i)->getPhysicalGeometry()->getLocalNodeCoordinates(node->getNodeNumber(i));
            logger.assert_always((pFirst == pOther), "node positioning");
        }
    }
    logger.assert_always((test->getNumberOfNodes() == test->getNumberOfNodeCoordinates()), "total amount of grid points (%) is not equal to the number of nodes (%)",
                         test->getNumberOfNodeCoordinates(), test->getNumberOfNodes());
}


//Also test the copy constructor of Mesh/MeshManipulator here after this test has been fixed.
int main(int argc, char** argv)
{
    Base::parse_options(argc, argv);
    
    
    
    Base::ConfigurationData config2(1);
    
    config2.numberOfUnknowns_       = 1;
    config2.numberOfTimeLevels_     = 1;
    config2.numberOfBasisFunctions_ = 1;
    
    Base::MeshManipulator<2> myTwoDDemoMesh(&config2);
    
    std::stringstream filename;
    
    filename << Base::getCMAKE_hpGEM_SOURCE_DIR() << "/tests/files/square_Test.hyb";
    
    myTwoDDemoMesh.readCentaurMesh(filename.str());
    
    std::cout << myTwoDDemoMesh << std::endl;
    testMesh(&myTwoDDemoMesh);
    
    std::ofstream file2D;
    file2D.open ("Savedsquare_15March16.dat");
    
    
    Output::TecplotDiscontinuousSolutionWriter<2> out1(file2D,"QuadMinimum Test Mesh","01","x, y");
    Dummy2 d2;
    out1.write(&myTwoDDemoMesh,"holi",false, std::function<void(const Base::Element* el, const Geometry::PointReference<2>& p, ostream& os)>{d2});
    
    file2D.close();
    
    
    
    
    
    
    //Next we test it for a 3D mesh generated using Centaur with SOLID WALL boundary conditions in all directions
    Base::ConfigurationData config3(1);
    
    config3.numberOfUnknowns_       = 1;
    config3.numberOfTimeLevels_     = 1;
    config3.numberOfBasisFunctions_ = 1;
    
    Base::MeshManipulator<3> myThreeDDemoMesh(&config3);
    
    std::stringstream filename3;
    
    filename3 << Base::getCMAKE_hpGEM_SOURCE_DIR() << "/tests/files/newCube.hyb";
    
    myThreeDDemoMesh.readCentaurMesh(filename3.str());
    
    std::cout << myThreeDDemoMesh << std::endl;
    
    testMesh(&myThreeDDemoMesh);
    std::ofstream file3D;
    file3D.open ("SavedCubeMesh.dat");
    
    
    Output::TecplotDiscontinuousSolutionWriter<3> out3(file3D,"QuadMinimum Test Mesh","012","x, y, z");
    Dummy3 d3;
    out3.write(&myThreeDDemoMesh,"holi",false, std::function<void(const Base::Element* el, const Geometry::PointReference<3>& p, ostream& os)>{d3});
    
    file3D.close();
    
    
    // Next we test it for 3D mesh generated using Centaur with Periodic Boundary Conditions and containing hexahedrons as well
    
    Base::ConfigurationData config4(1);
    
    config4.numberOfUnknowns_       = 1;
    config4.numberOfTimeLevels_     = 1;
    config4.numberOfBasisFunctions_ = 1;
    
    Base::MeshManipulator<3> myThreeDDemoMesh2(&config4);
    
    std::stringstream filename4;
    
    filename4 << Base::getCMAKE_hpGEM_SOURCE_DIR() << "/tests/files/Bragg_Test.hyb";
    
    myThreeDDemoMesh2.readCentaurMesh(filename4.str());
    
    std::cout << myThreeDDemoMesh2 << std::endl;
    
    testMesh(&myThreeDDemoMesh2);
    std::ofstream file3D2;
    file3D2.open ("SavedCubeMesh2.dat");
    
    
    Output::TecplotDiscontinuousSolutionWriter<3> out4(file3D2,"QuadMinimum Test Mesh","012","x, y, z");
    Dummy3 d4;
    out4.write(&myThreeDDemoMesh2,"holi",false, std::function<void(const Base::Element* el, const Geometry::PointReference<3>& p, ostream& os)>{d4});
    
    file3D2.close();
    

}

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

//This test validates that, for meshes generated by the rectangular mesh generator:
//-all IDs assigned to elements, faces or edges are unique (but there may be an edge, a face and an element that share the same ID)
//-all elements connected to a face or an edge agree on the location of this face/edge and on the identity of the face/edge
//-all faces and all edges of an element are defined
//added all of the above for vertices
//periodic meshes need a different test
#include "Base/MeshManipulator.h"
#include "Base/Edge.h"
#include "Base/Element.h"
#include "Base/Face.h"
#include "Base/ConfigurationData.h"
#include "Base/RectangularMeshDescriptor.h"
#include "Geometry/PhysicalGeometry.h"
#include "Geometry/ReferenceGeometry.h"
#include "Base/ElementCacheData.h"
#include "Base/FaceCacheData.h"
#include "Geometry/PointReference.h"
#include "Base/CommandLineOptions.h"

#include "unordered_set"
#include "Logger.h"
#include "Base/CommandLineOptions.h"

void testMesh(Base::MeshManipulator* test)
{
    std::unordered_set<std::size_t> elementIDs, faceIDs, edgeIDs, vertexIDs;
    std::cout << test->getElementsList(Base::IteratorType::GLOBAL).size() << std::endl;
    std::cout << test->getElementsList(Base::IteratorType::LOCAL).size() << std::endl;
    for (Base::Element* element : test->getElementsList())
    {
        logger.assert_always((elementIDs.find(element->getID()) == elementIDs.end()), "duplicate element ID");
        elementIDs.insert(element->getID());
        logger.assert_always((element->getNrOfFaces() == element->getReferenceGeometry()->getNrOfCodim1Entities()), "confusion about the number of faces");
        if (test->dimension() == 2)
        {
            logger.assert_always((element->getNrOfEdges() == 0), "confusion about the number of edges");
        }
        else
        {
            logger.assert_always((element->getNrOfEdges() == element->getReferenceGeometry()->getNrOfCodim2Entities()), "confusion about the number of edges");
        }
        logger.assert_always((element->getNrOfNodes() == element->getReferenceGeometry()->getNumberOfNodes()), "confusion about the number of vertices");
        for (std::size_t i = 0; i < element->getNrOfFaces(); ++i)
        {
            logger.assert_always((element->getFace(i) != nullptr), "missing Face no. % in element %", i, element->getID());
        }
        for (std::size_t i = 0; i < element->getNrOfEdges(); ++i)
        {
            logger.assert_always((element->getEdge(i) != nullptr), "missing Face");
        }
        for (std::size_t i = 0; i < element->getNrOfNodes(); ++i)
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
                logger.assert_always((found), "face positioning");
            }
            logger.assert_always((leftNodes.size() == rightNodes.size()), "face positioning");
            logger.assert_always((face->getPtrElementRight()->getFace(face->localFaceNumberRight()) == face), "element<->face matching");
        }
    }
    for (Base::Edge* edge : test->getEdgesList())
    {
        logger.assert_always((edgeIDs.find(edge->getID()) == edgeIDs.end()), "duplicate edge ID");
        edgeIDs.insert(edge->getID());
        logger.assert_always((edge->getElement(0)->getEdge(edge->getEdgeNr(0)) == edge), "element<->edge matching");
        std::vector<std::size_t> firstNodes(edge->getElement(0)->getReferenceGeometry()->getCodim2ReferenceGeometry(edge->getEdgeNr(0))->getNumberOfNodes()), otherNodes(firstNodes);
        firstNodes = edge->getElement(0)->getReferenceGeometry()->getCodim2EntityLocalIndices(edge->getEdgeNr(0));
        for (std::size_t i = 0; i < firstNodes.size(); ++i)
        {
            firstNodes[i] = edge->getElement(0)->getPhysicalGeometry()->getNodeIndex(firstNodes[i]);
        }
        for (std::size_t i = 1; i < edge->getNrOfElements(); ++i)
        {
            logger.assert_always((edge->getElement(i)->getEdge(edge->getEdgeNr(i)) == edge), "element<->edge matching");
            otherNodes = edge->getElement(i)->getReferenceGeometry()->getCodim2EntityLocalIndices(edge->getEdgeNr(i));
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
    for (Base::Node* vertex : test->getNodesList())
    {
        logger.assert_always((vertexIDs.find(vertex->getID()) == vertexIDs.end()), "duplicate vertex ID");
        vertexIDs.insert(vertex->getID());
        logger.assert_always((vertex->getElement(0)->getNode(vertex->getNodeNr(0)) == vertex), "element<->vertex matching");
        Geometry::PointPhysical pFirst(test->dimension()), pOther(pFirst);
        pFirst = vertex->getElement(0)->getPhysicalGeometry()->getLocalNodeCoordinates(vertex->getNodeNr(0));
        for (std::size_t i = 1; i < vertex->getNrOfElements(); ++i)
        {
            logger.assert_always((vertex->getElement(i)->getNode(vertex->getNodeNr(i)) == vertex), "element<->vertex matching");
            pOther = vertex->getElement(i)->getPhysicalGeometry()->getLocalNodeCoordinates(vertex->getNodeNr(i));
            logger.assert_always((pFirst == pOther), "vertex positioning");
        }
    }
    logger.assert_always((test->getNumberOfNodes() == test->getNumberOfNodeCoordinates()), "total amount of grid points (%) is not equal to the number of nodes (%)",
                         test->getNumberOfNodeCoordinates(), test->getNumberOfNodes());
}

void testPointPhysicalsOfElementsOfCopiedMesh(Base::MeshManipulator* mesh, Base::MeshManipulator* meshCopy)
{
    std::vector<Base::Element*> eltsMesh = mesh->getElementsList();
    std::vector<Base::Element*> eltsMeshCopy = meshCopy->getElementsList();
    logger.assert_always(eltsMesh.size() == eltsMeshCopy.size(), "The copy does not have the same number of elements as the original MeshManipulator.");
    for (std::size_t i = 0; i < eltsMesh.size(); ++i)
    {
        logger.assert_always(eltsMesh[i]->getPhysicalGeometry()->getNodes() == eltsMeshCopy[i]->getPhysicalGeometry()->getNodes(), "The points of Element % are different.", i);        
    }
}

int main(int argc, char** argv)
{
    Base::parse_options(argc, argv);
    // dim 1
    Base::RectangularMeshDescriptor description1D(1), description2D(2), description3D(3);
    description1D.bottomLeft_[0] = 0;
    description2D.bottomLeft_[0] = 0;
    description2D.bottomLeft_[1] = 0;
    description3D.bottomLeft_[0] = 0;
    description3D.bottomLeft_[1] = 0;
    description3D.bottomLeft_[2] = 0;
    description1D.topRight_[0] = 1;
    description2D.topRight_[0] = 1;
    description2D.topRight_[1] = 1;
    description3D.topRight_[0] = 1;
    description3D.topRight_[1] = 1;
    description3D.topRight_[2] = 1;
    description1D.boundaryConditions_[0] = Base::BoundaryType::SOLID_WALL;
    description2D.boundaryConditions_[0] = Base::BoundaryType::SOLID_WALL;
    description2D.boundaryConditions_[1] = Base::BoundaryType::SOLID_WALL;
    description3D.boundaryConditions_[0] = Base::BoundaryType::SOLID_WALL;
    description3D.boundaryConditions_[1] = Base::BoundaryType::SOLID_WALL;
    description3D.boundaryConditions_[2] = Base::BoundaryType::SOLID_WALL;
    
    description1D.numElementsInDIM_[0] = 2;
    
    Base::MeshManipulator *test = new Base::MeshManipulator(new Base::ConfigurationData(1, 1, 2, 0), Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, 2, 0);
    test->createRectangularMesh(description1D.bottomLeft_, description1D.topRight_, description1D.numElementsInDIM_);
    
    testMesh(test);
    logger.assert_always((test->getNumberOfElements() == 2), "number of elements");
    
    
    delete test;
    description1D.numElementsInDIM_[0] = 3;
    
    test = new Base::MeshManipulator(new Base::ConfigurationData(1, 1, 2, 0), Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, 2, 0);
    test->createRectangularMesh(description1D.bottomLeft_, description1D.topRight_, description1D.numElementsInDIM_);    
    
    testMesh(test);
    logger.assert_always((test->getNumberOfElements() == 3), "number of elements");
    
    // dim 2
    
    delete test;
    description2D.numElementsInDIM_[0] = 2;
    description2D.numElementsInDIM_[1] = 3;
    
    test = new Base::MeshManipulator(new Base::ConfigurationData(2, 1, 2, 0), Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, 2, 0);
    test->createRectangularMesh(description2D.bottomLeft_, description2D.topRight_, description2D.numElementsInDIM_);
    
    testMesh(test);
    logger.assert_always((test->getNumberOfElements() == 6), "number of elements");
    
    delete test;
    description2D.numElementsInDIM_[0] = 3;
    description2D.numElementsInDIM_[1] = 2;
    
    test = new Base::MeshManipulator(new Base::ConfigurationData(2, 1, 2, 0), Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, 2, 0);
    test->createRectangularMesh(description2D.bottomLeft_, description2D.topRight_, description2D.numElementsInDIM_);
    
    testMesh(test);
    logger.assert_always((test->getNumberOfElements() == 6), "number of elements");
    
    // dim 3
    
    delete test;
    description3D.numElementsInDIM_[0] = 2;
    description3D.numElementsInDIM_[1] = 2;
    description3D.numElementsInDIM_[2] = 3;
    
    test = new Base::MeshManipulator(new Base::ConfigurationData(3, 1, 2, 0), Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, 2, 0);
    test->createRectangularMesh(description3D.bottomLeft_, description3D.topRight_, description3D.numElementsInDIM_);
    
    testMesh(test);
    logger.assert_always((test->getNumberOfElements() == 12), "number of elements");
    
    delete test;
    description3D.numElementsInDIM_[0] = 2;
    description3D.numElementsInDIM_[1] = 3;
    description3D.numElementsInDIM_[2] = 2;
    
    test = new Base::MeshManipulator(new Base::ConfigurationData(3, 1, 2, 0), Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, 2, 0);
    test->createRectangularMesh(description3D.bottomLeft_, description3D.topRight_, description3D.numElementsInDIM_);
    
    testMesh(test);
    logger.assert_always((test->getNumberOfElements() == 12), "number of elements");
    
    delete test;
    description3D.numElementsInDIM_[0] = 3;
    description3D.numElementsInDIM_[1] = 2;
    description3D.numElementsInDIM_[2] = 2;
    
    test = new Base::MeshManipulator(new Base::ConfigurationData(3, 1, 2, 0), Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, 2, 0);
    test->createRectangularMesh(description3D.bottomLeft_, description3D.topRight_, description3D.numElementsInDIM_);
    
    testMesh(test);
    logger.assert_always((test->getNumberOfElements() == 12), "number of elements");    
    
    //test copy constructor of MeshManipulator, only most difficult case
    Base::MeshManipulator* test2 = new Base::MeshManipulator(*test);
    testMesh(test2);    
    testPointPhysicalsOfElementsOfCopiedMesh(test, test2);
    
    delete test2;
    delete test;
    
    return 0;
}

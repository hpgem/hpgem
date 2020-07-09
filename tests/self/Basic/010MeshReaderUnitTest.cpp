/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2017, University of Twente
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors
 may be used to endorse or promote products derived from this software without
 specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2014, University of Twente
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors
 may be used to endorse or promote products derived from this software without
 specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
using namespace hpgem;
template <std::size_t DIM>
void testMesh(Base::MeshManipulator<DIM>* test, bool isPeriodic) {
    std::unordered_set<std::size_t> elementIDs, faceIDs, edgeIDs, nodeIDs;
    logger.suppressWarnings([&]() {
        std::cout << test->getElementsList(Base::IteratorType::GLOBAL).size()
                  << std::endl;
    });
    std::cout << test->getElementsList(Base::IteratorType::LOCAL).size()
              << std::endl;
    for (Base::Element* element : test->getElementsList()) {
        logger.assert_always(
            (elementIDs.find(element->getID()) == elementIDs.end()),
            "duplicate element ID");
        elementIDs.insert(element->getID());
        logger.assert_always(
            (element->getNumberOfFaces() ==
             element->getReferenceGeometry()->getNumberOfCodim1Entities()),
            "confusion about the number of faces");
        if (test->dimension() == 2) {
            logger.assert_always((element->getNumberOfEdges() == 0),
                                 "confusion about the number of edges");
        } else {
            logger.assert_always(
                (element->getNumberOfEdges() ==
                 element->getReferenceGeometry()->getNumberOfCodim2Entities()),
                "confusion about the number of edges");
        }
        logger.assert_always(
            (element->getNumberOfNodes() ==
             element->getReferenceGeometry()->getNumberOfNodes()),
            "confusion about the number of vertices");
        for (std::size_t i = 0; i < element->getNumberOfFaces(); ++i) {
            logger.assert_always((element->getFace(i) != nullptr),
                                 "missing Face no. % in element %", i,
                                 element->getID());
        }
        for (std::size_t i = 0; i < element->getNumberOfEdges(); ++i) {
            logger.assert_always((element->getEdge(i) != nullptr),
                                 "missing Face");
        }
        for (std::size_t i = 0; i < element->getNumberOfNodes(); ++i) {
            logger.assert_always((element->getNode(i) != nullptr),
                                 "missing Face");
        }
    }
    for (Base::Face* face : test->getFacesList()) {
        logger.assert_always((faceIDs.find(face->getID()) == faceIDs.end()),
                             "duplicate face ID");
        faceIDs.insert(face->getID());
        logger.assert_always((face->getPtrElementLeft()->getFace(
                                  face->localFaceNumberLeft()) == face),
                             "element<->face matching");
        if (!isPeriodic && face->isInternal()) {
            std::vector<std::size_t> leftNodes(
                face->getReferenceGeometry()->getNumberOfNodes()),
                rightNodes(leftNodes);
            leftNodes =
                face->getPtrElementLeft()
                    ->getPhysicalGeometry()
                    ->getGlobalFaceNodeIndices(face->localFaceNumberLeft());
            rightNodes =
                face->getPtrElementRight()
                    ->getPhysicalGeometry()
                    ->getGlobalFaceNodeIndices(face->localFaceNumberRight());
            logger.assert_always(
                (leftNodes.size() == rightNodes.size()),
                "face positioning: inconsistent face size, % vs. %",
                leftNodes.size(), rightNodes.size());
            for (unsigned long leftNode : leftNodes) {
                std::cout << leftNode << " ";
            }
            std::cout << endl;
            for (unsigned long rightNode : rightNodes) {
                std::cout << rightNode << " ";
            }
            std::cout << endl;
            for (unsigned long leftNode : leftNodes) {
                bool found = false;
                for (unsigned long rightNode : rightNodes) {
                    found |= leftNode == rightNode;
                }
                logger.assert_always((found), "face positioning");
            }
            logger.assert_always((face->getPtrElementRight()->getFace(
                                      face->localFaceNumberRight()) == face),
                                 "element<->face matching");
        }
        if (isPeriodic) {
            logger.assert_always(face->isInternal(),
                                 "Boundary face detected in a periodic mesh");
        }
    }
    for (Base::Edge* edge : test->getEdgesList()) {
        logger.assert_always((edgeIDs.find(edge->getID()) == edgeIDs.end()),
                             "duplicate edge ID");
        edgeIDs.insert(edge->getID());
        logger.assert_always(
            (edge->getElement(0)->getEdge(edge->getEdgeNumber(0)) == edge),
            "element<->edge matching");
        std::vector<std::size_t> firstNodes(
            edge->getElement(0)
                ->getReferenceGeometry()
                ->getCodim2ReferenceGeometry(edge->getEdgeNumber(0))
                ->getNumberOfNodes()),
            otherNodes(firstNodes);
        firstNodes = edge->getElement(0)
                         ->getReferenceGeometry()
                         ->getCodim2EntityLocalIndices(edge->getEdgeNumber(0));
        for (unsigned long& firstNode : firstNodes) {
            firstNode =
                edge->getElement(0)->getPhysicalGeometry()->getNodeIndex(
                    firstNode);
        }
        for (std::size_t i = 1; i < edge->getNumberOfElements(); ++i) {
            logger.assert_always(
                (edge->getElement(i)->getEdge(edge->getEdgeNumber(i)) == edge),
                "element<->edge matching");
            otherNodes =
                edge->getElement(i)
                    ->getReferenceGeometry()
                    ->getCodim2EntityLocalIndices(edge->getEdgeNumber(i));
            for (unsigned long& otherNode : otherNodes) {
                otherNode =
                    edge->getElement(i)->getPhysicalGeometry()->getNodeIndex(
                        otherNode);
            }
            for (unsigned long firstNode : firstNodes) {
                bool found = false;
                for (unsigned long otherNode : otherNodes) {
                    found |= firstNode == otherNode;
                }
                logger.assert_always(found || isPeriodic, "edge positioning");
            }
            logger.assert_always((firstNodes.size() == otherNodes.size()),
                                 "edge positioning");
        }
    }
    for (Base::Node* node : test->getNodesList()) {
        logger.assert_always((nodeIDs.find(node->getID()) == nodeIDs.end()),
                             "duplicate node ID");
        nodeIDs.insert(node->getID());
        logger.assert_always(
            (node->getElement(0)->getNode(node->getNodeNumber(0)) == node),
            "element<->node matching");
        Geometry::PointPhysical<DIM> pFirst, pOther;
        pFirst =
            node->getElement(0)->getPhysicalGeometry()->getLocalNodeCoordinates(
                node->getNodeNumber(0));
        for (std::size_t i = 1; i < node->getNumberOfElements(); ++i) {
            logger.assert_always(
                (node->getElement(i)->getNode(node->getNodeNumber(i)) == node),
                "element<->node matching");
            pOther = node->getElement(i)
                         ->getPhysicalGeometry()
                         ->getLocalNodeCoordinates(node->getNodeNumber(i));
            logger.assert_always((pFirst == pOther) || isPeriodic,
                                 "node positioning");
        }
    }
    logger.assert_always(
        (test->getNumberOfNodes() == test->getNumberOfNodeCoordinates()) ||
            isPeriodic,
        "total amount of grid points (%) is not equal to the number of nodes "
        "(%)",
        test->getNumberOfNodeCoordinates(), test->getNumberOfNodes());
}

// Also test the copy constructor of Mesh/MeshManipulator here after this test
// has been fixed.
int main(int argc, char** argv) {
    using namespace std::string_literals;
    Base::parse_options(argc, argv);

    // 1d

    for (const std::string& baseFileName :
         {"1Drectangular1mesh"s, "1Drectangular2mesh"s}) {
        Base::ConfigurationData config(1);

        config.numberOfUnknowns_ = 1;
        config.numberOfTimeLevels_ = 1;
        config.numberOfBasisFunctions_ = 1;

        Base::MeshManipulator<1> mesh(&config);

        std::stringstream filename;

        filename << Base::getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/"s
                 << baseFileName << ".hpgem";

        mesh.readMesh(filename.str());

        std::cout << mesh << std::endl;
        testMesh(&mesh, false);

        std::ofstream outFile;
        outFile.open(baseFileName + "_output.dat"s);

        Output::TecplotDiscontinuousSolutionWriter<1> out(
            outFile, "1D Test Mesh", "0", "x");
        out.write(&mesh, "holi", false,
                  [](const Base::Element*, const Geometry::PointReference<1>&,
                     ostream&) {});

        outFile.close();
    }

    // 2d

    for (const std::string& baseFileName :
         {"2Drectangular1mesh"s, "2Drectangular2mesh"s, "2Dtriangular1mesh"s,
          "2Dtriangular2mesh"s, "square_Test"s}) {
        Base::ConfigurationData config(1);

        config.numberOfUnknowns_ = 1;
        config.numberOfTimeLevels_ = 1;
        config.numberOfBasisFunctions_ = 1;

        Base::MeshManipulator<2> mesh(&config);

        std::stringstream filename;

        filename << Base::getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/"s
                 << baseFileName << ".hpgem";

        mesh.readMesh(filename.str());

        std::cout << mesh << std::endl;
        testMesh(&mesh, false);

        std::ofstream outFile;
        outFile.open(baseFileName + "_output.dat"s);

        Output::TecplotDiscontinuousSolutionWriter<2> out(
            outFile, "2D Test Mesh", "01", "x, y");
        out.write(&mesh, "holi", false,
                  [](const Base::Element*, const Geometry::PointReference<2>&,
                     ostream&) {});

        outFile.close();
    }

    // 3d

    for (const std::string& baseFileName :
         {"3Drectangular1mesh"s, "3Drectangular2mesh"s, "3Drectangular3mesh"s,
          "3Dtriangular1mesh"s, "3Dtriangular2mesh"s, "3Dtriangular3mesh"s}) {
        Base::ConfigurationData config(1);

        config.numberOfUnknowns_ = 1;
        config.numberOfTimeLevels_ = 1;
        config.numberOfBasisFunctions_ = 1;

        Base::MeshManipulator<3> mesh(&config);

        std::stringstream filename;

        filename << Base::getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/"s
                 << baseFileName << ".hpgem";

        mesh.readMesh(filename.str());

        std::cout << mesh << std::endl;
        testMesh(&mesh, false);

        std::ofstream outFile;
        outFile.open(baseFileName + "_output.dat"s);

        Output::TecplotDiscontinuousSolutionWriter<3> out(
            outFile, "3D Test Mesh", "012", "x, y, z");
        out.write(&mesh, "holi", false,
                  [](const Base::Element*, const Geometry::PointReference<3>&,
                     ostream&) {});

        outFile.close();
    }

    for (const std::string& baseFileName : {"Bragg_Test"s, "newCube"s}) {
        Base::ConfigurationData config(1);

        config.numberOfUnknowns_ = 1;
        config.numberOfTimeLevels_ = 1;
        config.numberOfBasisFunctions_ = 1;

        Base::MeshManipulator<3> mesh(&config);

        std::stringstream filename;

        filename << Base::getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/"s
                 << baseFileName << ".hpgem";

        mesh.readMesh(filename.str());

        std::cout << mesh << std::endl;
        testMesh(&mesh, true);

        std::ofstream outFile;
        outFile.open(baseFileName + "_output.dat"s);

        Output::TecplotDiscontinuousSolutionWriter<3> out(
            outFile, "3D Test Mesh", "012", "x, y, z");
        out.write(&mesh, "holi", false,
                  [](const Base::Element*, const Geometry::PointReference<3>&,
                     ostream&) {});

        outFile.close();
    }
}

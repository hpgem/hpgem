/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2020, University of Twente
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

#include "CMakeDefinitions.h"

#include "Base/MeshManipulator.h"
#include "Base/CommandLineOptions.h"

#include "Utilities/GlobalIndexing.h"
#include "Utilities/SparsityEstimator.h"
using namespace hpgem;
using Layout = Utilities::GlobalIndexing::Layout;

/// Test configuration with GlobalIndexing layout and whether to include face
/// coupling.
struct TestConfiguration {
    TestConfiguration(Layout layout, bool faceCoupling)
        : layout_(layout), faceCoupling_(faceCoupling) {}

    Layout layout_;
    bool faceCoupling_;
};

Layout TEST_LAYOUTS[] = {Layout::SEQUENTIAL, Layout::BLOCKED_PROCESSOR,
                         Layout::BLOCKED_GLOBAL};

TestConfiguration TEST_CONFIGURATIONS[] = {
    TestConfiguration(Layout::SEQUENTIAL, false),
    TestConfiguration(Layout::SEQUENTIAL, true),
    TestConfiguration(Layout::BLOCKED_PROCESSOR, false),
    TestConfiguration(Layout::BLOCKED_PROCESSOR, true),
    TestConfiguration(Layout::BLOCKED_GLOBAL, false),
    TestConfiguration(Layout::BLOCKED_GLOBAL, true),
};

enum class BasisFunctionType {
    CG4,  // Fourth order CG, which in 3D has DoFs with support on every piece
    // of geometry
    DG1  // First order DG.
};

template <std::size_t DIM>
void applyBasisFunctionsType(
    Base::MeshManipulator<DIM>& mesh,
    const std::vector<BasisFunctionType>& basisFunctions) {

    if (!basisFunctions.empty()) {
        switch (basisFunctions[0]) {
            case BasisFunctionType::CG4:
                mesh.useDefaultConformingBasisFunctions(4);
                break;
            case BasisFunctionType::DG1:
                mesh.useDefaultDGBasisFunctions(1);
                break;
            default:
                logger.assert_always(false, "Unknown basis function type");
        }
        for (std::size_t i = 1; i < basisFunctions.size(); ++i) {
            switch (basisFunctions[i]) {
                case BasisFunctionType::CG4:
                    mesh.useDefaultConformingBasisFunctions(4, i);
                    break;
                case BasisFunctionType::DG1:
                    mesh.useDefaultDGBasisFunctions(1, i);
                    break;
                default:
                    logger.assert_always(false, "Unknown basis function type");
            }
        }
    }
}

/// Check that the estimate holds for the DoFs related to part of the geometry
/// \tparam GEOM The type of geometry part (Element, Face, etc.)
/// \param geom The geometry part
/// \param indexing The index to find the local indices for the geometry part
/// \param owned The sparsity estimate (number of non zero DoFs in each row) for
/// the owned DoFs. \param nonOwned The sparsity estimate (number of non zero
/// DoFs in each row) for the non owned DoFs. \param expectedOwned The expected
/// sparsity estimate for the owned DoFs \param expectedNonOnwed The expected
/// sparsity estimate for the non owned DoFs.
template <typename GEOM>
void check(GEOM geom, const Utilities::GlobalIndexing& indexing,
           std::vector<int>& owned, std::vector<int>& nonOwned,
           std::size_t expectedOwned, std::size_t expectedNonOnwed) {
    for (std::size_t unknown : indexing.getIncludedUnknowns()) {
        // Check for each of the basis functions
        std::size_t dof = indexing.getProcessorLocalIndex(geom, unknown);
        for (std::size_t dofOff = 0;
             dofOff < geom->getLocalNumberOfBasisFunctions(unknown); ++dofOff) {
            logger.assert_always(owned[dof + dofOff] == expectedOwned,
                                 "Expected % owned non zero entries but got %",
                                 expectedOwned, owned[dof + dofOff]);
            logger.assert_always(
                nonOwned[dof + dofOff] == expectedNonOnwed,
                "Expected % non owned non zero entries but got %",
                expectedNonOnwed, nonOwned[dof + dofOff]);
        }
    }
}

/// Select an object based on whether the geometrical object is owned or not
/// \tparam GEOM The type of geometrical object (Element, Face, etc.)
/// \tparam T The type of object to be selected
/// \param geom The geometrical object that is either owned or not
/// \param owned The value to return when owned
/// \param nonOwned The value to return when not owned.
/// \return Either owned or nonOwned
template <typename GEOM, typename T>
T& selectByOwner(const GEOM* geom, T& owned, T& nonOwned) {
    if (geom->isOwnedByCurrentProcessor()) return owned;

    return nonOwned;
}

/// Storage for geometrical objects.
struct GeomStorage {
    std::set<const Base::Element*> elements;
    std::set<const Base::Face*> faces;
    std::set<const Base::Edge*> edges;
    std::set<const Base::Node*> nodes;

    /// \brief Add the Element and all its faces, edges and nodes
    ///
    /// \param element The element to add
    /// \param neighbours Whether to include the neighbouring elements and their
    /// faces, edges and nodes
    void addGeometryAround(const Base::Element* element,
                           bool neighbours = false) {
        elements.insert(element);
        for (auto* face : element->getFacesList()) {
            faces.insert(face);
            if (face->isInternal() && neighbours) {
                addGeometryAround(face->getPtrOtherElement(element), false);
            }
        }
        for (auto* edge : element->getEdgesList()) edges.insert(edge);
        for (auto* node : element->getNodesList()) nodes.insert(node);
    }

    void addGeometryAround(const Base::Face* face,
                           bool secondNeighbours = false) {
        addGeometryAround(face->getPtrElementLeft(), secondNeighbours);
        if (face->isInternal())
            addGeometryAround(face->getPtrElementRight(), secondNeighbours);
    }
    void addGeometryAround(const Base::Edge* edge,
                           bool secondNeighbours = false) {
        for (const Base::Element* element : edge->getElements())
            addGeometryAround(element, secondNeighbours);
    }
    void addGeometryAround(const Base::Node* node,
                           bool secondNeighbours = false) {
        for (const Base::Element* element : node->getElements())
            addGeometryAround(element, secondNeighbours);
    }

    std::size_t countBasisFunctions(
        std::vector<std::size_t>* unknowns = nullptr) {
        // Assumes everything is locally owned
        if (unknowns == nullptr) {
            std::size_t basisFunctions = 0;
            for (const Base::Element* element : elements)
                basisFunctions +=
                    element->getTotalLocalNumberOfBasisFunctions();
            for (const Base::Face* face : faces)
                basisFunctions += face->getTotalLocalNumberOfBasisFunctions();
            for (const Base::Edge* edge : edges)
                basisFunctions += edge->getTotalLocalNumberOfBasisFunctions();
            for (const Base::Node* node : nodes)
                basisFunctions += node->getTotalLocalNumberOfBasisFunctions();
            return basisFunctions;
        }
        std::size_t basisFunctions = 0;
        for (std::size_t unknown : *unknowns) {
            for (const Base::Element* element : elements)
                basisFunctions +=
                    element->getLocalNumberOfBasisFunctions(unknown);
            for (const Base::Face* face : faces)
                basisFunctions += face->getLocalNumberOfBasisFunctions(unknown);
            for (const Base::Edge* edge : edges)
                basisFunctions += edge->getLocalNumberOfBasisFunctions(unknown);
            for (const Base::Node* node : nodes)
                basisFunctions += node->getLocalNumberOfBasisFunctions(unknown);
        }
        return basisFunctions;
    }
};

/// Test with DG basis, which has the (dis)advantage that all basis functions
/// are confined to an element
void testWithDGBasis(std::size_t unknowns, std::string meshFile) {
    Base::ConfigurationData config(unknowns);
    Base::MeshManipulator<3> mesh(&config);

    using namespace std::string_literals;
    std::stringstream filename;
    filename << Base::getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/"s << meshFile
             << ".hpgem";
    mesh.readMesh(filename.str());

    mesh.useDefaultDGBasisFunctions(0);
    for (std::size_t unknown = 1; unknown < unknowns; ++unknown) {
        // Use more than zeroth order DG
        mesh.useDefaultDGBasisFunctions(unknown, unknown);
    }

    for (TestConfiguration configuration : TEST_CONFIGURATIONS) {
        Utilities::GlobalIndexing indexing(&mesh, configuration.layout_);
        Utilities::SparsityEstimator estimator(indexing);

        std::vector<int> owned, nonOwned;
        estimator.computeSparsityEstimate(owned, nonOwned,
                                          configuration.faceCoupling_);
        logger.assert_always(
            owned.size() == indexing.getNumberOfLocalBasisFunctions(),
            "Wrong size owned");
        logger.assert_always(
            nonOwned.size() == indexing.getNumberOfLocalBasisFunctions(),
            "Wrong size owned");

        for (const Base::Element* element : mesh.getElementsList()) {
            std::size_t numberOfOwnDoFs =
                element->getTotalLocalNumberOfBasisFunctions();
            std::size_t numberOfNonOwnDoFs = 0;
            // Skip face loop if there is no face coupling
            if (configuration.faceCoupling_) {
                for (const Base::Face* face : element->getFacesList()) {
                    if (face->getFaceType() ==
                            Geometry::FaceType::SUBDOMAIN_BOUNDARY ||
                        face->getFaceType() ==
                            Geometry::FaceType::PERIODIC_SUBDOMAIN_BC) {
                        numberOfNonOwnDoFs +=
                            face->getPtrOtherElement(element)
                                ->getTotalLocalNumberOfBasisFunctions();
                    } else if (face->isInternal()) {
                        numberOfOwnDoFs +=
                            face->getPtrOtherElement(element)
                                ->getTotalLocalNumberOfBasisFunctions();
                    }
                }
            }
            check(element, indexing, owned, nonOwned, numberOfOwnDoFs,
                  numberOfNonOwnDoFs);
        }
    }
}

/**
 * Test with conforming basis functions on a 1D mesh. Here the sparisty pattern
 * is far more predictable than in a higher dimensions.
 */
void testConformingWith1DMesh() {
    Base::ConfigurationData config(1);
    Base::MeshManipulator<1> mesh(&config);

    using namespace std::string_literals;
    std::stringstream filename;
    filename << Base::getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/"s
             << "1Drectangular2mesh"s
             << ".hpgem";
    mesh.readMesh(filename.str());

    mesh.useDefaultConformingBasisFunctions(3);

    for (TestConfiguration configuration : TEST_CONFIGURATIONS) {
        Utilities::GlobalIndexing indexing(&mesh, configuration.layout_);
        Utilities::SparsityEstimator estimator(indexing);

        std::vector<int> owned, nonOwned;
        estimator.computeSparsityEstimate(owned, nonOwned,
                                          configuration.faceCoupling_);
        logger.assert_always(
            owned.size() == indexing.getNumberOfLocalBasisFunctions(),
            "Wrong size owned");
        logger.assert_always(
            nonOwned.size() == indexing.getNumberOfLocalBasisFunctions(),
            "Wrong size owned");

        // Check for Element based DoF
        for (const Base::Element* element : mesh.getElementsList()) {
            // The following DoFs are contributing
            // The DoFs from the element and its neighbours
            // The DoFs from the faces to element
            // The DoFs from the faces of the neighbouring elements
            std::size_t numberOfOwnDoFs =
                element->getTotalLocalNumberOfBasisFunctions();
            std::size_t numberOfNonOwnDoFs = 0;
            for (const Base::Face* face : element->getFacesList()) {
                selectByOwner(face, numberOfOwnDoFs, numberOfNonOwnDoFs) +=
                    face->getTotalLocalNumberOfBasisFunctions();
                // If there is no face at the other side or we do not couple
                // through the face matrix skip the computation of the DoFs on
                // the other side
                if (!face->isInternal() || !configuration.faceCoupling_)
                    continue;
                // The neighbouring element
                const Base::Element* otherElement =
                    face->getPtrOtherElement(element);
                selectByOwner(otherElement, numberOfOwnDoFs,
                              numberOfNonOwnDoFs) +=
                    element->getTotalLocalNumberOfBasisFunctions();
                // Other face of the neighbour
                const Base::Face* otherFace = otherElement->getFace(
                    otherElement->getFace(0) == face ? 1 : 0);
                selectByOwner(otherFace, numberOfOwnDoFs, numberOfNonOwnDoFs) +=
                    otherFace->getTotalLocalNumberOfBasisFunctions();
            }
            check(element, indexing, owned, nonOwned, numberOfOwnDoFs,
                  numberOfNonOwnDoFs);
        }
        // Check for Face based DoFs
        for (const Base::Face* face : mesh.getFacesList()) {
            // The contributions are from
            // The DoFs on the face
            // The DoFs from the element adjacent to the face it its faces
            // The DoFs from the element adjecent to those and their faces.
            std::vector<const Base::Element*> elements(
                {face->getPtrElementLeft()});
            if (face->isInternal())
                elements.emplace_back(face->getPtrElementRight());
            std::size_t numberOfOwnDoFs =
                face->getTotalLocalNumberOfBasisFunctions();
            std::size_t numberOfNonOwnDoFs = 0;
            for (const Base::Element* element : elements) {
                selectByOwner(element, numberOfOwnDoFs, numberOfNonOwnDoFs) +=
                    element->getTotalLocalNumberOfBasisFunctions();
                const Base::Face* otherFace =
                    element->getFace(element->getFace(0) == face ? 1 : 0);
                selectByOwner(otherFace, numberOfOwnDoFs, numberOfNonOwnDoFs) +=
                    otherFace->getTotalLocalNumberOfBasisFunctions();
                // Next element
                if (!otherFace->isInternal() || !configuration.faceCoupling_)
                    continue;
                const Base::Element* nextElement =
                    otherFace->getPtrOtherElement(element);
                selectByOwner(nextElement, numberOfOwnDoFs,
                              numberOfNonOwnDoFs) +=
                    nextElement->getTotalLocalNumberOfBasisFunctions();
                const Base::Face* nextFace = nextElement->getFace(
                    nextElement->getFace(0) == otherFace ? 1 : 0);
                selectByOwner(nextFace, numberOfOwnDoFs, numberOfNonOwnDoFs) +=
                    nextFace->getTotalLocalNumberOfBasisFunctions();
            }
            check(face, indexing, owned, nonOwned, numberOfOwnDoFs,
                  numberOfNonOwnDoFs);
        }
    }
}

/**
 * Test the sparsity pattern of a standard CG mass/stiffness matrix, thus
 * without any face coupling.
 * @param basisFunctions The type of basis functions to use
 */
void testMassOnly(std::vector<BasisFunctionType> basisFunctions) {
    Base::ConfigurationData config(basisFunctions.size());
    Base::MeshManipulator<3> mesh(&config);

    using namespace std::string_literals;
    std::stringstream filename;
    filename << Base::getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/"s
             << "3Drectangular2mesh"s
             << ".hpgem";
    mesh.readMesh(filename.str());

    applyBasisFunctionsType(mesh, basisFunctions);

    // Sparsity estimate
    Utilities::GlobalIndexing indexing(&mesh, Layout::SEQUENTIAL);
    Utilities::SparsityEstimator estimator(indexing);

    std::vector<int> owned, nonOwned;
    // No face-face coupling
    estimator.computeSparsityEstimate(owned, nonOwned, false);
    logger.assert_always(
        owned.size() == indexing.getNumberOfLocalBasisFunctions(),
        "Wrong size owned");
    logger.assert_always(
        nonOwned.size() == indexing.getNumberOfLocalBasisFunctions(),
        "Wrong size owned");

    for (const Base::Element* element : mesh.getElementsList()) {
        GeomStorage storage;
        storage.addGeometryAround(element);
        std::size_t numDofs = storage.countBasisFunctions();
        check(element, indexing, owned, nonOwned, numDofs, 0);
    }
    for (const Base::Face* face : mesh.getFacesList()) {
        GeomStorage storage;
        storage.addGeometryAround(face);
        std::size_t numDofs = storage.countBasisFunctions();
        check(face, indexing, owned, nonOwned, numDofs, 0);
    }
    for (const Base::Edge* edge : mesh.getEdgesList()) {
        GeomStorage storage;
        storage.addGeometryAround(edge);
        std::size_t numDofs = storage.countBasisFunctions();
        check(edge, indexing, owned, nonOwned, numDofs, 0);
    }
    for (const Base::Node* node : mesh.getNodesList()) {
        GeomStorage storage;
        storage.addGeometryAround(node);
        std::size_t numDofs = storage.countBasisFunctions();
        check(node, indexing, owned, nonOwned, numDofs, 0);
    }
}

/**
 * Test sparisty estimator with for a matrix where the rows and columns use a
 * different subset of the unknowns.
 * @param meshFile The mesh to test on.
 */
void testRowColumnDifference(std::string meshFile) {
    // Test where the GlobalIndexing for the rows and columns differ.
    Base::ConfigurationData config(2);
    Base::MeshManipulator<3> mesh(&config);

    using namespace std::string_literals;
    std::stringstream filename;
    filename << Base::getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/"s << meshFile
             << ".hpgem";
    mesh.readMesh(filename.str());

    // Note different types of basis functions for both DoFs.
    mesh.useDefaultDGBasisFunctions(2);
    mesh.useDefaultConformingBasisFunctions(4, 1);

    for (Layout layout : TEST_LAYOUTS) {
        std::vector<std::size_t> dgUnknown({0}), cgUnknown({1});
        Utilities::GlobalIndexing indexing0(&mesh, layout, &dgUnknown);
        Utilities::GlobalIndexing indexing1(&mesh, layout, &cgUnknown);
        Utilities::GlobalIndexing indexing01(&mesh, layout);
        {
            // Rows are DG, columns are CG, no face coupling.
            // The row for each element consists of the number of CG basis
            // functions with support on the element.

            Utilities::SparsityEstimator estimator(indexing0, indexing1);

            std::vector<int> owned, nonOwned;
            estimator.computeSparsityEstimate(owned, nonOwned, false);
            logger.assert_always(
                owned.size() == indexing0.getNumberOfLocalBasisFunctions(),
                "Wrong size owned");
            logger.assert_always(
                nonOwned.size() == indexing0.getNumberOfLocalBasisFunctions(),
                "Wrong size owned");
            for (const Base::Element* element : mesh.getElementsList()) {
                check(element, indexing0, owned, nonOwned,
                      element->getNumberOfBasisFunctions(1), 0);
            }
        }
        {
            // Rows are CG, Columns are DG, no face coupling.
            // The row for each CG basis functions should contain the number of
            // DG basis functions on the elements on which the CG function has
            // support.
            Utilities::SparsityEstimator estimator(indexing1, indexing0);
            std::vector<int> owned, nonOwned;
            estimator.computeSparsityEstimate(owned, nonOwned, false);
            logger.assert_always(
                owned.size() == indexing1.getNumberOfLocalBasisFunctions(),
                "Wrong size owned");
            logger.assert_always(
                nonOwned.size() == indexing1.getNumberOfLocalBasisFunctions(),
                "Wrong size owned");

            std::size_t dgDoFsperElement =
                (*mesh.elementColBegin())->getNumberOfBasisFunctions(0);
            for (const Base::Element* element : mesh.getElementsList()) {
                check(element, indexing1, owned, nonOwned, dgDoFsperElement, 0);
            }
            for (const Base::Face* face : mesh.getFacesList()) {
                std::size_t expectedDoFs = face->isInternal()
                                               ? 2 * dgDoFsperElement
                                               : dgDoFsperElement;
                check(face, indexing1, owned, nonOwned, expectedDoFs, 0);
            }
            for (const Base::Edge* edge : mesh.getEdgesList()) {
                check(edge, indexing1, owned, nonOwned,
                      edge->getNumberOfElements() * dgDoFsperElement, 0);
            }
            for (const Base::Node* node : mesh.getNodesList()) {
                check(node, indexing1, owned, nonOwned,
                      node->getNumberOfElements() * dgDoFsperElement, 0);
            }
        }

        // Testing with face-connections
        {
            // Rows are DG, columns are DG+CG, face coupling.
            // The row count is the sum of functions with support on the element
            // or on the adjacent elements.
            Utilities::SparsityEstimator estimator(indexing0, indexing01);
            std::vector<int> owned, nonOwned;
            estimator.computeSparsityEstimate(owned, nonOwned, true);
            logger.assert_always(
                owned.size() == indexing0.getNumberOfLocalBasisFunctions(),
                "Wrong size owned");
            logger.assert_always(
                nonOwned.size() == indexing0.getNumberOfLocalBasisFunctions(),
                "Wrong size owned");
            for (const Base::Element* element : mesh.getElementsList()) {
                GeomStorage storage;
                storage.addGeometryAround(element, true);
                check(element, indexing0, owned, nonOwned,
                      storage.countBasisFunctions(), 0);
            }
        }
        {
            // Rows all, columns DG, face coupling
            // Basically counts the number of elements on which the element has
            // support
            Utilities::SparsityEstimator estimator(indexing01, indexing0);
            std::vector<int> owned, nonOwned;
            estimator.computeSparsityEstimate(owned, nonOwned, true);
            logger.assert_always(
                owned.size() == indexing01.getNumberOfLocalBasisFunctions(),
                "Wrong size owned");
            logger.assert_always(
                nonOwned.size() == indexing01.getNumberOfLocalBasisFunctions(),
                "Wrong size owned");
            for (const Base::Element* element : mesh.getElementsList()) {
                GeomStorage storage;
                storage.addGeometryAround(element, true);
                check(element, indexing01, owned, nonOwned,
                      storage.countBasisFunctions(&dgUnknown), 0);
            }
            for (const Base::Face* face : mesh.getFacesList()) {
                GeomStorage storage;
                storage.addGeometryAround(face, true);
                check(face, indexing01, owned, nonOwned,
                      storage.countBasisFunctions(&dgUnknown), 0);
            }
            for (const Base::Edge* edge : mesh.getEdgesList()) {
                GeomStorage storage;
                storage.addGeometryAround(edge, true);
                check(edge, indexing01, owned, nonOwned,
                      storage.countBasisFunctions(&dgUnknown), 0);
            }
            for (const Base::Node* node : mesh.getNodesList()) {
                GeomStorage storage;
                storage.addGeometryAround(node, true);
                check(node, indexing01, owned, nonOwned,
                      storage.countBasisFunctions(&dgUnknown), 0);
            }
        }
    }
}

/**
 * Test with an empty GlobalIndex, where there are no basis functions and thus
 * an empty sparsity pattern would be the result.
 */
void testEmptyIndex() {
    // Test with an empty GlobalIndex
    Utilities::GlobalIndexing emptyIndex;
    Utilities::SparsityEstimator estimator(emptyIndex);
    std::vector<int> owned, nonOwned;
    estimator.computeSparsityEstimate(owned, nonOwned);
    logger.assert_always(owned.empty(),
                         "Non empty owned estimate with empty index");
    logger.assert_always(nonOwned.empty(),
                         "Non empty non owned estimate with empty index");
}

int main(int argc, char** argv) {
    using namespace std::string_literals;
    Base::parse_options(argc, argv);

    testWithDGBasis(1, "3Drectangular1mesh"s);
    testWithDGBasis(3, "3Drectangular1mesh"s);
    testWithDGBasis(1, "3Dtriangular1mesh"s);

    testConformingWith1DMesh();

    testMassOnly({BasisFunctionType::CG4});
    testMassOnly({BasisFunctionType::CG4, BasisFunctionType::CG4,
                  BasisFunctionType::CG4});
    // This configuration actually caused a bug.
    testMassOnly({BasisFunctionType::DG1, BasisFunctionType::DG1,
                  BasisFunctionType::CG4});

    testRowColumnDifference("3Drectangular1mesh"s);
    testRowColumnDifference("3Dtriangular1mesh"s);

    testEmptyIndex();
}

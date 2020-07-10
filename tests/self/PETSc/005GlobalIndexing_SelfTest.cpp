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
using namespace hpgem;
struct IndexStore {
    // The used indices
    std::vector<bool> usedIndices;
    // For blocked based storage
    // The minimum index encountered for an unknown.
    std::vector<std::size_t> minIndices;
    // The maximum index encountered for an unknown
    std::vector<std::size_t> maxIndices;
};

// Check on a single geometrical unit
using Layout = Utilities::GlobalIndexing::Layout;
template <std::size_t DIM, typename GEOM>
void checkLocal(const GEOM* geom, Base::MeshManipulator<DIM>& mesh,
                Utilities::GlobalIndexing& index, Layout layout,
                IndexStore& indexStore) {
    std::size_t numberOfUnknowns = index.getTotalNumberOfUnknowns();
    std::size_t numberOfBasisFunctions = index.getNumberOfLocalBasisFunctions();

    // Last index of previous unknown with basis functions.
    std::size_t prevIndex = 0;
    // Whether prevIndex is set by a previous unknown with basis functions.
    bool prevIndexExists = false;
    for (std::size_t unknown : index.getIncludedUnknowns()) {
        std::size_t bindex = index.getGlobalIndex(geom, unknown);
        // Outside of MPI, so local must be global.
        logger.assert_always(
            bindex == index.getProcessorLocalIndex(geom, unknown),
            "Local == Global");
        std::size_t numLocalBasis =
            geom->getLocalNumberOfBasisFunctions(unknown);
        if (numLocalBasis == 0) {
            // No basis functions == nothing to check
            continue;
        }
        // Last index of the block with local basis functions
        std::size_t bindexLast = bindex + numLocalBasis - 1;
        logger.assert_always(bindexLast < numberOfBasisFunctions,
                             "Index out of range");
        // Update global max and minimum indices for the unknown
        indexStore.minIndices[unknown] =
            std::min(indexStore.minIndices[unknown], bindex);
        indexStore.maxIndices[unknown] =
            std::max(indexStore.maxIndices[unknown], bindexLast);
        // All the intermediate indices
        for (std::size_t ind = bindex; ind <= bindexLast; ++ind) {
            // No index may be used more than once
            logger.assert_always(!indexStore.usedIndices[ind],
                                 "Index already in use");
            indexStore.usedIndices[ind] = true;
            // Test conversion global -> local index, whithout MPI this should
            // be the identity
            logger.assert_always(
                index.globalToProcessorLocalIndex(ind) == ind,
                "Global index -> Local index should be identity");
        }
        if (layout == Layout::SEQUENTIAL && prevIndexExists) {
            // Test whether the indices follow sequentially from the indices
            // from the previous basis function.
            logger.assert_always(bindex == prevIndex + 1, "Non sequential");
        }
        // Update previous index (for checking sequentiality)
        prevIndex = bindexLast;
        prevIndexExists = true;
    }
}

template <std::size_t DIM>
void checkIndex(Base::MeshManipulator<DIM>& mesh,
                Utilities::GlobalIndexing& index, Layout layout) {
    index.verifyCompleteIndex();
    std::size_t numberOfUnknowns = index.getTotalNumberOfUnknowns();
    std::size_t numberOfBasisFunctions = index.getNumberOfLocalBasisFunctions();
    // Remove check when numberOfUnknowns is not in configData
    logger.assert_always(
        numberOfUnknowns == mesh.getConfigData()->numberOfUnknowns_,
        "Matching unknown count");

    IndexStore indexStore;
    indexStore.usedIndices.resize(numberOfBasisFunctions, false);
    indexStore.minIndices.resize(numberOfUnknowns,
                                 std::numeric_limits<std::size_t>::max());
    indexStore.maxIndices.resize(numberOfUnknowns, 0);

    for (Base::Element* element : mesh.getElementsList()) {
        checkLocal(element, mesh, index, layout, indexStore);
    }
    for (Base::Face* face : mesh.getFacesList()) {
        checkLocal(face, mesh, index, layout, indexStore);
    }
    for (Base::Edge* edge : mesh.getEdgesList()) {
        checkLocal(edge, mesh, index, layout, indexStore);
    }
    if (DIM > 1) {
        for (Base::Node* node : mesh.getNodesList()) {
            checkLocal(node, mesh, index, layout, indexStore);
        }
    }

    // Check that all indices are used.
    for (const bool b : indexStore.usedIndices) {
        logger.assert_always(b, "Unused index");
    }
    // Ensure no resize
    logger.assert_always(
        indexStore.usedIndices.size() == numberOfBasisFunctions,
        "Resized index");
    // Check blocked layouts
    if (layout != Layout::SEQUENTIAL) {
        const std::vector<std::size_t>& unknowns = index.getIncludedUnknowns();
        for (std::size_t i = 1; i < unknowns.size(); ++i) {
            // The indices of unknowns[i] should sequentially follow those of
            // unknowns[i-1]
            logger.assert_always(indexStore.minIndices[unknowns[i]] ==
                                     indexStore.maxIndices[unknowns[i - 1]] + 1,
                                 "Non blocked");
        }
    }
}

int main(int argc, char** argv) {
    // TODO: Untested getGlobalIndex(element), getGlobalIndex(face)
    using namespace std::string_literals;
    Base::parse_options(argc, argv);

    // Create a mesh to work with
    std::size_t numberOfUnknowns = 3;
    Base::ConfigurationData config(numberOfUnknowns);
    Base::MeshManipulator<3> mesh(&config);

    std::stringstream filename;
    filename << Base::getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/"s
             << "3Drectangular1mesh"s
             << ".hpgem";
    mesh.readMesh(filename.str());

    // Basis functions for the first tests. Using 0-th order DG basis functions
    // we have a single basis function per element for each unknown.
    mesh.useDefaultDGBasisFunctions(0);
    mesh.useDefaultDGBasisFunctions(0, 1);
    mesh.useDefaultDGBasisFunctions(0, 2);

    for (auto layout : {Layout::BLOCKED_PROCESSOR, Layout::BLOCKED_GLOBAL,
                        Layout::SEQUENTIAL}) {
        Utilities::GlobalIndexing index(&mesh, layout);
        std::size_t numberOfElements = mesh.getNumberOfElements();
        logger.assert_always(index.getNumberOfLocalBasisFunctions() ==
                                 numberOfUnknowns * numberOfElements,
                             "Number of basis functions");
        checkIndex(mesh, index, layout);
    }

    // Same tests, but with higher order basis functions. These have multiple
    // basis
    // functions per element and support on each geometric object (node, edge,
    // face).
    mesh.useDefaultDGBasisFunctions(2);
    // Fourth order to have basis functions with support based on nodes, edges,
    // faces and elements.
    mesh.useDefaultConformingBasisFunctions(4, 1);
    mesh.useDefaultDGBasisFunctions(0, 2);

    for (auto layout : {Layout::BLOCKED_PROCESSOR, Layout::BLOCKED_GLOBAL,
                        Layout::SEQUENTIAL}) {
        Utilities::GlobalIndexing index(&mesh, layout);
        checkIndex(mesh, index, layout);
    }

    // Test with subset of the unknowns.
    for (auto layout : {Layout::BLOCKED_PROCESSOR, Layout::BLOCKED_GLOBAL,
                        Layout::SEQUENTIAL}) {
        std::vector<std::size_t> usedBasisUnknowns = {0, 1};
        Utilities::GlobalIndexing index(&mesh, layout, &usedBasisUnknowns);
        logger.assert_always(index.getIncludedUnknowns() == usedBasisUnknowns,
                             "Different included unknowns");
        checkIndex(mesh, index, layout);
    }
}

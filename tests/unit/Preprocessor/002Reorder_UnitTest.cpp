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

#include <CMakeDefinitions.h>
#include "Base/MeshManipulator.h"
#include "Base/CommandLineOptions.h"

#include "Utilities/ElementLocalIndexing.h"

#define CATCH_CONFIG_RUNNER
#include "../catch.hpp"

using namespace hpgem;

void validate(const Utilities::ElementLocalIndexing& indexing,
              std::size_t numberOfUnknowns) {
    indexing.validateInternalState();

    // Check the values for all the included and not included unknowns
    std::vector<bool> included(numberOfUnknowns, false);
    for (auto& unknown : indexing.getIncludedUnknowns()) {
        included[unknown] = true;
    }
    std::size_t totalNumberOfDoFs = indexing.getNumberOfDoFs();
    for (std::size_t i = 0; i < numberOfUnknowns; ++i) {
        if (included[i]) {
            INFO("Unknown DoF index does not go over max")
            CHECK(indexing.getDoFOffset(i) + indexing.getNumberOfDoFs(i) <=
                  totalNumberOfDoFs);
        } else {
            INFO("Non included unknown has size 0")
            CHECK(indexing.getNumberOfDoFs(i) == 0);
            INFO("Non included unknown has correct offset")
            CHECK(indexing.getDoFOffset(i) ==
                  Utilities::ElementLocalIndexing::UNKNOWN_NOT_INCLUDED);
        }
    }
}

TEST_CASE("Basic indexing", "[LocalIndexing]") {
    // Create a simple 1D mesh to work with
    const std::size_t numberOfUnknowns = 3;
    Base::ConfigurationData config(numberOfUnknowns);
    Base::MeshManipulator<1> mesh(&config);

    {
        using namespace std::string_literals;
        std::stringstream filename;
        filename << Base::getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/"s
                 << "1Drectangular2mesh"s
                 << ".hpgem";
        mesh.readMesh(filename.str());
    }

    // Random mix of DoF types.
    mesh.useDefaultDGBasisFunctions(2);
    mesh.useDefaultConformingBasisFunctions(3, 1);
    mesh.useDefaultDGBasisFunctions(0, 2);

    // Validate with empty unknowns
    Utilities::ElementLocalIndexing indexing;
    INFO("Initial index should be empty");
    REQUIRE(indexing.getIncludedUnknowns().empty() == true);
    for (const Base::Element* element : mesh.getElementsList()) {
        indexing.reinit(element);
        INFO("Element should match the reinit");
        REQUIRE(indexing.getElement() == element);
        validate(indexing, numberOfUnknowns);
    }

    // Switch unknowns and revalidate
    std::vector<std::size_t> unknowns(1, 1);
    indexing.reinit(unknowns);
    indexing.validateInternalState();
    for (const Base::Element* element : mesh.getElementsList()) {
        indexing.reinit(element);
        validate(indexing, numberOfUnknowns);
    }

    // Verify with local data
    unknowns.resize(3);
    unknowns[0] = 0;
    unknowns[1] = 1;
    unknowns[2] = 2;
    indexing.reinit(unknowns);
    INFO("Included unknowns should match reinit");
    REQUIRE(indexing.getIncludedUnknowns() == unknowns);
    for (const Base::Element* element : mesh.getElementsList()) {
        indexing.reinit(element);
        INFO("Check total basis function count");
        REQUIRE(indexing.getNumberOfDoFs() ==
                element->getTotalNumberOfBasisFunctions());
        INFO("Offset for first unknown should be zero")
        REQUIRE(indexing.getDoFOffset(0) == 0);
        std::size_t dofZeroSize = element->getNumberOfBasisFunctions(0);

        INFO("Correct number of DoFs for unknown 0");
        REQUIRE(indexing.getNumberOfDoFs(0) == dofZeroSize);

        INFO("Only the first unknown DoFs are before those of the second one");
        REQUIRE(indexing.getDoFOffset(1) == dofZeroSize);
        // More can be added if desired
    }

    // Check with partial unknowns, dropping the first unknown
    unknowns.resize(2);
    unknowns[0] = 1;
    unknowns[1] = 2;
    indexing.reinit(unknowns);
    for (const Base::Element* element : mesh.getElementsList()) {
        indexing.reinit(element);
        INFO("First included unknown should have offset 0");
        REQUIRE(indexing.getDoFOffset(1) == 0);
        std::size_t dof1Size = element->getNumberOfBasisFunctions(1);
        INFO("Correct number of DoFs for first unknown");
        REQUIRE(indexing.getNumberOfDoFs(1) == dof1Size);
        INFO("Correct offset of second unknown");
        REQUIRE(indexing.getDoFOffset(2) == dof1Size);
    }

    // Deassociate with an element and check
    indexing.reinit(nullptr);
    INFO("Without element there are no DoFs");
    REQUIRE(indexing.getNumberOfDoFs() == 0);
}

int main(int argc, char* argv[]) {

    // The mesh reader needs an initialized MPI and possibly more. To achieve
    // this initialize the system using only the program name, as the arguments
    // to the program are most probably used to direct Catch.
    Base::parse_options(1, argv);

    int result = Catch::Session().run(argc, argv);

    return result;
}
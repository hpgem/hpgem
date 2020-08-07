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

#include "Utilities/LocalIndexing.h"

namespace hpgem {

void run() {
    // Create a simple 1D mesh to work with
    std::size_t numberOfUnknowns = 3;
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
    Utilities::ElementLocalIndexing indexing(numberOfUnknowns);
    logger.assert_always(indexing.getIncludedUnknowns().empty(), "Non empty starting index");
    for (const Base::Element* element : mesh.getElementsList()) {
        indexing.reinit(element);
        logger.assert_always(element == indexing.getElement(), "Element not matching after reinit");
        indexing.validate();
    }

    // Switch unknowns and revalidate
    std::vector<std::size_t> unknowns(1, 1);
    indexing.reinit(unknowns);
    indexing.validate();
    for (const Base::Element* element : mesh.getElementsList()) {
        indexing.reinit(element);
        indexing.validate();
    }

    // Verify with local data
    unknowns.resize(3);
    unknowns[0] = 0;
    unknowns[1] = 1;
    unknowns[2] = 2;
    indexing.reinit(unknowns);
    logger.assert_always(unknowns == indexing.getIncludedUnknowns(), "Not all unknowns are included.");
    for (const Base::Element* element : mesh.getElementsList()) {
        indexing.reinit(element);
        logger.assert_always(indexing.getNumberOfDoFs() ==
                                 element->getTotalNumberOfBasisFunctions(),
                             "Incorrect total basis function count");
        logger.assert_always(indexing.getDoFOffset(0) == 0, "Incorrect offset for DoF 0");
        std::size_t dofZeroSize = element->getNumberOfBasisFunctions(0);
        logger.assert_always(indexing.getDoFOffset(1) == dofZeroSize, "Incorrect offset for DoF 1");
        logger.assert_always(indexing.getNumberOfDoFs(0) == dofZeroSize, "Incorrect DoFSize");
        // More can be added if desired
    }

    // Check with partial unknowns
    unknowns.resize(2);
    unknowns[0] = 1;
    unknowns[1] = 2;
    indexing.reinit(unknowns);
    for (const Base::Element* element : mesh.getElementsList()) {
        indexing.reinit(element);
        logger.assert_always(indexing.getDoFOffset(1) == 0, "Incorrect offset of DoF1 (partial index)");
        std::size_t dof1Size = element->getNumberOfBasisFunctions(1);
        logger.assert_always(indexing.getNumberOfDoFs(1) == dof1Size, "Incorrect size for DoF1 (partial index)");
        logger.assert_always(indexing.getDoFOffset(2) == dof1Size, "Incorrect offset for DoF 2 (partial index)");
    }

    // Deassociate with an element and check
    indexing.reinit(nullptr);
    logger.assert_always(indexing.getNumberOfDoFs() == 0, "No element and positive DoFs");
}

}  // namespace hpgem

int main(int argc, char** argv) {
    hpgem::Base::parse_options(argc, argv);
    hpgem::run();
    return 0;
}

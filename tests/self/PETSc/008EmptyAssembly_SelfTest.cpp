/*
This file forms part of hpGEM. This package has been developed over a number of
years by various people at the University of Twente and a full list of
contributors can be found at http://hpgem.org/about-the-code/team

This code is distributed using BSD 3-Clause License. A copy of which can found
below.


Copyright (c) 2020, Univesity of Twenete
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <Base/CommandLineOptions.h>
#include <Base/ConfigurationData.h>
#include <Base/MeshManipulator.h>
#include <CMakeDefinitions.h>
#include "Utilities/GlobalMatrix.h"
#include "Utilities/GlobalVector.h"

void testWithEmptyIndex() {
    // Test that creating a GlobalMatrix and GlobalVector with an empty index is
    // possible
    using namespace Utilities;

    GlobalIndexing emptyIndex;
    GlobalPetscMatrix emptyMatrix(emptyIndex, -1);
    {
        PetscErrorCode ierr;
        PetscInt rows, cols;
        ierr = MatGetSize(emptyMatrix, &rows, &cols);
        CHKERRV(ierr);
        logger.assert_always(rows == 0,
                             "Non zero rows for matrix without index");
        logger.assert_always(cols == 0,
                             "Non zero cols for matrix without index");
    }

    // Should already be done, but reassemble here
    emptyMatrix.assemble();

    GlobalPetscVector emptyVector(emptyIndex);
    {
        PetscErrorCode ierr;
        PetscInt rows;
        ierr = VecGetSize(emptyVector, &rows);
        CHKERRV(ierr);
        logger.assert_always(rows == 0,
                             "Non zero rows for vector without index");
    }
    emptyVector.assemble();
}

void testReinit() {
    using namespace Utilities;
    GlobalIndexing indexing;
    GlobalPetscMatrix matrix(indexing, -1);
    GlobalPetscVector vector(indexing, -1, -1);

    // Load any mesh.
    Base::ConfigurationData data(1);  // 1 unknown
    Base::MeshManipulator<1>* mesh = new Base::MeshManipulator<1>(&data);
    using namespace std::string_literals;
    mesh->readMesh(Base::getCMAKE_hpGEM_SOURCE_DIR() +
                   "/tests/files/poissonMesh1.hpgem"s);

    // Associate some basis functions with the elements so that there will be
    // DoFs
    mesh->useDefaultDGBasisFunctions(0);

    // Reinit everything
    indexing.reset(
        mesh, GlobalIndexing::Layout::SEQUENTIAL);  // Layout does not matter
    matrix.reinit();
    vector.reinit();

    // Now the matrix and vector should have rows & columns
    {
        PetscErrorCode ierr;
        PetscInt rows, cols;
        ierr = MatGetSize(matrix, &rows, &cols);
        CHKERRV(ierr);
        logger.assert_always(rows > 0, "No rows after matrix reinit");
        logger.assert_always(cols > 0, "No columns after matrix reinit");

        ierr = VecGetSize(vector, &rows);
        CHKERRV(ierr);
        logger.assert_always(rows > 0, "No rows after vector reinit");
    }
}

int main(int argc, char** argv) {
    Base::parse_options(argc, argv);

    testWithEmptyIndex();
    testReinit();
    return 0;
}
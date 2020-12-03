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

#include "Base/CommandLineOptions.h"

#include "Utils/MatrixBlocks.h"

#include <petsc.h>
#define CATCH_CONFIG_RUNNER
#include "../catch.hpp"

namespace DGMax {

TEST_CASE("Matrix blocks init", "[MatrixBlocksUnitTest]") {
    std::size_t rows = 2;
    std::size_t columns = 3;
    LinearAlgebra::MiddleSizeMatrix block1(rows, columns, 3.0);

    std::vector<PetscInt> rowIndices({0, 1});
    std::vector<PetscInt> colIndices({2, 3, 4});

    SECTION("Single (non Hermitian) block") {
        // Single block
        MatrixBlocks blocks(rowIndices, colIndices, block1);

        INFO("Matrix block is not a pair");
        CHECK_FALSE(blocks.isPair());

        INFO("Block has the right size");
        REQUIRE(blocks.getBlockSize() == rows * columns);
    }
    SECTION("Symmetrically placed pair of blocks") {
        // Paired block
        LinearAlgebra::MiddleSizeMatrix block2(columns, rows, 3.0);
        MatrixBlocks blocks(rowIndices, colIndices, block1, block2);

        INFO("Block is a pair");
        CHECK(blocks.isPair());
        INFO("Block has the right size");
        REQUIRE(blocks.getBlockSize() == rows * columns);
    }
}

TEST_CASE("Block insertion", "[MatrixBlocksUnitTest]") {

    LinearAlgebra::MiddleSizeMatrix block(2, 2);
    block(0, 0) = 1.0;
    block(0, 1) = 2.0;
    block(1, 0) = 3.0;
    block(1, 1) = 4.0;
    std::vector<PetscInt> rows({0, 1});
    std::vector<PetscInt> columns({2, 3});

    SECTION("Single block") {
        MatrixBlocks blocks(rows, columns, block);
        // Test with a dense matrix, as that is the simplest to create
        Mat mat;
        MatCreateSeqDense(MPI_COMM_SELF, 5, 5, nullptr, &mat);
        MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);

        // Set to column orientation to match hpgem.
        MatSetOption(mat, MAT_ROW_ORIENTED, PETSC_FALSE);

        // Insert blocks
        std::vector<PetscScalar> tempStorage;
        blocks.loadBlocks(tempStorage);
        INFO("Sufficient space allocated in storage");
        REQUIRE(tempStorage.size() >= block.size());
        blocks.insertBlocks(tempStorage, mat);

        MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);

        PetscScalar value;
        PetscInt ri, ci;
        INFO("Entries placed correctly");
        for (std::size_t i = 0; i < rows.size(); ++i) {
            for (std::size_t j = 0; j < columns.size(); ++j) {
                ri = rows[i];
                ci = columns[j];
                MatGetValues(mat, 1, &ri, 1, &ci, &value);
                REQUIRE(value == block(i, j));
            }
        }
        // Scaling test
        blocks.loadBlocks(tempStorage);
        INFO("Sufficient space allocated in storage (2)");
        REQUIRE(tempStorage.size() >= block.size());
        for (std::size_t i = 0; i < block.size(); ++i) {
            tempStorage[i] *= 3.0;
        }

        blocks.insertBlocks(tempStorage, mat);

        MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
        INFO("Scaled entries placed correctly");
        for (std::size_t i = 0; i < rows.size(); ++i) {
            for (std::size_t j = 0; j < columns.size(); ++j) {
                ri = rows[i];
                ci = columns[j];
                MatGetValues(mat, 1, &ri, 1, &ci, &value);
                REQUIRE(value == 3.0 * block(i, j));
            }
        }
        // Cleanup
        MatDestroy(&mat);
    }

    SECTION("Symmetrically placed pair of blocks") {

        // Test with a symmetric pair of blocks
        LinearAlgebra::MiddleSizeMatrix block2(2, 2);
        block2(0, 0) = 5.0;
        block2(0, 1) = 6.0;
        block2(1, 0) = 7.0;
        block2(1, 1) = 8.0;
        MatrixBlocks blocks(rows, columns, block, block2);
        // Test with a dense matrix, as that is the simplest to create
        Mat mat;
        MatCreateSeqDense(MPI_COMM_SELF, 5, 5, nullptr, &mat);
        MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);

        // Set to column orientation to match hpgem.
        MatSetOption(mat, MAT_ROW_ORIENTED, PETSC_FALSE);

        // Insert blocks
        std::vector<PetscScalar> tempStorage;
        blocks.loadBlocks(tempStorage);
        INFO("Sufficient space allocated in storage");
        REQUIRE(tempStorage.size() >= 2 * block.size());
        blocks.insertBlocks(tempStorage, mat);

        MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);

        PetscScalar value;
        PetscInt ri, ci;
        for (std::size_t i = 0; i < rows.size(); ++i) {
            for (std::size_t j = 0; j < columns.size(); ++j) {
                ri = rows[i];
                ci = columns[j];
                MatGetValues(mat, 1, &ri, 1, &ci, &value);
                INFO("Correctly placed value in symmetric case");
                REQUIRE(value == block(i, j));
                MatGetValues(mat, 1, &ci, 1, &ri, &value);
                INFO("Correctly placed symmetric value in symmetric case");
                REQUIRE(value == block2(j, i));
            }
        }
    }
}

}  // namespace DGMax

int main(int argc, char* argv[]) {
    PetscErrorCode ierr = PetscInitialize(&argc, &argv, nullptr, nullptr);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);

    int result = Catch::Session().run(argc, argv);

    ierr = PetscFinalize();
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    return result;
}
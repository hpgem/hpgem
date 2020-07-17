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

namespace DGMax {

void basicInitTest() {
    std::size_t rows = 2;
    std::size_t columns = 3;
    LinearAlgebra::MiddleSizeMatrix block1(rows, columns, 3.0);

    std::vector<PetscInt> rowIndices({0, 1});
    std::vector<PetscInt> colIndices({2, 3, 4});

    {
        // Single block
        MatrixBlocks blocks(rowIndices, colIndices, block1);
        logger.assert_always(blocks.getBlockSize() == rows * columns,
                             "Incorrect block size");
        logger.assert_always(!blocks.isPair(),
                             "Incorrect pair value for single block");
    }
    {
        // Paired block
        LinearAlgebra::MiddleSizeMatrix block2(columns, rows, 3.0);
        MatrixBlocks blocks(rowIndices, colIndices, block1, block2);
        logger.assert_always(blocks.getBlockSize() == rows * columns,
                             "Incorrect block size");
        logger.assert_always(blocks.isPair(),
                             "Incorrect pair value for pair block");
    }
}

void matrixInsertionTest() {

    LinearAlgebra::MiddleSizeMatrix block(2, 2);
    block(0, 0) = 1.0;
    block(0, 1) = 2.0;
    block(1, 0) = 3.0;
    block(1, 1) = 4.0;
    std::vector<PetscInt> rows({0, 1});
    std::vector<PetscInt> columns({2, 3});

    {
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
        logger.assert_always(tempStorage.size() >= block.size(),
                             "Too small storage");
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
                logger.assert_always(value == block(i, j),
                                     "Incorrectly placed value");
            }
        }
        // Scaling test

        blocks.loadBlocks(tempStorage);
        logger.assert_always(tempStorage.size() >= block.size(),
                             "Too small storage");
        for (std::size_t i = 0; i < block.size(); ++i) {
            tempStorage[i] *= 3.0;
        }

        blocks.insertBlocks(tempStorage, mat);

        MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
        for (std::size_t i = 0; i < rows.size(); ++i) {
            for (std::size_t j = 0; j < columns.size(); ++j) {
                ri = rows[i];
                ci = columns[j];
                MatGetValues(mat, 1, &ri, 1, &ci, &value);
                logger.assert_always(value == 3.0 * block(i, j),
                                     "Incorrect scaled value");
            }
        }
        // Cleanup
        MatDestroy(&mat);
    }
    {
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
        logger.assert_always(tempStorage.size() >= 2 * block.size(),
                             "Too small storage");
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
                logger.assert_always(
                    value == block(i, j),
                    "Incorrectly placed value in symmetric case");
                MatGetValues(mat, 1, &ci, 1, &ri, &value);
                logger.assert_always(
                    value == block2(j, i),
                    "Incorrectly placed symmetric value in symmetric case");
            }
        }
    }
}

}  // namespace DGMax

int main(int argc, char** argv) {
    Base::parse_options(argc, argv);  // To init petsc

    DGMax::basicInitTest();
    DGMax::matrixInsertionTest();

    return 0;
}
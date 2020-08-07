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
#include "MatrixBlocks.h"

namespace DGMax {

void MatrixBlocks::insertBlock(std::vector<PetscScalar>& storage,
                               bool hermitianPart, Mat mat) const {
    PetscErrorCode err;
#ifdef HPGEM_ASSERTS
    // Check if the matrix is column oriented, as MiddleSizeMatrix stores it in
    // column order.
    PetscBool isRowOriented = PETSC_FALSE;
    // TODO: Does not seem to do anything due to what I expect is a bug in PETSC
    err = MatGetOption(mat, MAT_ROW_ORIENTED, &isRowOriented);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    logger.assert_debug(isRowOriented == PETSC_FALSE,
                        "Requires column orientation");
#endif
    for (std::size_t i = 0; i < rowIndices_.size(); ++i) {
        for (std::size_t j = 0; j < columnIndices_.size(); ++j) {
            insertValue(storage, hermitianPart, mat, i, j);
        }
    }
    //
    //    if (!hermitianPart) {
    //        // Normal block (block1)
    //        err =
    //            MatSetValues(mat, (PetscInt)rowIndices_.size(),
    //            rowIndices_.data(),
    //                         (PetscInt)columnIndices_.size(),
    //                         columnIndices_.data(), storage.data(),
    //                         INSERT_VALUES);
    //    } else {
    //        // Hermitian block (block2)
    //        PetscScalar* data = storage.data() + getBlockSize();
    //        err = MatSetValues(mat, (PetscInt)columnIndices_.size(),
    //                           columnIndices_.data(),
    //                           (PetscInt)rowIndices_.size(),
    //                           rowIndices_.data(), data, INSERT_VALUES);
    //    }
    CHKERRABORT(PETSC_COMM_WORLD, err);
}

void MatrixBlocks::insertValue(std::vector<PetscScalar>& storage,
                               bool hermitianPart, Mat mat, std::size_t r,
                               std::size_t c) const {
    PetscScalar val;
    if (!hermitianPart) {
        val = storage[r + c * rowIndices_.size()];
        r = rowIndices_[r];
        c = columnIndices_[c];
    } else {
        val = storage[c + columnIndices_.size() * r + getBlockSize()];
        std::swap(r, c);
        r = columnIndices_[r];
        c = rowIndices_[c];
    }
    if (val == 0.0) {
        return;
    }
    PetscErrorCode err = MatSetValue(mat, r, c, val, INSERT_VALUES);
    CHKERRABORT(PETSC_COMM_WORLD, err);
}

void MatrixBlocks::validate() const {
    logger.assert_debug(rowIndices_.size() == block1_.getNumberOfRows(),
                        "Wrong rows size for block1");
    logger.assert_debug(columnIndices_.size() == block1_.getNumberOfColumns(),
                        "Wrong column size for block1");

    if (pair_) {
        logger.assert_debug(rowIndices_.size() == block2_.getNumberOfColumns(),
                            "Wrong column size for block2");
        logger.assert_debug(columnIndices_.size() == block2_.getNumberOfRows(),
                            "Wrong rows size for block2");
    }
}

}  // namespace DGMax

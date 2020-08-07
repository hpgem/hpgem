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
#ifndef HPGEM_MATRIXBLOCKS_H
#define HPGEM_MATRIXBLOCKS_H

#include <vector>
#include <petsc.h>

#include "LinearAlgebra/MiddleSizeMatrix.h"

namespace DGMax {

using namespace hpgem;

/// A block from the global matrix, or a pair of blocks that are placed
/// symmetrically.
class MatrixBlocks {
   public:
    MatrixBlocks(std::vector<PetscInt> rowIndices,
                 std::vector<PetscInt> columnIndices,
                 LinearAlgebra::MiddleSizeMatrix block)
        : rowIndices_(std::move(rowIndices)),
          columnIndices_(std::move(columnIndices)),
          block1_(std::move(block)),
          pair_(false){};

    MatrixBlocks(std::vector<PetscInt> rowIndices,
                 std::vector<PetscInt> columnIndices,
                 LinearAlgebra::MiddleSizeMatrix block1,
                 LinearAlgebra::MiddleSizeMatrix block2)
        : rowIndices_(std::move(rowIndices)),
          columnIndices_(std::move(columnIndices)),
          block1_(std::move(block1)),
          block2_(std::move(block2)),
          pair_(true){};

    /// Number of entries in the matrix block. If this is pair, then this gives
    /// the size of a single block.
    std::size_t getBlockSize() const {
        return rowIndices_.size() * columnIndices_.size();
    }

    /// Whether this is a pair of blocks or a single block.
    bool isPair() const { return pair_; }

    /// \brief Load blocks into temporary storage for insertion.
    ///
    /// This stores the raw entries of the block in the storage space for use
    /// with insertBlocks(). The first getBlockSize() entries of the storage are
    /// used for storing the first matrix block. If this object is a pair then
    /// the subsequent getBlockSize() entries will be used for the second block.
    ///
    /// \param storage The storage to use, will be resized if necessary.
    void loadBlocks(std::vector<PetscScalar>& storage) const {
        // Make sure we have enough space
        std::size_t storageSize = getBlockSize();
        std::size_t storageMultiplier = pair_ ? 2 : 1;
        if (storage.size() < storageMultiplier * storageSize) {
            storage.resize(storageMultiplier * storageSize);
        }
        // Copy data
        for (std::size_t i = 0; i < storageSize; ++i) {
            storage[i] = block1_[i];
        }
        if (pair_) {
            for (std::size_t i = 0; i < storageSize; ++i) {
                storage[i + storageSize] = block2_[i];
            }
        }
    }

    /// \brief Insert the block or blocks in a PetscMatrix.
    void insertBlocks(std::vector<PetscScalar>& storage, Mat mat) const {
        insertBlock(storage, false, mat);
        if (pair_) {
            insertBlock(storage, true, mat);
        }
    }

   private:
    /// The global indices for the rows
    std::vector<PetscInt> rowIndices_;
    /// The global indices for the columns
    std::vector<PetscInt> columnIndices_;

    /// The primary block
    LinearAlgebra::MiddleSizeMatrix block1_;
    /// The symmetrically placed block
    LinearAlgebra::MiddleSizeMatrix block2_;

    /// Whether this is a single block or a symmetrically placed pair.
    bool pair_;

    /// Helper function inserting a single block.
    void insertBlock(std::vector<PetscScalar>& storage, bool hermitianPart,
                     Mat mat) const;

    void insertValue(std::vector<PetscScalar>& storage, bool hermitianPart,
                     Mat mat, std::size_t r, std::size_t c) const;

    /// Validate the correctnees
    void validate() const;
};

}  // namespace DGMax

#endif  // HPGEM_MATRIXBLOCKS_H

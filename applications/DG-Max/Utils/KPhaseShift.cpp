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
#include "KPhaseShift.h"

#include "LinearAlgebra/MiddleSizeMatrix.h"
#include "Utilities/GlobalIndexing.h"

namespace DGMax {

template <std::size_t DIM>
void KPhaseShiftBlock<DIM>::apply(LinearAlgebra::SmallVector<DIM> k,
                                  std::vector<PetscScalar>& storage,
                                  Mat mat) const {
    // Note we could optimize and only apply this if needed, that is if
    // exp(ixk_old) = exp(ix k_new). However, that may give a minimal
    // improvement in performance and would add the need to keep track of k_old.

    // Load values in storage
    blocks_.loadBlocks(storage);
    // Apply phase shift
    std::size_t blockSize = blocks_.getBlockSize();
    const std::complex<double> phase =
        std::exp(std::complex<double>(0, k * dx_));

    for (std::size_t i = 0; i < blockSize; ++i) {
        storage[i] *= phase;
    }
    if (blocks_.isPair()) {
        const std::complex<double> antiPhase =
            std::exp(std::complex<double>(0, -(k * dx_)));
        for (std::size_t i = 0; i < blockSize; ++i) {
            storage[i + blockSize] *= antiPhase;
        }
    }
    // Insert the now phase shifted blocks into the matrix
    blocks_.insertBlocks(storage, mat);
}

template <std::size_t DIM>
void KPhaseShiftBlock<DIM>::applyDerivative(LinearAlgebra::SmallVector<DIM> k,
                                            LinearAlgebra::SmallVector<DIM> dk,
                                            std::vector<PetscScalar>& storage,
                                            Mat mat) const {
    // Load value in storage
    blocks_.loadBlocks(storage);
    std::size_t blockSize = blocks_.getBlockSize();
    const std::complex<double> phaseDeriv =
        std::complex<double>(0, dk * dx_) *  // From taking dk . grad_k
        std::exp(std::complex<double>(0, k * dx_));
    std::cout << "Phase deriv" << phaseDeriv << std::endl;
    for (std::size_t i = 0; i < blockSize; ++i) {
        storage[i] *= phaseDeriv;
    }
    if (blocks_.isPair()) {
        const std::complex<double> antiPhaseDeriv =
            std::complex<double>(0, -dk * dx_) *
            std::exp(std::complex<double>(0, -k * dx_));
        std::cout << "Phase deriv" << antiPhaseDeriv << std::endl;
        for (std::size_t i = blockSize; i < 2 * blockSize; ++i) {
            storage[i] *= antiPhaseDeriv;
        }
        blocks_.insertBlocks(storage, mat);
    }
}

/// Phase Shifts
template <std::size_t DIM>
void KPhaseShifts<DIM>::apply(LinearAlgebra::SmallVector<DIM> k,
                              Mat mat) const {
    PetscErrorCode error = MatSetOption(mat, MAT_ROW_ORIENTED, PETSC_FALSE);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    std::vector<PetscScalar> storage;
    for (const KPhaseShiftBlock<DIM>& block : blocks_) {
        block.apply(k, storage, mat);
    }

    // Go back to the default row orientation
    error = MatSetOption(mat, MAT_ROW_ORIENTED, PETSC_TRUE);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    error = MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
    CHKERRABORT(PETSC_COMM_WORLD, error);
}

template <std::size_t DIM>
void KPhaseShifts<DIM>::applyDeriv(LinearAlgebra::SmallVector<DIM> k,
                                   LinearAlgebra::SmallVector<DIM> dk,
                                   Mat mat) const {
    // Normalize
    if (dk.l2NormSquared() > 1e-8) {
        dk /= dk.l2Norm();
    } else {
        logger(
            WARN,
            "Computing directional derivative with zero vector as direction");
    }
    PetscErrorCode error = MatSetOption(mat, MAT_ROW_ORIENTED, PETSC_FALSE);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    std::vector<PetscScalar> storage;
    for (const KPhaseShiftBlock<DIM>& block : blocks_) {
        block.applyDerivative(k, dk, storage, mat);
    }

    // Go back to the default row orientation
    error = MatSetOption(mat, MAT_ROW_ORIENTED, PETSC_TRUE);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    error = MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
    CHKERRABORT(PETSC_COMM_WORLD, error);
}

template class KPhaseShifts<2>;
template class KPhaseShifts<3>;

}  // namespace DGMax
/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2021, University of Twente
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
#ifdef HPGEM_USE_SLEPC

#include <numeric>
#include "Eigenpairs.h"




namespace hpgem {
namespace Utilities {

Eigenpairs::Eigenpairs()
    : numberOfEigenvectors_(0),
      eigenvectors_(nullptr),
      eigenvalues_(0),
      ordering_(0){};

Eigenpairs::Eigenpairs(Eigenpairs &&other) noexcept
    : numberOfEigenvectors_(0),
      eigenvectors_(nullptr),
      eigenvalues_(0),
      ordering_(0) {
    // Essential, as the other vector now does not have any eigenvalues
    std::swap(numberOfEigenvectors_, other.numberOfEigenvectors_);
    std::swap(eigenvectors_, other.eigenvectors_);
    std::swap(eigenvalues_, other.eigenvalues_);
    std::swap(ordering_, other.ordering_);
}

Eigenpairs::~Eigenpairs() {
    if (numberOfEigenvectors_ > 0) {
        PetscErrorCode err =
            VecDestroyVecs(numberOfEigenvectors_, &eigenvectors_);
        CHKERRABORT(PETSC_COMM_WORLD, err);
    }
}

Eigenpairs &Eigenpairs::operator=(Eigenpairs &&other) noexcept {
    std::swap(numberOfEigenvectors_, other.numberOfEigenvectors_);
    std::swap(eigenvectors_, other.eigenvectors_);
    std::swap(eigenvalues_, other.eigenvalues_);
    std::swap(ordering_, other.ordering_);
    return *this;
}

void Eigenpairs::reserve(std::size_t newSize, Vec sample) {
    if (newSize <= numberOfEigenvectors_) {
        return;
    }
    PetscErrorCode err = VecDestroyVecs(numberOfEigenvectors_, &eigenvectors_);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    err = VecDuplicateVecs(sample, newSize, &eigenvectors_);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    numberOfEigenvectors_ = newSize;
}

/**
 * \brief Load eigenpairs from a solver
 *
 * @param eps The solver to load the eigenpairs from
 * @param sample A sample vector to duplicate for storing eigenvectors
 */
void Eigenpairs::loadEigenpairs(EPS eps, Vec sample) {
    PetscInt converged;
    PetscErrorCode err;
    err = EPSGetConverged(eps, &converged);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    // Ensure that we have enough space
    reserve(converged, sample);
    eigenvalues_.resize(converged);
    ordering_.resize(converged);
    // Load the eigenpairs
    for (PetscInt i = 0; i < converged; ++i) {
        // nullptrs for the imaginary parts
        err = EPSGetEigenpair(eps, i, &eigenvalues_[i], nullptr,
                              eigenvectors_[i], nullptr);
        CHKERRABORT(PETSC_COMM_WORLD, err);
    }
    // Reset the ordering
    std::iota(ordering_.begin(), ordering_.end(), 0);
}

void Eigenpairs::loadEigenpairs(LinearAlgebra::JacobiDavidsonMaxwellSolver eps, Vec sample) {

}

void Eigenpairs::reorder(std::vector<std::size_t> ordering) {
    logger.assert_always(ordering.size() == size(),
                         "Ordering is of incorrect size");
    // New ordering is with respect to the current ordering. So compute the
    // actual indices of the new ordering.
    for (std::size_t i = 0; i < size(); ++i) {
        ordering[i] = ordering_[ordering[i]];
    }
    // Set the ordering
    ordering_ = ordering;
}

}  // namespace Utilities
}  // namespace hpgem

#endif
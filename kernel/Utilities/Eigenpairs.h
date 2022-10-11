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
#ifndef HPGEM_EIGENPAIRS_H
#define HPGEM_EIGENPAIRS_H

#ifdef HPGEM_USE_SLEPC

#include <memory>
#include <vector>
#include <slepc.h>
#include "Logger.h"
#include "../EigenSolvers/JacobiDavidsonMaxwell.h"

namespace hpgem {
namespace Utilities {

/**
 * A list of eigenpairs from solving an eigenvalue problem.
 *
 * Storage and handling of a list of N eigenpairs (lambda_i, x_i) with i=1..N
 * and lambda_i the i-th eigenvalue and x_i the corresponding eigenvector. The
 * ordering of the eigenpairs in the list can be adjusted.
 */
class Eigenpairs final {
   public:
    Eigenpairs();
    Eigenpairs(const Eigenpairs&) = delete;
    Eigenpairs(Eigenpairs&& other) noexcept;
    ~Eigenpairs();

    Eigenpairs& operator=(const Eigenpairs&) = delete;
    Eigenpairs& operator=(Eigenpairs&& other) noexcept;

    /**
     * Reserve space so that at least a certain number of eigenpairs can be
     * stored.
     * @param newSize The desired minimum size
     * @param sample A sample vector to clone as eigenvector
     */
    void reserve(std::size_t newSize, Vec sample);

    /**
     * \brief Load eigenpairs from a solver
     *
     * @param eps The solver to load the eigenpairs from
     * @param sample A sample vector to duplicate for storing eigenvectors
     */
    void loadEigenpairs(EPS& eps, Vec sample);
    void loadEigenpairs(EigenSolvers::JacobiDavidsonMaxwellSolver& eps,
                        Vec sample);

    /**
     * Reorder the currently stored eigenpairs
     *
     * @param ordering A permutation on (0, size()-1). Being the ordering of the
     * current eigenpairs in the new ordering.
     */
    void reorder(std::vector<std::size_t> ordering);

    /**
     * @return Number of eigenpairs stored (== getEigenvalues().size())
     */
    std::size_t size() const { return eigenvalues_.size(); }

    /**
     * Get a specific eigenvector
     * @param i The index of the eigenvector (<size())
     * @return The eigenvector
     */
    Vec getEigenvector(std::size_t i) {
        logger.assert_always(i < eigenvalues_.size(),
                             "Asking for eigenvector % with only % eigenpairs",
                             i, eigenvalues_.size());
        return eigenvectors_[ordering_[i]];
    }

    /**
     * Access to the raw eigenvector storage. No guarantee is made about the
     * ordering of the vectors.
     *
     * @return A pointer to the raw eigenvector storage, having at least size()
     * vectors.
     */
    Vec* getRawEigenvectors() { return eigenvectors_; }

    /**
     * Get a specific eigenvalue
     * @param i  The index of the eigenvalues (<size())
     * @return The eigenvalue
     */
    PetscScalar getEigenvalue(std::size_t i) const {
        logger.assert_always(i < eigenvalues_.size(),
                             "Asking for eigenvalue % with only % eigenpairs",
                             i, eigenvalues_.size());
        // return eigenvalues_[ordering_[i]];
        return eigenvalues_[i];
    }

   private:
    /// Number of vectors available in eigenvectors_, this may be larger than
    /// the number of eigenvalues. However, then only the first
    /// eigenvalue.size() vectors contain valid eigenvectors.
    PetscInt numberOfEigenvectors_;
    /// Storage for eigenvectors
    Vec* eigenvectors_;
    /// Storage for the most recent eigenvalues
    std::vector<PetscScalar> eigenvalues_;
    /// Reordering of the eigenvalues. The i-th entry gives the index of the
    /// i-th eigenpair in eigenvalues_ and eigenvectors_
    std::vector<std::size_t> ordering_;
};
}  // namespace Utilities
}  // namespace hpgem

#endif

#endif  // HPGEM_EIGENPAIRS_H

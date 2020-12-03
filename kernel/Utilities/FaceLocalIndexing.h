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
#ifndef HPGEM_FACELOCALINDEXING_H
#define HPGEM_FACELOCALINDEXING_H

#include "ElementLocalIndexing.h"

namespace hpgem {

namespace Base {
class Face;
}

namespace Utilities {

/**
 * Index of the local DoFs on a Face.
 *
 * For face matrices and vectors the rows and columns correspond to DoFs of a
 * particular (sub)set of unknowns. These are the DoFs which have support on the
 * adjacent element(s). This calss provides easy access to information about the
 * local ordering of these DoFs in matrices and vectors.
 */
class FaceLocalIndexing {
   public:
    explicit FaceLocalIndexing(std::size_t numberOfUnknowns);
    void reinit(const std::vector<std::size_t>& includedUnknowns);
    /// Reinit the values for a different face
    /// \param element The new element (may be nullptr)
    void reinit(const hpgem::Base::Face* face);

    /// Get the offset of the DoFs for a particular unknown at one side of the
    /// face.
    /// \param unknown The unknown
    /// \param side The side of the face
    /// \return The offset of the first DoF for that unknown at that side.
    std::size_t getDoFOffset(std::size_t unknown,
                             hpgem::Base::Side side) const {
        if (side == hpgem::Base::Side::LEFT) {
            return left_.getDoFOffset(unknown);
        } else {
            hpgem::logger.assert_debug(
                right_.getElement() != nullptr,
                "Asking for the right side of an internal face");
            return left_.getNumberOfDoFs() + right_.getNumberOfDoFs(unknown);
        }
    }

    /// Get the number of DoFs for a particular unknown at one side of the face.
    /// \param unknown The unknown
    /// \param side The side of the face.
    /// \return The number of DoFs for the unknown with support on the
    /// particular side of the face.
    std::size_t getNumberOfDoFs(std::size_t unknown,
                                hpgem::Base::Side side) const {
        if (side == hpgem::Base::Side::LEFT) {
            return left_.getNumberOfDoFs(unknown);
        } else {
            hpgem::logger.assert_debug(
                right_.getElement() != nullptr,
                "Asking for the right side of an internal face");
            return right_.getNumberOfDoFs(unknown);
        }
    }

    /// \return The total number of DoFs for all the included unknowns.
    std::size_t getNumberOfDoFs() const {
        return left_.getNumberOfDoFs() + right_.getNumberOfDoFs();
    }

    /// \return The unknowns that are included.
    const std::vector<std::size_t>& getIncludedUnknowns() {
        return left_.getIncludedUnknowns();
    }

    /// \return The face for which this indexing is currently.
    const hpgem::Base::Face* getFace() { return face_; }

   private:
    /// The current face
    const hpgem::Base::Face* face_;
    /// Indexing on the left
    hpgem::Utilities::ElementLocalIndexing left_;
    /// Indexing on the right
    hpgem::Utilities::ElementLocalIndexing right_;
};
}  // namespace Utilities
}  // namespace hpgem

#endif  // HPGEM_FACELOCALINDEXING_H

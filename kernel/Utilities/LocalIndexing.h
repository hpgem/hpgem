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
#ifndef HPGEM_LOCALINDEXING_H
#define HPGEM_LOCALINDEXING_H

#include <vector>
#include "Logger.h"
#include "Base/Side.h"

namespace hpgem {

namespace Base {
class Element;
class Face;
}

namespace Utilities {

class ElementLocalIndexing {
   public:
    /// Create a ElementLocalIndexing without any includedUnknowns.
    /// \param numberOfUnknowns The largest unknown that can be used.
    explicit ElementLocalIndexing(std::size_t numberOfUnknowns);
    /// Update the included unknowns. The unknowns should be unique and less
    /// than the maximum.
    ///
    /// \param includedUnknowns The new set of included unknowns.
    void reinit(const std::vector<std::size_t>& includedUnknowns);
    /// Reinit the values for a different element
    /// \param element The new element (may be nullptr)
    void reinit(const Base::Element* element);

    /// Get the Offset in an ElementMatrix/Vector for the unknown. Will be an
    /// invalid value for non included unknowns or when no element is available.
    std::size_t getDoFOffset(std::size_t unknown) const {
        logger.assert_debug(unknown < numberOfUnknowns_,
                            "Unknown % larger than maximum %", unknown,
                            numberOfUnknowns_);
        return offsets_[unknown];
    }
    /// Get the number DoFs (rows/columns) for the unknown in an element
    /// matrix/vector. Will be zero for non included unknowns or when no element is available.
    std::size_t getNumberOfDoFs(std::size_t unknown) const {
        logger.assert_debug(unknown < numberOfUnknowns_,
                            "Unknown % larger than maximum %", unknown,
                            numberOfUnknowns_);
        return sizes_[unknown];
    }
    /// Total number of DoFs of the included unknowns and thus rows/columns of
    /// an element matrix or vector. Without element this will be zero.
    std::size_t getNumberOfDoFs() const { return totalNumberOfDofs_; }

    /// Get the current element.
    const Base::Element* getElement() const { return element_; }
    /// The unknowns that are included.
    const std::vector<std::size_t>& getIncludedUnknowns() const {
        return includedUnknowns_;
    }

    /// Check the invariants of the datastructure to see if the internal
    /// constraints are satisfied. To be used for testing and debugging.
    void validate() const;

    /// Constructor alternative for an index where all unknowns are included.
    /// \param numberOfUnknowns The total number of unknowns
    /// \return A new instance with all unknowns included
    static ElementLocalIndexing createFullIndexing(std::size_t numberOfUnknowns) {
        std::vector<std::size_t> unknowns(numberOfUnknowns);
        std::iota(unknowns.begin(), unknowns.end(), 0);
        ElementLocalIndexing result(numberOfUnknowns);
        result.reinit(unknowns);
        return result;
    }

   private:
    /// The current element
    const Base::Element* element_;

    /// The maximum number of unknowns.
    std::size_t numberOfUnknowns_;

    /// The total number of DoFs of the included unknowns.
    std::size_t totalNumberOfDofs_;
    /// Vector of size numberOfUnknowns_ with the offsets of each unknown. Will
    /// be -1 if the unknown is not included.
    std::vector<std::size_t> offsets_;
    /// Vector of size numberOfUnknowns_ with the number of DoFs for each
    /// unknown. Will be 0 if the unknown is not included.
    std::vector<std::size_t> sizes_;
    /// The unknowns that are included.
    std::vector<std::size_t> includedUnknowns_;
};

class FaceLocalIndexing {
   public:
    explicit FaceLocalIndexing(std::size_t numberOfUnknowns);
    void reinit(const std::vector<std::size_t>& includedUnknowns);
    /// Reinit the values for a different face
    /// \param element The new element (may be nullptr)
    void reinit(const Base::Face* face);

    std::size_t getDoFOffset(std::size_t unknown, Base::Side side) const {
        if (side == Base::Side::LEFT) {
            return left_.getDoFOffset(unknown);
        } else {
            logger.assert_debug(
                right_.getElement() != nullptr, "Asking for the right side of an internal face");
            return left_.getNumberOfDoFs() + right_.getNumberOfDoFs(unknown);
        }
    }

    std::size_t getNumberOfDoFs(std::size_t unknown, Base::Side side) const {
        if (side == Base::Side::LEFT) {
            return left_.getNumberOfDoFs(unknown);
        } else {
            logger.assert_debug(
                right_.getElement() != nullptr, "Asking for the right side of an internal face");
            return right_.getNumberOfDoFs(unknown);
        }
    }

    std::size_t getNumberOfDoFs() const { return left_.getNumberOfDoFs() + right_.getNumberOfDoFs(); }

    const std::vector<std::size_t>& getIncludedUnknowns() {
        return left_.getIncludedUnknowns();
    }

    const Base::Face* getFace() {
        return face_;
    }

   private:
    const Base::Face* face_;
    ElementLocalIndexing left_;
    ElementLocalIndexing right_;
};

}  // namespace Utilities
}  // namespace hpgem

#endif  // HPGEM_LOCALINDEXING_H

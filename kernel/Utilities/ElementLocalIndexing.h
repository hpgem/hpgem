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
#ifndef HPGEM_ELEMENTLOCALINDEXING_H
#define HPGEM_ELEMENTLOCALINDEXING_H

#include <numeric>
#include <vector>
#include <limits>
#include "Logger.h"
#include "Base/Side.h"

namespace hpgem {

namespace Base {
class Element;
}

namespace Utilities {

/**
 * Index of the local DoFs on an Element.
 *
 * For element matrices and vectors the rows and columns correspond to DoFs of
 * a particular (sub)set of unknowns. Which are the DoFs for these unknowns that
 * have support on the element. This class provides easy access to information
 * about the local ordering of the DoFs in matrices and vectors.
 */
class ElementLocalIndexing {
   public:
    static const std::size_t UNKNOWN_NOT_INCLUDED;

    ElementLocalIndexing();
    /// Update the included unknowns. The unknowns should be unique and less
    /// than the maximum.
    ///
    /// \param includedUnknowns The new set of included unknowns.
    void reinit(const std::vector<std::size_t>& includedUnknowns);
    /// Reinit the values for a different element
    /// \param element The new element (may be nullptr)
    void reinit(const Base::Element* element);

    /// Combined reinit for both the element and included unknowns.
    ///
    /// \param element The new element (may be nullptr)
    /// \param includedUnknowns The new set of included unknowns.
    void reinit(const Base::Element* element,
                const std::vector<std::size_t>& includedUnknowns) {
        reinit(includedUnknowns);
        reinit(element);
    }

    /// Get the Offset in an ElementMatrix/Vector for the unknown.
    ///
    /// Will be UNKNOWN_NOT_INCLUDED for non included unknowns or when no
    /// element is available.
    std::size_t getDoFOffset(std::size_t unknown) const {
        if (unknown >= offsets_.size()) {
            return UNKNOWN_NOT_INCLUDED;
        }
        return offsets_[unknown];
    }
    /// Get the number DoFs (rows/columns) for the unknown in an element
    /// matrix/vector. Will be zero for non included unknowns or when no
    /// element is available.
    std::size_t getNumberOfDoFs(std::size_t unknown) const {

        if (unknown >= sizes_.size()) {
            return 0;
        }
        return sizes_[unknown];
    }
    /// Total number of DoFs of the included unknowns and thus rows/columns of
    /// an element matrix or vector. Without element this will be zero.
    std::size_t getNumberOfDoFs() const {
        std::size_t totalNumberOfDoFs = 0;
        for (std::size_t includedUnknown : includedUnknowns_) {
            totalNumberOfDoFs += sizes_[includedUnknown];
        }
        return totalNumberOfDoFs;
    }

    /// Get the current element.
    const Base::Element* getElement() const { return element_; }
    /// The unknowns that are included.
    const std::vector<std::size_t>& getIncludedUnknowns() const {
        return includedUnknowns_;
    }

    /// Check the invariants of the datastructure to see if the internal
    /// constraints are satisfied. To be used for testing and debugging.
    void validateInternalState() const;

    /// Constructor alternative for an index where all unknowns are included.
    /// \param numberOfUnknowns The total number of unknowns
    /// \return A new instance with all unknowns included
    static ElementLocalIndexing createFullIndexing(
        std::size_t numberOfUnknowns) {
        std::vector<std::size_t> unknowns(numberOfUnknowns);
        std::iota(unknowns.begin(), unknowns.end(), 0);
        ElementLocalIndexing result;
        result.reinit(unknowns);
        return result;
    }

   private:
    /// Clear the state of the vectors relating to the specific unknowns (i.e.
    /// offsets_ and sizes_)
    void clearUnknownVectors();

    /// The current element
    const Base::Element* element_;

    /// Vector of size numberOfUnknowns_ with the offsets of each unknown. Will
    /// be -1 if the unknown is not included.
    std::vector<std::size_t> offsets_;
    /// Vector of size numberOfUnknowns_ with the number of DoFs for each
    /// unknown. Will be 0 if the unknown is not included.
    std::vector<std::size_t> sizes_;
    /// The unknowns that are included.
    std::vector<std::size_t> includedUnknowns_;
};

}  // namespace Utilities
}  // namespace hpgem

#endif  // HPGEM_ELEMENTLOCALINDEXING_H

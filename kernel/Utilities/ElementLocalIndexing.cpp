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
#include "ElementLocalIndexing.h"

#include "Element.h"

namespace hpgem {
namespace Utilities {

ElementLocalIndexing::ElementLocalIndexing()
    : element_(nullptr), offsets_(), sizes_(), includedUnknowns_(0) {}

void ElementLocalIndexing::reinit(
    const std::vector<std::size_t> &includedUnknowns) {
    includedUnknowns_ = includedUnknowns;
    std::sort(includedUnknowns_.begin(), includedUnknowns_.end());
#ifdef HPGEM_ASSERTS
    // Check for duplicates, easy after sorting
    for (std::size_t i = 0; i < includedUnknowns_.size() - 1; ++i) {
        logger.assert_debug(includedUnknowns_[i] != includedUnknowns_[i + 1],
                            "Duplicate unknown");
    }
    // Check maximum and extend if necessary
    if (!includedUnknowns_.empty()) {
        std::size_t maxUnknown =
            includedUnknowns_[includedUnknowns_.size() - 1];
        if (maxUnknown >= sizes_.size()) {
            sizes_.resize(maxUnknown + 1, 0);
            offsets_.resize(maxUnknown + 1, -1);
        }
    }
#endif
    // Update the content
    if (element_ != nullptr) {
        // Clear the storage to the default values as reinit does only overwrite
        // the included unknowns.
        std::fill(offsets_.begin(), offsets_.end(), -1);
        std::fill(sizes_.begin(), sizes_.end(), 0);
        reinit(element_);
    }
}

void ElementLocalIndexing::reinit(const Base::Element *element) {
    element_ = element;
    if (element == nullptr) {
        std::fill(offsets_.begin(), offsets_.end(), -1);
        std::fill(sizes_.begin(), sizes_.end(), 0);
    } else {
        std::size_t totalNumberOfDoFs = 0;
        for (std::size_t unknown : includedUnknowns_) {
            std::size_t numDoFs = element->getNumberOfBasisFunctions(unknown);
            offsets_[unknown] = totalNumberOfDoFs;
            sizes_[unknown] = numDoFs;
            totalNumberOfDoFs += numDoFs;
        }
    }
}

void ElementLocalIndexing::validateInternalState() const {
    logger.assert_always(offsets_.size() == sizes_.size(),
                         "Size and offset tables of different length.");
    logger.assert_always(
        std::is_sorted(includedUnknowns_.begin(), includedUnknowns_.end()),
        "Unsorted included unknowns");
}

}  // namespace Utilities

}  // namespace hpgem
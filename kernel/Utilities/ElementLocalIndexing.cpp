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

ElementLocalIndexing::ElementLocalIndexing(std::size_t numberOfUnknowns)
    : element_(nullptr),
      numberOfUnknowns_(numberOfUnknowns),
      totalNumberOfDofs_(0),
      offsets_(numberOfUnknowns_, -1),
      sizes_(numberOfUnknowns_, 0),
      includedUnknowns_(0) {}

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
    // Check maximum
    if (!includedUnknowns_.empty()) {
        std::size_t maxUnknown =
            includedUnknowns_[includedUnknowns_.size() - 1];
        logger.assert_debug(
            maxUnknown < numberOfUnknowns_,
            "Largest unknown % is larger than maximum allowable unknown %",
            maxUnknown, numberOfUnknowns_);
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
        totalNumberOfDofs_ = 0.0;
    } else {
        totalNumberOfDofs_ = 0;
        for (std::size_t unknown : includedUnknowns_) {
            std::size_t numDoFs = element->getNumberOfBasisFunctions(unknown);
            offsets_[unknown] = totalNumberOfDofs_;
            sizes_[unknown] = numDoFs;
            totalNumberOfDofs_ += numDoFs;
        }
    }
}

void ElementLocalIndexing::validate() const {
    logger.assert_always(offsets_.size() == numberOfUnknowns_,
                         "Incorrectly sized offsets table");
    logger.assert_always(sizes_.size() == numberOfUnknowns_,
                         "Incorrectly sized sizes table");
    logger.assert_always(
        std::is_sorted(includedUnknowns_.begin(), includedUnknowns_.end()),
        "Unsorted included unknowns");
    std::vector<bool> included(numberOfUnknowns_, false);
    for (const auto &unknown : includedUnknowns_) included[unknown] = true;
    for (std::size_t i = 0; i < numberOfUnknowns_; ++i) {
        if (included[i]) {
            logger.assert_always(
                offsets_[i] + sizes_[i] <= totalNumberOfDofs_,
                "Too large offset/size for included dof %: (%,%)", i,
                offsets_[i], sizes_[i]);
        } else {
            logger.assert_always(offsets_[i] == -1,
                                 "Incorrect offset for non included unknown: %",
                                 offsets_[i]);
            logger.assert_always(sizes_[i] == 0,
                                 "Incorrect size for non included unknown: %",
                                 sizes_[i]);
        }
    }
}

}  // namespace Utilities

}  // namespace hpgem
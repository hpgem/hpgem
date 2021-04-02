/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2014, University of Twente
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

#include "ElementReorder.h"
#include <algorithm>
#include <numeric>
#include "Logger.h"

namespace Preprocessor {

const ElementReorder::Element& ElementReorder::FindElement(
    size_t dimension, size_t indices_size) const {
    auto iterator = std::find_if(
        orderPerElement_.begin(), orderPerElement_.end(),
        [&](const Element& e) {
            return e.dimension_ == dimension && e.order_.size() == indices_size;
        });

    return *iterator;
}

void ElementReorder::addElementType(size_t dimension, const std::string& name,
                                    const std::vector<size_t>& order) {
    auto iterator = std::find_if(
        orderPerElement_.begin(), orderPerElement_.end(),
        [&](const Element& e) {
            return e.dimension_ == dimension && e.order_.size() == order.size();
        });

    bool found = (iterator != orderPerElement_.end());

    hpgem::logger.assert_always(
        !found,
        "Element with dimension % and same number of nodes % "
        "already exists. Clear disambiguation not possible.",
        dimension, order.size());

    orderPerElement_.push_back(Element(dimension, name, order));
}

void ElementReorder::reorderToHpGem(size_t dimension,
                                    std::vector<size_t>& indeces) const {
    const Element& e = FindElement(dimension, indeces.size());

    std::vector<size_t> copy = indeces;

    for (size_t i = 0; i < indeces.size(); i++) {
        indeces[i] = copy[e.order_[i]];
    }
}

void ElementReorder::reorderFromHpGem(size_t dimension,
                                      std::vector<size_t>& indeces) const {
    const Element& e = FindElement(dimension, indeces.size());

    std::vector<size_t> copy = indeces;

    for (size_t i = 0; i < indeces.size(); i++) {
        indeces[e.order_[i]] = copy[i];
    }
}

void ElementReorder::Element::checkOrder(std::vector<size_t> order) const {
    std::sort(order.begin(), order.end());
    std::adjacent_difference(order.begin(), order.end(), order.begin());

    bool valid = (order[0] == 0) && (std::accumulate(order.begin(), order.end(),
                                                     0) == order.size() - 1);
    hpgem::logger.assert_always(valid,
                                "Order has to contain 0 as an index and then "
                                "contain all indices contiguously");
}

}  // namespace Preprocessor

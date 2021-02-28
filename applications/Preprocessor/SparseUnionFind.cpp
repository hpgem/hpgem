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
#include "SparseUnionFind.h"

namespace Preprocessor {

std::size_t SparseUnionFind::findSet(std::size_t elem) {
    auto optTarget = parents_.find(elem);
    if (optTarget == parents_.end()) {
        // Singleton set
        return elem;
    } else {
        std::size_t target = findSetPresent(optTarget->second);
        // Same update as in findSetPresent()
        parents_[elem] = target;
        return target;
    }
}

std::size_t SparseUnionFind::findSetPresent(std::size_t elem) {
    // Note reference
    std::size_t& target = parents_[elem];
    if (target != elem) {
        // Note: we update target here (the reason why target is a reference),
        // so that future lookups don't have to recurse down to the node
        target = findSetPresent(target);
    }
    return target;
}

void SparseUnionFind::unionSets(std::size_t elem1, std::size_t elem2) {
    if (elem1 == elem2) {
        return;
    }
    // Make sure that both sets are allocated
    if (parents_.find(elem1) == parents_.end()) {
        parents_[elem1] = elem1;
        ranks_[elem1] = 0;
    }
    if (parents_.find(elem2) == parents_.end()) {
        parents_[elem2] = elem2;
        ranks_[elem2] = 0;
    }
    // Safe to use findSetPresent by previous statements
    std::size_t target1 = findSetPresent(elem1);
    std::size_t target2 = findSetPresent(elem2);
    if (target1 == target2) {
        // Already in the same set
        return;
    }
    std::size_t rank1 = ranks_[target1];
    std::size_t& rank2 = ranks_[target2]; // Reference to allow for update

    if (rank1 > rank2) {
        parents_[target2] = target1;
    } else {
        parents_[target1] = target2;
        if (rank1 == rank2) {
            rank2++;
        }
    }
}

}  // namespace Preprocessor
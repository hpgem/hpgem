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
#ifndef HPGEM_SPARSEUNIONFIND_H
#define HPGEM_SPARSEUNIONFIND_H

#include <map>

namespace Preprocessor {

/**
 * Sparse version of a  union-find data structure
 *
 * A Union-Find data structure that is optimized for a large (unknown) index
 * range where almost all indices correspond to single sets. Thus only very few
 * indices are grouped into larger sets.
 *
 * The Union-Find data-structure represents N items, represented by the indices
 * 0 to N-1, grouped into sets. Each item is part of exactly one set. The
 * datastructure supports two primitive operations:
 *  - findSet(i): Give the representative entry of the set to which i belongs.
 *    Each entry in the same set as i will return the same representative, but
 *    which entry it is unspecified.
 *  - unionSet(i, j): Merge the sets to which i and j belong.
 * and the derived operations:
 *  - isSameSet(i, j): Determine if i and j are in the same set.
 *
 * See for more information:
 * https://en.wikipedia.org/wiki/Disjoint-set_data_structure
 * Implementation inspired by: Competitive Programming 3 (https://cpbook.net/)
 */
class SparseUnionFind {
   public:
    /**
     * Const iterator over the indices for non-singleton-set entries.
     */
    struct iterator {
       public:
        using innerIter_t =
            std::map<std::size_t,
                     std::pair<std::size_t, std::size_t>>::const_iterator;

        // Be a good LegacyInputIterator
        using reference = const std::size_t&;
        using value_type = std::size_t;
        using difference_type = std::size_t;
        using pointer = std::size_t const*;
        using iterator_category = std::forward_iterator_tag;

        iterator() = default;
        iterator(innerIter_t iter) : innerIter(iter){};

        std::size_t operator*() const { return innerIter->first; }
        std::size_t operator++() { return (++innerIter)->first; }
        std::size_t operator++(int) { return (innerIter++)->first; }
        bool operator==(const iterator& other) const {
            return innerIter == other.innerIter;
        }
        bool operator!=(const iterator& other) const {
            return innerIter != other.innerIter;
        }

       private:
        innerIter_t innerIter;
    };

    /**
     * Finds a representative for the set to which the given entity belongs
     * @param elem The sample member of the set.
     * @return The representative of that set.
     */
    std::size_t findSet(std::size_t elem);
    /**
     * Checks whether two elements belong to the same set.
     * @param elem1 First element
     * @param elem2 Second element
     * @return Whether elem1 and elem2 belong to the same set
     */
    bool inSameSet(std::size_t elem1, std::size_t elem2) {
        return findSet(elem1) == findSet(elem2);
    }
    /**
     * Merge the sets to which the two given elements belong.
     *
     * @param elem1 The element of the first set
     * @param elem2 The elemetn of the second set.
     */
    void unionSets(std::size_t elem1, std::size_t elem2);

    /**
     * Begin iterator to the indices of elements that are part of a set with at
     * least two elements.
     * @return The iterator
     */
    iterator begin() const { return iterator(storage_.cbegin()); }
    /**
     * End iterator to the indices of elements that are part of a set with at
     * least two elements.
     * @return The end iterator.
     */
    iterator end() const { return iterator(storage_.cend()); }

   private:
    /**
     * Same as findSet(std::size_t), but with the added requirement that
     * parents_[elem] MUST be present.
     * @param elem The element to lookup
     * @return The representative
     */
    std::size_t findSetPresent(std::size_t elem);

    /**
     * Spare storage of the data. For each element it stores two values
     *  - A parent entry, forming a tree structure for each set. The top entry
     *    of the set/tree points to itself and is the representative for that
     *    set.
     *  - A rank, estimating the depth of the tree under the entry can be. This
     *    allows for automatic balancing of the tree during merges.
     *
     * Sparsity means that only entries are stored that are part of a set of at
     * least two entries (i.e. not a singleton set). Entries that are not
     * included would be of the form: (i, 0), where i is the index.
     */
    std::map<std::size_t, std::pair<std::size_t, std::size_t>> storage_;
};

}  // namespace Preprocessor

#endif  // HPGEM_SPARSEUNIONFIND_H

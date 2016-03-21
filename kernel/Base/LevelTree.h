/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, University of Twente
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
//
#ifndef LevelTree_h
#define LevelTree_h
#include <deque>
#include <vector>
#include <iostream>
#include "Logger.h"
#include "TreeEntry.h"
#include "TreeIterator.h"

namespace Base
{
    
    template<typename V>
    class TreeEntry;

    /**
     * A tree structure that allows iteration over a single level (distance from the root entry) or over all levels
     * in addition to normal pre-order and post-order traversal. Entries may reside in multiple levels. This allows
     * for example locally refining a mesh on the next level without having to duplicate the elements that were not
     * refined. The structure supports having multiple entries on level 0.
     */
    //actually a forest
    template<typename V>
    class LevelTree
    {
    public:
        
        //! construct a new empty level tree
        LevelTree();

        //!clean up the entries (does not do memory management on V)
        ~LevelTree();

        //! false if there are entries in the tree
        bool empty() const;

        //! Number of trees in the LevelTree
        //!\note could also recursively descend the trees to return the size including children if this is more intuitive
        std::size_t size() const;

        //! highest level that contains an element
        std::size_t maxLevel() const;

        //! iterators created after this call will traverse only the indicated level. Does not affect existing iterators
        void setSingleLevelTraversal(std::size_t level);

        //! iterators created after this call will traverse all entries, visiting all entries on lower levels before visiting any entries on higher levels. Does not affect existing iterators
        void setAllLevelTraversal();

        //! iterators created after this call will traverse in pre-order (visiting parents before their children, but visiting children before visiting unvisited siblings). Does not affect existing iterators
        void setPreOrderTraversal();

        //! iterators created after this call will traverse in post-order (visiting parents after their children, but visiting parents before visiting unvisited siblings or their children). Does not affect existing iterators
        void setPostOrderTraversal();

        //! iterators created after this call will use the indicated traversal method. Does not affect existing iterators
        void setTraversalMethod(TreeTraversalMethod method);

        //! returns the level currently being visited for single level traversal, or an arbitrary number if any other traversal method is currently active
        std::size_t getActiveLevel() const;

        //! Getting the beginning of traversal
        TreeIterator<V> begin();

        //! Getting the end of traversal
        TreeIterator<V> end();

        //! Getting the beginning of traversal
        TreeIteratorConst<V> begin() const;

        //! Getting the end of traversal
        TreeIteratorConst<V> end() const;

        //! Getting the reverse iterator to the reverse beginning
        std::reverse_iterator<TreeIterator<V>> rbegin();

        //! Getting the reverse iterator to the reverse end
        std::reverse_iterator<TreeIterator<V>> rend();

        //! Getting the reverse iterator to the reverse beginning
        std::reverse_iterator<TreeIteratorConst<V>> rbegin() const;

        //! Getting the reverse iterator to the reverse end
        std::reverse_iterator<TreeIteratorConst<V>> rend() const;

        //! Add an additional tree to the forest
        TreeIterator<V> addRootEntry(const V& newEl);

        //! Add children to an entry. Provided iterator must point into this tree.
        void addChildren(TreeIteratorConst<V> parentEl, const std::vector<V>& subEntries);

        //! Add children of an entry. Provided iterator must point into this tree.
        TreeIterator<V> addChild(TreeIteratorConst<V> parentEl, const V& subEntry);

        //! Add children of an entry, making the children occupy a specific level. The parent will stretch to fill empty lower levels if necessary. Provided iterator must point into this tree.
        void addChildren(TreeIteratorConst<V> parentEl, const std::vector<V>& subEntries, std::size_t level);

        //! Add children of an entry, making the children occupy a specific level. The parent will stretch to fill empty lower levels if necessary. Provided iterator must point into this tree.
        TreeIterator<V> addChild(TreeIteratorConst<V> parentEl, const V& subEntry, std::size_t level);

        ConstIterableWrapper<TreeEntry<V>> getRootEntries() const;

        //! Make all leaves occupy a specified level, stretching leaves to fill empty lower levels if necessary.
        void fillToLevel(std::size_t level);

        //! Getting the beginning of single level traversal
        TreeIterator<V> beginLevel(std::size_t level);

        //! Getting the beginning of all level traversal
        TreeIterator<V> beginAllLevel();

        //! Getting the beginning of pre-order traversal
        TreeIterator<V> beginPreOrder();

        //! Getting the beginning of post-order traversal
        TreeIterator<V> beginPostOrder();

        //! Erase all descendants of an entry (and reset its depth to 1)
        void eraseChilds(TreeIterator<V> parentEl);

        //! Erase all entries (including descendants)
        void clear();

        //! Describe the LevelTree
        friend std::ostream& operator<<(std::ostream& os, const LevelTree<V>& e)
        {
            os << "LevelTree: ";
            os << "entries_.size()= " << e.entries_.size() << " ";
            os << "minLevel_= " << 0 << " ";
            os << "maxLevel_= " << e.maxLevel_ << " ";
            os << "activeLevel_= " << e.activeLevel_ << " ";
            return os;
        }

    private:
        std::vector<TreeEntry<V>*> entries_;
        std::size_t maxLevel_;
        std::size_t activeLevel_;
        TreeTraversalMethod traversalMethod_;
    };

} // close namespace Base

// merge the implementation file here.
#include "LevelTree_Impl.h"

#endif //LevelTree_h

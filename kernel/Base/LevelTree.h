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

namespace Base
{
    
    template<class V>
    class TreeEntry;
    
    template<class V, class T>
    class TreeIterator;
    
    //actually a forest
    template<class V>
    class LevelTree
    {
        
    public:
        using valueT = V;
        using treeEntryT = TreeEntry<V>;
        //  this is our (rather special) iterator 
        using iterator = TreeIterator<V,TreeEntry<V> >;

        LevelTree();

        ~LevelTree();

        bool empty() const;

        //! Number of entries in the LevelTree
        std::size_t size() const;

        std::size_t maxLevel() const;

        void setSingleLevelTraversal(std::size_t level);

        void setAllLevelTraversal();

        void setPreOrderTraversal();

        void setPostOrderTraversal();

        std::size_t getActiveLevel() const;

        //! Getting the beginning of traversal
        iterator begin();

        //! Getting the end of traversal
        iterator end();

        //! Getting the reverse iterator to the reverse beginning
        std::reverse_iterator<iterator> rbegin();

        //! Getting the reverse iterator to the reverse end
        std::reverse_iterator<iterator> rend();

        //! Add an additional tree to the forest
        iterator addRootEntry(const valueT& newEl);

        //! Add children of an entry.
        iterator addChildren(iterator parentEl, const std::vector<valueT>& subEntries);

        //! Add children of an entry.
        iterator addChild(iterator parentEl, const valueT& subEntries);

        //! Add children of an entry.
        iterator addChildren(iterator parentEl, const std::vector<valueT>& subEntries, std::size_t level);

        //! Add children of an entry.
        iterator addChild(iterator parentEl, const valueT& subEntries, std::size_t level);

        //! Getting the beginning of tree-level traversal
        iterator beginLevel(const int level);

        //! Getting the beginning of leaves traversal
        iterator beginAllLevel();

        //! Getting the beginning of pre-order traversal
        iterator beginPreOrder();

        //! Getting the beginning of post-order traversal
        iterator beginPostOrder();

        //! Erase all descendants of an entry
        void eraseChilds(iterator parentEl);

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
        std::vector<treeEntryT*> entries_;
        int maxLevel_;
        int activeLevel_;
        enum class TraversalMethod
        {
            SINGLELEVEL, ALLLEVEL, PREORDER, POSTORDER
        } traversalMethod_;
    };

} // close namespace Base

// merge the implementation file here.
#include "Base/LevelTree_Impl.h"

#endif //LevelTree_h

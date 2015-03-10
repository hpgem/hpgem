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
//------------------------------------------------------------------------------
// System includes and names imported from them:
#include <deque>
#include <vector>
#include <iostream>
#include "Logger.h"
//------------------------------------------------------------------------------
// Package includes:

//------------------------------------------------------------------------------

namespace Base
{
    
    template<class V>
    class TreeEntry;
    
    template<class V, class T>
    class TreeIterator;
    
    template<class V>
    class LevelTree
    {
        
    public:
        using valueT = V;
        using treeEntryT = TreeEntry<V>;
        //  this is our (rather special) iterator 
        using iterator = TreeIterator<V,TreeEntry<V> >;
        using DimT = unsigned int;

        LevelTree();

        ~LevelTree();

        bool empty() const;

        //! Number of entries in the LevelTree
        int size() const;

        int maxLevel() const;

        //protected:
        void setActiveLevel(const DimT level);

        void resetActiveLevel();

        int getActiveLevel() const;

    public:
        //! Getting the beginning of traversal
        iterator begin();

        //! Getting the end of traversal
        iterator end();

        //! Add new entry
        iterator addEntry(const valueT& newEl, const bool preserveLinks = false);

        //! Add new dummy TreeEntry
        iterator addDummyEntry(iterator it);

        //! Add children of an entry.
        iterator addChildren(iterator parentEl, const std::vector<valueT>& subEntries);

        //! Set the current entries as part of the coarsest mesh.
        void setAsTheCoarsestEntries();

        //! Getting the beginning of tree-level traversal
        iterator beginLevel(const int level);

        //! Getting the beginning of leaves traversal
        iterator beginLeaf();

        //! Getting the beginning of pre-order traversal
        iterator beginPreOrder();

        //! Getting the beginning of post-order traversal
        iterator beginPostOrder();

        //! Erase all descendants of an entry
        void eraseChildsOf(iterator fci);

        //! Describe the LevelTree
        friend std::ostream& operator<<(std::ostream& os, const LevelTree<V>& e)
        {
            os << "LevelTree: ";
            os << "noRootEntries_= " << e.noRootEntries_ << " ";
            os << "entries_.size()= " << e.entries_.size() << " ";
            os << "minLevel_= " << e.minLevel_ << " ";
            os << "maxLevel_= " << e.maxLevel_ << " ";
            os << "activeLevel_= " << e.activeLevel_ << " ";
            os << "coarsestEntriesSet_= " << e.coarsestEntriesSet_ << " ";
            return os;
        }
        
    private:
        //! Add new TreeEntry
        iterator addTreeEntry(treeEntryT* const newEnt, const bool preserveLinks = false);

        //! Add tree-children of an entry.
        iterator addTreeChildren(iterator parentEl, const std::vector<treeEntryT>& subEntries);

        //! Erase a leaf entry
        iterator eraseLastLeaf(iterator fci);

    private:
        int noRootEntries_;
        std::deque<treeEntryT*> entries_;
        int minLevel_;
        int maxLevel_;
        int activeLevel_;
        bool coarsestEntriesSet_;
    };

} // close namespace Base

// merge the implementation file here.
#include "Base/LevelTree_Impl.h"

#endif //LevelTree_h

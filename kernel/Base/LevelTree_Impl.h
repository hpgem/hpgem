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

#ifndef LevelTree_Impl_h
#define LevelTree_Impl_h

#include "Base/LevelTree.h"

namespace Base
{
    template<class V>
    LevelTree<V>::LevelTree()
            : noRootEntries_(0), minLevel_(0), maxLevel_(0), activeLevel_(-1), coarsestEntriesSet_(false)
    {
    }
    
    template<class V>
    LevelTree<V>::~LevelTree()
    {
        logger(VERBOSE, "-------------Deleting LevelTree....");
        while (!entries_.empty())
        {
            delete entries_.back();
            entries_.pop_back();
        }
        logger(VERBOSE, "-------------LevelTree is deleted.");
    }
    
    template<class V>
    bool LevelTree<V>::empty() const
    {
        return entries_.empty();
    }
    
    //! Number of entries in the LevelTree
    template<class V>
    int LevelTree<V>::size() const
    {
        return entries_.size();
    }
    
    template<class V>
    int LevelTree<V>::maxLevel() const
    {
        return maxLevel_;
    }
    
    template<class V>
    void LevelTree<V>::setActiveLevel(const DimT level)
    {
        activeLevel_ = level;
    }
    
    template<class V>
    void LevelTree<V>::resetActiveLevel()
    {
        activeLevel_ = -1;
    }
    
    template<class V>
    int LevelTree<V>::getActiveLevel() const
    {
        return activeLevel_;
    }
    
    //! Setting sequential list traversal
    template<class V>
    typename LevelTree<V>::iterator LevelTree<V>::begin()
    {
        if (activeLevel_ >= 0)
        {
            return beginLevel(activeLevel_);
        }
        else
        {
            iterator fci;
            fci.ptr_ = entries_.begin();
            
            return fci;
        }
    }
    
    //! Setting end of traversal
    template<class V>
    typename LevelTree<V>::iterator LevelTree<V>::end()
    {
        iterator fci;
        fci.ptr_ = entries_.end();
        
        return fci;
    }
    
    //! Add new entry
    template<class V>
    typename LevelTree<V>::iterator LevelTree<V>::addEntry(const valueT& newEl, const bool preserveLinks)
    {
        return addTreeEntry(new treeEntryT(newEl));
    }
    
    template<class V>
    typename LevelTree<V>::iterator LevelTree<V>::addTreeEntry(treeEntryT* const newEnt, const bool preserveLinks)
    {
        // put it on the back of the list
        entries_.push_back(newEnt);
        ++noRootEntries_;
        
        // update some statistics of the LevelTree
        const int level = newEnt->getLevel();
        if (minLevel_ > level)
            minLevel_ = level;
        else if (maxLevel_ < level)
            maxLevel_ = level;
        
        iterator fci;
        fci.ptr_ = --entries_.end(); // it's on the back, dear!
        if (!preserveLinks)
        {
            // set the default links of the new tree node (new entry)
            fci->setParent(fci); // self loop-back for parent
            fci->setSelf(fci); // self iterator
            fci->setChild(fci); // self loop-back for child
        }
        
        return fci;
    }
    
    template<class V>
    typename LevelTree<V>::iterator LevelTree<V>::addDummyEntry(iterator it)
    {
        // update some statistics of the LevelTree
        const int level = (*it)->getBasicLevel() + (*it)->getDepth() - 1;
        if (minLevel_ > level)
            minLevel_ = level;
        else if (maxLevel_ < level)
            maxLevel_ = level;
        
        return it;
    }
    
    //! Set the current entries as the coarsest mesh.  This means that all
    //! existing entries are the roots of the tree.  So, we just need to
    //! update sibling related numbers of the roots.
    template<class V>
    void LevelTree<V>::setAsTheCoarsestEntries()
    {
        if (coarsestEntriesSet_)
            return;
        coarsestEntriesSet_ = true;
        
        if (entries_.empty())
        {
            throw "LevelTree<V>::setAsTheCoarsestEntries() error: the mesh is empty!\n";
        }
        
        int numRoots = 0;
        for (iterator it = begin(); it != end(); ++it)
        {
            if ((*it)->getLevel() == 0)
            {
                (*it)->setSiblingIndex(numRoots);
                ++numRoots;
            }
            else
                break;
        }
        
        for (iterator it = begin(); it != end(); ++it)
        {
            if ((*it)->getLevel() == 0)
            {
                (*it)->setNumSiblings(numRoots);
            }
            else
                break;
        }
        
    }
    
    //! Setting tree-level traversal
    template<class V>
    typename LevelTree<V>::iterator LevelTree<V>::beginLevel(const int level)
    {
        if (!coarsestEntriesSet_)
        {
            throw "LevelTree<V>::beginLevel() error: the coarsest entries has not been set!";
        }
        
        iterator fci;
        
        if ((level < minLevel_) || (level > maxLevel_))
        {
            throw "LevelTree<V>::beginLevel() error: the level must in the valid range!";
            
            fci.end_ = entries_.end();
            fci.ptr_ = fci.end_;
            fci.last_ = fci.end_;
            
            return fci;
        }
        
        // set the traversal level
        fci.setBasicLevelTraversal(level);
        
        // mark the last on the level
        fci.ptr_ = entries_.begin();
        fci->resetDepthCounter();
        fci.lastOnLevel(level);
        fci.last_ = fci.ptr_;
        
        // move to the first on the level
        fci.ptr_ = entries_.begin();
        fci.firstOnLevel(level);
        
        // mark the end of the LevelTree
        fci.end_ = entries_.end();
        
        return fci;
    }
    
    //! Setting up leaves traversal
    template<class V>
    typename LevelTree<V>::iterator LevelTree<V>::beginLeaf()
    {
        if (!coarsestEntriesSet_)
        {
            throw "LevelTree<V>::beginLeaf() error: the coarsest entries has not been set!";
        }
        
        iterator fci;
        
        // set the traversal method to be visiting leaves
        fci.setLeavesTraversal();
        
        // mark the last leaf
        fci.ptr_ = entries_.begin();
        fci.LastLeaf();
        fci.last_ = fci.ptr_;
        
        // move to the first leaf
        fci.ptr_ = entries_.begin();
        fci.first_Leaf();
        fci.first_ = fci.ptr_;
        
        // mark the end of the LevelTree
        fci.end_ = entries_.end();
        
        return fci;
    }
    
    //! Setting up pre-order traversal
    template<class V>
    typename LevelTree<V>::iterator LevelTree<V>::beginPreOrder()
    {
        if (!coarsestEntriesSet_)
        {
            throw "LevelTree<V>::beginPreOrder() error: the coarsest entries has not been set!";
        }
        
        iterator fci;
        
        // set the traversal method to be pre-order
        fci.setPreOrderTraversal();
        
        // move to the first node
        fci.ptr_ = entries_.begin();
        fci.first_ = fci.ptr_;
        
        // mark the end of the LevelTree
        fci.end_ = entries_.end();
        
        // mark the last node
        fci.ptr_ = entries_.begin();
        fci.LastLeaf();
        fci.last_ = fci.ptr_;
        
        return fci;
    }
    
    //! Setting up post-order traversal
    template<class V>
    typename LevelTree<V>::iterator LevelTree<V>::beginPostOrder()
    {
        if (!coarsestEntriesSet_)
        {
            throw "LevelTree<V>::beginPostOrder() error: the coarsest entries has not been set!";
        }
        
        iterator fci;
        
        // set the traversal method to be post-order
        fci.setPostOrderTraversal();
        
        // mark the last node, the top rightmost node
        fci.ptr_ = entries_.begin();
        for (int i = 0; i < fci->getNumSiblings() - 1; ++i, ++(fci.ptr_))
            ; // the last brother
        fci.last_ = fci.ptr_;
        
        // move to the first leaf node
        fci.ptr_ = entries_.begin();
        fci.first_Leaf();
        fci.first_ = fci.ptr_;
        
        // mark the end of the LevelTree
        fci.end_ = entries_.end();
        
        return fci;
    }
    
    template<class V>
    typename LevelTree<V>::iterator LevelTree<V>::addChildren(iterator parentIt, const std::vector<valueT>& subEntries)
    {
        std::vector<treeEntryT> subTreeEntries;
        
        for (typename std::vector<valueT>::const_iterator it = subEntries.begin(); it != subTreeEntries.end(); ++it)
            subEntries.push_back(*(new treeEntryT(*it)));
        
        addTreeChildren(parentIt, subTreeEntries);
    }
    
    //! Add children of an entry.  Note that you are not allowed 
    //! to add more children to the entry later.
    //! \return iterator to the first child is returned
    template<class V>
    typename LevelTree<V>::iterator LevelTree<V>::addTreeChildren(iterator parentIt, const std::vector<treeEntryT>& subEntries)
    {
        if (!coarsestEntriesSet_)
        {
            throw "LevelTree<V>::addTreeChildren() error: the coarsest entries has not been set!";
        }
        
        DimT numSubEntries = subEntries.size();
        if (numSubEntries == 0)
        {
            parentIt->incrDepth();
            addDummyEntry(parentIt);
            return parentIt;
        }
        
        iterator firstChild;
        for (DimT i = 0; i < numSubEntries; ++i)
        {
            treeEntryT te = subEntries[i];
            te->setBasicLevel(parentIt->getLevel() + 1);
            te->setSiblingIndex(i);
            te->setNumSiblings(numSubEntries);
            iterator it = addTreeEntry(te);
            
            (*it)->setParent(parentIt);
            if (i == 0)
            {
                parentIt->setChild(it);
                firstChild = it;
            }
        } // end for
        
        return firstChild;
    }
    
    //! Erase all descendants of an entry
    template<class V>
    void LevelTree<V>::eraseChildsOf(iterator fci)
    {
        if (!coarsestEntriesSet_)
        {
            throw "LevelTree<V>::eraseChildsOf() error: the coarsest entries has not been set!";
        }
        
        if (fci->hasChild())
        {
            iterator it(fci->getChild().ptr_);
            it.LastLeaf();
            it.setPreOrderTraversal();
            do
            {
                eraseLastLeaf((*it)--);
            } while (!((*it)->id() == fci->id()));
        }
    }
    
    //! Erase a leaf entry
    template<class V>
    typename LevelTree<V>::iterator LevelTree<V>::eraseLastLeaf(iterator fci)
    {
        if ((!fci->isLastSibling()) || !(fci->isLeaf()))
        {
            throw "LevelTree<V>::eraseLastLeaf() error: eraseLastLeaf not the last leaf node!";
        }
        
        if (fci->canDecreaseCounter())
        {
            fci.parent();
            fci->decrDepth();
            return fci;
        }
        
        // correct the link from parent to its child
        if (fci->isFirstSibling())
            fci->getParent()->setChild(fci->getParent());
        else
        {
            // update numSiblings info on the undeleted brothers
            iterator bro;
            bro.ptr_ = fci->getParent()->getChild().ptr_;
            int numBro = fci->getNumSiblings() - 1;
            for (int i = 0; i < numBro; ++i, ++bro)
                bro->setNumSiblings(numBro);
        }
        
        iterator it(entries_.erase(fci.ptr_));
        
        return it;
    }

} // close namespace Base

#endif //LevelTree_Impl_h

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

namespace Base
{
    template<typename V>
    LevelTree<V>::LevelTree()
            : maxLevel_(0), activeLevel_(-1), traversalMethod_(TreeTraversalMethod::ALLLEVEL)
    {
    }
    
    template<typename V>
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
    
    template<typename V>
    bool LevelTree<V>::empty() const
    {
        return entries_.empty();
    }
    
    //! Number of trees in the LevelTree
    template<typename V>
    std::size_t LevelTree<V>::size() const
    {
        return entries_.size();
    }
    
    template<typename V>
    std::size_t LevelTree<V>::maxLevel() const
    {
        return maxLevel_;
    }
    
    template<typename V>
    void LevelTree<V>::setSingleLevelTraversal(const std::size_t level)
    {
        logger.assert(maxLevel_ > level, "Trying to iterate over level %, but there are only % levels", level, maxLevel_);
        activeLevel_ = level;
        traversalMethod_ = TreeTraversalMethod::SINGLELEVEL;
    }
    
    template<typename V>
    void LevelTree<V>::setAllLevelTraversal()
    {
        activeLevel_ = -1;
        traversalMethod_ = TreeTraversalMethod::ALLLEVEL;
    }

    template<typename V>
    void LevelTree<V>::setPreOrderTraversal()
    {
        traversalMethod_ = TreeTraversalMethod::PREORDER;
    }

    template<typename V>
    void LevelTree<V>::setPostOrderTraversal()
    {
        traversalMethod_ = TreeTraversalMethod::POSTORDER;
    }

    template<typename V>
    void LevelTree<V>::setTraversalMethod(TreeTraversalMethod method)
    {
        logger.assert(method != TreeTraversalMethod::SINGLELEVEL, "need to provide a level to traverse");
        traversalMethod_ = method;
    }
    
    template<typename V>
    std::size_t LevelTree<V>::getActiveLevel() const
    {
        return activeLevel_;
    }
    
    //! Setting sequential list traversal
    template<typename V>
    TreeIterator<V> LevelTree<V>::begin()
    {
        if(empty()) {
            return TreeIterator<V>();
        }
        TreeIterator<V> result = entries_.front()->getIterator(traversalMethod_, activeLevel_);
        if(traversalMethod_ == TreeTraversalMethod::POSTORDER)
        {
            while((*result.ptr_)->hasChild())
            {
                result.moveToChild(0);
            }
        }
        return result;
    }
    
    //! Setting end of traversal
    template<typename V>
    TreeIterator<V> LevelTree<V>::end()
    {
        if(empty()) {
            return TreeIterator<V>();
        }
        //it is easier to first get the last valid position
        TreeIterator<V> result = entries_.back()->getIterator(traversalMethod_, activeLevel_);
        if(traversalMethod_ == TreeTraversalMethod::PREORDER || traversalMethod_ == TreeTraversalMethod::ALLLEVEL)
        {
            while((*result.ptr_)->hasChild())
            {
                result.moveToChild((*result.ptr_)->getNumberOfChildren() - 1);
            }
        }
        if(traversalMethod_ == TreeTraversalMethod::SINGLELEVEL)
        {
            result.moveToLastOnLevel(activeLevel_);
        }
        //and then increment once
        ++result;
        return result;
    }

    //! Setting sequential list traversal
    template<typename V>
    TreeIteratorConst<V> LevelTree<V>::begin() const
    {
        if(empty()) {
            return TreeIterator<V>();
        }
        TreeIterator<V> result = entries_.front()->getIterator(traversalMethod_, activeLevel_);
        if(traversalMethod_ == TreeTraversalMethod::POSTORDER)
        {
            while((*result.ptr_)->hasChild())
            {
                result.moveToChild(0);
            }
        }
        return result;
    }

    //! Setting end of traversal
    template<typename V>
    TreeIteratorConst<V> LevelTree<V>::end() const
    {
        if(empty()) {
            return TreeIterator<V>();
        }
        //it is easier to first get the last valid position
        TreeIterator<V> result = entries_.back()->getIterator(traversalMethod_, activeLevel_);
        if(traversalMethod_ == TreeTraversalMethod::PREORDER || traversalMethod_ == TreeTraversalMethod::ALLLEVEL)
        {
            while((*result.ptr_)->hasChild())
            {
                result.moveToChild((*result.ptr_)->getNumberOfChildren() - 1);
            }
        }
        if(traversalMethod_ == TreeTraversalMethod::SINGLELEVEL)
        {
            result.moveToLastOnLevel(activeLevel_);
        }
        //and then increment once
        ++result;
        return result;
    }

    template<typename V>
    std::reverse_iterator<TreeIterator<V>> LevelTree<V>::rbegin()
    {
        return std::reverse_iterator<TreeIterator<V>>(end());
    }

    template<typename V>
    std::reverse_iterator<TreeIterator<V>> LevelTree<V>::rend()
    {
        return std::reverse_iterator<TreeIterator<V>>(begin());
    }

    template<typename V>
    std::reverse_iterator<TreeIteratorConst<V>> LevelTree<V>::rbegin() const
    {
        return std::reverse_iterator<TreeIteratorConst<V>>(end());
    }

    template<typename V>
    std::reverse_iterator<TreeIteratorConst<V>> LevelTree<V>::rend() const
    {
        return std::reverse_iterator<TreeIteratorConst<V>>(begin());
    }
    
    //! Add new entry
    template<typename V>
    TreeIterator<V> LevelTree<V>::addRootEntry(const V& newEl)
    {
        maxLevel_ = std::max(1UL, maxLevel_);
        entries_.push_back(new TreeEntry<V>(newEl));
        entries_.back()->setSiblings(entries_.size() - 1, &entries_);
        return entries_.back()->getIterator(traversalMethod_);
    }

    template<typename V>
    void LevelTree<V>::addChildren(TreeIterator<V> parentEl, const std::vector<V>& subEntries)
    {
        (*parentEl.ptr_)->addChildren(subEntries);
        maxLevel_ = std::max((*parentEl.ptr_)->getLastChild()->getDepth() + (*parentEl.ptr_)->getLastChild()->getLevel(), maxLevel_);
    }

    template<typename V>
    TreeIterator<V> LevelTree<V>::addChild(TreeIterator<V> parentEl, const V& subEntry)
    {
        (*parentEl.ptr_)->addChild(subEntry);
        maxLevel_ = std::max((*parentEl.ptr_)->getLastChild()->getDepth() + (*parentEl.ptr_)->getLastChild()->getLevel(), maxLevel_);
        return (*parentEl.ptr_)->getLastChild()->getIterator(traversalMethod_);
    }

    template<typename V>
    void LevelTree<V>::addChildren(TreeIterator<V> parentEl, const std::vector<V>& subEntries, std::size_t level)
    {
        logger.assert((*parentEl.ptr_)->getLevel() < level, "trying to at children at level %, but the parent lives at level %", level, (*parentEl.ptr_)->getLevel());
        (*parentEl.ptr_)->setDepth(level - (*parentEl.ptr_)->getLevel());
        (*parentEl.ptr_)->addChildren(subEntries);
        maxLevel_ = std::max((*parentEl.ptr_)->getLastChild()->getDepth() + (*parentEl.ptr_)->getLastChild()->getLevel(), maxLevel_);
    }

    template<typename V>
    TreeIterator<V> LevelTree<V>::addChild(TreeIterator<V> parentEl, const V& subEntry, std::size_t level)
    {
        logger.assert((*parentEl.ptr_)->getLevel() < level, "trying to at children at level %, but the parent lives at level %", level, (*parentEl.ptr_)->getLevel());
        (*parentEl.ptr_)->setDepth(level - (*parentEl.ptr_)->getLevel());
        (*parentEl.ptr_)->addChild(subEntry);
        maxLevel_ = std::max((*parentEl.ptr_)->getLastChild()->getDepth() + (*parentEl.ptr_)->getLastChild()->getLevel(), maxLevel_);
        return (*parentEl.ptr_)->getLastChild()->getIterator(traversalMethod_);
    }

    template<typename V>
    ConstIterableWrapper<TreeEntry<V>> LevelTree<V>::getRootEntries() const
    {
        return ConstIterableWrapper<TreeEntry<V>>{entries_};
    }

    template<typename V>
    void LevelTree<V>::fillToLevel(std::size_t level)
    {
        maxLevel_ = std::max(level + 1, maxLevel_);
        TreeTraversalMethod active = traversalMethod_;
        traversalMethod_ = TreeTraversalMethod::POSTORDER;
        //no range based loop because I have to modify 'inside' the iterator
        for(TreeIterator<V> iterator = begin(); iterator != end(); ++iterator)
        {
            if((*iterator.ptr_)->isLeaf() && (*iterator.ptr_)->getLevel() + (*iterator.ptr_)->getDepth() - 1 < level)
            {
                (*iterator.ptr_)->setDepth(level - (*iterator.ptr_)->getLevel() + 1);
            }
        }
        traversalMethod_ = active;
    }
    
    //! Setting tree-level traversal
    template<typename V>
    TreeIterator<V> LevelTree<V>::beginLevel(std::size_t level)
    {
        logger.assert(maxLevel_ > level, "Trying to iterate over level %, but there are only % levels", level, maxLevel_);
        return entries_.front()->getIterator(TreeTraversalMethod::SINGLELEVEL, level);
    }
    
    //! Setting up leaves traversal
    template<typename V>
    TreeIterator<V> LevelTree<V>::beginAllLevel()
    {
        return entries_.front()->getIterator(TreeTraversalMethod::ALLLEVEL);
    }
    
    //! Setting up pre-order traversal
    template<typename V>
    TreeIterator<V> LevelTree<V>::beginPreOrder()
    {
        return entries_.front()->getIterator(TreeTraversalMethod::PREORDER);
    }
    
    //! Setting up post-order traversal
    template<typename V>
    TreeIterator<V> LevelTree<V>::beginPostOrder()
    {
        TreeIterator<V> result = entries_.front()->getIterator(TreeTraversalMethod::POSTORDER);
        while((*result.ptr_)->hasChild())
        {
            result.moveToChild(0);
        }
        return result;
    }
    
    //! Erase all descendants of an entry
    template<typename V>
    void LevelTree<V>::eraseChilds(TreeIterator<V> parent)
    {
        (*parent.ptr_)->removeChildren();
        TreeIterator<V> testIterator = entries_.front()->getIterator(TreeTraversalMethod::ALLLEVEL);
        //can probably become more efficient
        while(!testIterator.moveDownToLevelBegin(maxLevel_ - 1))
        {
            testIterator = entries_.front()->getIterator(TreeTraversalMethod::ALLLEVEL);
            --maxLevel_;
        }
    }

    template<typename V>
    void Base::LevelTree<V>::clear()
    {
        while (!entries_.empty())
        {
            delete entries_.back();
            entries_.pop_back();
        }
        maxLevel_ = 0;
    }

} // close namespace Base

#endif //LevelTree_Impl_h

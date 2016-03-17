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
#ifndef TreeEntry_h
#define TreeEntry_h
#include <typeinfo>

#include <iostream>
#include "Logger.h"
#include "ConstIterableWrapper.h"

namespace Base
{
    enum class TreeTraversalMethod
    {
        SINGLELEVEL, ALLLEVEL, PREORDER, POSTORDER
    };

    template<typename T>
    struct PointerToConst {
        using type = const T;
    };

    template<typename T>
    struct PointerToConst<T*> {
        using type = const typename PointerToConst<T>::type*;
    };

    template<typename V, bool isConst>
    class TreeIterator;
    
    //! data structure that allows LevelTree to provide the tree behavior it promises. It is responsible for its children, but not for its siblings or its parent
    template<typename V>
    class TreeEntry
    {
    public:
        //! create a new entry without siblings
        TreeEntry(const V& data, std::size_t level = 0, std::size_t depth = 1)
            : data_(data), level_(level), depth_(depth), parent_(nullptr), siblingIndex_(0), siblings_(nullptr)
        {
        }

        TreeEntry(const TreeEntry<V>& other) = delete;

        //! clean up, don't assume all siblings still exist
        ~TreeEntry()
        {
            for(TreeEntry* child : children_)
            {
                delete child;
            }
        }
        
        //! get whatever is stored in this tree
        V& getData()
        {
            return data_;
        }
        
        typename PointerToConst<V>::type getData() const
        {
            return data_;
        }

        operator V()
        {
            return data_;
        }

        //! create an iterator with desired traversal method that points to this entry. Second argument is used only when single level traversal is requested
        TreeIterator<V, true> getIterator(TreeTraversalMethod type, std::size_t singleLevelLevel) const
        {
            return TreeIterator<V, true>(siblings_->begin() + siblingIndex_, type, singleLevelLevel);
        }

        //! create an iterator with desired traversal method that points to this entry. Second argument is used only when single level traversal is requested
        TreeIterator<V, true> getIterator(TreeTraversalMethod type) const
        {
            return TreeIterator<V, true>(siblings_->begin() + siblingIndex_, type);
        }

        //! create an iterator with desired traversal method that points to this entry. Second argument is used only when single level traversal is requested
        TreeIterator<V, false> getIterator(TreeTraversalMethod type, std::size_t singleLevelLevel)
        {
            return TreeIterator<V, false>(siblings_->begin() + siblingIndex_, type, singleLevelLevel);
        }

        //! create an iterator with desired traversal method that points to this entry. Second argument is used only when single level traversal is requested
        TreeIterator<V, false> getIterator(TreeTraversalMethod type)
        {
            return TreeIterator<V, false>(siblings_->begin() + siblingIndex_, type);
        }

        //! get the lowest level that this entry resides in
        std::size_t getLevel() const
        {
            return level_;
        }

        //! true if this entry resides in level 0
        bool isRoot() const
        {
            return level_ == 0;
        }

        //! set the lowest level this entry resides in. Will not alter the number of levels this entry resides in
        void setBasicLevel(std::size_t level)
        {
            level_ = level;
        }

        //! get the number of levels this entry resides in.
        std::size_t getDepth() const
        {
            return depth_;
        }

        //! set the number of levels this entry resides in.
        void setDepth(std::size_t depth)
        {
            logger.assert(depth > 0, "a node cannot occupy no levels");
            depth_ = depth;
        }

        //! modify the number of levels this entry resides in, increasing or decreasing by the specified amount
        void shiftDepth(int increment)
        {
            logger.assert(-increment < depth_, "a node cannot occupy no levels");
            depth_ += increment;
        }

        //! get the position of this entry in the list of siblings
        std::size_t getSiblingIndex() const
        {
            return siblingIndex_;
        }

        //! get the current number of siblings
        std::size_t getNumberOfSiblings() const
        {
            return siblings_->size();
        }

        //! get the current number of siblings
        std::vector<TreeEntry*>& getSiblings()
        {
            return *siblings_;
        }

        //! get the current number of siblings
        ConstIterableWrapper<TreeEntry> getSiblings() const
        {
            return ConstIterableWrapper<TreeEntry>(*siblings_);
        }

        //! true if this element is at the beginning of the list of siblings
        bool isFirstSibling() const
        {
            return siblingIndex_ == 0;
        }

        //! true if this element is at the end of the list of siblings
        bool isLastSibling() const
        {
            return siblingIndex_ + 1 == siblings_->size();
        }

        //! give this entry a new set of siblings (and a new position in the new set)
        void setSiblings(std::size_t siblingIndex, std::vector<TreeEntry*>* siblings)
        {
            logger.assert((*siblings)[siblingIndex] == this, "this sibling is not in the indicated position in the vector of siblings");
            siblingIndex_ = siblingIndex;
            siblings_ = siblings;
        }

        //! get the entry that has this entry as one of its children
        TreeEntry* getParent()
        {
            logger.assert(!this->isRoot(), "root has no parents");
            logger.assert(parent_ != nullptr, "no parent has been set for this entry");
            return parent_;
        }

        //! get the entry that has this entry as one of its children
        const TreeEntry* getParent() const
        {
            logger.assert(!this->isRoot(), "root has no parents");
            logger.assert(parent_ != nullptr, "no parent has been set for this entry");
            return parent_;
        }

        //! set the entry that has this entry as one of its children
        void setParent(TreeEntry* parent)
        {
            parent_ = parent;
        }

        //! get the child at the specified position in the list
        TreeEntry* getChild(std::size_t childIndex)
        {
            logger.assert(childIndex < children_.size(), "asked for the %th child, but there are only % children", childIndex, children_.size());
            return children_[childIndex];
        }

        //! get the child at the specified position in the list
        const TreeEntry* getChild(std::size_t childIndex) const
        {
            logger.assert(childIndex < children_.size(), "asked for the %th child, but there are only % children", childIndex, children_.size());
            return children_[childIndex];
        }

        std::vector<TreeEntry*>& getChildren()
        {
            return children_;
        }

        ConstIterableWrapper<TreeEntry> getChildren() const
        {
            return ConstIterableWrapper<TreeEntry>(children_);
        }

        //! get the child at index 0
        TreeEntry* getFirstChild()
        {
            logger.assert(children_.size() > 0, "cannot get to the children of a leaf");
            return children_.front();
        }

        //! get the child at index 0
        const TreeEntry* getFirstChild() const
        {
            logger.assert(children_.size() > 0, "cannot get to the children of a leaf");
            return children_.front();
        }

        //! get the child at the end of the list
        TreeEntry* getLastChild()
        {
            logger.assert(children_.size() > 0, "cannot get to the children of a leaf");
            return children_.back();
        }

        //! get the child at the end of the list
        const TreeEntry* getLastChild() const
        {
            logger.assert(children_.size() > 0, "cannot get to the children of a leaf");
            return children_.back();
        }

        //! false if this is a leaf
        bool hasChild() const
        {
            return children_.size() > 0;
        }

        std::size_t getNumberOfChildren() const
        {
            return children_.size();
        }

        //! false if this entry has children
        bool isLeaf() const
        {
            return !hasChild();
        }

        //! add a new child to the list
        void addChild(const V& newChild)
        {
            children_.push_back(new TreeEntry(newChild, level_ + depth_));
            getLastChild()->setSiblings(getNumberOfChildren() - 1, &children_);
            getLastChild()->setParent(this);
        }

        //! add a lot of new children to the list
        void addChildren(const std::vector<V>& newChildren)
        {
            children_.reserve(newChildren.size() + children_.size());
            for(const V& child : newChildren)
            {
                addChild(child);
            }
        }

        //! remove all children, cleaning up their memory
        void removeChildren()
        {
            for(TreeEntry<V>* child : children_)
            {
                delete child;
            }
            children_.resize(0);
            depth_ = 1;
        }

    private:
        V data_;

        std::size_t level_;
        std::size_t depth_;
        std::size_t siblingIndex_;
        std::vector<TreeEntry*>* siblings_;
        std::vector<TreeEntry*> children_;
        TreeEntry* parent_;
    };
} // close namespace Base

#endif // TreeEntry_h

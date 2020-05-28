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
#ifndef TreeIterator_h
#define TreeIterator_h
#include <vector>
#include <iostream>

#include "Logger.h"
#include "TreeEntry.h"
//------------------------------------------------------------------------------
namespace Base {

namespace Detail {
template <bool doAdd, typename V>
struct maybeAddConst {
    using type = const V;
    using iterator = typename std::vector<TreeEntry<V>*>::const_iterator;
};

template <typename V>
struct maybeAddConst<false, V> {
    using type = V;
    using iterator = typename std::vector<TreeEntry<V>*>::iterator;
};
}  // namespace Detail

template <typename V>
class LevelTree;

// for some routines a faster, more elaborate implementation might be available.
// Recommended to first test if the current implementation is fast enough for
// low overhead in most simulations
template <typename V, bool isConst = false>
class TreeIterator {
   public:
    // typedefs expected by the standard c++ library (do not rename to conform
    // to code standards)
    using difference_type = std::ptrdiff_t;
    using value_type = V;
    using pointer = typename Detail::maybeAddConst<isConst, V>::type*;
    using reference = typename Detail::maybeAddConst<isConst, V>::type&;
    // can go forward or back, but cant skip
    using iterator_category = std::bidirectional_iterator_tag;

    // levelTree is able to step iterators to their expected positions while
    // ignoring the ordering
    friend LevelTree<V>;

    // TreeIteratorConst may construct by directly reading variables from
    // non-const TreeIterator
    friend TreeIterator<V, true>;

    //! the standard c++ library expects iterators to be default constructible
    TreeIterator() {
        // by default just do everything
        ptr_ = std::vector<TreeEntry<V>*>().begin();
        end_ = std::vector<TreeEntry<V>*>().end();
        traversalMethod_ = TreeTraversalMethod::ALLLEVEL;
        depthCounter_ = 0;
    }

    //! actual constructor; may not directly be used to construct a past-the-end
    //! iterator
    TreeIterator(typename std::vector<TreeEntry<V>*>::iterator position,
                 TreeTraversalMethod method, std::size_t singleLevelLevel) {
        ptr_ = position;
        end_ = (*position)->getSiblings().end();
        traversalMethod_ = method;
        depthCounter_ = 0;
        // shift to the correct level if we only iterate part of the tree
        if (method == TreeTraversalMethod::SINGLELEVEL) {
            setSingleLevelTraversal(singleLevelLevel);
        }
    }

    //! allow to construct an iterator on a specific entry with single level
    //! traversal, even when it is on the wrong level
    TreeIterator(typename std::vector<TreeEntry<V>*>::iterator position,
                 TreeTraversalMethod method) {
        ptr_ = position;
        end_ = (*position)->getSiblings().end();
        traversalMethod_ = method;
        depthCounter_ = 0;
    }

    //! Copy constructor
    TreeIterator(const TreeIterator& i) {
        ptr_ = i.ptr_;  // current position
        end_ = i.end_;
        traversalMethod_ = i.traversalMethod_;
        depthCounter_ = i.depthCounter_;
    }

    //! allow implicit non-const->const conversion (template magic to make sure
    //! the copy constructor is used for the non-const tree-iterator)
    template <bool otherConst,
              typename = typename std::enable_if<!otherConst>::type>
    TreeIterator(const TreeIterator<V, otherConst>& i) {
        ptr_ = i.ptr_;  // current position
        end_ = i.end_;
        traversalMethod_ = i.traversalMethod_;
        depthCounter_ = i.depthCounter_;
    }

    //! Move constructor
    TreeIterator(TreeIterator&& i) {
        ptr_ = i.ptr_;  // current position
        end_ = i.end_;
        traversalMethod_ = i.traversalMethod_;
        depthCounter_ = i.depthCounter_;
    }

    //! Copy Assignment operator
    TreeIterator& operator=(const TreeIterator& i) {
        ptr_ = i.ptr_;  // current position
        end_ = i.end_;
        traversalMethod_ = i.traversalMethod_;
        depthCounter_ = i.depthCounter_;
        return *this;
    }

    //! Move Assignment operator
    TreeIterator& operator=(TreeIterator&& i) {
        ptr_ = i.ptr_;  // current position
        end_ = i.end_;
        traversalMethod_ = i.traversalMethod_;
        depthCounter_ = i.depthCounter_;
        return *this;
    }

    // NEVER write a non-const version of this routine. (non-const) Tree-entry
    // allows you to do nasty things with your tree, such as create cycles or
    // break parent/child symmetry
    const TreeEntry<V>* getTreeEntry() const {
        logger.assert_debug(ptr_ != end_,
                            "cannot read from past-the-end position");
        return *ptr_;
    }

    reference operator*() {
        logger.assert_debug(ptr_ != end_,
                            "cannot read from past-the-end position");
        return (*ptr_)->getData();
    }

    pointer operator->() const {
        logger.assert_debug(ptr_ != end_,
                            "cannot read from past-the-end position");
        return &(*ptr_)->getData();
    }

    TreeIterator& operator++() {
        logger.assert_debug(ptr_ != end_,
                            "illegal to increment an iterator when it is "
                            "already equal to end()");
        switch (traversalMethod_) {
            case TreeTraversalMethod::ALLLEVEL:
                moveToNextMultiLevel();
                break;

            case TreeTraversalMethod::PREORDER:
                moveToNextPreOrder();
                break;

            case TreeTraversalMethod::POSTORDER:
                moveToNextPostOrder();
                break;

            default:  // SINGLELEVEL
                moveToNextOnLevel();
                break;
        }

        return (*this);
    }

    TreeIterator& operator--() {
        // no convenient check for begin available validity will be checked in
        // called routines
        switch (traversalMethod_) {
            case TreeTraversalMethod::ALLLEVEL:
                moveToPreviousMultiLevel();
                break;

            case TreeTraversalMethod::PREORDER:
                moveToPreviousPreOrder();
                break;

            case TreeTraversalMethod::POSTORDER:
                moveToPreviousPostOrder();
                break;

            default:  // SINGLELEVEL
                moveToPreviousOnLevel();
                break;
        }

        return (*this);
    }

    TreeIterator operator++(int)  // postinc
    {
        logger.assert_debug(ptr_ != end_,
                            "illegal to increment an iterator when it is "
                            "already equal to end()");
        TreeIterator tmp(*this);
        ++*this;
        return (tmp);
    }

    TreeIterator operator--(int)  // postdec
    {
        // no convenient check for begin available validity will be checked in
        // called routines
        TreeIterator tmp(*this);
        --*this;
        return (tmp);
    }

    //! Are they the same iterators? (ignore const-ness of the other iterator)
    template <bool otherConst>
    bool operator==(const TreeIterator<V, otherConst>& i) const {
        // technically comparing two iterators to two different vectors is not
        // guaranteed to give meaningful results (if any) so has to be slightly
        // more elaborate
        if (ptr_ != end_ && traversalMethod_ == i.traversalMethod_ &&
            i.ptr_ != i.end_) {
            return *ptr_ == *i.ptr_;
        } else {
            // a == b implies ++a == ++b, and --a == --b implies a == b, but the
            // reverses are not required and having one end() for all
            // treeTraversalMethods is easier
            return ptr_ == end_ && i.ptr_ == i.end_;
        }
    }

    //! Are they different iterators?
    template <bool otherConst>
    bool operator!=(const TreeIterator<V, otherConst>& i) const {
        return !(*this == i);
    }

    //! set traversal method to pre order (also skips to the indicated level)
    void setSingleLevelTraversal(std::size_t level) {
        traversalMethod_ = TreeTraversalMethod::SINGLELEVEL;
        if (ptr_ == end_) {
            --ptr_;
            moveToLastOnLevel(level);
            moveToNextOnLevel();
        } else if (level < (*ptr_)->getLevel() + depthCounter_) {
            moveUpToLevel(level);
        } else {
            moveDownToLevelBegin(level);
        }
    }

    //! set traversal method to leaves traversal
    void setAllLevelTraversal() {
        traversalMethod_ = TreeTraversalMethod::ALLLEVEL;
        if (ptr_ == end_) {
            // if the iterator is in the past-the-end position, update to the
            // new past-the-end position
            --ptr_;
            while ((*ptr_)->hasChild()) {
                moveToChild((*ptr_)->getNumberOfChildren() - 1);
            }
            ++ptr_;
        }
    }

    //! set traversal method to pre order
    void setPreOrderTraversal() {
        traversalMethod_ = TreeTraversalMethod::PREORDER;
        if (ptr_ == end_) {
            // if the iterator is in the past-the-end position, update to the
            // new past-the-end position
            --ptr_;
            while ((*ptr_)->hasChild()) {
                moveToChild((*ptr_)->getNumberOfChildren() - 1);
            }
            ++ptr_;
        }
    }

    //! set traversal method to post order
    void setPostOrderTraversal() {
        traversalMethod_ = TreeTraversalMethod::POSTORDER;
        if (ptr_ == end_) {
            // if the iterator is in the past-the-end position, update to the
            // new past-the-end position
            --ptr_;
            while (!(*ptr_)->isRoot()) {
                moveToParent();
            }
            ++ptr_;
        }
    }

   private:
    bool canDecreaseCounter() const { return depthCounter_ > 0; }

    bool canIncreaseCounter() const {
        return !(ptr_ == end_) && depthCounter_ < (*ptr_)->getDepth() - 1;
    }

    void increaseCounter() {
        logger.assert_debug(canIncreaseCounter(),
                            "Tree iterator was asked to increase the depth "
                            "counter, while this is currently impossible");
        ++depthCounter_;
    }

    void decreaseCounter() {
        logger.assert_debug(canDecreaseCounter(),
                            "Tree iterator was asked to decrease the depth "
                            "counter, while this is currently impossible");
        --depthCounter_;
    }

    void setMaxDepthCounter() {
        depthCounter_ = (ptr_ == end_ ? 0 : (*ptr_)->getDepth() - 1);
    }

    //! Move the iterator to the parent
    void moveToParent() {
        logger.assert_debug(
            ptr_ != end_,
            "cannot move to the parent of the past the end iterator");
        if (canDecreaseCounter()) {
            decreaseCounter();
        } else {
            logger.assert_debug(!(*ptr_)->isRoot(),
                                "cannot move to the parent of the root node");
            *this = (*ptr_)->getParent()->getIterator(traversalMethod_);
            setMaxDepthCounter();
        }
        end_ = (*ptr_)->getSiblings().end();
    }

    //! Move the iterator to a child
    void moveToChild(std::size_t iChild = 0) {
        logger.assert_debug(
            ptr_ != end_,
            "cannot move to the children of the past the end iterator");
        if (canIncreaseCounter()) {
            increaseCounter();
        } else {
            logger.assert_debug((*ptr_)->hasChild(),
                                "cannot move to the child of a leaf node");
            *this = (*ptr_)->getChild(iChild)->getIterator(traversalMethod_);
            // depthcounter is automatically 0
        }
        end_ = (*ptr_)->getSiblings().end();
    }

    //! find the first node at a specified level
    void moveToFirstOnLevel(std::size_t level) {
        if (ptr_ == end_) {
            --ptr_;
        }
        moveUpToLevel(0);
        // there might be a forest
        *this = (*ptr_)->getSiblings()[0]->getIterator(traversalMethod_);
        moveDownToLevelBegin(level);
    }

    //! find the last node at a specified level
    void moveToLastOnLevel(std::size_t level) {
        if (ptr_ == end_) {
            --ptr_;
        }
        moveUpToLevel(0);
        *this = (*ptr_)
                    ->getSiblings()[(*ptr_)->getNumberOfSiblings() - 1]
                    ->getIterator(traversalMethod_);
        moveDownToLevelEnd(level);
    }

    void moveUpToLevel(std::size_t level) {
        logger.assert_debug(ptr_ != end_,
                            "cannot move up from the past-the-end iterator");
        while ((*ptr_)->getLevel() + depthCounter_ > level) {
            moveToParent();
        }
    }

    // might move to siblings if it can't find the appropriate level among
    // descendants
    bool moveDownToLevelBegin(std::size_t level) {
        while (!(ptr_ == end_) &&
               ((*ptr_)->getLevel() + depthCounter_) < level) {
            moveToNextPreOrder();
        }
        return !(ptr_ == end_);
    }

    // might move to siblings if it can't find the appropriate level among
    // descendants
    bool moveDownToLevelEnd(std::size_t level) {
        // don't rely on the existence of --begin()
        while (ptr_ == end_ || ((*ptr_)->getLevel() + depthCounter_ < level &&
                                hasPreviousPostOrder())) {
            moveToPreviousPostOrder();
        }
        return (*ptr_)->getLevel() + depthCounter_ == level;
    }

    //! move to next node on the same level
    void moveToNextOnLevel() {
        logger.assert_debug(ptr_ != end_,
                            "can only increment dereferenceable iterators");
        // if we can't move to the next, this should be the last one
        TreeIterator backup = *this;
        std::size_t level = (*ptr_)->getLevel() + depthCounter_;
        while ((*ptr_)->isLastSibling() && !(*ptr_)->isRoot()) {
            // move up until we can move to the next
            moveToParent();
        }
        if ((*ptr_)->isLastSibling()) {
            // if we still can't move to the next entry, restore the old pointer
            // so we can move back to the last on the level with --
            *this = backup;
            ++ptr_;
        } else {
            ++ptr_;
            depthCounter_ = 0;
            if (!moveDownToLevelBegin(level)) {
                // moveDownToLevelBegin will automatically try to move to the
                // next sibling and/or up to another parent if it fails initially
                // so if it fails this was the last on the level
                *this = backup;
                ++ptr_;
            } else {
                end_ = (*ptr_)->getSiblings().end();
            }
        }
    }

    //! move to previous node on the same level
    void moveToPreviousOnLevel() {
        if (ptr_ == end_) {
            // the increment operator made sure we end up at the back
            --ptr_;
        } else {
            //--begin() is undefined so we don't have to care about restoring
            //the current position
            std::size_t level = (*ptr_)->getLevel() + depthCounter_;
            while ((*ptr_)->isFirstSibling()) {
                // move up until we can go to the previous
                logger.assert_debug(
                    !(*ptr_)->isRoot(),
                    "cannot move to previous from the first element");
                moveToParent();
            }
            --ptr_;
            depthCounter_ = 0;
            // please keep as assert_always, routine needs to run and must
            // return true
            logger.assert_always(
                moveDownToLevelEnd(level),
                "cannot move to previous from the first element");
            end_ = (*ptr_)->getSiblings().end();
        }
    }

    bool hasPreviousOnLevel() {
        if (ptr_ == end_ || !(*ptr_)->isFirstSibling()) {
            return true;
        } else {
            std::size_t level = (*ptr_)->getLevel() + depthCounter_;
            TreeIterator duplicate = *this;
            while (!(*duplicate.ptr_)->isRoot() &&
                   (*duplicate.ptr_)->isFirstSibling()) {
                duplicate.moveToParent();
            }
            return !(*duplicate.ptr_)->isFirstSibling() &&
                   (--duplicate).moveDownToLevelEnd(level);
        }
    }

    //! move to next node in the pre-order traversal
    void moveToNextPreOrder() {
        logger.assert_debug(ptr_ != end_,
                            "can only increment dereferenceable iterators");
        if ((*ptr_)->hasChild() || canIncreaseCounter()) {
            moveToChild(0);
        } else {
            while (!(*ptr_)->isRoot() && (*ptr_)->isLastSibling()) {
                moveToParent();
            }
            if ((*ptr_)->isLastSibling()) {
                // move back to the final node before moving to end_
                while ((*ptr_)->hasChild()) {
                    moveToChild((*ptr_)->getNumberOfChildren() - 1);
                }
            }
            ++ptr_;
            depthCounter_ = 0;
        }
    }

    //! move to previous node in the pre-order traversal
    void moveToPreviousPreOrder() {
        if (ptr_ == end_ ||
            (!(*ptr_)->isFirstSibling() && depthCounter_ == 0)) {
            --ptr_;
            while ((*ptr_)->hasChild()) {
                moveToChild((*ptr_)->getNumberOfChildren() - 1);
            }
        } else {
            logger.assert_debug(
                !(*ptr_)->isRoot() || depthCounter_ > 0,
                "may not decrement the first entrty of a container");
            moveToParent();
        }
    }

    //! move to next node in the post-order traversal
    void moveToNextPostOrder() {
        logger.assert_debug(ptr_ != end_,
                            "can only increment dereferenceable iterators");
        if (!(*ptr_)->isLastSibling() && depthCounter_ == 0) {
            ++ptr_;
            while ((*ptr_)->hasChild()) {
                moveToChild(0);
            }
        } else if ((*ptr_)->isRoot() && depthCounter_ == 0) {
            ++ptr_;
        } else {
            moveToParent();
        }
    }

    //! move to previous node in the post-order traversal
    void moveToPreviousPostOrder() {
        // we have the function anyway, might as well use it
        logger.assert_debug(
            hasPreviousPostOrder(),
            "decrementing an iterator to the beginning of a range is undefined "
            "behaviour, please don't do that.");
        if (ptr_ != end_ && ((*ptr_)->hasChild() || canIncreaseCounter())) {
            moveToChild((*ptr_)->getNumberOfChildren() - 1);
        } else if (ptr_ == end_) {
            --ptr_;
        } else {
            while ((*ptr_)->isFirstSibling()) {
                moveToParent();
            }
            --ptr_;
        }
    }

    bool hasPreviousPostOrder() {
        if (ptr_ == end_ || (*ptr_)->hasChild() || !(*ptr_)->isFirstSibling() ||
            canIncreaseCounter()) {
            return true;
        } else {
            TreeIterator duplicate = *this;
            while (!(*duplicate.ptr_)->isRoot() &&
                   (*duplicate.ptr_)->isFirstSibling()) {
                duplicate.moveToParent();
            }
            return !(*duplicate.ptr_)->isFirstSibling();
        }
    }

    void moveToNextMultiLevel() {
        logger.assert_debug(
            ptr_ != end_,
            "iterator has to be dereferenceable to be incrementable");
        moveToNextOnLevel();
        if (ptr_ == end_) {
            --ptr_;
            std::size_t level = (*ptr_)->getLevel() + depthCounter_ + 1;
            // will move to end_ again if this was the last level
            moveToFirstOnLevel(level);
        }
    }

    void moveToPreviousMultiLevel() {
        if (hasPreviousOnLevel()) {
            moveToPreviousOnLevel();
        } else {
            std::size_t level = (*ptr_)->getLevel() + depthCounter_;
            logger.assert_debug(level > 0,
                                "Cannot move to the previous from the first "
                                "entry in the range");
            moveToLastOnLevel(level - 1);
        }
    }

   private:
    typename Detail::maybeAddConst<isConst, V>::iterator ptr_;  // current
                                                                // position
    typename Detail::maybeAddConst<isConst, V>::iterator end_;  // past the end
                                                                // iterator for
                                                                // the current
                                                                // node

    TreeTraversalMethod traversalMethod_;

    std::size_t depthCounter_;
};

template <typename V>
using TreeIteratorConst = TreeIterator<V, true>;
}  // namespace Base

#endif  // TreeIterator_h

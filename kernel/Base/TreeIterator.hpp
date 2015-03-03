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
#ifndef TreeIterator_hpp
#define TreeIterator_hpp
#include <deque>
#include <iostream>

//------------------------------------------------------------------------------
namespace Base
{
    // some forward declaration
    template <class V>
    class LevelTree;

    template <class V, class T>
    class TreeIterator
    {
        friend class LevelTree<V>;
        
    public:
        using valueType = T*;
        using reference = T&;
        using valueVType = V;
        using referenceV = V&;
        using iterator_ref = TreeIterator&;
        using difference_type = int;

        TreeIterator()
            {
                // on default, the iterator behaves as a std::deque::iterator
                traversalMethod_ = traversalList_;
            }

        //! Copy constructor
        TreeIterator(const TreeIterator& i)
            {
                if (this!=&i) 
                {
                    ptr_ = i.ptr_;     // current position
                    end_ = i.end_;     // end position
                    first_ = i.first_; // the last valid position
                    last_ = i.last_;   // the last valid position
                    traversalMethod_ = i.traversalMethod_;
                }
                else
                  std::cout << "Gotcha!  TreeIterator trap.\n";
            }

        //! Copy constructor
        TreeIterator(const typename std::deque<valueType>::iterator& li)
            {      
                ptr_=li;

                // on default, the iterator behaves as a std::deque::iterator
                traversalMethod_ = traversalList_;
            }

        //! Assignment operator
        TreeIterator& operator=(const TreeIterator& i)
            {
                if (this!=&i) 
                {
                    ptr_ = i.ptr_;     // current position
                    end_ = i.end_;     // end position
                    first_ = i.first_; // the last valid position
                    last_ = i.last_;   // the last valid position
                    traversalMethod_ = i.traversalMethod_;
                }
                
                return *this;
            }
        
        //! Return the object reference
//         reference operator*() const
//             {      
//                 return (**ptr_);
//             }
        referenceV operator*() const
            {      
                return ((**ptr_).getData());
            }


        //! Return the object instance
        valueType operator->() const
            {      
                return (*ptr_);
            }

        //! Preincrement iterator
        iterator_ref operator++()
            {
                switch (traversalMethod_)
                {
                    case traversalList_:
                        ++ptr_;
                        break;

                    case traversalTreeLeaves_:
                        if (ptr_ == last_)
                            ptr_ = end_;
                        else if (ptr_ == end_)
                            ptr_ = first_;
                        else
                            NextLeaf();
                        break;

                    case traversalTreePreOrder_:
//                         if ((ptr_ == last_) && (!(*this)->hasChild()))
                        if ((ptr_ == last_) && (!(*this)->hasChild()))
                            ptr_ = end_;
                        else if (ptr_ == end_)
                            ptr_ = first_;
                        else
                            NextPreOrder();
                        break;

                    case traversalTreePostOrder_:
                        if (ptr_ == last_)
                            ptr_ = end_;
                        else if (ptr_ == end_)
                            ptr_ = first_;
                        else
                            NextPostOrder();
                        break;

                    default: // TraversalTreeOnLevel
                        if (ptr_ == last_)
                            ptr_ = end_;
                        else
                            NextOnLevel();
                        break;
                }
                
                return (*this);
            }

        //! Predecrement iterator
        iterator_ref operator--()
            {      
                switch (traversalMethod_)
                {
                    case traversalList_:
                        --ptr_;
                        break;

                    case traversalTreeLeaves_:
                        break;  // not implemented yet

                    case traversalTreePreOrder_:
                        if (ptr_ == first_)
                            ptr_ = end_;
                        else if (ptr_ == end_)
                            ptr_ = last_;
                        else
                            PreviousPreOrder();
                        break;

                    case traversalTreePostOrder_:
                        break;  // not implemented yet

                    default: // TraversalTreeOnLevel
                        break;  // not implemented yet
                }

                return (*this);
            }

        TreeIterator operator++(int) // postinc
            {
                TreeIterator tmp(*this);
                ++*this;
                return (tmp);
            }

        TreeIterator operator--(int) // postdec
            {
                TreeIterator tmp(*this);
                --*this;
                return (tmp);
            }

        //! Move the iterator to the parent
        iterator_ref parent()
            {      
                if ((*this)->canDecreaseCounter())
                    (*this)->decreaseCounter();
                else
                {
                    ptr_ = (*this)->getParent().ptr_;
                    (*this)->setCounterMax();
                }

                return (*this);
            }

        //! Move the iterator to a child
        iterator_ref child(const int iChild = 0)
            {
                if ((*this)->canIncreaseCounter())
                {
                    (*this)->increaseCounter();
                }
                else
                {
                    ptr_ = (*this)->getChild().ptr_;
                    (*this)->setCounterMin();
                }

                if (iChild < (*this)->getNumSiblings())
                   for (int i=0; i < iChild-1; ++i, ++ptr_);
                else std::cerr << "iterator::child() exceed NumSiblings\n";  // error message!

                return (*this);
            }

        //! Are they the same iterators?
        bool operator==(const TreeIterator &i) const
            {      
                return (ptr_ == i.ptr_);
            }
        
        //! Are they different iterators?
        bool operator!=(const TreeIterator &i) const
            { 
                return (ptr_ != i.ptr_);
            }

        //! set traversal method to pre order
        void setBasicLevelTraversal(const int level)
            {
                traversalMethod_ = level;
            }

        //! set traversal method to leaves traversal
        void setLeavesTraversal()
            {
                traversalMethod_ = traversalTreeLeaves_;
            }

        //! set traversal method to pre order
        void setPreOrderTraversal()
            {
                traversalMethod_ = traversalTreePreOrder_;
            }

        //! set traversal method to post order
        void setPostOrderTraversal()
            {
                traversalMethod_ = traversalTreePostOrder_;
            }


    private:

        //! find the first node at a specified level
        iterator_ref firstOnLevel(const int level)
            {
                for (int i=0; i<(*this)->getNumSiblings(); ++i, ++ptr_)
                {
                   if ((*this)->getLevel() >= level) break;  // success

                   if (!(*this)->hasChild()) // no childs
                   {
                        // the last brother, failed to find the level
                        if (i == (*this)->getNumSiblings()-1) break; 
                   }
                   else
                   {
                        child();  // try to find it below
                        firstOnLevel(level);  // recursive call
                        if ((*this)->getLevel() >= level) break;  // success
                        parent(); // can't find it below

                        // the last brother, failed to find the level
                        if (i == (*this)->getNumSiblings()-1) break;
                   }
                }

                return (*this);
            }

        //! find the last node at a specified level
        iterator_ref lastOnLevel(const int level)
            {
                // move to the last brother
                int numBro = (*this)->getNumSiblings();
                for (int i=(*this)->getSiblingIndex(); i<numBro-1; ++i, ++ptr_);
                
                // loop over all brothers
                for (int i=(*this)->getNumSiblings()-1; i>=0; --i, --ptr_)  
                {
                   if ((*this)->getLevel() >= level) break;  // success

                   if (!(*this)->hasChild()) // no childs
                   {
                        // the first brother, failed to find the level
                        if (i == 0) break;
                   }
                   else
                   {
                        child();  // try to find it below
                        lastOnLevel(level);  // recursive call
                        if ((*this)->getLevel() >= level) break;  // success
                        parent(); // can't find if below

                        // the first brother, failed to find the level
                        if (i == 0) break;
                   }
                }

                return (*this);
            }

        //! move to next node on the same level
        iterator_ref NextOnLevel()
            {
                if ((*this)->getLevel() == traversalMethod_)
                {
                    if (!(*this)->isLastSibling())
                    {
                        // move right if I can
                        ++ptr_;
                        (*this)->setDepthCounter((*this)->getLevel() - traversalMethod_ + 1);
                    }
                    else
                    {
                        // can't move right, then move upward until the top level
                        // or until I can move right at some level
                        do {
                            parent();
                        } while ((!(*this)->isRoot()) && (*this)->isLastSibling());

                        if (!(*this)->isLastSibling()) 
                        // move right if I can
                        {
                            ++ptr_;
                            if ((*this)->getLevel() != traversalMethod_)
                                NextOnLevel();
                        }
                    }
                }
                else if ((*this)->getLevel() < traversalMethod_)  // I am above the specified level)
                {
                    // move downward until reaching the specified level or until I can't move downward
                    do {
                        child();
                    } while (((*this)->getLevel() < traversalMethod_) && (*this)->hasChild());

                    // if I reach the specified level
                    if ((*this)->getLevel() == traversalMethod_)
                    {
                        return (*this);
                    }

                    if (!(*this)->isLastSibling())
                    // move right if I can
                    {
                        ++ptr_;
                        if ((*this)->getLevel() != traversalMethod_)
                            NextOnLevel();
                    }
                    else
                    {
                        // can't move right, then move upward until the top level
                        // or until I can move right at some level
                        do {
                            parent();
                        } while ((!(*this)->isRoot()) && (*this)->isLastSibling());

                        if (!(*this)->isLastSibling()) 
                        // move right if I can
                        {
                            ++ptr_;
                            if ((*this)->getLevel() != traversalMethod_)
                                NextOnLevel();
                        }
                    }
                }
                else // I am below the specified level
                {
                    // move upward until reaching the specified level
                    do {
                        parent();
                    } while ((*this)->getLevel() > traversalMethod_);
                }

                return (*this);
            }

        //! move to previous node on the same level
        iterator_ref PreviousOnLevel()
            {
                if ((*this)->getLevel() == traversalMethod_)
                {
                    if (!(*this)->isfirst_Sibling())
                    {
                        // move left if I can
                        --ptr_;
                    }
                    else
                    {
                        // can't move left, then move upward until the top level
                        // or until I can move right at some level
                        do {
                            parent();
                        } while ((!(*this)->isRoot()) && (*this)->isfirst_Sibling());

                        if (!(*this)->isfirst_Sibling()) 
                        {
                            // move left if I can
                            --ptr_;
                            PreviousOnLevel();
                        }
                    }
                }
                else if ((*this)->getLevel() < traversalMethod_)  // I am above the specified level)
                {
                    // not implemented yet
                }
                else // I am below the specified level
                {
                    // not implemented yet
                }

                return (*this);
            }

        //! find the first leaf node
        iterator_ref firstLeaf()
            {
                if ((*this)->hasChild()) 
                {   // find it below
                    child();  
                    firstLeaf();  // recursive call
                }

                return (*this);
            }


        //! find the last leaf node
        iterator_ref lastLeaf()
            {
                // move to the last brother
                for (int i=(*this)->getSiblingIndex(); i<(*this)->getNumSiblings()-1; ++i, ++ptr_);

                if ((*this)->hasChild()) 
                {   // find it below
                    child();
                    lastLeaf();  // recursive call
                }

                return (*this);
            }

        //! move to next leaf node
        iterator_ref NextLeaf()
            {
                // find it on the right, if I can move right
                if (!(*this)->isLastSibling())
                {
                    ++ptr_;  // move right
                    firstLeaf();
                }
                else
                {
                    // can't move right, then move upward until the top level
                    // or until I can move right at some point
                    do {
                        parent();
                    } while ((!(*this)->isRoot()) && (*this)->isLastSibling());

                    // find it on the right, if I can move right
                    if (!(*this)->isLastSibling()) 
                    {
                        ++ptr_;  // move right
                        firstLeaf();
                    }
                }

                return (*this);
            }

        //! move to next node in the pre-order traversal
        iterator_ref NextPreOrder()
            {
                if ((*this)->hasChild())
                    child();   // visit the first child
                else
                {
                    // visit the right node, if I can move right
                    if (!(*this)->isLastSibling()) 
                        ++ptr_;  // move right
                    else
                    {
                        // can't move right, then move upward until the top level
                        // or until I can move right at some point
                        do {
                            parent();
                        } while ((!(*this)->isRoot()) && (*this)->isLastSibling());
                        
                        // visit the right node, if I can move right
                        if (!(*this)->isLastSibling())
                            ++ptr_;  // move right
                    }
                }

                return (*this);
            }

        //! move to previous node in the pre-order traversal
        iterator_ref PreviousPreOrder()
            {
                if (!(*this)->isfirst_Sibling())
                {
                    --ptr_;  // move left
                    if ((*this)->hasChild())
                    {
                        child();
                        lastLeaf();
                    }
                }
                else
                    parent();   // visit the parent

                return (*this);
            }

        //! move to next node in the post-order traversal
        iterator_ref NextPostOrder()
            {
                // visit the right node, if I can move right
                if (!(*this)->isLastSibling()) 
                {
                    ++ptr_;  // visit the right
                    firstLeaf();
                }
                else
                    if (!(*this)->isRoot()) parent();

                return (*this);
            }

        //! move to previous node in the post-order traversal
        iterator_ref PreviousPostOrder()
            {
                // not implemented yet, sorry!
                return (*this);
            }

    private:
        typename std::deque<valueType>::iterator ptr_;   // current position
        typename std::deque<valueType>::iterator end_;   // end position
        typename std::deque<valueType>::iterator first_; // the last valid position
        typename std::deque<valueType>::iterator last_;  // the last valid position

        int traversalMethod_;
        enum {  traversalList_            = -1,
                traversalTreeLeaves_      = -2,
                traversalTreePreOrder_    = -3,
                traversalTreePostOrder_   = -4   };

    };
} // close namespace Base

#endif  // TreeIterator_hpp

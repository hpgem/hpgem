//--------------------------------------
// Copyright 2010 M.T. Julianto
//
#ifndef TreeBase_hpp
#define TreeBase_hpp

#include "TreeIterator.hpp"

namespace Base
{
    //! This class gives important tree properties to its derived class.
    /*!
        Some important properties of a node tree are defined:
        - All nodes of the tree always have one link, which is an iterator, to 
          their \b parent.  Parent of a root node is the node itself.
        - A node may have a \b child, children, or none.  All nodes keep a link 
          to its first child, if any.  If a node has no child, the child link 
          points to itself.
        - A node may have \b siblings or none.  A node is connected to its 
          siblings by C++ STL \c std::list connections.
        - A node may reside on one \b level or on multi level.  We define a
          rather special here, a node reside on multi level, to save memory.
          To enable this special property, we introduce \e basic \e level, 
          \e depth and \e depth \e counter.
        - \b Depth means the number of levels where the node resides on.
        - The value of \b depth \b counter together with the value of \b basic 
          \b level denote at which level the node actually resides.
          
        
     */
    template <class Iterator>
    class TreeBase
    {
    public:
        //! The constructor. 
        TreeBase() 
                : level_(0), siblingIndex_(0), numSiblings_(1), 
                  depth_(1), depthCounter_(1), dummy_(false)
            { }

        //! The constructor. 
        TreeBase(const int level, const int iBro, const int nBro) 
                : level_(level), siblingIndex_(iBro), numSiblings_(nBro), 
                  depth_(1), depthCounter_(1), dummy_(false) 
            { }

        Iterator getIterator() const
            {
                return self_;
            }

        //! Set self iterator.
        /*!
            \param
                parent an iterator to the parent
         */
        void setSelf(Iterator self) 
            { 
                self_ = self;
            }

        //! Get the actual level. 
        /*! 
            \return the actual level.

            The actual level is the basic level + depth counter - 1. 
            Its value is started from 0.
         */
        unsigned int getLevel() const
            { 
                return level_ + depthCounter_ - 1; 
            }

        //! Get the actual level. 
        /*! 
            \return the actual level.

            The actual level is the basic level + depth counter - 1. 
            Its value is started from 0.
         */
        unsigned int getHLevel() const
            { 
                return getLevel(); 
            }

        //! Get basic level. 
        /*!
            \return the basic level.
         */
        unsigned int getBasicLevel() const
            { 
                return level_; 
            }

        //! Set basic level. 
        /*!
            \param level a int value
            Please note the actual level is the basic level + depth counter - 1.
         */
        void setBasicLevel(const unsigned int level) 
            { 
                level_ = level; 
            }

        //! Get sibling index. 
        /*!
            The value is in range 0..numSiblings-1
         */
        unsigned int getSiblingIndex() const
            { 
                return depthCounter_==1 ? siblingIndex_ : 1; 
            }

        //! Set sibling index.  
        /*!
           The value must be in range 0..numSiblings-1
            
         */
        void setSiblingIndex(const unsigned int SiblingIndex) 
            { 
                siblingIndex_ = SiblingIndex; 
            }

        //! Get number of siblings. 
        /*!
            \returns 
                the number of siblings on the current level.

            For multilevel entries, the number of siblings is 1 for deeper level, i.e. 
                
         */
        unsigned int getNumSiblings() const
            { 
                return depthCounter_==1 ? numSiblings_ : 1; 
            }

        //! Set number of siblings. 
        void setNumSiblings(const unsigned int NumSiblings) 
            { 
                numSiblings_ = NumSiblings; 
            }

        //! Increase number of siblings. 
        /*!
            
         */
        void incrNumSiblings(const unsigned int incr = 1) 
            { 
                numSiblings_ += incr;
            }

        //! Decrease number of siblings. 
        /*!
            
         */
        void decrNumSiblings(const unsigned int decr = 1) 
            { 
                numSiblings_ -= decr;
                numSiblings_ = (numSiblings_ < 0 ? 0 : numSiblings_);
            }

        //! Is this the first sibling? 
        /*!
            \returns TRUE if this is the first sibling, 
            otherwise FALSE will be returned.
         */
        bool isFirstSibling() const
            { 
                return depthCounter_==1 ? (siblingIndex_==0) : true;
            }

        //! Is this the last sibling? 
        /*!
            \returns TRUE if this is the last sibling, 
            otherwise FALSE will be returned.
         */
        bool isLastSibling() const
            { 
                return depthCounter_==1 ? siblingIndex_==numSiblings_-1 : true;
            }

        //! Set number of multiple levels of this entry.
        /*!
            \param Depth a positive integer value

            Set the depth of this entry.
         */
        void setDepth(const unsigned int Depth) 
            { 
                depth_ = Depth; 
            }

        //! Get depth. 
        /*!
            
         */
        unsigned int getDepth() const
            { 
                return depth_; 
            }

        //! Increase depth value. 
        /*!
            
         */
        unsigned int incrDepth(const unsigned int incr = 1)
            { 
                depth_+=incr; 
                return depth_;
            }

        //! Decrease depth value. 
        /*!
            
         */
        unsigned int decrDepth(const unsigned int decr = 1)
            { 
                depth_-=decr;
                return depth_;
            }

        //! Set depth counter. 
        /*!
           Depth counter contributes to the actual level and is used for 
           moving up-down or parent-child. 
           Its value must be in range 1..depth
         */
        void setDepthCounter(const unsigned int DepthCounter) 
            { 
                if (DepthCounter <= depth_)
                    depthCounter_ = DepthCounter; 
                else
                    depthCounter_ = depth_;
            }

        //! Set depth counter to 1 which is its minimal value. 
        /*!
           Depth counter contributes to the actual level and is used for 
           moving up-down or parent-child. 
           Its value must be in range 1..depth
            
         */
        void resetDepthCounter() 
            { 
                depthCounter_ = 1; 
            }

        //! Get depth counter.  
        /*!
            Depth counter is used for moving up-down or parent-child.  
            Its value is in range 1..depth.
         */
        unsigned int getDepthCounter() const
            { 
                return depthCounter_; 
            }

        //! Increase depth counter value by \e incr.
        /*!
            \param incr a positif integer value.  Default value of param incr is 1.

            \returns the increased depth counter value.
         */
        unsigned int increaseCounter(const unsigned int incr = 1)
            { 
                if (depthCounter_+incr <= depth_) depthCounter_+=incr; 
                return depthCounter_;
            }

        //! Decrease depth counter value by decr.
        /*!
            \param decr a positif integer value.  Default value of decr is 1.

            \returns the decreased depth counter value.
         */
        unsigned int decreaseCounter(const unsigned int decr = 1)
            { 
                if (depthCounter_-decr >= 0) depthCounter_-=decr; 
                return depthCounter_;
            }

        //! Set depth counter value to its minimum value
        void setCounterMin()
            { 
                depthCounter_ = 1; 
            }

        //! Set depth counter value to its maximum value
        void setCounterMax()
            { 
                depthCounter_ = depth_; 
            }

        //! Can depth counter be increased?
        /*!
            \returns TRUE if depth counter can be increased, i.e. still smaller
                than depth, otherwise FALSE will be returned.
         */
        bool canIncreaseCounter() const
            { 
                return depthCounter_ < depth_;
            }

        //! Can depth counter be decreased?
        /*!
            \returns TRUE if depth counter can be decreased, i.e. still bigger 
                than zero, otherwise FALSE will be returned.
         */
        bool canDecreaseCounter() const
            { 
                return depthCounter_ > 1; 
            }

        //! Set parent.
        /*!
            \param
                parent an iterator to the parent
         */
        void setParent(Iterator parent) 
            { 
                parent_ = parent;
            }

        //! Get parent.
        /*!
            \returns
                iterator to the parent
         */
        Iterator getParent() const
            { 
                return parent_; 
            }

        //! Set first child.
        /*!
            \param child an iterator to the first child
         */
        void setChild(Iterator child) 
            { 
                child_ = child; 
            }

        //! Get first child.
        /*!
            \returns iterator to the first child
         */
        Iterator getChild(const unsigned int iChild = 0) const
            {
                if (depthCounter_ < depth_)
                    {  // return myself
                        if (getLevel()==0)
                            return parent_;
                        else
                            return this->child_;
                    }
                else
                {
                    Iterator it(child_);
                    for (unsigned int i=0; i<iChild; ++i, ++it);
                    return it;
                }
            }

        //! Is this on the top of the tree?
        /*!
            \returns TRUE if this is the root of the tree, 
                otherwise FALSE will be returned.
         */
        bool isRoot() const
            { 
                return (level_==0) && (depthCounter_==1);
            }

        //! Does this has a child?
        /*!
            \returns TRUE if this has a child, 
                otherwise FALSE will be returned.
         */
        bool hasChild() const
            {
                if (canIncreaseCounter())
                    return true;
                else
                    return (self_ != child_);
            }

        //! Is this a leaf of the tree?
        /*!
            \returns TRUE if this is a leaf of the tree, 
                otherwise FALSE will be returned.
         */
        bool isLeaf() const
            { 
                return ((self_ == child_) && 
                        (!canIncreaseCounter())); 
            }

        //! Is this a dummy node?
        /*!
            \returns TRUE if this is a dummy, 
                otherwise FALSE will be returned.
         */
        bool isDummy() const
            { 
                return dummy_; 
            }
            
        //! Set this as a dummy node.
        void setDummy()
            { 
                dummy_ = true;
            }
            
        //! Set this as a dummy node.
        void setDummy(const bool dummyVal)
            { 
                dummy_ = dummyVal;
            }
            
        //! Unset this as a dummy node.
        void unsetDummy() 
            { 
                dummy_ = false;
            }
            
        void describe()
        {
            std::cout << "\t" << level_ << "\t" << siblingIndex_ << "\t" << numSiblings_;
            std::cout << "\t" << depth_ << "\t" << depthCounter_ << std::endl;
        }
        

    private:
        unsigned int level_;
        unsigned int siblingIndex_;
        unsigned int numSiblings_;
        unsigned int depth_;
        unsigned int depthCounter_;
        bool dummy_;
        Iterator self_;
        Iterator parent_;
        Iterator child_;
    };
} // close namespace Base

#endif // TreeBase_hpp


//--------------------------------------
// Copyright 2010 M.T. Julianto
//
#ifndef TreeEntry_hpp
#define TreeEntry_hpp
#include <typeinfo>

#include <iostream>
#include "TreeBase.hpp"
#include "TreeIterator.hpp"

namespace Base
{
    template <class V>
    class TreeEntry : public Base::TreeBase<Base::TreeIterator<V,TreeEntry<V> > >
    {
    public:
        typedef TreeEntry<V> Type;
        typedef TreeEntry<V>* PtrType;
        typedef Base::TreeIterator<V,TreeEntry<V> > Iterator;

        //! Constructor
        TreeEntry(const V& data)
              : data_(data), 
                refinementType_(-1), 
                justRefined_(false),
                referenceCounter_(0)
            { 
//               std::cout << "TreeEntry is created\n"; 
            }
            
        //! Copy constructor.
        //! We need this shallow copy for duplicating a mesh.  As the 
        //! compensation, we implement reference counter (referenceCounter_) to avoid 
        //! double freeing the memory
        TreeEntry(const TreeEntry<V>& entry)
            : data_(entry.data_),
              refinementType_(entry.refinementType_), 
              justRefined_(entry.justRefined_),
              referenceCounter_(entry.referenceCounter_+1)
            {
                Base::TreeBase<Iterator>::setBasicLevel(entry.getLevel());
                Base::TreeBase<Iterator>::setDepth(entry.getDepth());
                Base::TreeBase<Iterator>::setDepthCounter(entry.getDepthCounter());
                Base::TreeBase<Iterator>::setNumSiblings(entry.getNumSiblings());
                Base::TreeBase<Iterator>::setSiblingIndex(entry.getSiblingIndex());
                Base::TreeBase<Iterator>::setDummy(entry.isDummy());
                
//                 std::cout << "TreeEntry is copied\n";
            }

        ~TreeEntry() 
        { 
//           std::cout << "TreeEntry is deleted\n"; 
        }

        //! Return the object reference
        V& operator*()
            {      
                return data_;
            }

    public:
        V data_;
        int  refinementType_;    // this element is a result of previous refinement of this type 
        bool justRefined_;       // a flag for marking unfinished refinement proccess
        int  toBeRefined_;       // type of refinement to be applied
        int  referenceCounter_;  // counter for reference an original element:
                                 //    0= the original
                                 //    >0 a shallow copy of an original element
    private:
      TreeEntry();
    };
} // close namespace Base

#endif // TreeEntry_hpp

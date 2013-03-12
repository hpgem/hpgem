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
        typedef TreeEntry<V>                                            Type;
        typedef TreeEntry<V>*                                           PtrType;
        typedef Base::TreeIterator<V,TreeEntry<V> >                     Iterator;
        typedef Base::TreeBase<Base::TreeIterator<V,TreeEntry<V> > >    BaseType;

        //! Constructor
        TreeEntry(const V& data)
              : data_(data), 
                referenceCounter_(0)
            { 
//               std::cout << "TreeEntry is created\n"; 
            }
            
        //! Copy constructor.
        //! We need this shallow copy for duplicating a mesh.  As the 
        //! compensation, we implement reference counter (referenceCounter_) to avoid 
        //! double freeing the memory
        TreeEntry(const TreeEntry<V>& other)
            : BaseType(other),
              data_(other.data_),
              referenceCounter_(other.referenceCounter_+1)
            {
//                 std::cout << "TreeEntry is copied\n";
            }

        ~TreeEntry() 
        { 
//           std::cout << "TreeEntry is deleted\n"; 
        }

        //! Return the data reference
        V& getData()
            {      
                return data_;
            }

    public:
        //! data stored at this TreeEntry
        V data_;
        
        //! counter for reference an original element: 0= the original, >0 a shallow copy of an original element
        int  referenceCounter_;

    private:
      TreeEntry();
    };
} // close namespace Base

#endif // TreeEntry_hpp

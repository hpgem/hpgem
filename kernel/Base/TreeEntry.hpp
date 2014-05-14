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

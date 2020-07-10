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

#ifndef HPGEM_KERNEL_CONSTITERABLEWRAPPER_H
#define HPGEM_KERNEL_CONSTITERABLEWRAPPER_H

#include <vector>

namespace hpgem {

template <typename T>
class IteratorOverConst;

/**
 * small wrapper to iterate over vectors of pointers where you shouldn't modify
 * whatever is in the pointer
 */
template <typename T>
class ConstIterableWrapper {
   public:
    ConstIterableWrapper(const std::vector<T*>& data) : data(data) {}

    IteratorOverConst<typename std::vector<T*>::const_iterator> begin() const {
        return data.begin();
    }

    IteratorOverConst<typename std::vector<T*>::const_iterator> end() const {
        return data.end();
    }

   private:
    const std::vector<T*>& data;
};

template <template <typename, typename...> class Iterator, typename T,
          typename... Q>
class IteratorOverConst<Iterator<T* const*, Q...>> {
   public:
    using difference_type = typename std::iterator_traits<
        Iterator<T* const*, Q...>>::difference_type;
    using value_type =
        typename std::iterator_traits<Iterator<T* const*, Q...>>::value_type;
    using pointer = const T* const*;  // pointer to const pointer to const T
    using reference =
        const T* const;  // const pointer to const T (no need to keep reference,
                         // since pointers are cheap to copy, also the implicit
                         // conversion from T*const to const T*const would cause
                         // the reference to be a reference to the wrong thing)
    using iterator_category = typename std::iterator_traits<
        Iterator<T* const*, Q...>>::iterator_category;

    //! default constructor
    IteratorOverConst() = default;

    //! Copy constructor
    IteratorOverConst(const IteratorOverConst& i) = default;
    //! Move constructor

    IteratorOverConst(IteratorOverConst&& i) = default;

    //! actual constructor;
    IteratorOverConst(Iterator<T* const*, Q...> nested) : inner(nested) {}

    //! Copy Assignment operator
    IteratorOverConst& operator=(const IteratorOverConst& i) = default;

    //! Move Assignment operator
    IteratorOverConst& operator=(IteratorOverConst&& i) = default;

    reference operator*() const { return *inner; }

    pointer operator->() const { return &(*inner); }

    IteratorOverConst& operator++() {
        ++inner;
        return *this;
    }

    IteratorOverConst& operator--() {
        --inner;
        return *this;
    }

    IteratorOverConst operator++(int)  // postinc
    {
        IteratorOverConst result = *this;
        ++inner;
        return result;
    }

    IteratorOverConst operator--(int)  // postdec
    {
        IteratorOverConst result = *this;
        --inner;
        return result;
    }

    //! Are they the same iterators?
    bool operator==(const IteratorOverConst& other) const {
        return inner == other.inner;
    }

    //! Are they different iterators?
    bool operator!=(const IteratorOverConst& i) const {
        return inner != i.inner;
    }

    IteratorOverConst& operator+=(difference_type n) {
        inner += n;
        return *this;
    }

    IteratorOverConst operator+(difference_type n) const {
        IteratorOverConst result{*this};
        return result += n;
    }

    IteratorOverConst& operator-=(difference_type n) {
        inner -= n;
        return *this;
    }

    IteratorOverConst operator-(difference_type n) const {
        IteratorOverConst result{*this};
        return result -= n;
    }

    difference_type operator-(const IteratorOverConst& other) const {
        return inner - other.inner;
    }

    reference operator[](difference_type n) const { return *(*this + n); }

    bool operator<(const IteratorOverConst& other) const {
        return inner < other.inner;
    }

    bool operator>(const IteratorOverConst& other) const {
        return inner > other.inner;
    }

    bool operator<=(const IteratorOverConst& other) const {
        return inner <= other.inner;
    }

    bool operator>=(const IteratorOverConst& other) const {
        return inner >= other.inner;
    }

   private:
    Iterator<T* const*, Q...> inner;
};

template <typename T>
IteratorOverConst<T> operator+(typename IteratorOverConst<T>::difference_type n,
                               const IteratorOverConst<T>& a) {
    IteratorOverConst<T> result{a};
    return result += n;
}

}  // namespace hpgem

#endif  // HPGEM_KERNEL_CONSTITERABLEWRAPPER_H

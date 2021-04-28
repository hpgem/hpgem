/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2021, University of Twente
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
#ifndef HPGEM_TEMPLATEARRAY_H
#define HPGEM_TEMPLATEARRAY_H

#include "tag.h"

namespace Preprocessor {

template <std::size_t size, template <std::size_t> class V>
class TemplateArray;

/**
 * Array with index-templates, that is an array of the form
 * [V<0>, V<1>, ..., V<size-1>]
 *
 * Note that this is not doable in a normal array, as the type of each entry is
 * different.
 * @tparam size The size of the array
 * @tparam V The entry template
 */
template <std::size_t size, template <std::size_t> class V>
class TemplateArray : public TemplateArray<size - 1, V> {
   public:
    // Default stuff
    TemplateArray() = default;
    ~TemplateArray() = default;
    TemplateArray(const TemplateArray&) = default;
    TemplateArray(TemplateArray&&) noexcept = default;
    TemplateArray& operator=(const TemplateArray&) = default;
    TemplateArray& operator=(TemplateArray&&) noexcept = default;

    ///
    template <typename... superArgs>
    TemplateArray(superArgs... args, V<size - 1> value)
        : TemplateArray<size - 1, V>(args...), value(value) {};

    /**
     * Access the values at a specific index
     * @tparam index The index to access
     * @return The value at the index
     */
    template <std::size_t index>
    typename std::enable_if<(index < size), V<index>&>::type get() {
        return TemplateArray<index + 1, V>::value;
    }

    template <std::size_t index>
    typename std::enable_if<(index < size), const V<index>&>::type get() const {
        return TemplateArray<index + 1, V>::value;
    }

    /**
     * Tag version of get()
     * @tparam index The index
     * @param tag The tag for the index
     * @return The value at the index
     */
    template <std::size_t index>
    V<index>& operator[](tag<index> tag) {
        return get<index>();
    }

    template <std::size_t index>
    const V<index>& operator[](tag<index> tag) const {
        return get<index>();
    }

   protected:  // Protected to allow access from other dimensions
    /// The value
    V<size - 1> value;
};

// Empty array
template <template <std::size_t> class V>
class TemplateArray<0, V> {
   public:
    // Default stuff
    TemplateArray() = default;
    ~TemplateArray() = default;
    TemplateArray(const TemplateArray&) = default;
    TemplateArray(TemplateArray&&) noexcept = default;
    TemplateArray& operator=(const TemplateArray&) = default;
    TemplateArray& operator=(TemplateArray&&) noexcept = default;
};

template <std::size_t>
class Sample;

}  // namespace Preprocessor

#endif  // HPGEM_TEMPLATEARRAY_H

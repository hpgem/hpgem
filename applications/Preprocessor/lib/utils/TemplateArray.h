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
#include "Logger.h"

namespace Preprocessor {

template <std::size_t size, template <std::size_t> class V>
class TemplateArray;

/**
 * Array for indexed-template, that is given a template V<std::size_t>
 * an array of the form  [V<0>, V<1>, ... V<size-1>].
 *
 * Note that each entry of the array has a different type, which would prevent a
 * regular array from being used.
 *
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

    /// Constructor providing the array items in reverse order
    ///
    /// This allows constructing a TemplateArray using:
    /// TemplateArray<size>(v<size-1>, .... v<1>, v<0>), where v<n> is the value
    /// at position n (and thus type V<n>). Note that contrary to regular array
    /// construction the values are in reverse order. This is needed to put the
    /// parameter packing as last argument.
    ///
    /// \tparam superArgs The types of the values
    /// \param value The value at index 'size-1'
    /// \param args The values at smaller indices (in reverse order)
    template <typename... superArgs>
    TemplateArray(V<size - 1> value, superArgs... args)
        : TemplateArray<size - 1, V>(args...), value(value){};

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

    template <int index>
    typename std::enable_if<index >= 0, V<index>&>::type operator[](
        itag<index> tag) {
        return get<index>();
    }

    template <int index>
    typename std::enable_if<index >= 0, const V<index>&>::type operator[](
        itag<index> tag) const {
        return get<index>();
    }

    // Dynamic access

    /**
     * @brief Dynamic mapping of the TemplateArray.
     *
     * This function is similar to functor(this[index]), but allows index to be
     * a runtime variable instead of a template argument. The functor must be
     * such that the return type is independent of the index that is accessed.
     *
     * Example use:
     *
     *     // The template array contains vectors of templated types
     *     template<std::size_t index>
     *     using VectorContent = std::vector<SomeTemplateType<index>>;
     *     TemplateArray<VectorContent> array;
     *     // Compute the size of the vector at 'parameter'.
     *     array.map(parameter, [](const auto& entry){return entry.size();}
     *
     * @param index The index to access
     * @param mapping The function to apply to the value identified by index.
     * Most likely a generic lambda, i.e. [](const auto&){....}
     * @tparam T The return type of the mapping
     * @tparam M The type of the mapping, should be equivalent to
     * template<i> T(V<i>&).
     * @return The return value of the function
     */
    template <typename T, typename M>
    T map(std::size_t index, M mapping) const {
        hpgem::logger.assert_debug(index < size, "Index % out of bounds [0,%)",
                                   index, size);
        // We can't index the values by a dynamic value. Hence, we need to do it
        // recursively via a helper function.
        return map(index, mapping, itag<size - 1>{});
    }

   protected:  // Protected to allow access from other dimensions
    /// The value
    V<size - 1> value;

   private:
    /// recursive helper for T map(std::size_t, F)
    template <typename T, typename F, int tindex>
    T map(std::size_t index, F functor, itag<tindex> tag) const {
        if (tindex == index) {
            return functor(value);
        } else {
            return map(index, functor, itag<tindex - 1>{});
        }
    }
    /// Base case
    template <typename T, typename F>
    T map(std::size_t, F, itag<-1>) const {
        return {};
    }
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

}  // namespace Preprocessor

#endif  // HPGEM_TEMPLATEARRAY_H

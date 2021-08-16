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
#ifndef HPGEM_LAZYCACHED_H
#define HPGEM_LAZYCACHED_H

#include <functional>
#include "Logger.h"

namespace hpgem {
namespace Base {

/**
 * A cache for a single value that is computed lazily.
 *
 * Note: As with most of hpgem, this is not safe for multithreading.
 *
 * @tparam T The type of the value that is cached
 */
template <typename T>
class LazyCached {
   public:
    /// Create a cache backed by a computation doing nothing
    LazyCached() : compute_([](T&) {}), hasValue_(false), value_(){};
    /**
     *  Create a cache backed by a computation to obtain the value.
     *
     *  Note the computation is only run when the value is first obtained. This
     *  should be taken into account for any side effects that it may cause.
     *
     *  Note: For the computation, the value is initialized to the default.
     *  \param computation The computation to update the value
     */
    explicit LazyCached(std::function<void(T&)> computation)
        : compute_(computation), hasValue_(false), value_() {}

    /**
     * Utility constructor using a non static member function as update
     * computation.
     * @tparam C The type of the class
     * @param instance The instance to call the member function on
     * @param func Pointer to the member function.
     */
    template <typename C>
    LazyCached(C* instance, void (C::*func)(T&))
        : LazyCached(std::bind(func, instance, std::placeholders::_1)){};

    /// Get the cached value
    const T& get() {
        if (!hasValue_) {
            compute_(value_);
            hasValue_ = true;
        }
        return value_;
    }

    /// Reset the cache to a state with no computed value
    void reset() { hasValue_ = false; }

   private:
    /// Function to update the value
    ///
    /// Design note: the signature is chosen so that the value can be updated
    /// in place. This can for example be used in conjunction with a vector
    /// as value to prevent allocating new vectors each time.
    std::function<void(T&)> compute_;
    /// Whether or not a value has been computed
    bool hasValue_{};
    T value_;
};

/**
 * The equivalent to a vector of LazyCached instances with a single computation.
 *
 * Note: Just like LazyCached this class is non thread-safe.
 *
 * @tparam T The type of value that is actually stored
 */
template <typename T>
class LazyVectorCached {
    using Entry = std::pair<bool, T>;

   public:
    /**
     * Default constructor, uses a comptuation that does not update the value.
     */
    LazyVectorCached() : compute_([](T&, std::size_t) {}), values_(){};

    /**
     * Create a cache vector backed by a computation that is used to compute the
     * entries of the vector on demand.
     *
     * Note the computation will be run only the first time the value for that
     * index is requested, and for each time after the first reset. This needs
     * to be taken into account for any side effects the computation may have.
     *
     * The first computation for a specific index after growing the cache to
     * include the index, will get a T() as value argument.
     *
     * @param computation The computation used to update the value.
     */
    explicit LazyVectorCached(std::function<void(T&, std::size_t)> computation)
        : compute_(computation), values_(){};

    /**
     * Utility constructor to use a non static member function as computation.
     * @tparam C The class
     * @param instance The instance to call the member function on
     * @param func The member function used to update the cache.
     */
    template <typename C>
    LazyVectorCached(C* instance, void (C::*func)(T&, std::size_t))
        : LazyVectorCached(std::bind(func, instance, std::placeholders::_1,
                                     std::placeholders::_2)){};

    const T& get(std::size_t index) {
        logger.assert_debug(index < values_.size(),
                            "Asking for index % with only % values", index,
                            values.size());
        Entry& entry = values_[index];
        if (!entry.first) {
            compute_(entry.second, index);
            entry.first = true;
        }
        return entry.second;
    }

    /**
     * Resets all the entries in the cache
     */
    inline void reset() {
        std::for_each(values_.begin(), values_.end(),
                      [](Entry& entry) { entry.first = false; });
    }
    /**
     * Reset and resize the cache.
     *
     * The parameterless constructor will be used to construct the new entries.
     *
     * @param newSize The new size of this cache
     */
    void reset(std::size_t newSize) {
        // A tiny performance gain could be had from comparing newSize to the
        // currentSize, but this is probably inconsequential to the cost of the
        // computation that is being cached.
        reset();
        values_.resize(newSize, Entry(false, T()));
    }

    std::size_t size() const { return values_.size(); }

    const T& operator[](std::size_t index) { return get(index); }

   private:
    std::function<void(T&, std::size_t)> compute_;
    std::vector<Entry> values_;
};

}  // namespace Base
}  // namespace hpgem

#endif  // HPGEM_LAZYCACHED_H

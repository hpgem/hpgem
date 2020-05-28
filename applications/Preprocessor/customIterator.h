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

#ifndef HPGEM_CUSTOMITERATOR_H
#define HPGEM_CUSTOMITERATOR_H

#include <cstddef>
#include <limits>
#include <memory>

#include <iterator>
#include <functional>
#include <utility>
#include <Logger.h>

namespace Preprocessor {

namespace Detail {
template <bool... predicates>
struct boolList {};

// just write it out with on ore two examples to validate
template <bool... predicates>
using OneOf = std::integral_constant<
    bool, !std::is_same<boolList<false, predicates...>,
                        boolList<predicates..., false>>::value>;

template <bool... predicates>
using AllOf =
    std::integral_constant<bool,
                           std::is_same<boolList<true, predicates...>,
                                        boolList<predicates..., true>>::value>;

template <typename T, typename... Ts>
using Contains = OneOf<std::is_same<T, Ts>::value...>;

template <typename T, typename... Ts>
struct Index_of;

template <typename T, typename U, typename... Ts>
struct Index_of<T, U, Ts...> {
    static_assert(Contains<T, Ts...>::value,
                  "Asking for the index of T, but it is not contained in Ts");
    static constexpr std::size_t value = 1 + Index_of<T, Ts...>::value;
};

template <typename T, typename... Ts>
struct Index_of<T, T, Ts...> {
    static_assert(!Contains<T, Ts...>::value,
                  "The type T appears more than once in the sequence");
    static constexpr std::size_t value = 0;
};

template <typename L, typename R>
std::enable_if_t<!std::is_convertible<R, L>::value &&
                     !std::is_convertible<L, R>::value,
                 bool>
    equals(L* left, R* right) {
    return false;
}

template <typename L, typename R>
std::enable_if_t<
    std::is_convertible<R, L>::value || std::is_convertible<L, R>::value, bool>
    equals(L* left, R* right) {
    return *left == *right;
}

template <typename L, typename R>
std::enable_if_t<!std::is_convertible<R, L>::value &&
                     !std::is_convertible<L, R>::value,
                 std::ptrdiff_t>
    difference(L* left, R* right) {
    logger(ERROR, "programmer error");
    return std::numeric_limits<std::ptrdiff_t>::quiet_NaN();
}

template <typename L, typename R>
std::enable_if_t<std::is_convertible<R, L>::value ||
                     std::is_convertible<L, R>::value,
                 std::ptrdiff_t>
    difference(L* left, R* right) {
    return *left - *right;
}

template <typename L, typename R>
std::enable_if_t<!std::is_convertible<R, L>::value> convert(L*, R*) {
    logger(ERROR, "programmer error");
};

template <typename L, typename R>
std::enable_if_t<std::is_convertible<R, L>::value> convert(L* to, R* from) {
    new (to) std::remove_cv_t<L>(*from);
};

}  // namespace Detail

// std::function cannot be moved; MoveFunction cannot be copied
template <typename T>
class MoveFunction;

template <typename Result, typename... Args>
class MoveFunction<Result(Args...)> {
    struct impl_base {
        impl_base() noexcept = default;
        virtual ~impl_base() = default;
        impl_base(const impl_base&) = delete;
        impl_base& operator=(const impl_base&) = delete;
        virtual Result operator()(Args...) = 0;

       protected:
        impl_base(impl_base&&) = default;
        impl_base& operator=(impl_base&&) = default;
    };
    template <typename F = Result (*)(Args...)>
    struct impl final : public impl_base {
        impl() noexcept = default;
        ~impl() override = default;
        impl(const impl&) = delete;
        impl(impl&&) = default;

        impl(F function) : function(std::move(function)) {}

        impl& operator=(const impl&) = delete;
        impl& operator=(impl&&) = default;

        Result operator()(Args... args) override {
            return function(std::forward<Args>(args)...);
        }
        F function;
    };
    std::unique_ptr<impl_base> function_;

   public:
    MoveFunction() noexcept = default;
    ~MoveFunction() = default;
    MoveFunction(const MoveFunction&) = delete;
    MoveFunction(MoveFunction&&) noexcept = default;

    template <typename F, typename = std::enable_if_t<!std::is_same<
                              std::decay_t<F>, MoveFunction>::value>>
    MoveFunction(F&& function)
        : function_(
              std::make_unique<impl<F>>(impl<F>(std::forward<F>(function)))) {}

    MoveFunction& operator=(const MoveFunction&) = delete;
    MoveFunction& operator=(MoveFunction&&) noexcept = default;

    Result operator()(Args... args) const {
        return (*function_)(std::forward<Args>(args)...);
    }
};

// void doesn't return
template <typename... Args>
class MoveFunction<void(Args...)> {
    struct impl_base {
        impl_base() noexcept = default;
        virtual ~impl_base() = default;
        impl_base(const impl_base&) = delete;
        impl_base& operator=(const impl_base&) = delete;
        virtual void operator()(Args...) = 0;

       protected:
        impl_base(impl_base&&) = default;
        impl_base& operator=(impl_base&&) = default;
    };
    template <typename F = void (*)(Args...)>
    struct impl final : public impl_base {
        impl() noexcept = default;
        ~impl() override = default;
        impl(const impl&) = delete;
        impl(impl&&) = default;

        impl(F function) : function(std::move(function)) {}

        impl& operator=(const impl&) = delete;
        impl& operator=(impl&&) = default;

        void operator()(Args... args) override {
            function(std::forward<Args>(args)...);
        }
        F function;
    };
    std::unique_ptr<impl_base> function_;

   public:
    MoveFunction() = default;
    ~MoveFunction() = default;
    MoveFunction(const MoveFunction&) = delete;
    MoveFunction(MoveFunction&&) noexcept = default;

    template <typename F, typename = std::enable_if_t<!std::is_same<
                              std::decay_t<F>, MoveFunction>::value>>
    MoveFunction(F&& function)
        : function_(
              std::make_unique<impl<F>>(impl<F>(std::forward<F>(function)))) {}

    MoveFunction& operator=(const MoveFunction&) = delete;
    MoveFunction& operator=(MoveFunction&&) noexcept = default;

    void operator()(Args... args) { (*function_)(std::forward<Args>(args)...); }
};

/**
 * @brief convert an action into an iterator
 * @tparam T the value type of the iterator
 * An action is a function that converts its argument into the next value in the
 * sequence. For ease of finding the end of a range, the default-constructed
 * ActionIterator will compare equal to all past-the-end ActionIterators.
 */
template <typename T>
class ActionIterator {
   public:
    using difference_type = std::size_t;
    using value_type = T;
    using pointer = T*;
    using reference = T&;
    using iterator_category = std::input_iterator_tag;

    ActionIterator() = default;

    ActionIterator(const ActionIterator&) = delete;

    ActionIterator(ActionIterator&& other) noexcept
        : sentinel(other.sentinel),
          increment(std::move(other.increment)),
          valid(std::move(other.valid)),
          value(std::move(other.value)),
          remaining(other.remaining) {
        if (remaining > 0) {
            valid = [this](const T&) { return remaining > 0; };
        }
    }

    ~ActionIterator() = default;

    ActionIterator& operator=(const ActionIterator&) = delete;

    ActionIterator& operator=(ActionIterator&& other) noexcept {
        sentinel = other.sentinel;
        value = std::move(other.value);
        increment = std::move(other.increment);
        remaining = other.remaining;
        if (remaining > 0) {
            valid = [this](const T&) { return remaining > 0; };
        } else {
            valid = std::move(other.valid);
        }
    }

    //! construct a general iterator from an initial value, a function that can
    //! increment the iterator (find a next value, given a value) and a
    //! predicate
    //! that is false when the end is reached
    ActionIterator(T first, MoveFunction<void(T&)>&& increment,
                   MoveFunction<bool(const T&)>&& valid)
        : sentinel(false),
          increment(std::move(increment)),
          valid(std::move(valid)),
          value(first){};

    //! construct a general iterator that will take n steps from an initial
    //! value and a function that can increment the iterator
    //
    ActionIterator(T first, MoveFunction<void(T&)>&& increment, std::size_t n)
        : sentinel(false),
          increment(std::move(increment)),
          valid([this](const T&) { return remaining > 0; }),
          value(first),
          remaining(n){};

    reference operator*() { return value; }

    const T& operator*() const { return value; }

    pointer operator->() { return &value; }

    const T* operator->() const { return &value; }

    ActionIterator& operator++() {
        if (remaining > 0) {
            remaining--;
        }
        if (valid(value)) {
            increment(value);
        }
        return *this;
    }

    ActionIterator& operator++(int) {
        ActionIterator result = *this;
        ++*this;
        return result;
    }

    // todo: compare predicates
    bool operator==(const ActionIterator& other) const {
        if (sentinel || !valid(value)) {
            return other.sentinel || !other.valid(value);
        }
        if (other.sentinel || !other.valid(value)) {
            return false;
        }
        if (remaining > 0 && other.remaining > 0) {
            return value == other.value && remaining == other.remaining;
        }
        return value == other.value;
    }

    bool operator!=(const ActionIterator& other) const {
        return !(*this == other);
    }

   private:
    bool sentinel = true;
    MoveFunction<void(T&)> increment;
    MoveFunction<bool(const T&)> valid;
    T value;
    std::size_t remaining = 0;
};

template <typename T>
class Range {
   public:
    Range() = default;

    Range(const Range&) = delete;

    Range(Range&&) noexcept = default;

    ~Range() = default;

    Range& operator=(const Range&) = delete;

    Range& operator=(Range&&) noexcept = default;

    explicit Range(ActionIterator<T> begin) : beginIterator(begin){};

    Range(T first, MoveFunction<void(T&)>&& increment,
          MoveFunction<bool(const T&)>&& valid)
        : beginIterator{first, std::move(increment), std::move(valid)} {};

    Range(T first, MoveFunction<void(T&)>&& increment, std::size_t n)
        : beginIterator{first, std::move(increment), n} {};

    ActionIterator<T>&& begin() { return std::move(beginIterator); }

    ActionIterator<T>&& end() { return std::move(endIterator); }

   private:
    ActionIterator<T> beginIterator;
    ActionIterator<T> endIterator{};
};

template <typename T>
ActionIterator<T> begin(Range<T> range) {
    return range.begin();
}

template <typename T>
ActionIterator<T> end(Range<T> range) {
    return range.end();
}

/// An union of iterators that have the same element type.
///
/// Unfortunately the move constructor and move-assignment operator give
/// errors in gcc, therefore it is not currently used.
/// \tparam T The element type
/// \tparam Iterators The types of the iterators
template <typename T, typename... Iterators>
class hybridIterator {
    template <std::size_t I, typename... Ts>
    using get_t = std::tuple_element_t<I, std::tuple<Ts...>>;

   public:
    using value_type = T;
    using difference_type = std::common_type_t<
        typename std::iterator_traits<Iterators>::difference_type...>;
    using reference =
        typename std::iterator_traits<get_t<0, Iterators...>>::reference;
    using pointer = std::common_type_t<
        typename std::iterator_traits<Iterators>::pointer...>;
    using iterator_category = std::common_type_t<
        typename std::iterator_traits<Iterators>::iterator_category...>;

    static_assert(Detail::AllOf<std::is_same<
                      reference, typename std::iterator_traits<
                                     Iterators>::reference>::value...>::value,
                  "Some iterators have a different reference type");

    hybridIterator() {
        static_assert(sizeof...(Iterators) > 0,
                      "hybridIterator needs to wrap at least one iterator");
        new (&data) get_t<0, Iterators...>{};
    }

    hybridIterator(const hybridIterator& other) : active(other.active) {
        process([&](auto* data) {
            new (data) std::remove_cv_t<std::remove_pointer_t<decltype(data)>>(
                reinterpret_cast<const std::remove_cv_t<
                    std::remove_pointer_t<decltype(data)>>&>(other.data));
        });
    }

    hybridIterator(hybridIterator&& other) noexcept : active(other.active) {
        process([&](auto* data) {
            new (data) std::remove_cv_t<std::remove_pointer_t<decltype(data)>>(
                reinterpret_cast<
                    std::remove_cv_t<std::remove_pointer_t<decltype(data)>>&&>(
                    std::move(other.data)));
        });
    }

    template <typename S,
              typename = std::enable_if_t<Detail::Contains<
                  std::decay_t<S>, std::decay_t<Iterators>...>::value>>
    hybridIterator(const S& data) {
        new (&this->data) std::decay_t<S>(data);
        active = Detail::Index_of<std::decay_t<S>,
                                  std::decay_t<Iterators>...>::value;
    };

    std::size_t getActive() const noexcept { return active; }

    template <typename... otherIterators,
              typename = std::enable_if_t<
                  Detail::AllOf<std::is_convertible<
                      otherIterators, Iterators>::value...>::value &&
                  !std::is_same<std::tuple<Iterators...>,
                                std::tuple<otherIterators...>>::value>>
    hybridIterator(const hybridIterator<T, otherIterators...>& other)
        : active(other.getActive()) {
        process([other](auto* data) {
            other.process(
                [data](auto* otherData) { Detail::convert(data, otherData); });
            // static_cast<const
            // std::remove_pointer_t<decltype(data)>&>(other));
        });
    };

    template <typename S,
              typename = std::enable_if_t<Detail::Contains<
                  std::decay_t<S>, std::decay_t<Iterators>...>::value>>
    operator S() {
        return forward([](auto* data) {
            logger.assert_debug(typeid(S) == typeid(*data),
                                "Can only convert to active iterator");
            return reinterpret_cast<S&>(*data);
        });
    };

    bool operator==(const hybridIterator& other) const {
        if (active != other.active) return false;
        return forward([other](auto* data) {
            return other.forward([data](auto* otherData) {
                return Detail::equals(data, otherData);
            });
        });
    }

    bool operator!=(const hybridIterator& other) const {
        return !(*this == other);
    }

    reference operator*() {
        return forward([](auto* data) -> reference { return *(*data); });
    }

    pointer operator->() { return &(*(*this)); }

    hybridIterator& operator++() {
        process([](auto* data) { ++(*data); });
        return *this;
    }

    hybridIterator operator++(int) {
        auto temp = *this;
        ++(*this);
        return temp;
    }

    hybridIterator& operator+=(std::ptrdiff_t n) {
        process([&](auto* data) { (*data) += n; });
        return *this;
    }

    reference operator[](std::ptrdiff_t n) { return *(*this + n); }

    hybridIterator& operator--() {
        process([](auto* data) { --(*data); });
        return *this;
    }

    hybridIterator operator--(int) {
        auto temp = *this;
        --(*this);
        return temp;
    }

    hybridIterator& operator-=(std::ptrdiff_t n) {
        process([&](auto* data) { (*data) -= n; });
        return *this;
    }

    std::ptrdiff_t operator-(const hybridIterator& other) const {
        logger.assert_debug(
            active == other.active,
            "You cannot compare two pointers into different ranges");
        return forward([other](auto* data) {
            return other.forward([data](auto* otherData) {
                return Detail::difference(data, otherData);
            });
        });
    }

    hybridIterator& operator=(const hybridIterator& other) {
        process([&](auto* data) { clear(data); });
        active = other.active;
        process([&](auto* data) {
            new (data) std::remove_cv_t<std::remove_pointer_t<decltype(data)>>(
                reinterpret_cast<const std::remove_cv_t<
                    std::remove_pointer_t<decltype(data)>>&>(other.data));
        });
        return *this;
    }

    hybridIterator& operator=(hybridIterator&& other) noexcept {
        process([&](auto* data) { clear(data); });
        active = other.active;
        process([&](auto* data) {
            new (data) std::remove_cv_t<std::remove_pointer_t<decltype(data)>>(
                reinterpret_cast<
                    std::remove_cv_t<std::remove_pointer_t<decltype(data)>>&&>(
                    std::move(other.data)));
        });
        return *this;
    }

    ~hybridIterator() {
        process([&](auto* data) { clear(data); });
    }

    template <typename F>
    decltype(auto) forward(const F& function) {
        return forward(function,
                       std::make_index_sequence<sizeof...(Iterators)>());
    }

    template <typename F>
    decltype(auto) forward(const F& function) const {
        return forward(function,
                       std::make_index_sequence<sizeof...(Iterators)>());
    }

    template <typename F>
    void process(const F& function) const {
        process(function, std::make_index_sequence<sizeof...(Iterators)>());
    }

    template <typename F>
    void process(const F& function) {
        process(function, std::make_index_sequence<sizeof...(Iterators)>());
    }

   private:
    template <typename S>
    void clear(S* data) {
        data->~S();
    }

    template <typename F, std::size_t N, std::size_t... Ns>
    void process(const F& function, std::index_sequence<N, Ns...>) {
        if (active == N) {
            function(reinterpret_cast<get_t<N, Iterators...>*>(&data));
        } else {
            process(function, std::index_sequence<Ns...>{});
        }
    };

    template <typename F>
    void process(const F& other, std::index_sequence<>) {
        logger(ERROR, "programmer error");
    }

    template <typename F, std::size_t N, std::size_t... Ns>
    void process(const F& function, std::index_sequence<N, Ns...>) const {
        if (active == N) {
            function(reinterpret_cast<const get_t<N, Iterators...>*>(&data));
        } else {
            process(function, std::index_sequence<Ns...>{});
        }
    };

    template <typename F>
    void process(const F& other, std::index_sequence<>) const {
        logger(ERROR, "programmer error");
    }

    template <typename F, std::size_t N, std::size_t... Ns>
    decltype(auto) forward(const F& function, std::index_sequence<N, Ns...>) {
        if (active == N) {
            return function(reinterpret_cast<get_t<N, Iterators...>*>(&data));
        } 
            return forward(function, std::index_sequence<Ns...>{});
        
    };

    template <typename F>
    decltype(auto) forward(const F& function, std::index_sequence<>) {
        logger(ERROR, "programmer error");
        auto temp = get_t<0, Iterators...>{};
        return function(&temp);
    }

    template <typename F, std::size_t N, std::size_t... Ns>
    decltype(auto) forward(const F& function,
                           std::index_sequence<N, Ns...>) const {
        if (active == N) {
            return function(
                reinterpret_cast<const get_t<N, Iterators...>*>(&data));
        } 
            return forward(function, std::index_sequence<Ns...>{});
        
    };

    template <typename F>
    decltype(auto) forward(const F& function, std::index_sequence<>) const {
        logger(ERROR, "programmer error");
        auto temp = get_t<0, Iterators...>{};
        return function(&temp);
    }

   private:
    std::size_t active = 0;
    std::aligned_union_t<0, Iterators...> data;
};

template <typename T, typename... Iterators>
hybridIterator<T, Iterators...> operator+(
    const hybridIterator<T, Iterators...>& it, std::size_t n) {
    auto temp = it;
    return temp += n;
}

template <typename T, typename... Iterators>
hybridIterator<T, Iterators...> operator+(
    std::size_t n, const hybridIterator<T, Iterators...>& it) {
    auto temp = it;
    return temp += n;
}

template <typename T, typename... Iterators>
hybridIterator<T, Iterators...> operator-(
    const hybridIterator<T, Iterators...>& it, std::size_t n) {
    auto temp = it;
    return temp -= n;
}

template <typename T, typename... Iterators>
bool operator<(const hybridIterator<T, Iterators...>& a,
               const hybridIterator<T, Iterators...>& b) {
    return (b - a) > 0;
}

template <typename T, typename... Iterators>
bool operator>(const hybridIterator<T, Iterators...>& a,
               const hybridIterator<T, Iterators...>& b) {
    return b < a;
}

template <typename T, typename... Iterators>
bool operator<=(const hybridIterator<T, Iterators...>& a,
                const hybridIterator<T, Iterators...>& b) {
    return !(b < a);
}

template <typename T, typename... Iterators>
bool operator>=(const hybridIterator<T, Iterators...>& a,
                const hybridIterator<T, Iterators...>& b) {
    return !(a < b);
}
}  // namespace Preprocessor

#endif  // HPGEM_CUSTOMITERATOR_H

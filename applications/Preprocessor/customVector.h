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

#ifndef HPGEM_APP_CUSTOMVECTOR_H
#define HPGEM_APP_CUSTOMVECTOR_H

#include <memory>
#include <vector>
#include <Logger.h>
#include "customIterator.h"

namespace Preprocessor {

/**
 * @brief A thing that behaves like vector but can also store a bit of data
 * locally
 * @tparam T value type to be stored. T should not manage a limited resource
 * @tparam N number of values that should be storable internally
 * if the size of the vector regularly crosses the in-place threshold
 * performance will degrade
 *
 * Note, currently unused (and replaced by std::vector) as the iterator does not
 * compile with gcc.
 */
template <typename T, std::size_t N = 8>
class stackVector {
   public:
    using value_type = T;
    using iterator = hybridIterator<T, typename std::vector<T>::iterator,
                                    typename std::array<T, N>::iterator>;
    using const_iterator =
        hybridIterator<T, typename std::vector<T>::const_iterator,
                       typename std::array<T, N>::const_iterator>;
    using reverse_iterator =
        hybridIterator<T, typename std::vector<T>::reverse_iterator,
                       typename std::array<T, N>::reverse_iterator>;
    using const_reverse_iterator =
        hybridIterator<T, typename std::vector<T>::const_reverse_iterator,
                       typename std::array<T, N>::const_reverse_iterator>;

    stackVector() : stackVector(0){};
    stackVector(std::size_t count, const T& value = T{}) : size_(count) {
        if (size_ <= N) {
            data_.array = std::array<T, N>();
            data_.array.fill(value);
        } else {
            data_.array.~array();
            new (&data_.vector) std::vector<T>(count, value);
        }
    }
    template <typename inputIt>
    stackVector(inputIt first, inputIt last) {
        size_ = std::distance(first, last);
        if (size_ <= N) {
            data_.array = std::array<T, N>();
            for (std::size_t i = 0; i < size_; ++i) {
                data_.array[i] = *first++;
            }
        } else {
            data_.array.~array();
            new (&data_.vector) std::vector<T>(first, last);
        }
    }
    stackVector(const stackVector& other) : size_(other.size_) {
        if (size_ <= N) {
            data_.array = other.data_.array;
        } else {
            data_.array.~array();
            new (&data_.vector) std::vector<T>(other.data_.vector);
        }
    }
    stackVector(stackVector&&) = default;
    stackVector(std::initializer_list<T> init) {
        size_ = init.size();
        if (size_ <= N) {
            data_.array = std::array<T, N>();
            std::copy(init.begin(), init.end(), data_.array.begin());
        } else {
            data_.array.~array();
            new (&data_.vector) std::vector<T>(init);
        }
    }

    ~stackVector() {
        if (size_ <= N) {
            data_.array.~array();
        } else {
            data_.vector.~vector();
        }
    };

    stackVector& operator=(const stackVector& other) {
        if (size_ <= N) {
            data_.array.~array();
        } else {
            data_.vector.~vector();
        }
        size_ = other.size_;
        if (other.size_ <= N) {
            data_.array = other.data_.array;
        } else {
            new (&data_.vector) std::vector<T>(other.data_.vector);
        }
        return *this;
    }
    stackVector& operator=(stackVector&&) = default;
    stackVector& operator=(std::initializer_list<T> init) {
        return *this = stackVector(init);
    }

    void assign(std::size_t count, const T& value) {
        *this = stackVector(count, value);
    }

    template <typename inputIt>
    void assign(inputIt first, inputIt last) {
        *this = stackVector(first, last);
    }

    void assign(std::initializer_list<T> init) { *this = stackVector(init); }

    T& at(std::size_t p) {
        if (p < size_) {
            return (*this)[p];
        }
        throw std::out_of_range{"vector access out of range"};
    }

    const T& at(std::size_t p) const {
        if (p < size_) {
            return (*this)[p];
        }
        throw std::out_of_range{"vector access out of range"};
    }

    T& operator[](std::size_t p) {
        if (size_ <= N) {
            return data_.array[p];
        }
        return data_.vector[p];
    }

    const T& operator[](std::size_t p) const {
        if (size_ <= N) {
            return data_.array[p];
        }
        return data_.vector[p];
    }

    T& front() { return (*this)[0]; }

    const T& front() const { return (*this)[0]; }

    T& back() { return (*this)[size_ - 1]; }

    const T& back() const { return (*this)[size_ - 1]; }

    T* data() noexcept {
        if (size_ <= N) {
            return data_.array.data();
        }
        return data_.vector.data();
    }

    const T* data() const noexcept {
        if (size_ <= N) {
            return data_.array.data();
        }
        return data_.vector.data();
    }

    iterator begin() noexcept {
        if (size_ <= N) {
            return iterator{data_.array.begin()};
        }
        return iterator{data_.vector.begin()};
    }

    const_iterator begin() const noexcept {
        if (size_ <= N) {
            return const_iterator{data_.array.begin()};
        }
        return const_iterator{data_.vector.begin()};
    }

    const_iterator cbegin() const noexcept {
        if (size_ <= N) {
            return const_iterator{data_.array.cbegin()};
        }
        return const_iterator{data_.vector.cbegin()};
    }

    iterator end() noexcept {
        if (size_ <= N) {
            return iterator{data_.array.begin() + size_};
        }
        return iterator{data_.vector.end()};
    }

    const_iterator end() const noexcept {
        if (size_ <= N) {
            return const_iterator{data_.array.begin() + size_};
        }
        return const_iterator{data_.vector.end()};
    }

    const_iterator cend() const noexcept {
        if (size_ <= N) {
            return const_iterator{data_.array.cbegin() + size_};
        }
        return const_iterator{data_.vector.cend()};
    }

    reverse_iterator rbegin() noexcept {
        if (size_ <= N) {
            return reverse_iterator{data_.array.rend() - size_};
        }
        return reverse_iterator{data_.vector.rbegin()};
    }

    const_reverse_iterator rbegin() const noexcept {
        if (size_ <= N) {
            return const_reverse_iterator{data_.array.rend() - size_};
        }
        return const_reverse_iterator{data_.vector.rbegin()};
    }

    const_reverse_iterator crbegin() const noexcept {
        if (size_ <= N) {
            return const_reverse_iterator{data_.array.crend() - size_};
        }
        return const_reverse_iterator{data_.vector.crbegin()};
    }

    reverse_iterator rend() noexcept {
        if (size_ <= N) {
            return reverse_iterator{data_.array.rend()};
        }
        return reverse_iterator{data_.vector.rend()};
    }

    const_reverse_iterator rend() const noexcept {
        if (size_ <= N) {
            return const_reverse_iterator{data_.array.rend()};
        }
        return const_reverse_iterator{data_.vector.rend()};
    }

    const_reverse_iterator crend() const noexcept {
        if (size_ <= N) {
            return const_reverse_iterator{data_.array.crend()};
        }
        return const_reverse_iterator{data_.vector.crend()};
    }

    bool empty() const noexcept { return size_ == 0; }

    std::size_t size() const noexcept { return size_; }

    std::size_t max_size() const noexcept {
        return std::max(N, std::vector<T>{}.max_size());
    }

    void reserve(std::size_t new_cap) {
        if (size_ <= N) {
            vector_limit = std::max(new_cap, vector_limit);
        } else {
            data_.vector.reserve(new_cap);
        }
    }

    std::size_t capacity() const noexcept {
        if (size_ <= N) {
            return vector_limit;
        }
        return data_.vector.capacity();
    }

    void shrink_to_fit() {
        if (size_ <= N) {
            vector_limit = N + 1;
        } else {
            data_.vector.shrink_to_fit();
        }
    }

    void clear() noexcept {
        if (size_ > N) {
            vector_limit = data_.vector.capacity();
            data_.vector.~vector();
            data_.array = std::array<T, N>();
        }
        data_.array.fill(T{});
        size_ = 0;
    }

    iterator insert(const_iterator position, const T& value) {
        return insert(position, 1, value);
    }

    iterator insert(const_iterator position, T&& value) {
        if (size_ < N) {
            std::copy_backward(position, data_.array.begin() + size_,
                               data_.array.begin() + size_ + 1);
            size_++;
            std::size_t distance = position - data_.array.begin();
            data_.array[distance] = std::move(value);
            return iterator{data_.array.begin() + distance};
        }
        if (size_ == N) {
            std::vector<T> new_data;
            new_data.reserve(vector_limit);
            new_data.insert(new_data.end(), data_.array.begin(), position);
            iterator result = data_.vector.end();
            std::size_t distance = std::distance(data_.vector.begin(), result);
            new_data.push_back(std::move(value));
            new_data.insert(new_data.end(), position,
                            data_.array.begin() + size_);
            data_.array.~array();
            new (&data_.vector) std::vector<T>(std::move(new_data));
            size_++;
            result = data_.vector.begin() + distance;
            return result;
        } else {
            size_++;
            return iterator{data_.vector.insert(position, std::move(value))};
        }
    }

    iterator insert(const_iterator position, std::size_t count,
                    const T& value) {
        if (size_ <= N) {
            if (size_ + count <= N) {
                std::copy_backward(position, cend(),
                                   data_.array.begin() + size_ + count);
                size_ += count;
                std::size_t distance = position - cbegin();
                for (std::size_t i = 0; i < count; ++i) {
                    data_.array[distance + i] = value;
                }
                return iterator{data_.array.begin() + distance};
            }
            std::vector<T> new_data;
            new_data.reserve(vector_limit);
            new_data.insert(new_data.end(), cbegin(), position);
            iterator result = new_data.insert(new_data.end(), count, value);
            std::size_t distance = std::distance({new_data.begin()}, result);
            new_data.insert(new_data.end(), position, cend());
            data_.array.~array();
            new (&data_.vector) std::vector<T>(std::move(new_data));
            size_ += count;
            result = data_.vector.begin() + distance;
            return result;

        } else {
            size_ += count;
            return iterator{data_.vector.insert(position, count, value)};
        }
    }

    template <typename InputIt>
    iterator insert(const_iterator position, InputIt first, InputIt last) {
        std::size_t count = std::distance(first, last);
        if (size_ <= N) {
            if (size_ + count <= N) {
                if (count > 0)
                    std::copy_backward(position, cend(), end() + count);
                size_ += count;
                std::size_t distance = position - cbegin();
                std::copy(first, last, data_.array.begin() + distance);
                return iterator{data_.array.begin() + distance};
            }
            std::vector<T> new_data;
            new_data.reserve(vector_limit);
            new_data.insert(new_data.end(), cbegin(), position);
            iterator result = new_data.insert(new_data.end(), first, last);
            std::size_t distance = std::distance({new_data.begin()}, result);
            new_data.insert(new_data.end(), position, cend());
            data_.array.~array();
            new (&data_.vector) std::vector<T>(std::move(new_data));
            size_ += count;
            result = data_.vector.begin() + distance;
            return result;

        } else {
            size_ += count;
            return iterator{data_.vector.insert(position, first, last)};
        }
    }

    iterator insert(const_iterator position, std::initializer_list<T> ilist) {
        std::size_t count = ilist.size();
        if (size_ <= N) {
            if (size_ + count <= N) {
                if (count > 0)
                    std::copy_backward(position, data_.array.begin() + size_,
                                       data_.array.begin() + size_ + count);
                size_ += count;
                std::size_t distance = position - cbegin();
                std::copy(ilist.begin(), ilist.end(),
                          data_.array.begin() + distance);
                return iterator{data_.array.begin() + distance};
            }
            std::vector<T> new_data;
            new_data.reserve(vector_limit);
            new_data.insert(new_data.end(), data_.array.begin(), position);
            iterator result = new_data.insert(data_.vector.end(), ilist);
            std::size_t distance = std::distance(data_.vector.begin(), result);
            new_data.insert(new_data.end(), position,
                            data_.array.begin() + size_);
            data_.array.~array();
            new (&data_.vector) std::vector<T>(std::move(new_data));
            size_ += count;
            result = data_.vector.begin() + distance;
            return result;

        } else {
            size_ += count;
            return iterator{data_.vector.insert(position, ilist)};
        }
    }

    template <typename... Args>
    iterator emplace(const_iterator position, Args&&... args) {
        return insert(position, T{std::forward<Args>(args)...});
    }

    iterator erase(const_iterator position) {
        return erase(position, position + 1);
    }

    iterator erase(const_iterator first, const_iterator last) {
        std::size_t count = std::distance(first, last);
        if (size_ <= N) {
            auto result = begin() + std::distance(cbegin(), first);
            if (count > 0) std::copy(last, cend(), result);
            size_ -= count;
            return result;
        }
        if (size_ - count <= N) {
            std::size_t first_offset = std::distance(cbegin(), first);
            std::size_t last_offset = std::distance(cbegin(), last);
            std::vector<T> old_data = std::move(data_.vector);
            data_.vector.~vector();
            data_.array = std::array<T, N>{};
            iterator result =
                std::copy(old_data.begin(), old_data.begin() + first_offset,
                          data_.array.begin());
            std::copy(old_data.begin() + last_offset, old_data.end(), result);
            size_ -= count;
            return result;
        } else {
            size_ -= count;
            return iterator{data_.vector.erase(first, last)};
        }
    }

    void push_back(const T& value) { insert(cend(), value); }

    void push_back(T&& value) { insert(cend(), std::move(value)); }

    template <typename... Args>
    void emplace_back(Args&&... args) {
        emplace(cend(), std::forward<Args>(args)...);
    }

    void pop_back() { erase(end() - 1); }

    void resize(std::size_t count) { resize(count, T{}); }

    void resize(std::size_t count, const T& elem) {
        if (size_ <= N && count <= N) {
            for (std::size_t i = size_; i < count; ++i) {
                data_.array[i] = elem;
            }
            for (std::size_t i = count; i < size_; ++i) {
                data_.array[i] = T{};
            }
        } else if (size_ > N && count > N) {
            data_.vector.resize(count, elem);
        } else if (size_ > N) {
            vector_limit = data_.vector.capacity();
            std::vector<T> old_data = std::move(data_.vector);
            data_.vector.~vector();
            data_.array = std::array<T, N>{};
            std::copy(old_data.begin(), old_data.begin() + count,
                      data_.array.begin());
        } else {
            std::vector<T> new_data;
            new_data.reserve(vector_limit);
            new_data.insert(new_data.end(), data_.array.begin(),
                            data_.array.end() + size_);
            new_data.resize(count, elem);
            data_.array.~array();
            new (&data_.vector) std::vector<T>(std::move(new_data));
        }
        size_ = count;
    }

    void swap(stackVector& other) {
        if (size_ <= N && other.size_ <= N) {
            std::swap(data_.array, other.data_.array);
        } else if (size_ > N && other.size_ > N) {
            std::swap(data_.vector, other.data_.vector);
        } else if (size_ > N) {
            std::vector<T> old_data = std::move(data_.vector);
            data_.vector.~vector();
            data_.array = other.data_.array;
            other.data_.array.~array();
            new (&other.data_.vector) std::vector<T>(std::move(old_data));
        } else {
            other.swap(*this);
        }
        std::swap(size_, other.size_);
        std::swap(vector_limit, other.vector_limit);
    }

   private:
    std::size_t size_;
    std::size_t vector_limit = N + 1;
    union U {
        std::vector<T> vector;
        std::array<T, N> array;
        U() : array(std::array<T, N>()){};
        ~U() = default;
        ;
    } data_;
};

template <typename T, std::size_t N>
bool operator==(const stackVector<T, N>& l, const stackVector<T, N>& r) {
    if (l.size() != r.size()) return false;
    for (auto lIt = l.begin(), rIt = r.begin(); lIt != l.end(); ++lIt, ++rIt) {
        if (*lIt != *rIt) return false;
    }
    return true;
};

template <typename T, std::size_t N>
bool operator!=(const stackVector<T, N>& l, const stackVector<T, N>& r) {
    return !(l == r);
};

template <typename T, std::size_t N>
bool operator<(const stackVector<T, N>& l, const stackVector<T, N>& r) {
    return std::lexicographical_compare(l.begin(), l.end(), r.begin(), r.end());
};

template <typename T, std::size_t N>
bool operator>(const stackVector<T, N>& l, const stackVector<T, N>& r) {
    return r < l;
};

template <typename T, std::size_t N>
bool operator<=(const stackVector<T, N>& l, const stackVector<T, N>& r) {
    return !(r < l);
};

template <typename T, std::size_t N>
bool operator>=(const stackVector<T, N>& l, const stackVector<T, N>& r) {
    return !(l < r);
};

template <typename T, std::size_t N>
void swap(stackVector<T, N>& l, stackVector<T, N>& r) {
    l.swap(r);
};
}  // namespace Preprocessor

#endif // HPGEM_APP_CUSTOMVECTOR_H

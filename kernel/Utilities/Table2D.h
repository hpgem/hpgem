/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2020, University of Twente
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
#ifndef HPGEM_TABLE2D_H
#define HPGEM_TABLE2D_H

#include <limits>
#include <memory>
#include <vector>

#include "Logger.h"

namespace hpgem {
namespace Utilities {

/// A safe 2D vector.
///
/// \tparam type The entry type
template <typename EntryType>
class Table2D {
   public:
    using type = bool;
    Table2D() : rows_(0), columns_(0){};

    Table2D(std::size_t rows, std::size_t columns)
        : rows_(rows), columns_(columns), entries_(rows * columns) {
        checkSize(rows, columns);
    };

    Table2D(std::size_t rows, std::size_t columns, const type& fill)
        : rows_(rows), columns_(columns), entries_(rows * columns, fill) {
        checkSize(rows, columns);
    };

    /// Fill the table with a constant value.
    void fill(type value) {
        std::fill(entries_.begin(), entries_.end(), value);
    }

    /// Resize the table. The content of the new table is undefined.
    void resize(std::size_t rows, std::size_t columns) {
        checkSize(rows, columns);
        rows_ = rows;
        columns_ = columns;
        entries_.resize(rows * columns);
    }

    std::size_t getNumberOfRows() const { return rows_; }

    std::size_t getNumberOfColumns() const { return columns_; }

    /// Total number of entries
    std::size_t getSize() const { return rows_ * columns_; }

    /// Get the value of a certain entry
    ///
    /// \param n The row number
    /// \param m The column number
    /// \return The value
    std::vector<type>::reference operator()(std::size_t n, std::size_t m) {
        logger.assert_debug(n < rows_,
                            "Requested row % for a table with only % rows", n,
                            rows_);
        logger.assert_debug(
            m < columns_, "Requested column % for a table with only % columns",
            m, columns_);

        return entries_[n + m*rows_];
    }

    std::vector<type>::const_reference operator()(std::size_t n, std::size_t m) const {
        logger.assert_debug(n < rows_,
                            "Requested row % for a table with only % rows", n,
                            rows_);
        logger.assert_debug(
            m < columns_, "Requested column % for a table with only % columns",
            m, columns_);
        const type& entry = entries_[n + m * rows_];
        return entry;
    }

    /// Linear index
    const type& operator[](std::size_t i) const {
        logger.assert_debug(i < rows_ * columns_, "Linear index out of range.");
        const type& entry = entries_[i];
        return entry;
    }

    /// Fill the table with a constant value
    /// \param value The value.
    /// \return This.
    Table2D& operator=(const type& value) {
        fill(value);
        return *this;
    }

   private:
    std::vector<type> entries_;
    std::size_t columns_;
    std::size_t rows_;

    static void checkSize(std::size_t rows, std::size_t cols) {
        logger.assert_debug(rows <= std::numeric_limits<int>::max(),
                            "Too large row count");
        logger.assert_debug(cols <= std::numeric_limits<int>::max(),
                            "Too large column count");
    }
};

}  // namespace Utilities
}  // namespace hpgem

#endif  // HPGEM_TABLE2D_H

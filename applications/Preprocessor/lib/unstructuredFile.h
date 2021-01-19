/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2017, University of Twente
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

#ifndef HPGEM_APP_UNSTRUCTUREDFILE_H
#define HPGEM_APP_UNSTRUCTUREDFILE_H

/**
 * @brief wrapper class allow streamed input (using >>) on unstructured classes
 * The standard library assumes a certain structure in files that are streamed
 * in using >>. For example, it skips over whitespace and binary data that looks
 * like whitespace
 */
#include <utility>
#include <ios>
#include "Logger.h"

using namespace hpgem;

namespace Preprocessor {
template <typename inStream>
class UnstructuredInputStream {
    // when needed this can be reimplemented to look like moveFunction (and not
    // need the template), but this might degrade performance a bit
   public:
    using char_type = typename inStream::char_type;
    using traits_type = typename inStream::traits_type;
    using int_type = typename inStream::int_type;
    using pos_type = typename inStream::pos_type;
    using off_type = typename inStream::off_type;

    UnstructuredInputStream() = default;
    UnstructuredInputStream(inStream&& stream) : stream(std::move(stream)){};
    UnstructuredInputStream(const UnstructuredInputStream&) = delete;
    UnstructuredInputStream(UnstructuredInputStream&&) = default;
    ~UnstructuredInputStream() = default;

    UnstructuredInputStream& operator=(const UnstructuredInputStream&) = delete;
    UnstructuredInputStream& operator=(UnstructuredInputStream&&) = default;

    template <typename T>
    UnstructuredInputStream& operator>>(T& data) {
        stream.read(reinterpret_cast<char_type*>(&data), sizeof(data));
        return *this;
    }

    UnstructuredInputStream& read(char_type* s, std::streamsize count) {
        stream.read(s, count);
        return *this;
    }

    pos_type tellg() { return stream.tellg(); }

    void seekg(pos_type pos) { stream.seekg(pos); }

    void open(const std::string& filename,
              std::ios_base::openmode mode = std::ios_base::in) {
        stream.open(filename, mode);
    }

    bool is_open() const { return stream.is_open(); }

    bool good() const { return stream.good(); }

    std::basic_string<char_type> str() const { return stream.str(); }

    explicit operator bool() { return !!stream; }

    bool eof() { return stream.eof(); }

    // fixme: complete as needed
   private:
    inStream stream;
};
}  // namespace Preprocessor

#endif  // HPGEM_APP_UNSTRUCTUREDFILE_H

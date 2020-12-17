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
#ifndef HPGEM_FORTRANUNFORMATTEDFILE_H
#define HPGEM_FORTRANUNFORMATTEDFILE_H

#include <fstream>
#include <streambuf>

class FortranUnformattedFile {
   public:
    FortranUnformattedFile(std::string path);

    /**
     * Read a single record of known size, if the size is wrong it will fail
     *
     * @param size The expected size of the record in bytes.
     * @param buffer Buffer of at least size bytes used for reading
     * @return The number of bytes read
     */
    std::uint32_t readRawRecord(std::uint32_t size, char* buffer);

    /**
     * Read the size of the record without affecting the position of the object
     * @return The record size.
     */
    std::uint32_t peekRecordSize();

    /**
     * Skip the next record.
     */
    void skipRecord();

    /**
     * Skip the next record, will fail if the size is wrong
     * @param size The known size of the next record
     */
    void skipRecord(std::uint32_t size);

   private:
    // Internal buffer
    std::filebuf file_;

    /**
     * Read a record size marker
     * @param size [out] The marker size
     * @return Whether successful
     */
    bool readRecordSize(std::uint32_t& size);
};

#endif  // HPGEM_FORTRANUNFORMATTEDFILE_H

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
#include "FortranUnformattedFile.h"
#include "Logger.h"

using namespace hpgem;

FortranUnformattedFile::FortranUnformattedFile(const std::string &path) {
    std::filebuf *buf =
        file_.open(path, std::ios_base::in | std::ios_base::binary);
    logger.assert_debug(buf != nullptr, "Opening failed");
}

void FortranUnformattedFile::readRawRecord(std::uint32_t size, char *buffer) {
    // Read start marker
    std::uint32_t recordSize;
    logger.assert_debug(readRecordSize(recordSize),
                        "record size reading failed.");
    logger.assert_debug(recordSize == size,
                        "Record size does not match expected size.");
    // Read content
    std::streamsize read = file_.sgetn(buffer, recordSize);
    logger.assert_debug(read == recordSize, "Incomplete record read");
    // Read end marker
    std::uint32_t endSize;
    readRecordSize(endSize);
    logger.assert_debug(recordSize == endSize, "Record markers do not match");
}

std::uint32_t FortranUnformattedFile::peekRecordSize() {
    // store current place
    std::uint32_t recordSize;
    std::filebuf::off_type read =
        file_.sgetn(reinterpret_cast<char *>(&recordSize), sizeof(recordSize));
    // Reset reading pointer
    file_.pubseekoff(-read, std::ios_base::cur, std::ios_base::in);
    // Check that we actually read enough
    logger.assert_debug(read == sizeof(recordSize),
                        "Not enough bytes for reading a record size.");
    return read;
}

void FortranUnformattedFile::skipRecord(std::uint32_t size) {
    std::uint32_t recordSize;
    logger.assert_debug(readRecordSize(recordSize),
                        "Record size reading failed");
    logger.assert_debug(recordSize == size,
                        "Record size does not match expected size");
    // Skip past record content
    file_.pubseekoff(recordSize, std::ios_base::cur, std::ios_base::in);
    // Check
    std::uint32_t endRecordSize;
    logger.assert_debug(readRecordSize(endRecordSize),
                        "End record size reading failed");
    logger.assert_debug(recordSize == endRecordSize,
                        "Record sizes don't match");
}

bool FortranUnformattedFile::readRecordSize(std::uint32_t &size) {
    std::filebuf::off_type readSize =
        file_.sgetn(reinterpret_cast<char *>(&size), sizeof(size));
    if (readSize == sizeof(size)) {
        return true;
    } else {
        // Undo in case of failure
        file_.pubseekoff(-readSize, std::ios_base::cur, std::ios_base::in);
        return false;
    }
}

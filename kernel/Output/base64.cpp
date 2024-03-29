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

#include "base64.h"
#include "Logger.h"

#include <cstdlib>
#include <cstdint>

namespace hpgem {

bool Output::Detail::isBigEndian() {
    std::uint32_t test = 0x01020304;
    // Do NOT let the compiler touch the underlying binary data
    std::uint8_t* pFirstByte = reinterpret_cast<std::uint8_t*>(&test);
    return *pFirstByte == 1;
}

std::string Output::Detail::toBase64(void* rawData, std::size_t len) {
    // Break out fast, even before checking the pointer
    if (len == 0) {
        return "";
    }
    logger.assert_debug(rawData != nullptr, "no raw data passed");
    const unsigned char* cRawData = (const unsigned char*)rawData;
    std::string result;
    result.reserve((len + 2) / 3 * 4);
    while (len > 2)  // 3 or larger
    {
        // ugly bitmasking to pick out groups of 6 bits
        result.push_back(base64Encode[(cRawData[0] >> 2) & 0x3F]);
        result.push_back(base64Encode[((cRawData[0] << 4) & 0x30) |
                                      ((cRawData[1] >> 4) & 0x0F)]);
        result.push_back(base64Encode[((cRawData[1] << 2) & 0x3C) |
                                      ((cRawData[2] >> 6) & 0x03)]);
        result.push_back(base64Encode[(cRawData[2]) & 0x3F]);
        len -= 3;
        cRawData += 3;
    }
    if (len == 2) {
        result.push_back(base64Encode[(cRawData[0] >> 2) & 0x3F]);
        result.push_back(base64Encode[((cRawData[0] << 4) & 0x30) |
                                      ((cRawData[1] >> 4) & 0x0F)]);
        result.push_back(base64Encode[((cRawData[1] << 2) & 0x3C)]);
        result.push_back('=');
    }
    if (len == 1) {
        result.push_back(base64Encode[(cRawData[0] >> 2) & 0x3F]);
        result.push_back(base64Encode[((cRawData[0] << 4) & 0x30)]);
        result.push_back('=');
        result.push_back('=');
    }
    return result;
}

}  // namespace hpgem

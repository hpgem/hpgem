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
#ifndef HPGEM_COMMANDLINEHELPERS_H
#define HPGEM_COMMANDLINEHELPERS_H

#include "LinearAlgebra/SmallVector.h"

#include <string>

namespace hpgem {
namespace Base {

/// Parse DIM comma separated numbers as the coordinates of a point.
/// \tparam DIM The dimension of the point
/// \param pointString The string containing the point coordinates
/// \param start The starting index in pointString
/// \param point The point (out)
/// \return The first index in pointString after the number.
template <std::size_t DIM>
std::size_t parsePoint(const std::string& pointString, std::size_t start,
                       LinearAlgebra::SmallVector<DIM>& point) {
    for (std::size_t i = 0; i < DIM; ++i) {
        logger.assert_always(start < pointString.size(),
                             "Last point has too few coordinates % expecting %",
                             i, DIM);
        std::size_t len = 0;
        try {
            point[i] = std::stod(pointString.substr(start), &len);
        } catch (const std::invalid_argument&) {
            // No parse, i.e. len == 0
            logger.fail("Number point parsing failed at '" +
                        pointString.substr(start) + "', expected a coordinate");
        }
        start += len;
        if (i < DIM - 1) {
            logger.assert_always(start < pointString.size(),
                                 "Reached the end of the string, expecting "
                                 "a comma and more coordinates");
            logger.assert_always(pointString[start] == ',',
                                 "Expected a comma to separate coordinates not "
                                 "'%' at position %",
                                 pointString[start], start);
            // Skip the comma, space, whatever that ended the point
            start++;
        }
    }
    return start;
}

}  // namespace Base
}  // namespace hpgem

#endif  // HPGEM_COMMANDLINEHELPERS_H

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
#ifndef HPGEM_ZONESTRUCTUREDESCRIPTION_H
#define HPGEM_ZONESTRUCTUREDESCRIPTION_H

#include <vector>
#include <regex>
#include <cmath>
#include "Logger.h"

#include "StructureDescription.h"

namespace DGMax {
/// Structure definition based on the zone information of the mesh
class ZoneInfoStructureDefinition : public StructureDescription {

   public:
    /// Create a structure definition based on the zone information in the mesh.
    ///
    /// This will determine the material properties of each element based on the
    /// zone it belongs to. The name of the zone will be matched in order to the
    /// given regexes. The material information (epsilon) of the first matching
    /// regex will be used. If no regex matches, it will use a default epsilon
    /// or generate an error when the default is nan.
    ///
    /// \param regexes vector of regexes used to match zone names
    /// \param epsilons vector with corresponding values of epsilon
    /// \param defaultEpsilon Optional default epsilon when no regex matches
    ///   (nan, default, will generate an error)
    ZoneInfoStructureDefinition(std::vector<std::regex> regexes,
                                std::vector<double> epsilons,
                                double defaultEpsilon = std::nan(""))
        : regexes_(std::move(regexes)),
          epsilons_(std::move(epsilons)),
          defaultEpsilon_(defaultEpsilon) {
        logger.assert_always(regexes_.size() == epsilons_.size(),
                             "Regexes not matching material information");
    };

    ElementInfos* createElementInfo(const Base::Element* element) final;

   private:
    std::vector<std::regex> regexes_;
    std::vector<double> epsilons_;
    double defaultEpsilon_;
};

}  // namespace DGMax

#endif  // HPGEM_ZONESTRUCTUREDESCRIPTION_H

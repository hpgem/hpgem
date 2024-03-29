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
#ifndef HPGEM_ZONEINFORMATION_H
#define HPGEM_ZONEINFORMATION_H

#include <string>
#include <regex>

namespace hpgem {
namespace Base {

/**
 * Description of a zone or region in the mesh
 */
class Zone {

   public:
    Zone(const std::string& name, std::size_t zoneId)
        : name_(name), zoneId_(zoneId){};

    /// The identifier for this zone in the mesh it belongs to.
    ///
    /// It should satisfy mesh.getZones()[i].getZoneId() == i
    /// \return Identifier for this zone in the mesh.
    std::size_t getZoneId() const { return zoneId_; }

    /// \return Name of the zone
    const std::string& getName() const { return name_; }

    /// Test whether the name matches the given regex
    bool matchesName(std::regex regex) const;

    /// Find the index of the first matching regex or -1 if none.
    int matchName(const std::vector<std::regex>& regexes) const;

   private:
    std::size_t zoneId_;
    std::string name_;
};
}  // namespace Base
}  // namespace hpgem

#endif  // HPGEM_ZONEINFORMATION_H

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

#include "../catch.hpp"

#include "Utils/ZoneStructureDescription.h"
#include "Base/Element.h"
#include "Base/Zone.h"

using namespace hpgem;

const std::vector<std::size_t> dummyNodesIndices = {0, 1};
std::vector<hpgem::Geometry::PointPhysical<1>> allNodes = {
    Geometry::PointPhysical<1>({0.0}), Geometry::PointPhysical<1>({1.0})};
const Base::Element::CollectionOfBasisFunctionSets bfs;

std::unique_ptr<Base::Element> createDummyElement(Base::Zone& zone) {
    static std::size_t id;
    return std::make_unique<Base::Element>(dummyNodesIndices, &bfs, allNodes, 0,
                                           0, id++, zone);
}

std::vector<std::regex> fromStrings(std::vector<std::string> input) {
    std::vector<std::regex> result;
    std::transform(input.begin(), input.end(), std::back_inserter(result),
                   [](const std::string& inp) { return std::regex(inp); });
    return result;
}

TEST_CASE("Two zones", "[ZoneStructureDescription]") {
    using namespace hpgem::Base;
    using namespace DGMax;

    Zone zone1("foo", 1);
    Zone zone2("fob", 2);
    auto element1 = createDummyElement(zone1);
    auto element2 = createDummyElement(zone2);

    // Intentionally overlapping patterns to see that the first matching pattern
    // is used.
    ZoneInfoStructureDefinition def(fromStrings({"fo*", "c.*", "f.*"}),
                                    {1.0, 2.0, 3.0});

    ElementInfos* infos1 = def.createElementInfo(element1.get());
    REQUIRE(infos1->getPermittivity() == 1.0);
    ElementInfos* infos2 = def.createElementInfo(element2.get());
    REQUIRE(infos2->getPermittivity() == 3.0);
    delete infos1;
    delete infos2;
}

TEST_CASE("Default value", "[ZoneStructureDescription]") {
    using namespace hpgem::Base;
    using namespace DGMax;

    Zone zone("foo", 1);
    auto element = createDummyElement(zone);

    ZoneInfoStructureDefinition def(fromStrings({"g.*"}), {1.0}, 13);
    ElementInfos* infos = def.createElementInfo(element.get());
    REQUIRE(infos->getPermittivity() == 13);
    delete infos;
}

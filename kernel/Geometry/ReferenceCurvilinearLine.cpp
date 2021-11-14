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
#include "ReferenceCurvilinearLine.h"

#include <map>

#include "ReferenceCurvilinearElement_Impl.h"
#include "ReferenceLine.h"
#include "ReferencePoint.h"

namespace hpgem {
namespace Geometry {

ReferenceCurvilinearLine& ReferenceCurvilinearLine::getReferenceLagrangeLine(
    std::size_t order) {
    logger.assert_always(order != 0, "A zero-th order line does not exist");
    logger.assert_always(order != 1,
                         "Use the regular ReferenceLine for first order");
    // Multiton
    static std::map<std::size_t, ReferenceCurvilinearLine*> lines;

    ReferenceCurvilinearLine*& line = lines[order];
    if (line == nullptr) {
        line = new ReferenceCurvilinearLine(order);
    }
    return *line;
}

std::string getReferenceLagrangeLineName(std::size_t order) {
    std::stringstream name;
    name << "Lagrange-ReferenceLine-" << order;
    return name.str();
}

std::vector<Geometry::PointReference<1>> createReferencePoints(
    std::size_t order) {
    // order + 1 equally spaced points from -1 to 1
    double h = 2.0 / order;
    std::vector<Geometry::PointReference<1>> result(order + 1);
    for (std::size_t i = 0; i < order + 1; ++i) {
        result[i] = {-1.0 + i * h};
    }
    return result;
}

ReferenceCurvilinearLine::ReferenceCurvilinearLine(std::size_t order)
    : ReferenceCurvilinearElement<1>(&ReferenceLine::Instance(),
                                  // codim 1 = points -> Not included
                                  std::vector<ReferenceGeometry*>(0),
                                  // No codim 2
                                  std::vector<ReferenceGeometry*>(0),
                                  // Lagrange points
                                  createReferencePoints(order),
                                  getReferenceLagrangeLineName(order), order) {}

}  // namespace Geometry
}  // namespace hpgem
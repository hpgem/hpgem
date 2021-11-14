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
#ifndef HPGEM_REFERENCECURVILINEARTRIANGLE_H
#define HPGEM_REFERENCECURVILINEARTRIANGLE_H

#include "ReferenceCurvilinearElement.h"

namespace hpgem {
namespace Geometry {
/// Lagrange triangle of arbitrary order based on the standard reference
/// triangle.
///
/// The vertices are generated such that if the coordinates are listed (y,x)
/// they are in lexicographical order. That way the first order triangle has
/// ordering (0,0) - (1,0) - (0,1) =(x,y), matching that of the regular triangle
class ReferenceCurvilinearTriangle : public ReferenceCurvilinearElement<2> {
   public:
    /// Reverse computation finding the order from the number of points.
    /// Negative values are used to denote that no such element exists
    static int getOrderFromPoints(std::size_t numberOfPoints);
    static ReferenceCurvilinearTriangle& getReferenceLagrangeTriangle(
        std::size_t order);

   private:
    explicit ReferenceCurvilinearTriangle(std::size_t order);
    static std::vector<Geometry::PointReference<2>> createPoints(
        std::size_t order);
};
}  // namespace Geometry
}  // namespace hpgem

#endif  // HPGEM_REFERENCECURVILINEARTRIANGLE_H

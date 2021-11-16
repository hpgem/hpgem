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
#ifndef HPGEM_REFERENCEGEOMETRYFACTORY_H
#define HPGEM_REFERENCEGEOMETRYFACTORY_H

#include "ReferenceGeometry.h"

namespace hpgem {
namespace Geometry {

class ReferenceGeometryFactory {
   public:
    static ReferenceGeometryFactory& Instance() {
        static ReferenceGeometryFactory instance;
        return instance;
    }

    /// Lookup the ReferenceGeometry by the dimension of the shape and the
    /// number of vertices.
    ///
    /// \param dimension The dimension of the geometry
    /// \param numberOfPoints The number of vertices
    /// \return A reference to the geometry. Will terminate the program if it
    /// does not exist.
    ReferenceGeometry& getGeometry(std::size_t dimension,
                                   std::size_t numberOfPoints);

   private:
    ReferenceGeometryFactory() = default;
    ReferenceGeometry& getGeometry0(std::size_t numberOfPoints);
    ReferenceGeometry& getGeometry1(std::size_t numberOfPoints);
    ReferenceGeometry& getGeometry2(std::size_t numberOfPoints);
    ReferenceGeometry& getGeometry3(std::size_t numberOfPoints);
    ReferenceGeometry& getGeometry4(std::size_t numberOfPoints);
    /// Cache for higher order 2D elements, indexed by numberOfPoints.
    std::map<std::size_t, ReferenceGeometry*> cached2DGeometries_;
};
}  // namespace Geometry
}  // namespace hpgem

#endif  // HPGEM_REFERENCEGEOMETRYFACTORY_H

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
#ifndef HPGEM_MESHCELL_H
#define HPGEM_MESHCELL_H

#include <vector>
#include <map>
#include "Base/Face.h"
#include "Base/MeshManipulator.h"
#include "Geometry/PointReference.h"
#include "Geometry/Mappings/MappingReferenceToReference.h"

namespace DGMax {
namespace Utils {

using namespace hpgem;

/**
 * A Polytope in a certain dimension
 * @tparam DIM The dimension
 */
template <std::size_t DIM>
class Polytope {
   public:
    virtual ~Polytope() = default;
    virtual const std::vector<Geometry::PointReference<DIM>>& getVertices()
        const = 0;
};

class TesselationElement {
   public:
    virtual ~TesselationElement() = default;
    /** Indices of the corner points with respect to the tesselation **/
    virtual const std::vector<std::size_t>& getCornerIds() const = 0;
    virtual const Geometry::ReferenceGeometry* getReferenceGeometry() const = 0;

    //    /**
    //     * Mapping from local coordiantes to coordinates on the teselation
    //     * @return
    //     */
    //    virtual const Geometry::MappingReferenceToReference<0>&
    //        toTesselationMapping() const;
};

template <std::size_t DIM>
class Tesselation : public Polytope<DIM> {
   public:
    virtual const std::vector<Geometry::PointReference<DIM>>& getSubVertices()
        const = 0;
};

template <std::size_t DIM>
class DomainFacet {
   public:
    const Tesselation<DIM - 1>& getFacetTesselation() const;

    /**
     * Correspondence between elements on the tesselation and faces. As the
     * ordering of the nodes may differ a reordering is needed.
     * @return
     */
    const std::map<const TesselationElement*,
                   std::pair<const Base::Face*, std::size_t>>&
        getFaceMapping() const;
    //    /**
    //     * Mapping from the tesselation coordinates to domain coordinates
    //     * @return
    //     */
    //    Geometry::MappingReferenceToReference<1>&
    //    getTesselationToDomainMapping()
    //        const;
};

template <std::size_t DIM>
class MeshDomain : public Polytope<DIM> {
   public:
    MeshDomain(std::shared_ptr<Base::MeshManipulator<DIM>> mesh) : mesh(mesh){};

    const Base::MeshManipulator<DIM>& getMesh() { return mesh; }
    const std::vector<Geometry::PointReference<DIM>>& getVertices()
        const final {
        return vertices;
    }
    const std::vector<DomainFacet<DIM>> getFacets() const final {
        return facets;
    }

   private:
    std::shared_ptr<Base::MeshManipulator<DIM>> mesh;
    std::vector<Geometry::PointReference<DIM>> vertices;
    std::vector<DomainFacet<DIM>> facets;
};

}  // namespace Utils
}  // namespace DGMax

#endif  // HPGEM_MESHCELL_H

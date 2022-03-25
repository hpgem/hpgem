/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2022, University of Twente
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
#ifndef HPGEM_FLUXFACETS_H
#define HPGEM_FLUXFACETS_H

#include <memory>
#include <vector>
#include <map>

#include <LinearAlgebra/SmallVector.h>
#include <Base/MeshManipulator.h>

#include <ProblemTypes/AbstractHarmonicResult.h>
#include <ProblemTypes/FieldPattern.h>

namespace DGMax {

/**
 * Facets of the mesh through which the flux needs to be computed
 */
class FluxFacets {

   public:
    FluxFacets(const Base::MeshManipulatorBase& mesh);

    template <std::size_t dim>
    std::vector<LinearAlgebra::SmallVector<4>> computeFluxes(
        DGMax::AbstractHarmonicResult<dim>& result, double wavenumber,
        const DGMax::FieldPattern<dim>* background) const;

    // Public to allow writing the header
    /**
     * Face for which an energy flux needs to be computed
     */
    struct FluxFace {
        /**
         * The face
         */
        Base::Face* face;
        /**
         * The side from which it needs to be computed
         */
        Base::Side side;
    };

    const std::vector<std::string>& getFacetNames() const {
        return facetNames_;
    }

   private:
    /**
     * Pair of indices for the zones adjacent to a face.
     *  - For external faces the second index is -1
     *  - For internal faces the first index should be smaller than the second
     */
    using FaceZoneIndexPair = std::pair<int, int>;

    /**
     * Generate the names of the facets
     */
    void generateFacetNames(const Base::MeshManipulatorBase& mesh);
    void generateZoneOrdering(std::size_t numberOfZones);
    FaceZoneIndexPair getZoneIndexPair(const Base::Face& face) const;

    /**
     * For a consistent output and communication we order the zone pairs of the
     * faces.
     */
    std::map<FaceZoneIndexPair, std::size_t> faceIndexOrdering_;
    /**
     * Mapping from zone pair index to the corresponding name
     */
    std::vector<std::string> facetNames_;

    std::map<FaceZoneIndexPair, std::vector<FluxFace>> facets;
};

}  // namespace DGMax

#endif  // HPGEM_FLUXFACETS_H

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
#include "FluxFacets.h"

namespace DGMax {

FluxFacets::FluxFacets(const Base::MeshManipulatorBase& mesh) {
    // When using more than 1 MPI rank, we would need to aggregate data. But the
    // first processor may not know all the regions and possible boundaries.
    logger.assert_always(Base::MPIContainer::Instance().getNumProcessors() == 1,
                         "FluxFacets not suitable for more than 1 MPI rank");
    for (Base::Face* face : mesh.getFacesList()) {
        using Base::Side;
        if (!face->isOwnedByCurrentProcessor()) {
            continue;
        }
        const Base::Zone& leftZone = face->getPtrElementLeft()->getZone();
        FluxFace fface;
        fface.face = face;
        std::string facetName;

        if (face->isInternal()) {
            const Base::Zone& rightZone = face->getPtrElementRight()->getZone();
            if (leftZone.getZoneId() == rightZone.getZoneId()) {
                // Only interested in zone boundaries
                continue;
            }
            // Stable direction
            fface.side = leftZone.getZoneId() < rightZone.getZoneId()
                             ? Side::LEFT
                             : Side::RIGHT;
            // Name
            std::stringstream name;
            name << "flux-"
                 << (fface.side == Side::LEFT ? leftZone : rightZone).getName()
                 << "--"
                 << (fface.side == Side::LEFT ? rightZone : leftZone).getName();
            facetName = name.str();
        } else {
            // Outward fluxes are always included
            fface.side = Side::LEFT;
            std::stringstream name;
            name << "flux-" << leftZone.getName();
            facetName = name.str();
        }
        facets[facetName].push_back(fface);
    }
}

template <std::size_t dim>
std::map<std::string, LinearAlgebra::SmallVector<4>> FluxFacets::computeFluxes(
    DGMax::AbstractHarmonicResult<dim>& result, double wavenumber,
    const DGMax::FieldPattern<dim>* background) const {
    std::map<std::string, LinearAlgebra::SmallVector<4>> fluxes;
    for (const auto& facet : facets) {
        LinearAlgebra::SmallVector<4> flux = {};
        for (const auto& face : facet.second) {
            flux += result.computeEnergyFlux(*face.face, face.side, wavenumber,
                                             background);
        }
        fluxes[facet.first] += flux;
    }
    // TODO: MPI
    return fluxes;
}

// Explicit instantiation

template
std::map<std::string, LinearAlgebra::SmallVector<4>>
    FluxFacets::computeFluxes(DGMax::AbstractHarmonicResult<2>& result, double wavenumber, const DGMax::FieldPattern<2>* background) const;
template
std::map<std::string, LinearAlgebra::SmallVector<4>>
    FluxFacets::computeFluxes(DGMax::AbstractHarmonicResult<3>& result, double wavenumber, const DGMax::FieldPattern<3>* background) const;

}
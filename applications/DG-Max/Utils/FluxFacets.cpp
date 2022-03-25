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
        FluxFace fface;
        fface.face = face;

        const Base::Zone& leftZone = face->getPtrElementLeft()->getZone();

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
        } else {
            // Outward fluxes are always included
            fface.side = Side::LEFT;
        }
        FaceZoneIndexPair zoneIndexPair = getZoneIndexPair(*face);
        facets[zoneIndexPair].push_back(fface);
    }
    generateZoneOrdering(mesh.getZones().size());
    generateFacetNames(mesh);
}

FluxFacets::FaceZoneIndexPair FluxFacets::getZoneIndexPair(
    const Base::Face& face) const {
    FaceZoneIndexPair result;
    result.first = face.getPtrElementLeft()->getZone().getZoneId();
    if (face.isInternal()) {
        result.second = face.getPtrElementRight()->getZone().getZoneId();
        if (result.first > result.second) {
            std::swap(result.first, result.second);
        }
    } else {
        result.second = -1;
    }
    return result;
}

void FluxFacets::generateZoneOrdering(std::size_t numberOfZones) {
    // Converts the index pair into a linear index using lexicographical
    // ordering.

    // For each lexicographical index, whether such facet is available.
    std::vector<int> ordering(numberOfZones * (numberOfZones + 1));
    for (const auto& facet : facets) {
        const FaceZoneIndexPair& zoneIndexPair = facet.first;
        // For each zone it can be connected to one of the N zones, or to the
        // outside.
        ordering[zoneIndexPair.first * (numberOfZones + 1) +
                 zoneIndexPair.second + 1] = 1;
    }
#ifdef HPGEM_USE_MPI
    // Share the information between processes
    MPI_Allreduce(MPI_IN_PLACE, ordering.data(), ordering.size(), MPI_INT,
                  MPI_MAX, MPI_COMM_WORLD);
#endif
    // Convert the 'is presence' into an ordering index.
    std::accumulate(ordering.begin(), ordering.end(), -1);
    // Extract the ordering for the local facets.
    for (const auto& facet : facets) {
        const FaceZoneIndexPair& zoneIndexPair = facet.first;
        std::size_t index = zoneIndexPair.first * (numberOfZones + 1) +
                            zoneIndexPair.second + 1;
        faceIndexOrdering_[zoneIndexPair] = ordering[index];
    }
}

void FluxFacets::generateFacetNames(const Base::MeshManipulatorBase& mesh) {
    facetNames_.resize(faceIndexOrdering_.size());
    for (const auto& facetPair : faceIndexOrdering_) {
        const FaceZoneIndexPair& zonePair = facetPair.first;
        std::stringstream name;
        name << "flux-";
        name << mesh.getZones()[zonePair.first]->getName();
        if (zonePair.second != -1) {
            name << "--" << mesh.getZones()[zonePair.second]->getName();
        }
        facetNames_[facetPair.second] = name.str();
    }
}

template <std::size_t dim>
std::vector<LinearAlgebra::SmallVector<4>> FluxFacets::computeFluxes(
    DGMax::AbstractHarmonicResult<dim>& result, double wavenumber,
    const DGMax::FieldPattern<dim>* background) const {
    const std::size_t numberOfFacets = faceIndexOrdering_.size();
    // Flat representation of the output for MPI communication
    std::vector<double> flatFluxes(4 * numberOfFacets);
    for (const auto& facet : facets) {
        LinearAlgebra::SmallVector<4> flux = {};
        for (const auto& face : facet.second) {
            flux += result.computeEnergyFlux(*face.face, face.side, wavenumber,
                                             background);
        }
        auto indexIter = faceIndexOrdering_.find(facet.first);
        logger.assert_always(indexIter != faceIndexOrdering_.end(),
                             "Facet index not found");
        std::size_t index = 4 * indexIter->second;
        for (std::size_t i = 0; i < 4; ++i) {
            flatFluxes[index + i] = flux[i];
        }
    }
#ifdef HPGEM_USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, flatFluxes.data(), flatFluxes.size(),
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    // Convert from flat to result
    std::vector<LinearAlgebra::SmallVector<4>> allFluxes(numberOfFacets);
    for (std::size_t i = 0; i < numberOfFacets; ++i) {
        LinearAlgebra::SmallVector<4> fluxes;
        for (std::size_t j = 0; j < 4; ++j) {
            fluxes[j] = flatFluxes[4 * i + j];
        }
        allFluxes[i] = fluxes;
    }
    return allFluxes;
}

// Explicit instantiation

template std::vector<LinearAlgebra::SmallVector<4>> FluxFacets::computeFluxes(
    DGMax::AbstractHarmonicResult<2>& result, double wavenumber,
    const DGMax::FieldPattern<2>* background) const;
template std::vector<LinearAlgebra::SmallVector<4>> FluxFacets::computeFluxes(
    DGMax::AbstractHarmonicResult<3>& result, double wavenumber,
    const DGMax::FieldPattern<3>* background) const;

}  // namespace DGMax
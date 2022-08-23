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
#include "PMLTransmission.h"

#include "PMLElementInfos.h"

namespace DGMax {

template <std::size_t dim>
PMLTransmission<dim>::PMLTransmission(
    const hpgem::Base::MeshManipulator<dim>& mesh) {
    // Relevant area by zoneId, for communicating by MPI
    const std::size_t numZones = mesh.getZones().size();
    std::vector<double> surfaceAreaByZone(numZones);
    std::vector<std::vector<Base::Face*>> localFacetsByZone(numZones);

    for (Base::Face* face : mesh.getFacesList()) {
        if (face->isInternal() || !face->isOwnedByCurrentProcessor()) {
            continue;
        }
        ElementInfos* elementInfos =
            ElementInfos::get(face->getPtrElementLeft());
        auto* pmlInfos = dynamic_cast<PMLElementInfos<dim>*>(elementInfos);
        if (pmlInfos == nullptr) {
            continue;
        }
        // Test if this face is at the far end of the PML
        auto directions = pmlInfos->getDampeningDirections();
        auto normal = face->getNormalVector(
            face->getReferenceGeometry()->getCenter().castDimension<dim - 1>());
        double normalization = normal.l2Norm();
        normal /= normalization;
        if (normal * directions < 1e-3) {
            // Not dampening in this direction
            continue;
        }
        double area = face->getReferenceGeometry()->measure() * normalization;
        auto zoneId = face->getPtrElementLeft()->getZone().getZoneId();
        surfaceAreaByZone[zoneId] += area;
        localFacetsByZone[zoneId].push_back(face);
    }
#ifdef HPGEM_USE_MPI
    // Sum all the area's
    MPI_Allreduce(MPI_IN_PLACE, surfaceAreaByZone.data(),
                  surfaceAreaByZone.size(), MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
#endif
    // surfaceAreaByZone will now have non-zero entries for all the zones
    // that we need to include in the output. Now we need to fill the structures
    // used for the actual computation.
    for (std::size_t i = 0; i < numZones; ++i) {
        if (surfaceAreaByZone[i] == 0.0) {
            // Hard comparison against starting value 0.0 to see if any face
            // contributed.
            continue;
        }
        // Facets included in the output
        surfaceArea_.push_back(surfaceAreaByZone[i]);
        resultFacets_.push_back(std::move(localFacetsByZone[i]));
        facetNames_.push_back(mesh.getZones()[i]->getName());
    }
}

template <std::size_t dim>
std::vector<double> PMLTransmission<dim>::pmlTransmission(
    AbstractHarmonicResult<dim>& result) const {
    std::vector<double> transmission(facetNames_.size());
    for (std::size_t i = 0; i < facetNames_.size(); ++i) {
        double totalField = 0.0;
        for (Base::Face* face : resultFacets_[i]) {
            totalField += result.computeFieldL2Integral(*face, Base::Side::LEFT);
        }
        transmission[i] = totalField / surfaceArea_[i];
    }
    return transmission;
}

template class PMLTransmission<2>;
template class PMLTransmission<3>;

}  // namespace DGMax
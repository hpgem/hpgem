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
#ifndef HPGEM_TRANSFORMATIONCONNECTOR_H
#define HPGEM_TRANSFORMATIONCONNECTOR_H

#include "mesh/Mesh.h"
#include "MergePlan.h"

#include <set>
#include <LinearAlgebra/SmallVector.h>
#include <LinearAlgebra/SmallMatrix.h>

namespace Preprocessor {

template <std::size_t dim>
void translationConnect(Mesh<dim>& mesh,
                        LinearAlgebra::SmallVector<dim> translation,
                        LinearAlgebra::SmallMatrix<dim, dim> transform = LinearAlgebra::SmallMatrix<dim, dim>::identity()) {
    // Stage 1: Find coordinates on the boundary
    std::set<CoordId> boundaryCoords;
    for (const MeshEntity<dim - 1, dim>& face : mesh.getFaces()) {
        if (face.getNumberOfElements() != 1) {
            // Not a boundary face, also excludes unused faces
            continue;
        }
        const Element<dim>& element = face.getElement(0);
        // Find the coord id for the nodes on the face
        for (EntityLId localNodeId :
             element.template getLocalIncidenceListAsIndices<0>(face)) {
            boundaryCoords.insert(element.getCoordinateIndex(localNodeId));
        }
    }

    // Stage 2: Find connected coordinates
    // We don't have any efficient datastructures for quickly matching the nodes
    // Hence perform an O(N^2) loop to match them
    std::map<CoordId, CoordId> pairing;
    for (CoordId coord1 : boundaryCoords) {
        LinearAlgebra::SmallVector<dim> startCoord = mesh.getCoordinate(coord1);
        auto translatedCoord = transform * startCoord + translation;
        for (CoordId coord2 : boundaryCoords) {
            auto diff = translatedCoord - mesh.getCoordinate(coord2);
            if (diff.l2NormSquared() < 1e-16) {
                // Matching coordinates
                pairing[coord1] = coord2;
            }
        }
    }
    // Stage 3: Actual merging
    MergePlan<dim> mergePlan = MergePlan<dim>::computeMergePlan(&mesh, pairing);
    mergePlan.executeMerge();
}

}  // namespace Preprocessor

#endif  // HPGEM_TRANSFORMATIONCONNECTOR_H

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
#ifndef HPGEM_PMLTRANSMISSION_H
#define HPGEM_PMLTRANSMISSION_H

#include <memory>
#include <vector>
#include <map>

#include <LinearAlgebra/SmallVector.h>
#include <Base/MeshManipulator.h>

#include <ProblemTypes/AbstractHarmonicResult.h>

namespace DGMax {

/**
 * Computes the transmission through the PML in an actual solution.
 *
 * As proxy for the actual transmission we consider the boundary faces where the
 * normal is in a dampening direction. On these faces we compute the L2-norm of
 * the solution and normalize this by the surface area.
 *
 * @tparam dim The dimension of the mesh
 */
 template<std::size_t dim>
class PMLTransmission {
   public:
    PMLTransmission(const hpgem::Base::MeshManipulator<dim>& mesh);

    /**
     * Compute the transmission through the PML for a given solution.
     *
     * @param result The computed result
     * @return For each of the included faces the integral of the square of the
     * solution.
     */
    std::vector<double> pmlTransmission(
        AbstractHarmonicResult<dim>& result) const;

    const std::vector<std::string>& getFacetNames() const {
        return facetNames_;
    }

   private:
    /**
     * For each of the facets in the output, the local faces that contribute to
     * that facet.
     */
    std::vector<std::vector<Base::Face*>> resultFacets_;
    /**
     * Names of the facets as used in the output.
     */
    std::vector<std::string> facetNames_;
    /**
     * Surface area of the facets (for normalization)
     */
    std::vector<double> surfaceArea_;
};

}  // namespace DGMax

#endif  // HPGEM_PMLTRANSMISSION_H

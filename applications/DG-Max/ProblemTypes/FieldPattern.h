/*
This file forms part of hpGEM. This package has been developed over a number of
years by various people at the University of Twente and a full list of
contributors can be found at http://hpgem.org/about-the-code/team

This code is distributed using BSD 3-Clause License. A copy of which can found
below.


Copyright (c) 2021, University of Twente
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef HPGEM_FIELDPATTERN_H
#define HPGEM_FIELDPATTERN_H

#include "LinearAlgebra/SmallVector.h"
#include "Geometry/PointPhysical.h"
#include "MaterialTensor.h"


namespace DGMax {

/**
 * Describes a field pattern for the electric field.
 *
 * This describes a field pattern E(x) for the electric field. It includes the
 * required Curl's to compute the harmonic source term.
 *
 * @tparam dim The dimension of the space
 */
template <std::size_t dim>
class FieldPattern {
   public:
    using VecC = hpgem::LinearAlgebra::SmallVectorC<dim>;
    using VecR = hpgem::LinearAlgebra::SmallVector<dim>;
    using PPhys = hpgem::Geometry::PointPhysical<dim>;

    virtual ~FieldPattern() = default;

    /**
     * Electric field vector
     * @param p The point to at which to compute it.
     */
    virtual VecC field(const PPhys& p) const = 0;
    /**
     * Curl of the electric field
     * @param p The point at which to compute it.
     */
    virtual VecC fieldCurl(const PPhys& p) const = 0;
    /**
     * Double curl of the field, with material constant: Curl (Mu^{-1} Curl E).
     *
     * The material constant is assumed to be locally piecewise constant, and
     * thus have zero derivatives for the outer curl.
     * @param p The point at which to compute it.
     */
    virtual VecC fieldDoubleCurl(const PPhys& p,
                                 const MaterialTensor& material) const = 0;

    /**
     * Compute local time averaged energy flux (Poynting vector) for a harmonic
     * field
     * @param p The point to compute at
     * @param material The material tensor for the permeability
     * @param omega The time harmonic frequency for scaling the result 1/omega
     */
    VecR localFlux(const PPhys& p, const MaterialTensor& material,
                   double omega) const {
        // Poynting vector =
        //    Re(E x H*)
        //  = Re(i/omega E x 1/mu* Curl E*)
        //  = 1/omega Im(E x 1/mu* Curl E*)
        auto Efield = field(p);
        auto EfieldCurlConj = fieldCurl(p).conj();
        return 1 / omega *
               hpgem::LinearAlgebra::leftDoubledCrossProduct(
                   Efield, material.adjoint().applyCurl(EfieldCurlConj))
                   .imag();
    }
};

}  // namespace DGMax

#endif  // HPGEM_FIELDPATTERN_H

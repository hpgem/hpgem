/*
This file forms part of hpGEM. This package has been developed over a number of
years by various people at the University of Twente and a full list of
contributors can be found at http://hpgem.org/about-the-code/team

This code is distributed using BSD 3-Clause License. A copy of which can found
below.


Copyright (c) 2021, Univesity of Twenete
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
#ifndef HPGEM_PLANEWAVE_H
#define HPGEM_PLANEWAVE_H

#include <complex>
#include <LinearAlgebra/SmallVector.h>
#include <Geometry/PointPhysical.h>

#include "../Material.h"
#include "ProblemTypes/FieldPattern.h"

namespace DGMax {

/// A plane wave E0 exp(ik.x) (with E0 . k == 0)
///
/// Note: This is just a field profile, it is not associated with a material or
/// frequency.
///
/// \tparam dim The dimension of the domain (2 or 3)
template <std::size_t dim>
class PlaneWave : public FieldPattern<dim> {
   public:
    PlaneWave() : k_(), E0_(), phase_(0.0) {}
    PlaneWave(LinearAlgebra::SmallVector<dim> k,
              LinearAlgebra::SmallVectorC<dim> E0, double phase)
        : k_(k), E0_(E0), phase_(phase) {

        logger.assert_debug(
            std::abs(E0 * k) <= 1e-8 * (k.l2Norm() * E0.l2Norm()),
            "Non orthogonal wavevector and field.");
    };

    static PlaneWave onDispersion(double omega,
                                  LinearAlgebra::SmallVector<dim> khat,
                                  LinearAlgebra::SmallVectorC<dim> E0,
                                  DGMax::Material material, double phase) {
        // Set the length of the wave vector to be on dispersion relation
        khat *= omega * material.getRefractiveIndex() / khat.l2Norm();
        return PlaneWave(khat, E0, phase);
    }

    LinearAlgebra::SmallVectorC<dim> field(
        const Geometry::PointPhysical<dim>& p) const final {
        return E0_ * phaseFactor(p);
    }

    LinearAlgebra::SmallVectorC<dim> fieldCurl(
        const Geometry::PointPhysical<dim>& p) const final {
        // To ensure the complex valued cross product is used
        LinearAlgebra::SmallVectorC<dim> kc = k_;
        using namespace std::complex_literals;

        return kc.crossProduct(E0_) * (1i * phaseFactor(p));
    }

    LinearAlgebra::SmallVectorC<dim> fieldDoubleCurl(const Geometry::PointPhysical<dim>& p,
                         const MaterialTensor& material) const final {
        // To ensure the complex valued cross product is used
        LinearAlgebra::SmallVectorC<dim> kc = k_;
        using namespace std::complex_literals;
        // -1.0 is from two derivatives, each giving a factor 1i
        return LinearAlgebra::leftDoubledCrossProduct(
                   kc, material.applyCurl(kc.crossProduct(E0_))) *
               (-1.0 * phaseFactor(p));
    }

    const LinearAlgebra::SmallVector<dim>& waveVector() const { return k_; }

   private:
    std::complex<double> phaseFactor(
        const Geometry::PointPhysical<dim>& p) const {
        using namespace std::complex_literals;
        return std::exp(1i * (k_ * p.getCoordinates() + phase_));
    }

    LinearAlgebra::SmallVector<dim> k_;
    LinearAlgebra::SmallVectorC<dim> E0_;
    double phase_;
};

}  // namespace DGMax

#endif  // HPGEM_PLANEWAVE_H

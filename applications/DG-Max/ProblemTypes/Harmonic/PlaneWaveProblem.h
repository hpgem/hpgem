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

#include "SampleHarmonicProblems.h"

#include <Utils/PlaneWave.h>

#ifndef HPGEM_PLANEWAVEPROBLEM_H
#define HPGEM_PLANEWAVEPROBLEM_H

namespace DGMax {

/// Harmonic test case for the plane wave
///     E0 exp[i (k.x + phi)]
/// where
///   - E0 is the direction of the field (may be complex valued)
///   - k is the wave vector
///   - phi is a phase
/// The wave is assumed to propagate in a medium of constant material
/// coefficients.
///
template <std::size_t dim>
class [[maybe_unused]] PlaneWaveProblem : public SampleHarmonicProblem<dim> {
   public:
    PlaneWaveProblem(LinearAlgebra::SmallVector<dim> k,
                     LinearAlgebra::SmallVectorC<dim> E0, double omega,
                     double phase, Material material)
        : planeWave_(k, E0, phase), material_(material), omega_(omega) {}

    static PlaneWaveProblem onDispersionPlaneWave(
        LinearAlgebra::SmallVector<dim> k, LinearAlgebra::SmallVectorC<dim> E0,
        Material material, double phase = 0.0) {
        return PlaneWaveProblem<dim>{
            k, E0, std::sqrt(k.l2NormSquared()) / material.getRefractiveIndex(),
            phase, material};
    }

    double omega() const override { return omega_; }

    LinearAlgebra::SmallVectorC<dim> exactSolution(
        const Geometry::PointPhysical<dim>& point) const override {
        return planeWave_.field(point);
    }
    LinearAlgebra::SmallVectorC<dim> exactSolutionCurl(
        const Geometry::PointPhysical<dim>& point) const override {
        return planeWave_.fieldCurl(point);
    }
    LinearAlgebra::SmallVectorC<dim> sourceTerm(
        const Geometry::PointPhysical<dim>& point) const override {
        return (planeWave_.waveVector().l2NormSquared() /
                    material_.getPermeability() -
                omega_ * omega_ * material_.getPermittivity()) *
               planeWave_.field(point);
    }

    LinearAlgebra::SmallVector<dim> localFlux(
        const Geometry::PointPhysical<dim>& point) {
        // Poynting vector =
        //    Re(E x H*)
        //  = Re(i/omega E x 1/mu Curl E*)
        //  = 1/omega Im(E x 1/mu Curl E*)
        auto field = exactSolution(point);
        auto curlFieldConj = exactSolutionCurl(point).conj();
        return 1 / (omega_ * material_.getPermeability()) *
               LinearAlgebra::leftDoubledCrossProduct(field, curlFieldConj)
                   .imag();
    }

   private:
    DGMax::PlaneWave<dim> planeWave_;
    Material material_;
    double omega_;
};

}  // namespace DGMax

#endif  // HPGEM_PLANEWAVEPROBLEM_H

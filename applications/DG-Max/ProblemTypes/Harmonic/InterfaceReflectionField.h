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
#ifndef HPGEM_INTERFACEREFLECTIONFIELD_H
#define HPGEM_INTERFACEREFLECTIONFIELD_H

#include "../FieldPattern.h"
#include "../../Utils/PlaneWave.h"
#include "../../DGMaxLogger.h"

namespace DGMax {

/// Reflection of a plane wave from an interface.
///
/// This field pattern assumes that an incoming wave
///   E(x,y) = <0,1> exp[i (k.x + phi)]
/// is reflected by an interface at x = x0, where the material changes
/// changes from mat1 (x < x0) to mat2 (x > x0). The wave is assumed to
/// be on the dispersion curve.
///
/// Natural boundary conditions would be:
///  - x-normal: Silver Muller with incoming plane wave imposed
///  - y-normal: PEC (as field is in y direction)
///  - z-normal: PMC (field is parallel to the surface)
template <std::size_t dim>
class InterfaceReflectionField : public FieldPattern<dim> {
   public:
    using typename FieldPattern<dim>::VecC;
    using typename FieldPattern<dim>::PPhys;

    InterfaceReflectionField(double omega, double phase, Material mat1,
                             Material mat2, double position)
        : interfacePosition_(position) {
        using namespace std::complex_literals;
        using namespace hpgem;

        // Wavevectors in x-direction with the right magnitude
        using VecR = LinearAlgebra::SmallVector<dim>;
        VecR kin, kref, ktrans, E0;
        kin[0] = omega * mat1.getRefractiveIndex();
        kref[0] = -kin[0];
        ktrans[0] = omega * mat2.getRefractiveIndex();
        E0[1] = 1.0;

        // The reflection and transmission coefficients are derived in several
        // standard textbooks.
        // beta = sqrt(epsilon_2 mu_1/epsilon_1 mu_2)
        // r=(1-beta)/(1+beta), t=2/(1+beta)
        double beta = mat2.getImpedance() / mat1.getImpedance();
        double r = (1 - beta) / (1 + beta);
        double t = 2 / (1 + beta);
        DGMaxLogger(INFO, "Interface reflection with r=%, t=%", r, t);

        // At the interface the phases of the fields should match. The
        // reflection and transmission coefficients would ensure this when using
        // three plane waves and an interface at x=0:
        // Ei = <0,1> exp(ik_i x + phi)
        // Er = <0,r> exp(ik_r x + phi)
        // Et = <0,t> exp(ik_t x + phi)
        // However, we have in interface at x = x_L, as such we need to
        // compensate the phase of each wave.
        auto rPhase = phase + kin[0] * position - kref[0] * position;
        auto tPhase = phase + kin[0] * position - ktrans[0] * position;
        incident_ = PlaneWave<dim>(kin, E0, phase);
        reflected_ = PlaneWave<dim>(kref, r * E0, rPhase);
        transmitted_ = PlaneWave<dim>(ktrans, t * E0, tPhase);
    };

    VecC field(const PPhys& p) const override {
        if (p[0] < interfacePosition_) {
            return incident_.field(p) + reflected_.field(p);
        } else {
            return transmitted_.field(p);
        }
    }
    VecC fieldCurl(const PPhys& p) const override {
        if (p[0] < interfacePosition_) {
            return incident_.fieldCurl(p) + reflected_.fieldCurl(p);
        } else {
            return transmitted_.fieldCurl(p);
        }
    }
    VecC fieldDoubleCurl(const PPhys& p,
                         const MaterialTensor& material) const override {
        if (p[0] < interfacePosition_) {
            return incident_.fieldDoubleCurl(p, material) +
                   reflected_.fieldDoubleCurl(p, material);
        } else {
            return transmitted_.fieldDoubleCurl(p, material);
        }
    }

   private:
    double interfacePosition_;
    PlaneWave<dim> incident_;
    PlaneWave<dim> reflected_;
    PlaneWave<dim> transmitted_;
};

}  // namespace DGMax

#endif  // HPGEM_INTERFACEREFLECTIONFIELD_H

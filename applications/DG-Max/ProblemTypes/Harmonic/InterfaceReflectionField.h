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
/// This field pattern assumes that a plane wave is normal incident on an
/// interface. This is then partially reflected and partially transmitted.
/// The interface is specified by the condition n . x = c, for some normal
/// vector n and scalar c, with n . x < c being the first material with the
/// incident wave and n . x > c the material with the transmitted wave.
///
/// For a logical interpretation one should have n . k > 0.
template <std::size_t dim>
class InterfaceReflectionField : public FieldPattern<dim> {
   public:
    using typename FieldPattern<dim>::VecC;
    using typename FieldPattern<dim>::PPhys;

    /// \param incidentWave The incident plane wave
    /// \param mat1 Material 1 (incident, reflected wave)
    /// \param mat2 Material 2 (transmitted wave)
    /// \param interfaceNormal Normal direction of the interface
    /// \param position Position c of the interface
    InterfaceReflectionField(PlaneWave<dim> incidentWave, Material mat1,
                             Material mat2,
                             LinearAlgebra::SmallVector<dim> interfaceNormal,
                             double position)
        : incident_(incidentWave),
          interfaceNormal_(interfaceNormal),
          interfacePosition_(position) {
        logger.assert_always(
            std::abs(std::abs(incidentWave.waveVector() * interfaceNormal_) -
                     incidentWave.waveVector().l2Norm()) <
                1e-8 * incidentWave.waveVector().l2Norm(),
            "Wave vector and interface normal are not parallel");

        using namespace std::complex_literals;
        using namespace hpgem;

        // Wavevectors in x-direction with the right magnitude
        using VecR = LinearAlgebra::SmallVector<dim>;
        VecR kin, kref, ktrans;
        kin = incidentWave.waveVector();
        kref = kin - 2 * (interfaceNormal_ * kin) *
                         interfaceNormal_;  // non-normal incident
        // Only for normal incidence
        ktrans = kin * mat2.getRefractiveIndex() / mat1.getRefractiveIndex();

        VecC E0 = incidentWave.field({});  // Field at the origin

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
        // three plane waves and an interface at x=0 (normal = e_x):
        // Ei = <0,1> exp(ik_i x + phi)
        // Er = <0,r> exp(ik_r x + phi)
        // Et = <0,t> exp(ik_t x + phi)
        // However, we have in interface away from the origin, so we have to
        // compensate with a phase factor.
        // Choose a point on the interface
        auto pos = interfaceNormal_ * position;
        auto inPhase = pos * kin;
        auto rPhase = inPhase - kref * pos;
        auto tPhase = inPhase - ktrans * pos;
        reflected_ = PlaneWave<dim>(kref, r * E0, rPhase);
        transmitted_ = PlaneWave<dim>(ktrans, t * E0, tPhase);
    }

    /// Lecagy constructor. Assumes normal in <0,1,0> direction and polarization
    /// <1, 0, 0>.
    /// \param omega Frequency of the wave
    /// \param phase Phase offset
    /// \param mat1 Material 1 (incident, reflected)
    /// \param mat2 Material 2 (transmitted)
    /// \param position position of interface (y coordinate).
    InterfaceReflectionField(double omega, double phase, Material mat1,
                             Material mat2, double position);

    VecC field(const PPhys& p) const override {
        if (interfaceNormal_ * p.getCoordinates() < interfacePosition_) {
            return incident_.field(p) + reflected_.field(p);
        } else {
            return transmitted_.field(p);
        }
    }
    VecC fieldCurl(const PPhys& p) const override {
        if (interfaceNormal_ * p.getCoordinates() < interfacePosition_) {
            return incident_.fieldCurl(p) + reflected_.fieldCurl(p);
        } else {
            return transmitted_.fieldCurl(p);
        }
    }
    VecC fieldDoubleCurl(const PPhys& p,
                         const MaterialTensor& material) const override {
        if (interfaceNormal_ * p.getCoordinates() < interfacePosition_) {
            return incident_.fieldDoubleCurl(p, material) +
                   reflected_.fieldDoubleCurl(p, material);
        } else {
            return transmitted_.fieldDoubleCurl(p, material);
        }
    }

   private:
    LinearAlgebra::SmallVector<dim> interfaceNormal_;
    double interfacePosition_;
    PlaneWave<dim> incident_;
    PlaneWave<dim> reflected_;
    PlaneWave<dim> transmitted_;
};

}  // namespace DGMax

#endif  // HPGEM_INTERFACEREFLECTIONFIELD_H

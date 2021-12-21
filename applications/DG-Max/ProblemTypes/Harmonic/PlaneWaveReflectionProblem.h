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

#ifndef HPGEM_PLANEWAVEREFLECTIONPROBLEM_H
#define HPGEM_PLANEWAVEREFLECTIONPROBLEM_H

#include "SampleHarmonicProblems.h"

namespace DGMax {
/// Reflection of a plane wave from an interface.
///
/// This problem assumes that an incoming wave
///   E(x,y) = <0,1> exp[i (k.x + phi)]
/// is reflected by an interface at x = x0, where the dielectric constant
/// changes from epsilon1 (x < x0) to epsilon2 (x > x0). The wave is assumed to
/// be on the dispersion curve omega = sqrt(epsilon1) k1 = sqrt(epsilon2) k2.
///
/// Boundary conditions are fixed:
///  - x-normal: Silver Muller with incomming plane wave imposed
///  - y-normal: PEC (as field is in y direction)
///  - z-normal: PMC (field is parallel to the surface)
/// \tparam dim The dimension of the problem, either 2 or 3
template <std::size_t dim>
class [[maybe_unused]] PlaneWaveReflectionProblem
    : public ExactHarmonicProblem<dim> {
   public:
    PlaneWaveReflectionProblem(double omega, double phase, Material material1,
                               Material material2, double interfacePosition)
        : omega_(omega),
          k1_(omega_ * material1.getRefractiveIndex()),
          k2_(omega_ * material2.getRefractiveIndex()),
          interfacePosition_(interfacePosition) {

        using namespace std::complex_literals;
        // Incident wave is expressed as
        // EI(x,y,z) = <0, phasor, 0> exp[i k1 x]
        incidentPhasor_ = std::exp(1i * phase);
        // The phasor corresponding to propagation from x=0 to
        // x=interfacePosition. The incident field at the interface is thus
        // incidentPhasor * transitionPhasor
        auto transitionPhasor = std::exp(1i * (k1_ * interfacePosition_));
        auto phasorAtInterface = incidentPhasor_ * transitionPhasor;
        // The reflection and transmission coefficients are derived in several
        // EM-textbooks.
        // beta = sqrt(epsilon_2 mu_1/epsilon_1 mu_2)
        // r=(1-beta)/(1+beta), t=2/(1+beta)
        double beta = material2.getImpedance() / material1.getImpedance();
        // The reflected wave is represented as
        // ER(x,y,z) = <0, phasor, 0> exp[-i k1 x]
        reflectionPhasor_ = phasorAtInterface *
                            // Reflection coefficient
                            (1 - beta) / (1 + beta) *
                            // The value so far is at x=interface, compensate to
                            // make it relative to x=0.
                            transitionPhasor;
        // The transmissted wave is represented as
        // ET(x,y,z) = <0, phasor, 0> exp[i k2 x]
        transmissionPhasor_ = phasorAtInterface *
                              // Transmission coefficient
                              2.0 / (1 + beta) *
                              // The value so far is at x=interface, compensate
                              // to make it relative to x=0.
                              std::exp(-1i * k2_ * interfacePosition_);
    }

    double omega() const override { return omega_; }

    LinearAlgebra::SmallVectorC<dim> exactSolution(
        const Geometry::PointPhysical<dim>& point) const override;
    LinearAlgebra::SmallVectorC<dim> exactSolutionCurl(
        const Geometry::PointPhysical<dim>& point) const override;
    LinearAlgebra::SmallVectorC<dim> sourceTerm(
        const Base::Element&,
        const Geometry::PointPhysical<dim>&) const override {
        // No source term, on dispersion curve
        return {};
    }

    LinearAlgebra::SmallVectorC<dim> boundaryCondition(Base::PhysicalFace<dim> &
                                                       face) const override;

    BoundaryConditionType getBoundaryConditionType(const Base::Face& face)
        const override;

   private:
    double omega_;
    double k1_;
    double k2_;
    double interfacePosition_;
    /// exp[i phi]
    std::complex<double> incidentPhasor_;
    /// Phasor for the reflection, so that it is described as
    /// <0,phasor> exp[-i k.x]
    std::complex<double> reflectionPhasor_;
    /// Phasor for the transmission, so that it is described as
    /// <0,phasor> exp[i kt.x] with kt = omega * sqrt[epsilon2]
    std::complex<double> transmissionPhasor_;
};

}  // namespace DGMax

#endif  // HPGEM_PLANEWAVEREFLECTIONPROBLEM_H

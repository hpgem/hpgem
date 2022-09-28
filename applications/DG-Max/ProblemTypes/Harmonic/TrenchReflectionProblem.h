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

#ifndef HPGEM_TRENCHREFLECTIONPROBLEM_H
#define HPGEM_TRENCHREFLECTIONPROBLEM_H

#include "SampleHarmonicProblems.h"

#include <Utils/PlaneWave.h>

namespace DGMax {

/// Problem describing plane waves propagating in a trench 0 < x < L.
///
/// The goal of this problem is to study the effect/correctness of transparent
/// boundary conditions, which have some reflection. For this we use a
/// rectangular trench/waveguide 0 < x < Lx, and truncate it to 0 < y < Ly. As
/// boundary conditions we use
///  - x = 0, x = Lx: PEC, mirroring the solution
///  - y = 0, y = Ly: Transparent boundary condition (of some form).
///
/// As incident field we impose two plane waves:
///  Ein(x,y) = [1, -kx/ky] exp(i[kx x + ky y] + [1, kx/ky] exp(i[ky y - kx x])
/// where the transverse wavevector kx = n pi / L (n = 1,2,...) is chosen to
/// match the boundary conditions at x=0, Lx. The wave vector ky follows from
/// the dispersion relation kx^2 + ky^2 = omega^2 mu epsilon.
///
/// When the boundary condition at y = Ly is perfectly transparent, the wave
/// will propagate through. However, most boundary conditions have some spurious
/// reflection with amplitude reflection coefficient r. As result the wave
/// reflects at the top boundary into
///  Er(x,y) = [1, kx/ky] exp(i[kx x - ky y] + [1, -kx/ky] exp(i[-ky y - kx x])
/// For the initial reflection of the incident wave this results in the
/// condition -r Ein(x,Ly)|_x = Er(x, Ly)|_x. This wave propagates to the
/// boundary at y=0, where it reflects into the same mode as Ein. Thus the
/// truncated trench forms a leaky Fabry Perot cavity.
///
/// For the 3D problem the solution is taken as constant in z-direction. The
/// corresponding boundary has PMC/Neumann boundary conditions
///
/// \tparam dim The dimension of the problem (2, 3)
// Implementation note: This is not using ExactFieldHarmonicProblem as base
// class. That class assumes that the 'source terms' at the boundary are to
// match the exact solution. Here we do the opposite, we set an incident field
// and the homogeneous+incident field boundary condition. Then we derive what
// the analytical reflection should be from the boundary.
template <std::size_t dim>
class TrenchReflectionProblem : public SampleHarmonicProblem<dim> {
   public:
    TrenchReflectionProblem(double omega, double width, std::size_t waveCount,
                            double length, Material material);

    double omega() const final { return omega_; }
    LinearAlgebra::SmallVectorC<dim> sourceTerm(
        const Base::Element&,
        const Geometry::PointPhysical<dim>& point) const final {
        // On dispersion by design
        return {};
    }

    /**
     * The boundary condition to use at y=Ly, i.e. the far end of the trench.
     * @param bct
     */
    void setFarEndBoundaryCondition(BoundaryConditionType bct) {
        farEndBCT_ = bct;
    }

    /**
     * Set that a region y > y_PML is ocupied by the rectilinear PML.
     *
     * This assumes the PML uses sigma = i/omega * a * (y - y_PML)^2/depth^3 as
     * scaling/conductivity.
     * @param ystart The start of the PML
     * @param scaling The scaling constant
     */
    void setPML(double ystart, double scaling, double depth) {
        pmlYstart_ = ystart;
        pmlScaling_ = scaling / (depth * depth * depth);
    }

    double farSideL2FieldIntegral() const;

    BoundaryConditionType getBoundaryConditionType(
        const Base::Face& face) const final;
    LinearAlgebra::SmallVectorC<dim> exactSolution(
        const Geometry::PointPhysical<dim>& point) const final;
    LinearAlgebra::SmallVectorC<dim> exactSolutionCurl(
        const Geometry::PointPhysical<dim>& point) const final;
    LinearAlgebra::SmallVectorC<dim> boundaryCondition(
        Base::PhysicalFace<dim>& face) const final;

   private:
    /**
     * Combine a quantity from the forward and (reflected) backward wave to
     * account for the Fabry-Perot effect of the boundaries and the dampening of
     * the PML.
     *
     * The forward quantity should correspond to the value from the incident
     * wave. The backward quantity should correspond the to value from the
     * backward wave, with zero phase at y=0.
     *
     * This is then combined to include the reflection at both the y=Ly and y=0,
     * boundaries and the attenuation of the PML.
     * @tparam T The value type
     * @param forward Quantity from the forward wave
     * @param backward Quantity from the backward wave
     * @param y
     * @return The combined value, including all reflections
     */
    template <typename T>
    T combineWaves(T forward, T backward, double y) const;

    double omega_;
    /**
     * Width of the trench (Lx)
     */
    double width_;
    /**
     * Length of the trench (Ly)
     */
    double length_;
    /**
     * Number of half waves in the x direction
     */
    std::size_t waveCount_;
    /**
     * Material in the channel
     */
    Material material_;
    /**
     * Boundary condition used at the y=Ly end.
     */
    DGMax::BoundaryConditionType farEndBCT_;

    double pmlYstart_;
    double pmlScaling_;

    /**
     * Solution functions
     */
    PlaneWave<dim> waveIn1_, waveIn2_, waveR1_, waveR2_;
};

}  // namespace DGMax

#endif  // HPGEM_TRENCHREFLECTIONPROBLEM_H

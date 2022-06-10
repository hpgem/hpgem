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

#include "TrenchReflectionProblem.h"

namespace DGMax {

template <std::size_t dim>
TrenchReflectionProblem<dim>::TrenchReflectionProblem(double omega,
                                                      double width,
                                                      std::size_t waveCount,
                                                      double length,
                                                      Material material)
    : omega_(omega),
      width_(width),
      waveCount_(waveCount),
      length_(length),
      material_(material),
      farEndBCT_(BoundaryConditionType::SILVER_MULLER),
      // Disable PML by setting scaling to 0
      pmlYstart_(0.0),
      pmlScaling_(0.0) {
    logger.assert_always(waveCount > 0, "Non zero wave count");
    double kx = waveCount_ * M_PI / width_;
    double ky2 = omega_ * omega_ * material.getPermittivity() *
                     material.getPermeability() -
                 kx * kx;
    logger.assert_always(ky2 > 0,
                         "Frequency too low to support the transverse wave");
    double ky = std::sqrt(ky2);
    waveIn1_ = PlaneWave<dim>({kx, ky}, {1, -kx / ky}, 0.0);
    waveIn2_ = PlaneWave<dim>({-kx, ky}, {1, kx / ky}, 0.0);
    waveR1_ = PlaneWave<dim>({kx, -ky}, {1, kx / ky}, 0.0);
    waveR2_ = PlaneWave<dim>({-kx, -ky}, {1, -kx / ky}, 0.0);
}

template <std::size_t dim>
BoundaryConditionType TrenchReflectionProblem<dim>::getBoundaryConditionType(
    const Base::Face& face) const {

    auto normal = face.getNormalVector(
        face.getReferenceGeometry()->getCenter().castDimension<1>());
    normal /= normal.l2Norm();
    if (std::abs(normal[0]) > 1e-8) {
        return BoundaryConditionType::DIRICHLET;
    } else if (normal[1] > 1e-1) {
        return farEndBCT_;
    } else if (dim > 2 && std::abs(normal[2]) > 1e-8) {
        return BoundaryConditionType::NEUMANN;
    } else {
        return BoundaryConditionType::SILVER_MULLER;
    }
}
template <std::size_t dim>
LinearAlgebra::SmallVectorC<dim> TrenchReflectionProblem<dim>::exactSolution(
    const Geometry::PointPhysical<dim>& point) const {
    auto forward = waveIn1_.field(point) + waveIn2_.field(point);
    auto backWard = waveR1_.field(point) + waveR2_.field(point);

    return combineWaves(forward, backWard, point[1]);
}

template <std::size_t dim>
LinearAlgebra::SmallVectorC<dim>
    TrenchReflectionProblem<dim>::exactSolutionCurl(
        const Geometry::PointPhysical<dim>& point) const {
    auto forward = waveIn1_.fieldCurl(point) + waveIn2_.fieldCurl(point);
    auto backWard = waveR1_.fieldCurl(point) + waveR2_.fieldCurl(point);

    return combineWaves(forward, backWard, point[1]);
}

template <std::size_t dim>
LinearAlgebra::SmallVectorC<dim>
    TrenchReflectionProblem<dim>::boundaryCondition(
        Base::PhysicalFace<dim>& face) const {
    BoundaryConditionType bct = getBoundaryConditionType(*face.getFace());

    if (bct == BoundaryConditionType::DIRICHLET ||
        bct == BoundaryConditionType::NEUMANN) {
        // By design, we have PEC boundaries
        return {};
    } else {
        using namespace std::complex_literals;
        const auto& p = face.getPointPhysical();
        const auto& normal = face.getUnitNormalVector();

        if (normal[1] > 0) {
            // Top surface has no incident field
            return {};
        }

        LinearAlgebra::SmallVectorC<dim> result;
        result =
            (waveIn1_.fieldCurl(p) + waveIn2_.fieldCurl(p)) /
                material_.getPermeability() +
            1i * omega_ * material_.getImpedance() *
                (waveIn1_.field(p) + waveIn2_.field(p)).crossProduct(normal);

        return result;
    }
}

template <std::size_t dim>
template <typename T>
T TrenchReflectionProblem<dim>::combineWaves(T forward, T backward,
                                             double y) const {
    using namespace std::complex_literals;

    double k2 = waveIn1_.waveVector().l2NormSquared();
    // Choosing wave in as ky > 0
    double ky = waveIn1_.waveVector()[1];

    // The effect of the PML is to attenuate the incoming wave and its
    // reflection. The incident and reflected wave are assumed to not include
    // the effects of the PML. Hence, we need to include several effects:
    //  - For a position inside the PML we need to attenuate the incident field
    //  - For a position inside the PML we need to exponentially increase the
    //    reflected wave.
    //  - The reflection coefficient needs to include the total attenuation of
    //    both the incident wave going into the PML and from the reflection
    //    traveling out of the PML.

    /// Attenuation factor from the PML
    double pmlFactor;
    {
        // PML rescaling from attenuation
        double pmly = (y - pmlYstart_);
        if (pmly > 0) {
            // Scales as exp(-1/omega int_(pmly, y) pmlscaling (y-pmly)^2)
            double exponent =
                -pmlScaling_ * ky / (3 * omega_) * pmly * pmly * pmly;
            double scaling = std::exp(exponent);
            forward *= scaling;
            // The backwards wave is from the reflection, so that grows
            // exponentially towards the far side.
            backward /= scaling;

        }
        double pmlDepth = (length_ - pmlYstart_);
        double exponent =
            -pmlScaling_ * ky / (3 * omega_) * pmlDepth * pmlDepth * pmlDepth;
        // 2 from both traveling into, and then out of it
        exponent *= 2;
        pmlFactor = std::exp(exponent);
    }

    // Compute reflection coefficient using the wavevector. This is equivalent
    // to -tan^2(theta/2), where ky = |k| cos(theta)
    // i.e. theta is the angle of the wave vector with respect to the normal.
    double r = 1.0 - 2 * k2 / (k2 + ky * std::sqrt(k2));
    // Reflection coefficient of the y=Ly side
    double rL;
    switch (farEndBCT_) {
        case BoundaryConditionType::DIRICHLET:
            rL = 1.0 * pmlFactor;
            break;
        case BoundaryConditionType::SILVER_MULLER:
            rL = r * pmlFactor;
            break;
        case BoundaryConditionType::NEUMANN:
            rL = -1.0 * pmlFactor;
            break;
        default:
            logger.fail(
                "Reflection coefficient not implemented for this boundary "
                "condition type");
    }

    // Phase gained from traversing to the back boundary and back.
    auto phase = std::exp(2.0i * (length_ * ky));

    return (forward - rL * phase * backward) / (1.0 - rL * r * phase);
}

template class TrenchReflectionProblem<2>;
template class TrenchReflectionProblem<3>;

}  // namespace DGMax

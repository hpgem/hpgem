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

#ifndef HPGEM_PMLELEMENTINFOS_H
#define HPGEM_PMLELEMENTINFOS_H

#include "ElementInfos.h"

#include "LinearAlgebra/SmallVector.h"
#include "Utils/Dispersive.h"

/**
 * PML ElementInfos, currently only supporting
 *  - Axis aligned PML regions
 *  - Quadratic conductance profile (cubic decay)
 *
 * The implementation is a standard PML. It can be derived from a coordinate
 * transformation into the complex plane. Here we implement the version from the
 * book by Monk, and only the rectilinear version. For each coordinate direction
 * xi with PML in the region xi > ai it stretches the coordinate using
 *  - xi + i/omega int_ai^xi sigma(s) ds   (xi > ai)
 *  - xi                                   (otherwise)
 * with sigma(s) = sigma0/di^3 (s-ai), sigma0 a scaling constant and di the
 * thickness.
 *
 * This results in an exponential decay of propagating waves in the xi direction
 * by a factor ki/omega int_ai^xi sigma(s) ds, where ki is the wave vector for
 * the xi direction. The scaling by di^3 is such that integral evaluates to
 * sigma0/3 (assuming integration from ai to ai+di)
 *
 * Note: The field in the PML region is not the attenuated electric field.
 * Instead, it is the attenuated electric field multiplied by a diagonal tensor
 * 1 + i/omega sigma(xi)
 *
 * @tparam dim The dimension of the problem
 */
template <std::size_t dim>
class PMLElementInfos : public ElementInfos, public DGMax::Dispersive {
    using VecC3 = hpgem::LinearAlgebra::SmallVectorC<3>;

   public:
    using VecR = hpgem::LinearAlgebra::SmallVector<dim>;

    PMLElementInfos() = default;
    /**
     * Constructor
     * @param material Background material (epsilon, sigma)
     * @param offset The edge of the PML for each coordinate.
     * @param directions The direction 1,0,-1 of the PML for each of the
     * coordinates. Use 0 for no PML boundary in that direction, 1 for PML for
     * x > offset, and -1 for x < offset
     * @param depths Depth of the layer in each direction, used to factor this
     * out.
     * @param scaling The scaling constant.
     */
    PMLElementInfos(const DGMax::Material& material, const VecR& offset,
                    const VecR& directions, const VecR& depths,
                    const VecR& scaling)
        : ElementInfos(material),
          offset_(offset),
          directions_(directions),
          scaling_(scaling) {
        // Factor out the scaling based on the material depth.
        for (std::size_t i = 0; i < dim; ++i) {
            scaling_[i] /= depths[i] * depths[i] * depths[i];
        }
    }

    PMLElementInfos(const DGMax::Material& material, const VecR& offset,
                    const VecR& directions, const VecR& depths, double scaling)
        : PMLElementInfos(material, offset, directions, depths,
                          VecR::constant(scaling)) {}

    DGMax::MaterialTensor getMaterialConstantDiv(
        const PointPhysicalBase& p) const override {
        auto diagTensor = diagonalTensor(p, getDispersionWavenumber()) *
                          ElementInfos::getMaterial().getPermittivity();
        return DGMax::MaterialTensor(diagTensor);
    }
    DGMax::MaterialTensor getMaterialConstantCurl(
        const PointPhysicalBase& p) const override {
        auto diagTensor = diagonalTensor(p, getDispersionWavenumber()) *
                          ElementInfos::getMaterial().getPermeability();
        return DGMax::MaterialTensor(diagTensor);
    }

    DGMax::MaterialTensor getFieldRescaling(
        const PointPhysicalBase& p) const override {
        VecC3 stretching =
            getCoordinateStretching(p, getDispersionWavenumber());
        for (std::size_t i = 0; i < 3; ++i) {
            stretching[i] = 1.0 / stretching[i];
        }
        return DGMax::MaterialTensor(stretching);
    }

    DGMax::MaterialTensor getCurlFieldRescaling(
        const PointPhysicalBase& p) const override {
        // The curl needs to be rescaled by the diagonal tensor
        // 1/d2d3, 1/d1d3, 1/d1d2.
        // applyCurl in MaterialConstant applies the inverse of the values so
        // compute d2d3 d1d3 d1d2.
        VecC3 stretching = getCoordinateStretching(p, getDispersionWavenumber());
        VecC3 result = VecC3 ::constant(1.0);
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                auto val = stretching[j];
                if (i != j) {
                    result[i] *= val;
                }
            }
        }

        return DGMax::MaterialTensor(result);
    }

    bool isDispersive() const final { return true; }

    /// \brief Compute the required scaling parameter for an attenuation.
    ///
    /// Computes the scaling parameter to get the desired attenuation. This
    /// computation is done independently for each PML direction
    /// (direction[j] == +- 1). This allows computing corner regions.
    ///
    /// The computation for a direction i assumes a plane wave propagating along
    /// the i-axis (direction depends on direction[i]). After this plane wave
    /// enters the PML it will travel to the interface at the far end, being
    /// attenuated on the way. It is reflected of the far interface and then
    /// attenuated again on its way back to the front interface (where it exits
    /// the PML). The scaling parameter is such that the attenuation of the PML
    /// (excluding the reflection at the far end) is the specified value.
    ///
    /// For inactive directions (i.e. direction[i] == 0) the scaling value will
    /// be 0.
    ///
    /// \param material The base material.
    /// \param direction The direction of PML
    /// \param attenuation The required attenuation in each direction
    /// \return The required scaling for each direction.
    static VecR computeScaling(const DGMax::Material& material, VecR direction,
                               VecR attenuation);

   private:
    VecC3 getCoordinateStretching(
        const hpgem::Geometry::PointPhysical<dim>& point, double omega) const {
        using namespace std::complex_literals;

        VecC3 stretching;
        stretching.set(1.0);
        for (std::size_t i = 0; i < dim; ++i) {
            double xi = (point[i] - offset_[i]) * directions_[i];
            if (xi > 0) {
                stretching[i] += (1i / omega) * scaling_[i] * xi * xi;
            }
        }
        return stretching;
    }

    VecC3 diagonalTensor(const hpgem::Geometry::PointPhysical<dim>& point,
                         double omega) const {
        VecC3 ds = getCoordinateStretching(point, omega);
        // Mix the di's
        VecC3 result;
        result.set(1.0);
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                auto val = ds[j];
                if (i == j) {
                    result[i] /= val;
                } else {
                    result[i] *= val;
                }
            }
        }
        return result;
    }

    VecR offset_;
    VecR directions_;
    VecR scaling_;
};

template <std::size_t dim>
hpgem::LinearAlgebra::SmallVector<dim> PMLElementInfos<dim>::computeScaling(
    const DGMax::Material& material, VecR direction, VecR attenuation) {
    VecR scaling;
    for (std::size_t d = 0; d < dim; ++d) {
        if (static_cast<int>(direction[d]) == 0) {
            scaling[d] = 0.0;
            continue;
        }
        // One way attenuation is given by:
        // exp(-scaling_i/3 n)
        // with
        //  - scaling_i the scaling factor, includes scaling due to thickness of
        //    the medium.
        //  - n the refractive index (needed for the ratio k_i / omega)
        // Two way attenuation is the square of one way attenuation

        // Divide by 2 for 2way -> 1way
        double att = -std::log(attenuation[d]) / 2.0;
        att /= material.getRefractiveIndex();
        att *= 3;  // From the integration
        scaling[d] = att;
    }
    return scaling;
}

#endif  // HPGEM_PMLELEMENTINFOS_H

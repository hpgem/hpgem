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
 * with sigma(s) = sigma0 (s-ai)^2 and sigma0 a scaling constant.
 *
 * This results in an exponential decay of propagating waves in the xi direction
 * by a factor ki/omega int_ai^xi sigma(s) ds, where ki is the wave vector for
 * the xi direction.
 *
 * Note: The field in the PML region is not the attenuated electric field.
 * Instead, it is the attenuated electric field multiplied by a diagonal tensor
 * 1 + i/omega sigma(xi)
 *
 * @tparam dim The dimension of the problem
 */
template <std::size_t dim>
class PMLElementInfos : public ElementInfos {
   public:
    /**
     * Constructor
     * @param material Background material (epsilon, sigma)
     * @param offset The edge of the PML for each coordinate.
     * @param directions The direction 1,0,-1 of the PML for each of the
     * coordinates. Use 0 for no PML boundary in that direction, 1 for PML for x
     * > offset, and -1 for x < offset
     * @param scaling The scaling constant.
     */
    PMLElementInfos(const DGMax::Material& material,
                    const LinearAlgebra::SmallVector<dim>& offset,
                    const LinearAlgebra::SmallVector<dim>& directions,
                    double scaling)
        : ElementInfos(material),
          offset_(offset),
          directions_(directions),
          scaling_(scaling) {}

    DGMax::MaterialTensor getMaterialConstantDiv(const PointPhysicalBase& p,
                                                 double omega) const override {
        auto diagTensor = diagonalTensor(p, omega) *
                          ElementInfos::getMaterial().getPermittivity();
        return DGMax::MaterialTensor(diagTensor);
    }
    DGMax::MaterialTensor getMaterialConstantCurl(const PointPhysicalBase& p,
                                                  double omega) const override {
        auto diagTensor = diagonalTensor(p, omega) *
                          ElementInfos::getMaterial().getPermeability();
        return DGMax::MaterialTensor(diagTensor);
    }

   private:
    LinearAlgebra::SmallVectorC<3> diagonalTensor(
        const Geometry::PointPhysical<dim>& point, double omega) const {
        using namespace std::complex_literals;

        LinearAlgebra::SmallVectorC<3> ds;
        ds.set(1.0);
        for (std::size_t i = 0; i < dim; ++i) {
            double xi = (point[i] - offset_[i]) * directions_[i];
            if (xi > 0) {
                ds[i] += 1i / omega * scaling_ * xi * xi;
            }
        }
        // Mix the di's
        LinearAlgebra::SmallVectorC<3> result;
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

    LinearAlgebra::SmallVector<dim> offset_;
    LinearAlgebra::SmallVector<dim> directions_;
    double scaling_;
};

#endif  // HPGEM_PMLELEMENTINFOS_H

/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2022, University of Twente
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors
 may be used to endorse or promote products derived from this software without
 specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef HPGEM_MAPPINGTESTHELPERS_H
#define HPGEM_MAPPINGTESTHELPERS_H

#include <Geometry/Mappings/MappingReferenceToReference.h>
#include <Geometry/Mappings/MappingReferenceToPhysical.h>
#include <Geometry/ReferenceGeometry.h>
#include "../../catch.hpp"

namespace hpgem {

template <std::size_t dim, std::size_t codim>
void testJacobian(const Geometry::PointReference<dim>& pref,
                  const Geometry::MappingReferenceToReference<codim>& mapping,
                  double h = 1e-8) {
    using PRef = Geometry::PointReference<dim>;
    auto pmapped = mapping.transform(pref);
    Geometry::Jacobian<dim, dim + codim> jac = mapping.calcJacobian(pref);

    for (std::size_t i = 0; i < dim; ++i) {
        PRef prefh = pref;
        prefh[i] += h;
        auto dmapped = mapping.transform(prefh) - pmapped;
        auto gradient = dmapped.getCoordinates() / h;
        // Factor 1000 is for rounding errors in the comptutation, could be
        // improved by taking into account the condition number of the Jacobian
        // and/or relative size of the columns.
        REQUIRE((jac.getColumn(i) - gradient).l2NormSquared() < 1e3 * h * h);
    }
}

template <std::size_t dim>
void testJacobian(const Geometry::PointReference<dim>& pref,
                  const Geometry::MappingReferenceToPhysical<dim>& mapping,
                  double h = 1e-8) {
    using PRef = Geometry::PointReference<dim>;
    auto pmapped = mapping.transform(pref);
    Geometry::Jacobian<dim, dim> jac = mapping.calcJacobian(pref);

    for (std::size_t i = 0; i < dim; ++i) {
        PRef prefh = pref;
        prefh[i] += h;
        auto dmapped = mapping.transform(prefh) - pmapped;
        auto gradient = dmapped.getCoordinates() / h;
        // Factor 1000 is for rounding errors in the comptutation, could be
        // improved by taking into account the condition number of the Jacobian
        // and/or relative size of the columns.
        REQUIRE((jac.getColumn(i) - gradient).l2NormSquared() < 1e3 * h * h);
    }
}

template <std::size_t dim>
void testLinearCodim2Mappings(const Geometry::ReferenceGeometry& geometry) {
    using namespace Geometry;
    for (std::size_t codimindex = 0;
         codimindex < geometry.getNumberOfCodim2Entities(); ++codimindex) {
        INFO("Codim index " << codimindex);
        const ReferenceGeometry* codimGeom =
            geometry.getCodim2ReferenceGeometry(codimindex);
        const MappingReferenceToReference<2>* codimMapping =
            geometry.getCodim2MappingPtr(codimindex);
        std::vector<std::size_t> nodeIndices =
            geometry.getCodim2EntityLocalIndices(codimindex);

        // Test mapping of each coordinate, including Jacobian
        for (std::size_t i = 0; i < nodeIndices.size(); ++i) {
            INFO("Checking mapping of corner point " << i);
            auto codimPoint = codimGeom->getReferenceNodeCoordinate(i)
                                  .castDimension<dim - 2>();
            auto actualPoint =
                geometry.getReferenceNodeCoordinate(nodeIndices[i])
                    .castDimension<dim>();
            REQUIRE((codimMapping->transform(codimPoint) - actualPoint)
                        .getCoordinates()
                        .l2NormSquared() < 1e-16);
            INFO("Checking Jacobian at corner point " << i);
            testJacobian<dim-2,2>(codimPoint, *codimMapping);
        }
        // Test Jacobian at the centre
        INFO("Checking Jacobian at the centre ");
        testJacobian<dim-2,2>(codimGeom->getCenter().castDimension<dim - 2>(),
                     *codimMapping);
    }
}

}  // namespace hpgem

#endif  // HPGEM_MAPPINGTESTHELPERS_H

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
#ifndef HPGEM_SCATTERINGPROBLEM_H
#define HPGEM_SCATTERINGPROBLEM_H

#include <memory>

#include "../FieldPattern.h"
#include "../HarmonicProblem.h"
#include "../../PMLElementInfos.h"

namespace DGMax {

template <std::size_t dim>
class ScatteringProblem : public HarmonicProblem<dim> {
   public:
    ScatteringProblem(double omega,
                      std::shared_ptr<FieldPattern<dim>> incidentField)
        : omega_(omega), incidentField_(std::move(incidentField)){};

    double omega() const override { return omega_; }

    LinearAlgebra::SmallVectorC<dim> sourceTerm(
        const Base::Element& element,
        const Geometry::PointPhysical<dim>& point) const override {

        const ElementInfos& material = ElementInfos::get(element);
        if (dynamic_cast<const PMLElementInfos<dim>*>(&material) != nullptr) {
            return {};
        }

        auto matCurl = material.getMaterialConstantCurl(point, omega_);
        auto matDiv = material.getMaterialConstantDiv(point, omega_);

        // Note the minus sign compared to the standard
        // Curl^2 E - omega^2 E = 0 of Maxwell's equations.
        return -incidentField_->fieldDoubleCurl(point, matCurl) +
               omega_ * omega_ * matDiv.applyDiv(incidentField_->field(point));
    }

    bool isScatterFieldProblem() const override {
        return true;
    }

    BoundaryConditionType getBoundaryConditionType(
        const Base::Face& face) const override {
        return BoundaryConditionType::DIRICHLET;
    }
    LinearAlgebra::SmallVectorC<dim> boundaryCondition(
        Base::PhysicalFace<dim>& face) const override {
        return {};
    }
    LinearAlgebra::SmallVectorC<dim> incidentField(
        Geometry::PointPhysical<dim>& p) const override {
        return incidentField_->field(p);
    }
    LinearAlgebra::SmallVectorC<dim> incidentFieldCurl(
        Geometry::PointPhysical<dim>& p) const override {
        return incidentField_->fieldCurl(p);
    }

   private:
    double omega_;
    std::shared_ptr<FieldPattern<dim>> incidentField_;
};

}  // namespace DGMax

#endif  // HPGEM_SCATTERINGPROBLEM_H

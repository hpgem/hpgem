/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2021, University of Twente
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
#ifndef HPGEM_STANDARDGAUSSQUADRATURERULE_H
#define HPGEM_STANDARDGAUSSQUADRATURERULE_H

#include "GaussQuadratureRule.h"
#include "Geometry/PointReference.h"
#include "Geometry/ReferenceGeometry.h"

namespace hpgem {
namespace QuadratureRules {

template <std::size_t dimension>
struct QuadraturePoint {
    QuadraturePoint(double weight, Geometry::PointReference<dimension>&& point)
        : weight_(weight), point_(point) {}

    Geometry::PointReference<dimension> point_;
    double weight_;
};

template <std::size_t dim>
class StandardGaussQuadratureRule : public GaussQuadratureRule {

   public:
    StandardGaussQuadratureRule(std::string name, std::size_t order,
                                Geometry::ReferenceGeometry* referenceGeometry,
                                std::vector<QuadraturePoint<dim>>&& points)
        : name_(name),
          order_(order),
          referenceGeometry_(referenceGeometry),
          points_(points) {
        logger.assert_always(referenceGeometry->getDimension() == dim,
                             "Reference geometry is for the wrong dimension");
    }

    std::string getName() const final { return name_; }

    std::size_t order() const final { return order_; }
    std::size_t dimension() const final { return dim; }
    std::size_t getNumberOfPoints() const final { return points_.size(); }

    double weight(std::size_t i) const final { return points_[i].weight_; }
    const Geometry::PointReferenceBase& getPoint(std::size_t i) const final {
        return points_[i].point_;
    }
    Geometry::ReferenceGeometry* forReferenceGeometry() const final {
        return referenceGeometry_;
    }

   private:
    std::string name_;
    std::size_t order_;
    std::vector<QuadraturePoint<dim>> points_;
    Geometry::ReferenceGeometry* referenceGeometry_;
};

}  // namespace QuadratureRules
}  // namespace hpgem

#endif  // HPGEM_STANDARDGAUSSQUADRATURERULE_H
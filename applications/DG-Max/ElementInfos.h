/*
This file forms part of hpGEM. This package has been developed over a number of
years by various people at the University of Twente and a full list of
contributors can be found at http://hpgem.org/about-the-code/team

This code is distributed using BSD 3-Clause License. A copy of which can found
below.


Copyright (c) 2014, University of Twente
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

// this file has a container data structure for everything you want to know on a
// per element basis
#ifndef HPGEM_APP_ELEMENTINFOS_H
#define HPGEM_APP_ELEMENTINFOS_H

#include <AbstractDimensionlessBase.h>
#include <Base/UserData.h>

#include <Logger.h>
#include <Base/Element.h>
#include <Geometry/PointPhysicalBase.h>

#include "Material.h"

class ElementInfos : public hpgem::Base::UserData {
   public:
    using Element = hpgem::Base::Element;
    using PointPhysicalBase = hpgem::Geometry::PointPhysicalBase;

    ElementInfos(double epsilon, double permeability = 1.0);
    ElementInfos(DGMax::Material material) : material_(material){};

    const double& getPermittivity() const {
        return material_.getPermittivity();
    }

    const double& getPermeability() const {
        return material_.getPermeability();
    }
    static ElementInfos* get(const Element* element) {
        using hpgem::logger;
        if (element == nullptr) {
            return nullptr;
        } else {
            return &get(*element);
        }
    }

    static ElementInfos& get(const Element& element) {
        using hpgem::logger;
        auto* data = element.getUserData();
        logger.assert_debug(data != nullptr,
                            "No material information available for element %",
                            element.getID());
        auto* elementInfo = dynamic_cast<ElementInfos*>(data);
        logger.assert_debug(elementInfo != nullptr,
                            "Element has incorrect material information");
        return *elementInfo;
    }

    /**
     * Material constant used in second order Maxwells equation.
     *
     * This is the material constant that is used in second order Maxwell's
     * equation at two points. It is the permittivity for the E-field
     * formulation, for H-field it would be the permeability. This constant is
     * used in two places:
     *
     *  - The divergence constrain Div(constant field) = 0 (Div(epsilon E) = 0)
     *  - The omega^2 term: -omega^2 constant field.
     *
     * For PMLs and similar materials it may vary inside the element.
     *
     * @param p The position
     * @return The constant
     */
    double getMaterialConstantDiv(const PointPhysicalBase& p) const {
        return material_.getPermittivity();
    }

    /**
     * Material constant used in second order Maxwells equation.
     *
     * This material constant that is used in second order Maxwells equation
     * between the two Curls. It is the inverse of the permeability for the
     * E-field formulation. For the H-field it would be the inverse of the
     * permittivity.
     *
     * For PMLs and similar materials it may vary inside the element.
     *
     * @param p The position
     * @return The constant
     */
    double getMaterialConstantCurl(const PointPhysicalBase& p) const {
        return 1.0 / material_.getPermeability();
    }

    const DGMax::Material& getMaterial() const { return material_; }

    double getImpedance() const { return material_.getImpedance(); }

    double getRefractiveIndex() const { return material_.getRefractiveIndex(); }

   private:
    DGMax::Material material_;
};
#endif  // HPGEM_APP_ELEMENTINFOS_H
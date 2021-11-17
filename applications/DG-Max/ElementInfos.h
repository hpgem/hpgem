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

#include "Base/UserData.h"
#include "Material.h"
#include "ProblemTypes/ProblemField.h"

class ElementInfos : public hpgem::Base::UserElementData {
   public:
    ElementInfos(double epsilon, double permeability = 1.0);
    ElementInfos(DGMax::Material material) : material_(material){};

    const double& getPermittivity() const {
        return material_.getPermittivity();
    }

    const double& getPermeability() const {
        return material_.getPermeability();
    }

    double getImpedance() const { return material_.getImpedance(); }

    double getRefractiveIndex() const { return material_.getRefractiveIndex(); }

    /// The material constant that is used between the two curls in the second
    /// order form. For E-field 1/permeability, for H-field 1/permittivity
    double getMaterialConstantCurl(DGMax::ProblemField field) const {
        if (field == DGMax::ProblemField::ELECTRIC_FIELD) {
            return 1.0 / getPermeability();
        } else {
            return 1.0 / getPermittivity();
        }
    }

    /// The material constant that is used for the divergence constraint.
    /// For E-field permittivity, for H-field permeability
    double getMaterialConstantDiv(DGMax::ProblemField field) const {
        return field == DGMax::ProblemField::ELECTRIC_FIELD ? getPermittivity()
                                                            : getPermeability();
    }

    /// Impedance for the given field, for E-field the standard impedance, for
    /// H-field the inverse.
    double getFieldImpedance(DGMax::ProblemField field) {
        if (field == DGMax::ProblemField::ELECTRIC_FIELD) {
            return std::sqrt(getPermittivity() / getPermeability());
        } else {
            return std::sqrt(getPermeability() / getPermittivity());
        }
    }

   private:
    DGMax::Material material_;
};
#endif  // HPGEM_APP_ELEMENTINFOS_H
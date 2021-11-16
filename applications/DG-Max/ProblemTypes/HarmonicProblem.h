/*
This file forms part of hpGEM. This package has been developed over a number of
years by various people at the University of Twente and a full list of
contributors can be found at http://hpgem.org/about-the-code/team

This code is distributed using BSD 3-Clause License. A copy of which can found
below.


Copyright (c) 2018, Univesity of Twenete
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

#ifndef HPGEM_APP_HARMONICPROBLEM_H
#define HPGEM_APP_HARMONICPROBLEM_H

#include "Geometry/PointPhysical.h"
#include "Base/PhysicalFace.h"
#include "LinearAlgebra/SmallVector.h"
#include "BoundaryConditionType.h"
#include "ElementInfos.h"

using namespace hpgem;

/// \brief Harmonic maxwell problem.
///
/// A problem description for a harmonic maxwell problem. Such a problem assumes
/// that the current density J(x,t) can be written as e^(i omega t) J(x) and
/// similarly decomposing both E(x,t) and H(x,t) (and no free charges).
///
/// The maxwell equations for such a harmonic problem reduce to finding E such
/// that
/// curl mu^{-1} curl E - omega^2 eps E = i omega J  in the domain
/// n x curl E = g (on the boundary)
/// For some known mu, eps, omega, J and g. Where it is required (but not
/// checked) that div J = 0.
///
/// For more information, see for example section 5.2.1 of Devashish's thesis.
///
template <std::size_t DIM>
class HarmonicProblem {
   public:
    virtual ~HarmonicProblem() = default;
    virtual double omega() const = 0;
    virtual LinearAlgebra::SmallVectorC<DIM> sourceTerm(
        const Geometry::PointPhysical<DIM>& point) const = 0;
    virtual DGMax::BoundaryConditionType getBoundaryConditionType(
        const Base::Face& face) const = 0;

    /**
     * Value of the boundary function. The use of this value depends on
     * getBoundaryConditionType(). For the different boundary conditions this
     * should be:
     *
     * - Dirichlet: n x g_D (TODO: Future refactor -> g_D)
     * - Neumann: g_N
     * - SilverMuller: g_N
     *
     * @param face The face and point to evaluate the BC on
     * @return The value.
     */
    virtual LinearAlgebra::SmallVectorC<DIM> boundaryCondition(
        Base::PhysicalFace<DIM>& face) const = 0;
};

/**
 * Bit field enumeration indicating which properties can be different between
 * two harmonic problems.
 *
 * Note: Not inside HarmonicProblem to not have template parameters
 */
enum class HarmonicProblemChanges : unsigned int {
    OMEGA = 0x1,
    CURRENT_SOURCE = 0x2,
    BOUNDARY_CONDITION_TYPE = 0x4,
    BOUNDARY_CONDITION_VALUE = 0x8,
    ANY = 0xf,
};

template <std::size_t DIM>
class ExactHarmonicProblem : public HarmonicProblem<DIM> {
   public:
    virtual LinearAlgebra::SmallVectorC<DIM> exactSolution(
        const Geometry::PointPhysical<DIM>& point) const = 0;
    virtual LinearAlgebra::SmallVectorC<DIM> exactSolutionCurl(
        const Geometry::PointPhysical<DIM>& point) const = 0;

    LinearAlgebra::SmallVectorC<DIM> boundaryCondition(
        Base::PhysicalFace<DIM>& face) const override {
        using BCT = DGMax::BoundaryConditionType;
        using Vec = LinearAlgebra::SmallVectorC<DIM>;
        BCT bct = this->getBoundaryConditionType(*face.getFace());
        switch (bct) {
            case BCT::DIRICHLET: {
                Vec efield = this->exactSolution(face.getPointPhysical());
                Vec normal = face.getUnitNormalVector();
                return normal.crossProduct(efield);
            }
            case BCT::NEUMANN: {
                auto* material = dynamic_cast<ElementInfos*>(
                    face.getFace()->getPtrElementLeft()->getUserData());
                return this->exactSolutionCurl(face.getPointPhysical()) /
                       material->getPermeability();
            }
            case BCT::SILVER_MULLER: {
                auto* material = dynamic_cast<ElementInfos*>(
                    face.getFace()->getPtrElementLeft()->getUserData());
                logger.assert_debug(material != nullptr, "No material.");
                Vec efield = this->exactSolution(face.getPointPhysical());
                Vec efieldCurl =
                    this->exactSolutionCurl(face.getPointPhysical());
                const Vec& normal = face.getUnitNormalVector();
                auto impedance = std::complex<double>(
                    0, this->omega() * material->impedance());
                // n x (Curl E + Z [E x n]) = n x g_N
                return efieldCurl / material->getPermeability() +
                       impedance * efield.crossProduct(normal);
            }
            default:
                logger(ERROR,
                       "Not implemented for this type of boundary condition.");
                return {};
        }
    }
};

#endif  // HPGEM_APP_HARMONICPROBLEM_H

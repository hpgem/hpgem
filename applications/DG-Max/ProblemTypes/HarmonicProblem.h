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
    virtual void sourceTerm(const Geometry::PointPhysical<DIM>& point,
                            LinearAlgebra::SmallVector<DIM>& result) const = 0;
    virtual void boundaryCondition(
        const Geometry::PointPhysical<DIM>& point,
        Base::PhysicalFace<DIM>& face,
        LinearAlgebra::SmallVector<DIM>& result) const = 0;
};

template <std::size_t DIM>
class ExactHarmonicProblem : public HarmonicProblem<DIM> {
   public:
    virtual void exactSolution(
        const Geometry::PointPhysical<DIM>& point,
        LinearAlgebra::SmallVector<DIM>& result) const = 0;
    virtual void exactSolutionCurl(
        const Geometry::PointPhysical<DIM>& point,
        LinearAlgebra::SmallVector<DIM>& result) const = 0;

    void boundaryCondition(
        const Geometry::PointPhysical<DIM>& point,
        Base::PhysicalFace<DIM>& face,
        LinearAlgebra::SmallVector<DIM>& result) const final {
        LinearAlgebra::SmallVector<DIM> sol;
        exactSolution(point, sol);
        const LinearAlgebra::SmallVector<DIM>& normal =
            face.getUnitNormalVector();
        normal.crossProduct(sol, result);
    }
};

#endif  // HPGEM_APP_HARMONICPROBLEM_H

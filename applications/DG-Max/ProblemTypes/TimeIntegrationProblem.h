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

#ifndef HPGEM_APP_TIMEINTEGRATIONPROBLEM_H
#define HPGEM_APP_TIMEINTEGRATIONPROBLEM_H

#include "Base/PhysicalFace.h"
#include "Geometry/PointPhysical.h"
#include "LinearAlgebra/SmallVector.h"

/// \brief Problem to intergrate a given starting field in time
///
/// This solves the standard time dependent problem
/// eps d^2E/dt^2 + sigma dE/dt + curl (mu^-1 curl E) = - dJ/dt.
/// For known eps, sigma and mu and J.
///
/// The code assumes that mu = 1
// TODO: Do we assume that div J = 0, mu = 1, use eps in the solver?
// Assume sigma = 0?
template <std::size_t DIM>
class TimeIntegrationProblem {
   public:
    virtual ~TimeIntegrationProblem() = default;

    /// \brief Initial value of the electric field E
    virtual void initialCondition(
        const Geometry::PointPhysical<DIM>& point,
        LinearAlgebra::SmallVector<DIM>& result) const = 0;
    /// \brief Initial value of the time derivative of the electric field
    virtual void initialConditionDerivative(
        const Geometry::PointPhysical<DIM>& point,
        LinearAlgebra::SmallVector<DIM>& result) const = 0;
    /// \brief The source term, -dJ/dt, where J is the electric current density.
    ///
    /// \param point The point at which to evaluate the source term.
    /// \param t  The time to evaluate it at
    /// \param result the value of -dJ/dt.
    virtual void sourceTerm(const Geometry::PointPhysical<DIM>& point, double t,
                            LinearAlgebra::SmallVector<DIM>& result) const = 0;
    /// \brief The boundary condition n x E
    ///
    /// \param point The point for which to evaluate the boundary condition.
    /// \param t The time to evaluate the boundary condition at
    /// \param face The face where it is evaluated (e.g. for the normal).
    /// \param result The resulting value of n x E
    virtual void boundaryCondition(
        const Geometry::PointPhysical<DIM>& point, double t,
        Base::PhysicalFace<DIM>& face,
        LinearAlgebra::SmallVector<DIM>& result) const = 0;

    /// \brief The conductivity of the medium.
    ///
    /// The conductivity of the medium, assumed to be a non negative scalar.
    /// See for example "Time-integration methods for finite element
    /// discretisations of the second-order Maxwell equation" by D Sarmany et.
    /// al.
    ///
    /// The default implementation assumes a non conduction medium, thus 0.
    /// Note: According to the older comments it is untested for non zero
    /// values. From the implementation it looks like it is only used with the
    /// CO2 method.
    virtual double conductivity() const { return 0; };
};

/// \brief Special case of a TimeIntgration problem where both source and
/// boundary terms can be separated in a spacial function and time dependent
/// function.
///
/// Separating the source term -dJ/dt (x,t) = Xj(x) Tj(t) and the boundary
/// condition n x E = G(x,t) = Xg(x) Tg(t) allows for a faster implementation
/// of the finite element methods. The idea is that the require inner products
/// with a test function 'v' (e.g. (G, v)) can be rewritten as Tg(t) (Xg(x), v).
/// The latter inner product is time independent and can thus be computed once
/// and then rescaled to obtain the value at different time.
template <std::size_t DIM>
class SeparableTimeIntegrationProblem
    : virtual public TimeIntegrationProblem<DIM> {
   public:
    /// \brief The time independent boundary field (Xg)
    virtual void boundaryConditionRef(
        const Geometry::PointPhysical<DIM>& point,
        Base::PhysicalFace<DIM>& face,
        LinearAlgebra::SmallVector<DIM>& result) const = 0;

    void boundaryCondition(
        const Geometry::PointPhysical<DIM>& point, double t,
        Base::PhysicalFace<DIM>& face,
        LinearAlgebra::SmallVector<DIM>& result) const override {
        boundaryConditionRef(point, face, result);
        result *= timeScalingBoundary(t);
    }

    /// \brief The time independent field of the source (Xj)
    virtual void sourceTermRef(
        const Geometry::PointPhysical<DIM>& point,
        LinearAlgebra::SmallVector<DIM>& result) const = 0;

    void sourceTerm(const Geometry::PointPhysical<DIM>& point, double t,
                    LinearAlgebra::SmallVector<DIM>& result) const override {
        sourceTermRef(point, result);
        result *= timeScalingSource(t);
    }

    // virtual double referenceTimeSource() const = 0;
    /// \brief The time dependent scaling of the boundary value (Tg)
    virtual double timeScalingBoundary(double t) const = 0;
    /// \brief The time depedent scaling of the source term (Tj)
    virtual double timeScalingSource(double t) const = 0;
};

/// \brief Special case of the TimeIntegrationProblem where we have an exact
/// and computable solution.
///
/// Using this solution we can directly implement the boundary and initial
/// condition (only the field, not its time derivative).
template <std::size_t DIM>
class ExactTimeIntegrationProblem : virtual public TimeIntegrationProblem<DIM> {
   public:
    /// \brief Analytical solution to the problem.
    virtual void exactSolution(
        const Geometry::PointPhysical<DIM>& point, double t,
        LinearAlgebra::SmallVector<DIM>& result) const = 0;
    /// \brief Curl of the analytical solution.
    virtual void exactSolutionCurl(
        const Geometry::PointPhysical<DIM>& point, double t,
        LinearAlgebra::SmallVector<DIM>& result) const = 0;

    void initialCondition(
        const Geometry::PointPhysical<DIM>& point,
        LinearAlgebra::SmallVector<DIM>& result) const override {
        exactSolution(point, 0, result);
    }
    void boundaryCondition(
        const Geometry::PointPhysical<DIM>& point, double t,
        Base::PhysicalFace<DIM>& face,
        LinearAlgebra::SmallVector<DIM>& result) const override {
        LinearAlgebra::SmallVector<DIM> eField;
        exactSolution(point, t, eField);
        LinearAlgebra::SmallVector<DIM> normal = face.getUnitNormalVector();
        normal.crossProduct(eField, result);
    }
};

/// \brief An time integration problem that has both separable source and
/// boundary terms and has an exact solution.
template <std::size_t DIM>
class ExactSeparableTimeIntegrationProblem
    : public ExactTimeIntegrationProblem<DIM>,
      public SeparableTimeIntegrationProblem<DIM> {

    void boundaryConditionRef(
        const Geometry::PointPhysical<DIM>& point,
        Base::PhysicalFace<DIM>& face,
        LinearAlgebra::SmallVector<DIM>& result) const final {
        LinearAlgebra::SmallVector<DIM> values;
        this->exactSolution(point, referenceTimeBoundary(), values);
        LinearAlgebra::SmallVector<DIM> normal = face.getUnitNormalVector();
        normal.crossProduct(values, result);
    }

    void boundaryCondition(
        const Geometry::PointPhysical<DIM>& point, double t,
        Base::PhysicalFace<DIM>& face,
        LinearAlgebra::SmallVector<DIM>& result) const final {
        // Note, both the Exact and Separable problems override
        // boundaryCondition, thus we have to tell the compiler which is the one
        // that we want to use for a problem that is both separable and exact.
        SeparableTimeIntegrationProblem<DIM>::boundaryCondition(point, t, face,
                                                                result);
    }

    /// \brief Time used to compute `boundaryConditionRef`.
    ///
    /// Note that timeScalingBoundary(referenceTimeBoundary()) should be 1.
    /// \return The time to sample the exact solution at for the reference
    /// boundary condition.
    virtual double referenceTimeBoundary() const = 0;
};

#endif  // HPGEM_APP_TIMEINTEGRATIONPROBLEM_H

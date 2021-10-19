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

#ifndef HPGEM_APP_SAMPLEHARMONICPROBLEMS_H
#define HPGEM_APP_SAMPLEHARMONICPROBLEMS_H

#include "../HarmonicProblem.h"

namespace DGMax {

using namespace hpgem;

/// Helper class to handle setting boundary conditions
template <std::size_t dim>
class SampleHarmonicProblem : public ExactHarmonicProblem<dim> {
   public:
    using BoundaryConditionIndicator =
        std::function<BoundaryConditionType(const Base::Face& face)>;

    BoundaryConditionType getBoundaryConditionType(
        const Base::Face& face) const override {
        if (boundaryConditionIndicator_) {
            return boundaryConditionIndicator_(face);
        } else {
            // Traditional default
            return BoundaryConditionType::DIRICHLET;
        }
    }

    void setBoundaryConditionIndicator(BoundaryConditionIndicator indicator) {
        boundaryConditionIndicator_ = indicator;
    }

   private:
    BoundaryConditionIndicator boundaryConditionIndicator_;
};

template <std::size_t dim>
class [[maybe_unused]] ConstantHarmonicProblem
    : public SampleHarmonicProblem<dim> {
   public:
    ConstantHarmonicProblem(double omega = 1) : omega_(omega) {
        for (std::size_t i = 0; i < dim; ++i) {
            field_[i] = static_cast<double>(i);
        }
    }
    ConstantHarmonicProblem(LinearAlgebra::SmallVectorC<dim> field,
                            double omega = 1)
        : field_(field), omega_(omega) {}

    double omega() const override { return omega_; }
    LinearAlgebra::SmallVectorC<dim> sourceTerm(
        const Geometry::PointPhysical<dim>& point) const override {
        return -(omega_ * omega_) * field_;
    }
    LinearAlgebra::SmallVectorC<dim> exactSolution(
        const Geometry::PointPhysical<dim>& point) const override {
        return field_;
    }
    LinearAlgebra::SmallVectorC<dim> exactSolutionCurl(
        const Geometry::PointPhysical<dim>& point) const override {
        return {};
    }

   private:
    double omega_;
    LinearAlgebra::SmallVectorC<dim> field_;
};

/// Sample problem in 3D using a product of sin-functions for each coordinate
///
/// Problem used in:
/// "Optimal Penalty Parameters for Symmetric Discontinous Galerkin
/// Discretizations of the Time-Harmonic Maxwell Equations"
/// D Sarmany, F Izsak and JJW van der Vegt
/// J. Sci. Comput. (2010) 44: 219-254
/// DOI: 10.1007/s10915-010-9366-1
class [[maybe_unused]] SarmanyHarmonicProblem
    : public SampleHarmonicProblem<3> {
   public:
    SarmanyHarmonicProblem(double omega) : omega_(omega) {}

    double omega() const override { return omega_; }

    LinearAlgebra::SmallVectorC<3> exactSolution(
        const Geometry::PointPhysical<3>& point) const override;
    LinearAlgebra::SmallVectorC<3> exactSolutionCurl(
        const Geometry::PointPhysical<3>& point) const override;
    LinearAlgebra::SmallVectorC<3> sourceTerm(
        const Geometry::PointPhysical<3>& point) const override;

   private:
    double omega_;
};

/// Harmonic test case for the plane wave
///     E0 exp[i (k.x + phi)]
/// where
///   - E0 is the direction of the field (may be complex valued)
///   - k is the wave vector
///   - phi is a phase
/// The wave is assumed to propagate in a medium of constant material
/// coefficients.
///
template <std::size_t dim>
class [[maybe_unused]] PlaneWaveProblem : public SampleHarmonicProblem<dim> {
   public:
    PlaneWaveProblem(LinearAlgebra::SmallVector<dim> k,
                     LinearAlgebra::SmallVectorC<dim> E0, double omega,
                     double phase, double epsilon)
        : k_(k), E0_(E0), omega_(omega), phase_(phase), epsilon_(epsilon) {
        logger.assert_debug(epsilon > 1e-3, "Negative or almost zero epsilon");
        logger.assert_debug(k.l2Norm() > 1e-9, "Zero wave vector k");
        logger.assert_debug(E0_.l2Norm() > 1e-9, "Zero field direction E0");
        logger.assert_debug(std::abs(k * E0) < 1e-9,
                            "Non orthogonal wave vector and field direction");
    }

    static PlaneWaveProblem onDispersionPlaneWave(
        LinearAlgebra::SmallVector<dim> k, LinearAlgebra::SmallVectorC<dim> E0,
        double epsilon, double phase = 0.0) {
        return PlaneWaveProblem<dim>{
            k, E0, std::sqrt(k.l2NormSquared() / epsilon), phase, epsilon};
    }

    double omega() const override { return omega_; }

    LinearAlgebra::SmallVectorC<dim> exactSolution(
        const Geometry::PointPhysical<dim>& point) const override {
        return E0_ * pointPhase(point);
    }
    LinearAlgebra::SmallVectorC<dim> exactSolutionCurl(
        const Geometry::PointPhysical<dim>& point) const override {
        using namespace std::complex_literals;

        return 1i * (k_.crossProduct(E0_)) * pointPhase(point);
    }
    LinearAlgebra::SmallVectorC<dim> sourceTerm(
        const Geometry::PointPhysical<dim>& point) const override {
        return (k_.l2NormSquared() - epsilon_ * omega_ * omega_) *
               exactSolution(point);
    }

   private:
    std::complex<double> pointPhase(
        const Geometry::PointPhysical<dim>& point) const {
        using namespace std::complex_literals;
        return std::exp(1i * (phase_ + k_ * point.getCoordinates()));
    }

    LinearAlgebra::SmallVector<dim> k_;
    LinearAlgebra::SmallVectorC<dim> E0_;
    double omega_;
    double phase_;
    double epsilon_;
};

}  // namespace DGMax
#endif  // HPGEM_APP_SAMPLEHARMONICPROBLEMS_H

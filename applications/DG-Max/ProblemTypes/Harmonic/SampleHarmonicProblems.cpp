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

#include "SampleHarmonicProblems.h"

using namespace hpgem;

namespace DGMax {

LinearAlgebra::SmallVectorC<3> SarmanyHarmonicProblem::exactSolution(
    const Geometry::PointPhysical<3> &point) const {
    LinearAlgebra::SmallVectorC<3> result;
    double sx = sin(M_PI * point[0]), sy = sin(M_PI * point[1]),
           sz = sin(M_PI * point[2]);
    result[0] = sy * sz;
    result[1] = sz * sx;
    result[2] = sx * sy;
    return result;
}
LinearAlgebra::SmallVectorC<3> SarmanyHarmonicProblem::exactSolutionCurl(
    const Geometry::PointPhysical<3> &point) const {
    LinearAlgebra::SmallVectorC<3> result;
    double x = point[0], y = point[1], z = point[2];
    result[0] = sin(M_PI * x) * (cos(M_PI * y) - cos(M_PI * z));
    result[1] = sin(M_PI * y) * (cos(M_PI * z) - cos(M_PI * x));
    result[2] = sin(M_PI * z) * (cos(M_PI * x) - cos(M_PI * y));
    result *= M_PI;
    return result;
}

LinearAlgebra::SmallVectorC<3> SarmanyHarmonicProblem::sourceTerm(
    const Geometry::PointPhysical<3> &point) const {
    return exactSolution(point) * (2 * M_PI * M_PI - omega_ * omega_);
}

template class ConstantHarmonicProblem<2>;
template class ConstantHarmonicProblem<3>;

template <std::size_t dim>
LinearAlgebra::SmallVectorC<dim> PlaneWaveReflectionProblem<dim>::exactSolution(
    const Geometry::PointPhysical<dim> &point) const {

    using namespace std::complex_literals;
    LinearAlgebra::SmallVectorC<dim> result;
    result[1] = 1.0;

    double x = point.getCoordinate(0);

    if (point.getCoordinates()[0] < interfacePosition_) {
        result *= (incidentPhasor_ * std::exp(1i * k1_ * x) +
                   reflectionPhasor_ * std::exp(-1i * k1_ * x));
    } else {
        result *= transmissionPhasor_ * std::exp(1i * k2_ * x);
    }
    return result;
}

template <std::size_t dim>
LinearAlgebra::SmallVectorC<dim>
    PlaneWaveReflectionProblem<dim>::exactSolutionCurl(
        const Geometry::PointPhysical<dim> &point) const {

    double x = point.getCoordinate(0);

    using namespace std::complex_literals;
    LinearAlgebra::SmallVectorC<dim> result;
    LinearAlgebra::SmallVectorC<dim> kdir;
    kdir[0] = 1.0;
    result[1] = 1.0;
    // Curl E0 exp(ikx) = ik x E0 exp(ikx)
    result = kdir.crossProduct(result);

    if (x < interfacePosition_) {
        result *= 1i * k1_ *
                  (incidentPhasor_ * std::exp(1i * k1_ * x) -
                   reflectionPhasor_ * std::exp(-1i * k1_ * x));
    } else {
        result *= 1i * k2_ * transmissionPhasor_ * std::exp(1i * k2_ * x);
    }
    return result;
}

template <std::size_t dim>
BoundaryConditionType PlaneWaveReflectionProblem<dim>::getBoundaryConditionType(
    const Base::Face &face) const {
    auto normal = face.getNormalVector(
        face.getReferenceGeometry()->getCenter().castDimension<dim - 1>());
    normal /= normal.l2Norm();

    if (std::abs(normal[0]) > 0.9) {
        // Normal in x direction
        return BoundaryConditionType::SILVER_MULLER;
    } else if (std::abs(normal[1]) > 0.9) {
        // Normal in y-direction, E-field parallel to the normal
        // use perfect electric conductor
        return BoundaryConditionType::DIRICHLET;
    } else if (dim > 2 && std::abs(normal[2]) > 0.9) {
        // Normal in z direction, field parallel to the face
        // use perfect magnetic conductor
        return BoundaryConditionType::NEUMANN;
    } else {
        logger.assert_always(
            false, "Problem not set up for a face in normal direction %",
            normal);
        return BoundaryConditionType::DIRICHLET;
    }
}
template <std::size_t dim>
LinearAlgebra::SmallVectorC<dim>
    PlaneWaveReflectionProblem<dim>::boundaryCondition(
        Base::PhysicalFace<dim> &face) const {
    return ExactHarmonicProblem<dim>::boundaryCondition(face);
}

template class PlaneWaveReflectionProblem<2>;
template class PlaneWaveReflectionProblem<3>;
}  // namespace DGMax
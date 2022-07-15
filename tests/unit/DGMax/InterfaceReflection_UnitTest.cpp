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

#include "../catch.hpp"
#include <ProblemTypes/Harmonic/InterfaceReflectionField.h>

using namespace DGMax;
using namespace hpgem;

template <std::size_t dim>
using VecR = hpgem::LinearAlgebra::SmallVector<dim>;

template <std::size_t dim>
void checkContinuity(const InterfaceReflectionField<dim>& field, VecR<dim> p) {
    VecR<dim> normal;
    normal[0] = 1.0;
    auto left = p, right = p;
    left -= 1e-8 * normal;
    right += 1e-8 * normal;
    auto leftTangentialField = field.field(left).crossProduct(normal);
    auto rightTangentialField = field.field(right).crossProduct(normal);
    auto diff = leftTangentialField - rightTangentialField;
    double l2Diff = diff.l2Norm();
    logger(INFO, "Difference % (=norm %)", l2Diff, diff);
    CHECK(l2Diff < 1e-4);
}

TEST_CASE("Continuity at zero", "[InterfaceReflectionField]") {
    Material mat1{12.1, 1.0}, mat2{1.0, 1.0};
    double omega = 2.0, phase = 0.0;
    InterfaceReflectionField<3> field(omega, phase, mat1, mat2, 0);

    INFO("At origin");
    checkContinuity(field, {0.0, 0.0, 0.0});
    INFO("At 0,1,1");
    checkContinuity(field, {0.0, 1.0, 1.0});
}

TEST_CASE("Continuity away from zero", "[InterfaceReflectionField]") {
    Material mat1{12.1, 1.1}, mat2{1.2, 5.3};
    double omega = 2.0, phase = 1.3;
    double xinterface = 5.3;
    InterfaceReflectionField<3> field(omega, phase, mat1, mat2, xinterface);

    INFO("At interface origin");
    checkContinuity(field, {xinterface, 0.0, 0.0});
    INFO("At xinterface,1,1");
    checkContinuity(field, {xinterface, 1.0, 1.0});
}

TEST_CASE("Check zero source", "InterfaceReflectionField") {
    Material mat1{7.1, 1.1}, mat2{1.1, 3.2};

    double omega = 3.0, phase = 1.5;
    double xinterface = 0.0;
    PlaneWave<3> incident({omega * mat1.getRefractiveIndex(), 0.0, 0.0},
                          {0, 1.0, 0}, phase);
    VecR<3> normal = {1.0, 0, 0};
    InterfaceReflectionField<3> field(incident, mat1, mat2, normal, xinterface);

    INFO("Material one")
    VecR<3> p = {-1.0, 32.1, 1.54};
    auto source =
        field.fieldDoubleCurl(
            p, MaterialTensor(mat1.getPermeability()).adjoint()) -
        omega * omega *
            MaterialTensor(mat1.getPermittivity()).applyDiv(field.field(p));
    CHECK(source.l2NormSquared() < 1e-8);
    INFO("Material 2")
    p[0] = 1.0;
    source =
        field.fieldDoubleCurl(
            p, MaterialTensor(mat2.getPermeability()).adjoint()) -
        omega * omega *
            MaterialTensor(mat2.getPermittivity()).applyDiv(field.field(p));
    CHECK(source.l2NormSquared() < 1e-8);
}
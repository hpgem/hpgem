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

#include "Base/CommandLineOptions.h"
#include "Base/MeshManipulator.h"
#include "hpgem-cmake.h"
#include "Integration/ElementIntegral.h"
#include "Integration/QuadratureRules/AllGaussQuadratureRules.h"

using namespace hpgem;

double runWithMesh(const std::string& fileName) {
    using namespace hpgem;
    Base::ConfigurationData configurationData(1);
    Base::MeshManipulator<2> mesh(&configurationData);
    mesh.readMesh(fileName);
    Integration::ElementIntegral<2> integral;
    double result = 0.0;
    auto* rule = QuadratureRules::AllGaussQuadratureRules::instance().getRule(
        &Geometry::ReferenceTriangle::Instance(), 2);
    for (const Base::Element* element : mesh.getElementsList()) {
        result += integral.integrate(
            element, [](Base::PhysicalElement<2>&) { return 1.0; }, rule);
    }
    return M_PI - result;
}

int main(int argc, char** argv) {
    using namespace std::string_literals;
    Base::parse_options(argc, argv);

    // Followed by [1-4].hpgem
    std::string meshPrefix =
        getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/circle-domain-quadratic-N";
    for (std::size_t i = 0; i < 4; ++i) {
        std::stringstream meshName;
        meshName << meshPrefix << (i + 1) << ".hpgem";
        double error = runWithMesh(meshName.str());
        std::cout << error << std::endl;
    }
}
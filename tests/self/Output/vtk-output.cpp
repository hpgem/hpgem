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
#include "Output/VTKSpecificTimeWriter.h"
#include "CMakeDefinitions.h"

/// TEST APPROACH ///
/////////////////////
//
// To test the output we generate the output and compare it to previous
// (manually checked) versions. When adding a new case beware that you need to
// add the reference, otherwise it won't be checked. Moreover, when adding/
// updating a reference file it MUST be checked manually to verify that it
// matches the intended output.
//
// The comparison of the output with the reference output is done via cmake.

namespace hpgem {

/**
 * Basic testing template, reads a mesh outputs a VTK file with a single data
 * set containing the function evaluation.
 * @tparam DIM The dimension of the mesh
 * @param meshFile The name of the mesh file
 * @param outputFile The base name for the VTK output file
 * @param valueFunction Function to evaluate and plot in the output. It will get
 *  the Physical position as input and should give a single double that will be
 *  used as output.
 */
template <std::size_t DIM>
void runBasicTest(
    std::string meshFile, std::string outputFile,
    std::function<double(Geometry::PointPhysical<DIM>)> valueFunction) {
    Base::ConfigurationData config(1);
    Base::MeshManipulator<DIM> mesh(&config);
    mesh.readMesh(Base::getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/" +
                  meshFile);
    Output::VTKSpecificTimeWriter<DIM> writer(outputFile, &mesh);
    writer.write(
        [&valueFunction](Base::Element* element,
                         const Geometry::PointReference<DIM>& pref,
                         std::size_t) {
            return valueFunction(element->referenceToPhysical(pref));
        },
        "value");
}

void test1D() {
    // Linear output 'x'
    runBasicTest<1>("unitLineN1.hpgem", "unitLineSegment-linear",
                    [](Geometry::PointPhysical<1> p) { return p[0]; });

    // Quadratic output: '1+x-x^2'
    runBasicTest<1>(
        "unitLineN1.hpgem", "unitLineSegment-quadratic",
        [](Geometry::PointPhysical<1> p) { return 1 + p[0] - p[0] * p[0]; });
}

void test2DSquare() {
    // Linear 1 + x + 2*y (=1+corner number on reference square)
    runBasicTest<2>(
        "unitSquareN1.hpgem", "unitSquare-linear",
        [](Geometry::PointPhysical<2> p) { return 1 + p[0] + 2 * p[1]; });
}

void test2DTriangles() {
    // Linear 1 + x + 2*y
    runBasicTest<2>(
        "unitSquareTrianglesN1.hpgem", "unitSquareTriangles-linear",
        [](Geometry::PointPhysical<2> p) { return 1 + p[0] + 2 * p[1]; });
}

void test3DCube() {
    // Linear 1 + x + 2*y + 4*z
    runBasicTest<3>("unitCubeN1.hpgem", "unitCube-linear",
                    [](Geometry::PointPhysical<3> p) {
                        return 1 + p[0] + 2 * p[1] + 4 * p[2];
                    });
}

void test3DTetrahedrons() {
    // Linear 1 + x + 2*y + 4*z
    runBasicTest<3>("unitCubeTetrahedronsN1.hpgem", "unitTetrahedrons-linear",
                    [](Geometry::PointPhysical<3> p) {
                        return 1 + p[0] + 2 * p[1] + 4 * p[2];
                    });
}

}  // namespace hpgem

int main(int argc, char** argv) {
    hpgem::Base::parse_options(argc, argv);
    hpgem::test1D();
    hpgem::test2DSquare();
    hpgem::test2DTriangles();
    hpgem::test3DCube();
    hpgem::test3DTetrahedrons();

    return 0;
}
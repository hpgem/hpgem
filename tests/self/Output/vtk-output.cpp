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
 * @param order The polynomial order of the plotting
 */
template <std::size_t DIM>
void runBasicTest(
    std::string meshFile, std::string outputFile,
    std::function<double(Geometry::PointPhysical<DIM>)> valueFunction,
    std::size_t order = 1) {
    Base::ConfigurationData config(1);
    Base::MeshManipulator<DIM> mesh(&config);
    mesh.readMesh(Base::getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/" +
                  meshFile);
    Output::VTKSpecificTimeWriter<DIM> writer(outputFile, &mesh, 0, order);
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

    // Quadratic output: '1+x-0.5*x^2'
    // Factor 0.5 to ensure that it has different values at 0, 0.5 and 1
    runBasicTest<1>(
        "unitLineN1.hpgem", "unitLineSegment-quadratic",
        [](Geometry::PointPhysical<1> p) {
            return 1 + p[0] - 0.5 * p[0] * p[0];
        },
        2);
}

void test2DSquare() {
    // Linear 1 + x + 2*y (=1+corner number on reference square)
    runBasicTest<2>(
        "unitSquareN1.hpgem", "unitSquare-linear",
        [](Geometry::PointPhysical<2> p) { return 1 + p[0] + 2 * p[1]; });
    // Quadratic 1 + 15x - 14x^2 + 18y - 20xy + 20x^2y - 16y^2 + 24xy^2 -
    // 24x^2y^2
    // This is constructed so that the corner points are 1-4 (in reference
    // order) Midpoints are 5-8 (x=0, y=0, y=1,x=1 edges) and midpoint is 9.
    runBasicTest<2>(
        "unitSquareN1.hpgem", "unitSquare-quadratic",
        [](Geometry::PointPhysical<2> p) {
            double x = p[0], x2 = x * x;
            double y = p[1], y2 = y * y;
            return 1 + 15 * x - 14 * x2 + 18 * y - 20 * x * y + 20 * x2 * y -
                   16 * y2 + 24 * x * y2 - 24 * x2 * y2;
        },
        2);

    runBasicTest<2>(
        "unitSquareN1.hpgem", "unitSquare-quartic",
        [](Geometry::PointPhysical<2> p) {
            double x = p[0], x2 = x * x;
            double y = p[1], y2 = y * y;
            return 1 + x2 * x2 - y2 * y2 + x2 * y2;
        },
        4);
}

void test2DTriangles() {
    // Linear 1 + x + 2*y
    runBasicTest<2>(
        "unitSquareTrianglesN1.hpgem", "unitSquareTriangles-linear",
        [](Geometry::PointPhysical<2> p) { return 1 + p[0] + 2 * p[1]; });

    // Quadratic 1 + 11x + 18y - 10x^2 - 16xy - 16y^2
    // Result: 1-2-3 at the corners, 4-5-6 at the midpoints (all in CCW order)
    runBasicTest<2>(
        "unitSquareTrianglesN1.hpgem", "unitSquareTriangles-quadratic",
        [](Geometry::PointPhysical<2> p) {
            double x = p[0];
            double y = p[1];
            return 1 + 11 * x + 18 * y - 10 * x * x - 16 * x * y - 16 * y * y;
        },
        2);

    // Higher order (4-th order)
    runBasicTest<2>(
        "unitSquareTrianglesN1.hpgem", "unitSquareTriangles-quartic",
        [](Geometry::PointPhysical<2> p) {
            double x = p[0];
            double y = p[1];
            // 1 + x^4 - y^4 + x^2y^2
            return 1 + x * x * x * x - y * y * y * y + x * x * y * y;
        },
        4);

    // Higher order (6-th order)
    runBasicTest<2>(
        "unitSquareTrianglesN1.hpgem", "unitSquareTriangles-sextic",
        [](Geometry::PointPhysical<2> p) {
            double x = p[0];
            double x3 = x * x * x;
            double y = p[1];
            double y2 = y * y;
            // 1 + x^6 - y^4 + x^3y^3
            return 1 + x3 * x3 - y2 * y2 + x3 * y2 * y;
        },
        6);
}

void test3DCube() {
    // Linear 1 + x + 2*y + 4*z
    runBasicTest<3>("unitCubeN1.hpgem", "unitCube-linear",
                    [](Geometry::PointPhysical<3> p) {
                        return 1 + p[0] + 2 * p[1] + 4 * p[2];
                    });
    // Quadratic function
    // Such that it has values 1-8 + (2z*8), where the values 1-8 are used
    // linearly for the points from x=0-1,y=0 then x=0-1,y=1/2 and then
    // x=0-1/2,y=1. The values 25-27 are used at x=1,y=1 with z=1 being 25 and
    // z=0 being 27.
    runBasicTest<3>(
        "unitCubeN1.hpgem", "unitCube-quadratic",
        [](Geometry::PointPhysical<3> p) {
            double x = p[0], x2 = x * x;
            double y = p[1], y2 = y * y;
            double z = p[2], z2 = z * z;
            return 1 + 2 * x + 6 * y + 17 * x * y - 34 * (x2 * y + x * y2) +
                   68 * x2 * y2 +
                   z * (16 - 15 * x * y + 30 * (x2 * y + x * y2) -
                        60 * x2 * y2) +
                   z2 * (-2 * x * y + 4 * (x2 * y + y2 * x) - 8 * x2 * y2);
        },
        2);

    runBasicTest<3>(
        "unitCubeN1.hpgem", "unitCube-quartic",
        [](Geometry::PointPhysical<3> p) {
            double x = p[0], x2 = x * x;
            double y = p[1], y2 = y * y;
            double z = p[2], z2 = z * z;
            return 1 + x2 + y * y2 + z2 * z2 - x * y * z;
        },
        4);
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
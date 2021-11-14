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

#include <chrono>

#include <Logger.h>
#include <Base/CommandLineOptions.h>
#include <Base/MeshManipulator.h>
#include <Integration/ElementIntegral.h>
#include <Integration/FaceIntegral.h>

#include "../ConvergenceTest.h"
#include "../TestMeshes.h"

using namespace hpgem;

double solveArea(const std::string& meshFileName, std::size_t) {
    // Compute the error in computing the area by integration
    Base::ConfigurationData config(0);
    Base::MeshManipulator<2> mesh(&config, 0, 0, 0, 0);
    mesh.readMesh(meshFileName);
    double area = 0.0;
    Integration::ElementIntegral<2> integral;
    for (Base::Element* element : mesh.getElementsList()) {
        area += integral.integrate(
            element, [](Base::PhysicalElement<2>&) { return 1.0; });
    }
    return area - M_PI;
}

double solvePerimeter(const std::string& meshFileName, std::size_t) {
    // Compute the error in computing the perimeter using a surface integral of
    // a constant.
    Base::ConfigurationData config(0);
    Base::MeshManipulator<2> mesh(&config, 0, 0, 0, 0);
    mesh.readMesh(meshFileName);
    double perimeter = 0.0;
    Integration::FaceIntegral<2> integrator;
    for (Base::Face* face : mesh.getFacesList()) {
        if (face->isInternal()) {
            continue;
        }
        perimeter += integrator.integrate(
            face, [](Base::PhysicalFace<2>&) { return 1.0; });
    }
    return perimeter - 2 * M_PI;
}

int main(int argc, char** argv) {
    using namespace std::string_literals;
    Base::parse_options(argc, argv);

    // The following tests are inspired by a deal.ii tutorial for the usage of
    // higher order mappings. At the time of writing this was the 'step 10'
    //

    // Define clocks for measuring simulation time.
    std::chrono::time_point<std::chrono::system_clock> startClock, endClock;
    startClock = std::chrono::system_clock::now();

    // For recomputing the error tables
    bool ignoreFailures = true;

    ConvergenceTestSet set = {
        getUnitCircleQuadraticTriangleMeshes(),
        {
            // Fourth order convergence = 2p, meaning super convergence
            -2.44508328e-03,  //------
            -1.54936886e-04,  // 15.78
            -9.71694789e-06,  // 15.95
            -6.07832133e-07,  // 15.99
        }};
    runConvergenceTest(set, ignoreFailures, solveArea);

    ConvergenceTestSet set2 = {
        getUnitCircleQuadraticTriangleMeshes(),
        {
            // Fourth order convergence = 2p, meaning super convergence
            -2.40538481e-03,  //------
            -1.54283939e-04,  // 15.59
            -9.70661212e-06,  // 15.89
            -6.07670078e-07,  // 15.97
        }};
    runConvergenceTest(set, ignoreFailures, solvePerimeter);

    endClock = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = endClock - startClock;
    logger(INFO, "Elapsed time %", elapsed_seconds.count());

    return 0;
}
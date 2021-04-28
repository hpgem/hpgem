/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2020, University of Twente
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

#include "DGMaxLogger.h"
#include "DGMaxProgramUtils.h"
#include "Utils/Verification/DGMaxEVConvergenceTest.h"
#include "Utils/Verification/EVTestPoint.h"

#include "Utils/HomogeneousBandStructure.h"

#include "testMeshes.h"

// Convergence test for DGMaxEigenvalue, based on computing the spectrum for
// Vacuum at a single k-point. The resulting frequencies are cached in the
// expected results below. If these are no longer correct the convergence test
// should be re-run and the expected results updated.

auto& runAsTestArg = Base::register_argument<bool>(
    't', "test",
    "Whether to run as test (true, default) or as convergence study (false)",
    false, true);

// clang-format off
// Leave the results as an easily readable table
DGMax::EVConvergenceResult expected ({
     {1.204336889,5.411086151,5.524294257,7.066481048,7.149443644,7.72780029},
     {1.204203888,5.43456121,5.548559589,7.1219632,7.208314184,7.695794404,8.910736644},
     {1.20417057,5.440372736,5.554561455,7.13560546,7.222790797,7.68700507,8.900260151,9.040610049},
     {1.204162236,5.441822184,5.556058059,7.139002763,7.226395987,7.684761581,8.897550726}

 });

DGMax::EVConvergenceResult expected2({
    {1.2043368891,5.4110861511,5.5242942569,7.0664810486,7.1494436437,7.7278002899,8.9463017271,9.0890054436,10.1789455152,11.3362127914,11.4352596910},
    {1.2042038880,5.4345612104,5.5485595891,7.1219632000,7.2083141835,7.6957944040,8.9107366441,9.0516696919,10.1143353023,11.6076760680,11.7128457455},
    {1.2041705699,5.4403727357,5.5545614554,7.1356054603,7.2227907973,7.6870050696,8.9002601508,9.0406100487,10.0948497388,11.6724130103},
    {1.2041622361,5.4418221838,5.5560580595,7.1390027626,7.2263959872,7.6847615808,8.8975507256,9.0377472103,10.0897957246,11.6884390538}
});

// clang-format on

int main(int argc, char** argv) {

    Base::parse_options(argc, argv);
    initDGMaxLogging();

    // Flag to switch between testing with known results (true) and checking the
    // results when they change (false).
    bool runAsTest = runAsTestArg.getValue();

    std::vector<std::string> meshes =
        DGMaxTest::singleProcessorRefinementMeshes2D();

    // Just a random point in vacuum
    DGMax::EVTestPoint<2> testPoint(LinearAlgebra::SmallVector<2>({0.8, 0.9}),
                                    0, 10);

    DGMaxEigenvalueBase::SolverConfig config;
    config.stab_ = 100;
    config.shiftFactor_ = 0;
    config.useHermitian_ = false;
    config.useProjector_ = DGMaxEigenvalueBase::NONE;
    DGMax::DGMaxEVConvergenceTest<2> testCase(testPoint, meshes, 1e-8, 1,
                                              config, &expected);
    DGMax::EVConvergenceResult result = testCase.run(runAsTest);

    // Second test of the same algorithm, but with completely different
    // settings, in an attempt to cover as many different paths with the least
    // number of tests.
    // TODO: Modify this to use the Braggstack (need: git fix theory to work in
    // 2D)
    config.shiftFactor_ = -1.0;
    config.useHermitian_ = true;
    config.useProjector_ = DGMaxEigenvalueBase::ALL;
    DGMax::DGMaxEVConvergenceTest<2> testCase2(testPoint, meshes, 1e-8, 1,
                                               config, &expected2);
    DGMax::EVConvergenceResult result2 = testCase2.run(runAsTest);

    // Code to check the results if they change
    // Expected convergence speed = 2^2p = 4
    if (!runAsTest) {
        std::array<LinearAlgebra::SmallVector<2>, 2> reciprocal;
        // Standard cubic lattice.
        reciprocal[0] = LinearAlgebra::SmallVector<2>({2 * M_PI, 0});
        reciprocal[1] = LinearAlgebra::SmallVector<2>({0, 2 * M_PI});
        HomogeneousBandStructure<2> analytical(reciprocal);
        std::vector<double> linear =
            analytical.computeLinearSpectrum(testPoint.getKPoint(), 20);

        std::cout << "Raw frequency table" << std::endl;
        result.printFrequencyTable(linear);

        result.printResultCode(10);

        std::cout << "Cleanup results" << std::endl;
        result.filterResults(
            0.01, true);  // Small eigenvalues are to be expected from DGMax
        result.printFrequencyTable(linear);
        result.printErrorTable(linear);

        // Some spacing between the results for the two seperate tests
        std::cout << "\n\n";

        std::cout << "Raw frequency table 2" << std::endl;
        result2.printFrequencyTable(linear);

        result2.printResultCode(10);

        std::cout << "Cleanup results" << std::endl;
        result2.filterResults(
            0.01, true);  // Small eigenvalues are to be expected from DGMax
        result2.printFrequencyTable(linear);
        result2.printErrorTable(linear);
    }
}

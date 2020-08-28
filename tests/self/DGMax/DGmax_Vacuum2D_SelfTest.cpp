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

// Note: Several nan results are from computing negative eigenvalues. As
// frequency = sqrt(eigenvalue), the resulting frequencies are NaN.
DGMax::EVConvergenceResult expected ({
     {2.268990566e-06,3.563100846e-06,1.204336889,5.411086151,5.524294257,7.066481048,7.149443644,7.72780029},
     {std::nan(""),2.812220581e-06,1.204203888,5.43456121,5.548559589,7.1219632,7.208314184,7.695794404,8.910736644},
     {std::nan(""),std::nan(""),1.641925264e-05,1.20417057,5.440372736,5.554561455,7.13560546,7.222790797,7.68700507,8.900260151,9.040610049},
     {std::nan(""),std::nan(""),0.0001017320878,1.204162236,5.441822184,5.556058059,7.139002763,7.226395987,7.684761581,8.897550726}
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
    }
}

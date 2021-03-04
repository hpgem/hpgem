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
#include "Utils/Verification/DivDGMaxEVConvergenceTest.h"
#include "Utils/Verification/EVTestPoint.h"

#include "Utils/HomogeneousBandStructure.h"

#include "testMeshes.h"

// Convergence test for DivDGMaxEigenvalue, based on computing the spectrum for
// Vacuum at a single k-point. The resulting frequencies are cached in the
// expected results below. If these are no longer correct the convergence test
//// should be re-run and the expected results updated.

auto& runAsTestArg = Base::register_argument<bool>(
    't', "test",
    "Whether to run as test (true, default) or as convergence study (false)",
    false, true);

// clang-format off
// Leave the results as an easily readable table
DGMax::EVConvergenceResult expected ({
    {1.2044351353,5.4116291469,5.5262052437,7.0674801624,7.1468756651,7.7213397212,8.9374986496,9.0778132093,10.1643930579,11.3111802309,11.4252815469},
    {1.2042284160,5.4346771252,5.5490171908,7.1221660674,7.2076240689,7.6940823890,8.9083377743,9.0486671389,10.1102877500,11.6012384097,11.7101182626},
    {1.2041766998,5.4404004795,5.5546745951,7.1356532352,7.2226152771,7.6865715589,8.8996495507,9.0398483410,10.0938166890,11.6707922206,11.7782466337},
    {1.2041637685,5.4418290426,5.5560862656,7.1390145224,7.2263519206,7.6846528664,8.8973974149,9.0375561136,10.0895361929,11.6880331460,11.7951266380}
});
// clang-format on

int main(int argc, char** argv) {
    Base::parse_options(argc, argv);
    initDGMaxLogging();

    // The build-in (default) Petsc LU solver has problems with zero pivots.
    // Other solvers (mumps, umfpack, etc.) do not have this problem, probably
    // due to a different (better) reordering. So we default the use of a
    // different solver.
    {
#if PETSC_VERSION_LE(3, 8, 0)
        // This was replaced. The corresponding function does not show up in the
        // documentation for petsc-3.9.
        const char option[] = "-st_pc_factor_mat_solver_package";
#else
        const char option[] = "-st_pc_factor_mat_solver_type";
#endif
        PetscBool present;
        PetscOptionsHasName(nullptr, nullptr, option, &present);
        if (!present) {
            DGMaxLogger(INFO, "Defaulting LU solver");
            // UMFPack is chosen because it is installed as dependency when
            // installing petsc on ubuntu (via apt) or mac (via homebrew).
            // Unfortunately this is the only compatible solver which is
            // installed via homebrew.
            PetscOptionsSetValue(nullptr, option, "umfpack");
        }
    }

    // Flag to switch between testing with known results (true) and checking the
    // results when they change (false).
    bool runAsTest = runAsTestArg.getValue();

    std::vector<std::string> meshes =
        DGMaxTest::singleProcessorRefinementMeshes2D();

    // Just a random point in vacuum
    DGMax::EVTestPoint<2> testPoint(LinearAlgebra::SmallVector<2>({0.8, 0.9}),
                                    0, 10);

    DivDGMaxDiscretization<2>::Stab stab;
    stab.stab1 = 5;
    stab.stab2 = 0;
    stab.stab3 = 5;
    stab.setAllFluxeTypes(DivDGMaxDiscretization<2>::FluxType::BREZZI);

    DGMax::DivDGMaxEVConvergenceTest<2> testCase(testPoint, meshes, 1e-8, 1,
                                                 stab, &expected);
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
        result.printFrequencyTable(linear);
        result.printErrorTable(linear);

        result.printResultCode(10);
    }
}

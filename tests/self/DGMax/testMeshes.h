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
#ifndef HPGEM_TESTMESHES_H
#define HPGEM_TESTMESHES_H

#include <string>
#include <vector>
#include "hpgem-cmake.h"

namespace DGMaxTest {
    bool isParallelRun(){
        // Check if code is running in parallel or serially
        bool isParallelRun = false;
#ifdef HPGEM_USE_MPI
        int numberOfProcessors = 1;
        int maxNumberOfProcessors = 2;  // the test is designed to run only for one (serial) or two (parallel) processors
        MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcessors);
        logger.assert_always((numberOfProcessors > 0) && (numberOfProcessors <= maxNumberOfProcessors),
                              "Parallel run only up to two processors");
        std::cout << "\nNumber of processors: " << numberOfProcessors;
        switch (numberOfProcessors){
            case 1:
                // Only one processor so we are running in serial mode
                std::cout << "\nSerial mode!\n";
                isParallelRun = false;
                break;
            case 2:
                // Two Processors so we are running in parallel mode (we can only run with 1 or 2 processors)
                std::cout << "\nParallel mode!\n";
                isParallelRun = true;
                break;
        }
#endif // HPGEM_USE_MPI

        return isParallelRun;
    }

    std::vector<std::string> singleProcessorRefinementMeshes2D() {

        using namespace std::string_literals;

        std::string prefix = hpgem::getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/"s;

        return std::vector<std::string>(
            {prefix + "unitPeriodicSimplexD2N8P1.hpgem",
            prefix + "unitPeriodicSimplexD2N16P1.hpgem",
            prefix + "unitPeriodicSimplexD2N32P1.hpgem",
            prefix + "unitPeriodicSimplexD2N64P1.hpgem"});
    }

    std::vector<std::string> dualProcessorRefinementMeshes2D() {

        using namespace std::string_literals;

        std::string prefix = hpgem::getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/parallel_meshes/"s;

        return std::vector<std::string>(
            {prefix + "unitPeriodicSimplexD2N8P1.hpgem",
            prefix + "unitPeriodicSimplexD2N16P1.hpgem",
            prefix + "unitPeriodicSimplexD2N32P1.hpgem",
            prefix + "unitPeriodicSimplexD2N64P1.hpgem"});
    }

    std::vector<std::string> refinementMeshes2D() {
        std::vector<std::string> meshFilenames;

        if (DGMaxTest::isParallelRun()){
          std::cout << "\nGenerating parallel meshes!\n";
          meshFilenames = dualProcessorRefinementMeshes2D();
          std::cout << "\n" << meshFilenames[0] << "\n";
        }
        else {
          std::cout << "\nGenerating serial meshes!\n";
          meshFilenames = singleProcessorRefinementMeshes2D();
          std::cout << "\n" << meshFilenames[0] << "\n";
        }

        return meshFilenames;
    }
}  // namespace DGMaxTest

#endif  // HPGEM_TESTMESHES_H

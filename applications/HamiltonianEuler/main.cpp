/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2014, University of Twente
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
#ifdef HPGEM_USE_MPI
// something breaks if this in included later on
#include <mpi.h>
#endif

#include <iostream>
#include <string>
using namespace std;

#include <petscmat.h>
#include <petscksp.h>

#include "HEuler.h"
// ============================================================
static char help[] = "Petsc Help";

int main(int argc, char **argv) {

    PetscInitialize(&argc, &argv, PETSC_NULL, help);

    const string fileName = "blabla.out";

    /// enum  SolutionType 		{INCOMPRESSIBLE_WALLS,
    /// INCOMPRESSIBLE_ONETHIRDPERIODIC, INCOMPRESSIBLE_PERIODIC,
    /// COMPRESSIBLE_WALLS, COMPRESSIBLE_PERIODIC};

    //    unsigned int numberOfUnknowns, unsigned int numberOfBasisFunctions,
    //    unsigned int  numberOfTimeLevels=1,

    HEulerConfigurationData config(4, 11, 1, fileName,
                                   HEulerConfigurationData::COMPRESSIBLE_WALLS);
    HEulerGlobalVariables globalVar;

    if (config.solutionType_ ==
            HEulerConfigurationData::INCOMPRESSIBLE_PERIODIC ||
        config.solutionType_ ==
            HEulerConfigurationData::COMPRESSIBLE_PERIODIC ||
        config.solutionType_ == HEulerConfigurationData::COMPRESSIBLE_WALLS) {
        config.lx_ = 1;
        config.ly_ = 1;
        config.lz_ = 1;
    } else if (config.solutionType_ ==
               HEulerConfigurationData::INCOMPRESSIBLE_ONETHIRDPERIODIC) {
        config.lx_ = 2 * 3.14159265359;
        config.ly_ = 2 * 3.14159265359;
        config.lz_ = 2 * 3.14159265359;
    } else if (config.solutionType_ ==
               HEulerConfigurationData::INCOMPRESSIBLE_WALLS) {
        config.lx_ = 2 * 3.14159265359;
        config.ly_ = 3.14159265359;
        config.lz_ = 3.14159265359;
    }
    config.nx_ = 10;
    config.ny_ = 10;
    config.nz_ = 10;

    HEuler eulerProblem(&globalVar, &config);

    eulerProblem.initialiseMesh();

    eulerProblem.initialConditions();

    eulerProblem.output(0);

    if (config.solutionType_ ==
            HEulerConfigurationData::COMPRESSIBLE_PERIODIC ||
        config.solutionType_ == HEulerConfigurationData::COMPRESSIBLE_WALLS) {
        eulerProblem.createCompressibleSystem();
    } else {
        eulerProblem.createIncompressibleSystem();
    }

    eulerProblem.output(0.0001);

    eulerProblem.solve();

    // eulerProblem.output(2.0);

    return 0;
}

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

#include "../TestMeshes.h"
#include "../ConvergenceTest.h"
#include "HarmonicErrorDriver.h"

#include <Base/CommandLineOptions.h>

#include <DGMaxLogger.h>
#include <DGMaxProgramUtils.h>
#include <Algorithms/HarmonicSolver.h>
#include <Algorithms/DGMaxDiscretization.h>
#include <Algorithms/DivDGMaxDiscretization.h>
#include <ProblemTypes/Harmonic/TrenchReflectionProblem.h>
#include <PMLElementInfos.h>

#include <petsc.h>

auto& runToDebug = Base::register_argument<bool>(
    'd', "debug", "Whether to run in debug mode", false, false);

using namespace DGMax;

struct ProblemData {

    static constexpr const double WIDTH = 1.0;
    static constexpr const double LENGTH = 1.0;
    static constexpr const double FREQUENCY = 2.0 * M_PI;
    static constexpr const std::size_t TRANSVERSE_HALF_WAVE_COUNT = 1;

    static constexpr const double PML_START = 0.5;
    static constexpr const double PML_DEPTH = LENGTH - PML_START;
    // This coefficient results in medium attenuation of
    // exp(-10/FREQUENCY * PML_DEPTH^3/3 * ky) ~= 0.66
    // for the incident wave to the back boundary, and the same factor for the
    // reflected wave. Assuming Dirichlet boundary conditions the reflection is
    // thus about 33%. Which is a significant dampening from no PML (100%
    // reflection) as well as large enough to affect the solution before the PML
    // and thus be detectable in the L2-error.
    static constexpr const double PML_SCALING =
        10.0 * (PML_DEPTH * PML_DEPTH * PML_DEPTH);

    // Non unity material parameters to check that those work.
    static constexpr const double permittivity = 1.5;
    static constexpr const double permeability = 1.3;

    ProblemData()
        : baseMaterial(permittivity, permeability),
          baseInfos(baseMaterial),
          pmlInfos(baseMaterial, {0, PML_START}, {0, 1}, {1, PML_DEPTH},
                   {0, PML_SCALING}),
          structureDescription(StructureDescription::fromFunction(
              [this](const Base::Element* element) -> ElementInfos* {
                  auto centre = element->referenceToPhysical(
                      element->getReferenceGeometry()
                          ->getCenter()
                          .castDimension<2>());
                  if (centre[1] > PML_START) {
                      return &pmlInfos;
                  } else {
                      return &baseInfos;
                  }
              })),
          problem(FREQUENCY, WIDTH, TRANSVERSE_HALF_WAVE_COUNT, LENGTH,
                  baseMaterial) {
        problem.setFarEndBoundaryCondition(BoundaryConditionType::DIRICHLET);
        problem.setPML(PML_START, PML_SCALING, PML_DEPTH);
    }

    Material baseMaterial;
    ElementInfos baseInfos;
    PMLElementInfos<2> pmlInfos;
    std::shared_ptr<StructureDescription> structureDescription;
    DGMax::TrenchReflectionProblem<2> problem;
};

double solve(std::string meshFile, std::size_t level, const std::string& prefix,
             std::shared_ptr<DGMax::AbstractDiscretization<2>> discretization) {
    ProblemData problemData;

    Base::ConfigurationData config(discretization->getNumberOfUnknowns());
    auto mesh =
        DGMax::readMesh<2>(meshFile, &config, *problemData.structureDescription,
                           discretization->getNumberOfElementMatrices());

    DGMax::HarmonicSolver<2> solver(discretization);

    std::stringstream outputfileName;
    outputfileName << prefix << "-" << level;
    Output::VTKSpecificTimeWriter<2> output(outputfileName.str(), mesh.get(), 0,
                                            discretization->getOrder());
    HarmonicErrorDriver<2> driver(problemData.problem);
    driver.setOutputPlotter(&output);
    solver.solve(*mesh, driver);
    return driver.getError();
}

int main(int argc, char** argv) {
    Base::parse_options(argc, argv);
    initDGMaxLogging();

    // Default the solver if not specified to a direct LU solver
    std::map<std::string, std::string> defaultOptions = {
        {"-ksp_type", "preonly"},
        {"-pc_type", "lu"},
        // The DivDGMax system can't be solved with the petsc solver
        {"-pc_factor_mat_solver_type", "umfpack"}};
    for (const auto& option : defaultOptions) {
        PetscBool present;
        PetscOptionsHasName(nullptr, nullptr, option.first.c_str(), &present);
        if (!present) {
            PetscOptionsSetValue(nullptr, option.first.c_str(),
                                 option.second.c_str());
        }
    }
    ConvergenceTestSet meshesDGMax{
        getUnitSquareTriangleMeshes(3),
    };
    meshesDGMax.addExpectedErrors("L2error",
                                  {
                                      4.71091751e-02,  //------
                                      1.17673957e-02,  //  4.00
                                      2.94181156e-03,  //  4.00
                                      7.35406176e-04,  //  4.00
                                  });
    auto dgmax = std::make_shared<DGMaxDiscretization<2>>(2, 100);
    runConvergenceTest(meshesDGMax, runToDebug.getValue(),
                       [&dgmax](std::string meshFile, std::size_t level) {
                           return solve(meshFile, level, "pml-reflection-dgmax",
                                        dgmax);
                       });
}
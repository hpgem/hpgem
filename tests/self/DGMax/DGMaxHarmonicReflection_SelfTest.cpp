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

#include "../ConvergenceTest.h"
#include "../TestMeshes.h"
#include <Base/CommandLineOptions.h>
#include <DGMaxProgramUtils.h>
#include <DGMaxLogger.h>
#include <Algorithms/DGMaxHarmonic.h>
#include <Algorithms/DivDGMaxHarmonic.h>
#include <ProblemTypes/Harmonic/SampleHarmonicProblems.h>

#include <petsc.h>

using namespace DGMax;

// Problem testing a 3D harmonic Maxwell problem with DGMax and DivDGMax.
// It uses first order basis functions to reduce the computational cost,
// additionally only a few coarse meshes are used for the same reason.

struct ProblemData {
    constexpr static const double epsilonLeft = 4.0;
    constexpr static const double epsilonRight = 1.0;
    constexpr static const double omega = M_PI;
    constexpr static const double phase = 1.0;
    constexpr static const double xInterface = 0.5;

    ProblemData()
        : infosLeft(epsilonLeft),
          infosRight(epsilonRight),
          structureDescription(StructureDescription::fromFunction(
              [&](const Base::Element* element) {
                  auto centerPos = element->referenceToPhysical(
                      element->getReferenceGeometry()
                          ->getCenter()
                          .castDimension<3>());
                  if (centerPos[0] < xInterface) {
                      return &infosLeft;
                  } else {
                      return &infosRight;
                  }
              })),
          problem(omega, phase, epsilonLeft, epsilonRight, xInterface) {}

    ElementInfos infosLeft, infosRight;
    std::shared_ptr<StructureDescription> structureDescription;
    PlaneWaveReflectionProblem<3> problem;
};

double solveDGMax(std::string meshFile, std::size_t level) {
    ProblemData problem;

    Base::ConfigurationData config(1);
    auto mesh =
        DGMax::readMesh<3>(meshFile, &config, *problem.structureDescription, 2);

    DGMaxHarmonic<3> solver(*mesh, 100, 1);
    solver.solve(problem.problem);

    std::stringstream fileName;
    fileName << "reflection-solution-dgmax-" << level;
    Output::VTKSpecificTimeWriter<3> output(fileName.str(), mesh.get(), 0, 1);
    solver.writeVTK(output);
    output.write(
        [&](Base::Element* element, const Geometry::PointReference<3>&,
            std::size_t) {
            return static_cast<ElementInfos*>(element->getUserData())->epsilon_;
        },
        "epsilon");

    return solver.computeL2Error(problem.problem);
}

double solveDivDGMax(std::string meshFile, std::size_t level) {
    ProblemData problem;

    Base::ConfigurationData config(2);
    auto mesh =
        DGMax::readMesh<3>(meshFile, &config, *problem.structureDescription, 2);

    DivDGMaxDiscretizationBase::Stab stab;
    stab.stab1 = 5;
    stab.stab2 = 0;
    stab.stab3 = 1;
    stab.setAllFluxeTypes(DivDGMaxDiscretizationBase::FluxType::BREZZI);

    DivDGMaxHarmonic<3> solver(*mesh, stab, 1);
    solver.solve(problem.problem);

    std::stringstream fileName;
    fileName << "reflection-solution-divdgmax-" << level;
    Output::VTKSpecificTimeWriter<3> output(fileName.str(), mesh.get(), 0, 1);
    solver.writeVTK(output);
    output.write(
        [&](Base::Element* element, const Geometry::PointReference<3>&,
            std::size_t) {
            return static_cast<ElementInfos*>(element->getUserData())->epsilon_;
        },
        "epsilon");

    return solver.computeL2Error(problem.problem);
}

int main(int argc, char** argv) {

    Base::parse_options(argc, argv);
    initDGMaxLogging();

    // For testing and updating => Should be false to actually use this test
    bool ignoreFailures = false;

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

    // Skip lowest level as that has too few elements for a left and right side
    // Expected convergence rate: 2 (=2^p with p=1)
    ConvergenceTestSet meshesDGMax = {getUnitCubeTetMeshes(1, 4),
                                      {
                                          6.42536944e-01,  //------
                                          2.47104771e-01,  //  2.60
                                          1.18920432e-01,  //  2.08
                                      }};
    runConvergenceTest(meshesDGMax, ignoreFailures, &solveDGMax);

    ConvergenceTestSet meshesDivDGMax = {getUnitCubeTetMeshes(1, 4),
                                         {
                                             1.02154508e+00,  //------
                                             4.76665783e-01,  //  2.14
                                             2.01204463e-01,  //  2.37
                                         }};
    runConvergenceTest(meshesDivDGMax, ignoreFailures, &solveDivDGMax);
}
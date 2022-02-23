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
#include "HarmonicErrorDriver.h"

#include <Base/CommandLineOptions.h>
#include <DGMaxProgramUtils.h>
#include <DGMaxLogger.h>
#include <Algorithms/HarmonicSolver.h>
#include <Algorithms/DGMaxDiscretization.h>
#include <Algorithms/DivDGMaxDiscretization.h>
#include <ProblemTypes/Harmonic/InterfaceReflectionField.h>
#include <ProblemTypes/Harmonic/ExactFieldHarmonicProblem.h>

#include <petsc.h>

using namespace DGMax;

// Problem testing a 3D harmonic Maxwell problem with DGMax and DivDGMax.
// It uses first order basis functions to reduce the computational cost,
// additionally only a few coarse meshes are used for the same reason.

struct ProblemData {
    constexpr static const Material leftMaterial = Material(4.0, 1.2);
    constexpr static const Material rightMaterial = Material(1.0, 1.1);
    constexpr static const double omega = M_PI;
    constexpr static const double phase = 1.0;
    constexpr static const double xInterface = -0.5;

    ProblemData()
        : infosLeft(leftMaterial),
          infosRight(rightMaterial),
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
          problem(omega,
                  std::make_shared<InterfaceReflectionField<3>>(
                      omega, phase, leftMaterial, rightMaterial, xInterface)) {
        problem.setBoundaryConditionIndicator([](const Base::Face& face) {
            auto normal = face.getNormalVector(
                face.getReferenceGeometry()->getCenter().castDimension<2>());
            normal /= normal.l2Norm();
            if (std::abs(normal[0]) > 1e-1) {
                return BoundaryConditionType::SILVER_MULLER;
            } else if (std::abs(normal[1]) > 1e-1) {
                // Field is in y-direction thus along the normal direction
                // = PEC boundary
                return BoundaryConditionType::DIRICHLET;
            } else if (std::abs(normal[2]) > 1e-1) {
                // Field is parallel to the face => PMC
                return BoundaryConditionType::DIRICHLET;
            } else {
                logger.fail("Problem with setting the boundary condition");
            }
        });
    }

    ElementInfos infosLeft, infosRight;
    std::shared_ptr<StructureDescription> structureDescription;
    ExactFieldHarmonicProblem<3> problem;
};

const Material ProblemData::leftMaterial;
const Material ProblemData::rightMaterial;
const double ProblemData::omega;
const double ProblemData::phase;
const double ProblemData::xInterface;

double solve(std::string meshFile, std::size_t level,
             std::shared_ptr<AbstractDiscretization<3>> discretization,
             std::string prefix) {
    ProblemData problemData;

    Base::ConfigurationData config(discretization->getNumberOfUnknowns());
    auto mesh =
        DGMax::readMesh<3>(meshFile, &config, *problemData.structureDescription,
                           discretization->getNumberOfElementMatrices());

    HarmonicSolver<3> solver(discretization);

    std::stringstream outputfileName;
    outputfileName << prefix << "-" << level;
    Output::VTKSpecificTimeWriter<3> output(outputfileName.str(), mesh.get(), 0,
                                            discretization->getOrder());
    HarmonicErrorDriver<3> driver(problemData.problem);
    driver.setOutputPlotter(&output);
    solver.solve(*mesh, driver);
    return driver.getError();
}

int main(int argc, char** argv) {

    Base::parse_options(argc, argv);
    initDGMaxLogging();

    // For testing and updating => Should be false to actually use this test
    bool ignoreFailures = true;

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
                                          7.39613968e-01,  //------
                                          3.81524577e-01,  //  1.94
                                          1.83728547e-01,  //  2.08
                                      }};
    auto dgmax = std::make_shared<DGMaxDiscretization<3>>(1, 100);
    runConvergenceTest(meshesDGMax, ignoreFailures,
                       [&dgmax](std::string meshFile, std::size_t order) {
                           return solve(meshFile, order, dgmax,
                                        "reflection-solution-dgmax");
                       });

    ConvergenceTestSet meshesDivDGMax = {getUnitCubeTetMeshes(1, 4),
                                         {
                                             1.09286633e+00,  //------
                                             5.36898591e-01,  //  2.04
                                             2.22872838e-01,  //  2.41
                                         }};
    DivDGMaxDiscretizationBase::Stab stab;
    stab.stab1 = 5;
    stab.stab2 = 0;
    stab.stab3 = 1;
    stab.setAllFluxeTypes(DivDGMaxDiscretizationBase::FluxType::BREZZI);
    auto divdgmax = std::make_shared<DivDGMaxDiscretization<3>>(1, stab);
    runConvergenceTest(meshesDivDGMax, ignoreFailures,
                       [&divdgmax](std::string meshFile, std::size_t order) {
                           return solve(meshFile, order, divdgmax,
                                        "reflection-solution-divdgmax");
                       });
}
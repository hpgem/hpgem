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
#include <ProblemTypes/Harmonic/SampleHarmonicProblems.h>

#include <petsc.h>

using namespace DGMax;

using Vec2 = LinearAlgebra::SmallVectorC<2>;
using Point = Geometry::PointPhysical<2>;

struct ProblemData {
   public:
    // Different from 1 to increase the likelihood of detecting a bug
    constexpr static const double epsilon = 2.0;
    // Not on the dispersion curve to require a source term
    // On dispersion would be k.l2Norm()/sqrt(epsilon)
    constexpr static const double omega = 1.5;
    constexpr static const double phase = 0.5;

    static const LinearAlgebra::SmallVector<2> k;
    static const LinearAlgebra::SmallVectorC<2> E0;

    ProblemData()
        : infos(epsilon),
          structureDescription(StructureDescription::fromFunction(
              [this](const Base::Element*) { return &this->infos; })),
          problem(k, E0, omega, phase, epsilon) {
        problem.setBoundaryConditionIndicator([](const Base::Face& face) {
            auto normal = face.getNormalVector(
                face.getReferenceGeometry()->getCenter().castDimension<1>());
            // Different boundary conditions on each boundary
            double xn = normal[0];
            if (std::abs(xn - 1) < 1e-9) {
                // normal = <1,0> => y = 1 boundary
                return BoundaryConditionType::DIRICHLET;
            } else if (std::abs(xn + 1) < 1e-9) {
                // normal = <-1,0> => y = 0 boundary
                return BoundaryConditionType::NEUMANN;
            } else {
                // normal = <0,+-1> => x=0,1 boundary
                return BoundaryConditionType::SILVER_MULLER;
            }
        });
    }

    ElementInfos infos;
    std::shared_ptr<StructureDescription> structureDescription;
    PlaneWaveProblem<2> problem;
};

class Driver : public AbstractHarmonicSolverDriver<2> {
   public:
    Driver(PlaneWaveProblem<2>& problem)
        : nextCalled_(false), problem_(&problem){};

    bool stop() const override { return nextCalled_; }
    void nextProblem() override { nextCalled_ = true; }
    size_t getExpectedNumberOfProblems() const override { return 1; }
    const HarmonicProblem<2>& currentProblem() const override {
        return *problem_;
    }

    void handleResult(AbstractHarmonicResult<2>& result) override {
        if (plotter_ != nullptr) {
            result.writeVTK(*plotter_);
        }
        errors_.resize(3);
        errors_[0] = result.computeL2Error(*problem_);
        double actualInFlux = 0.0;
        double expectedInFlux = 0.0;
        double actualOutFlux = 0.0;
        double expectedOutFlux = 0.0;

        for (Base::Face* face : result.getMesh().getFacesList()) {
            if (face->isInternal()) {
                continue;
            }
            double flux = result.computeEnergyFlux(*face, Base::Side::LEFT,
                                                   problem_->omega());
            Geometry::PointReference<1> midPoint =
                face->getReferenceGeometry()->getCenter();
            LinearAlgebra::SmallVector<2> normal =
                face->getNormalVector(midPoint);
            normal /= normal.l2Norm();
            // Point flux at the middle of the face
            double expectedFlux =
                normal *
                problem_->localFlux(face->referenceToPhysical(midPoint));
            // Assume flat flux profile
            expectedFlux *= face->getDiameter();

            if (expectedFlux < 0) {
                actualInFlux += flux;
                expectedInFlux += expectedFlux;
            } else {
                actualOutFlux += flux;
                expectedOutFlux += expectedFlux;
            }
        }
        errors_[1] = actualOutFlux - expectedOutFlux;
        errors_[2] = actualInFlux - expectedInFlux;
    }

    bool hasChanged(HarmonicProblemChanges change) const override {
        return true;
    }

    void setPlotter(Output::VTKSpecificTimeWriter<2>* plotter) {
        plotter_ = plotter;
    }

    std::vector<double> getErrors() { return errors_; }

   private:
    bool nextCalled_;
    PlaneWaveProblem<2>* problem_;
    Output::VTKSpecificTimeWriter<2>* plotter_;

    std::vector<double> errors_;
};

// k_ is choosen such that:
//  1. Not too small, as that would only show the linear part
//  2. Not too large, as multiple period will only converge on very
//     finer meshes.
//  3. Not in the cardinal directions or along the diagonals, as those
//     are the normal directions in the mesh.
// E0 follows from the requirement that it is orthogonal.
const LinearAlgebra::SmallVector<2> ProblemData::k{1 * M_PI, 0.5 * M_PI};
const LinearAlgebra::SmallVectorC<2> ProblemData::E0{0.5, -1};

std::vector<double> solve(
    std::string meshFile, std::size_t level,
    std::shared_ptr<AbstractDiscretization<2>> discretization,
    std::string prefix) {

    ProblemData problemData;
    Base::ConfigurationData config(discretization->getNumberOfUnknowns());

    auto mesh =
        DGMax::readMesh<2>(meshFile, &config, *problemData.structureDescription,
                           discretization->getNumberOfElementMatrices());

    HarmonicSolver<2> solver(discretization);

    std::stringstream fileName;
    fileName << "solution-dgmax-" << level;
    Output::VTKSpecificTimeWriter<2> output(fileName.str(), mesh.get(), 0, 2);

    Driver driver(problemData.problem);
    driver.setPlotter(&output);

    solver.solve(*mesh, driver);
    return driver.getErrors();
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

    // Using second order Nedelec basis functions. So convergence rate for the
    // L2 norm is expected at h^{p} => reduction of a factor 4 with each level
    auto dgmax = std::make_shared<DGMaxDiscretization<2>>(2, 100);

    ConvergenceTestSet meshes = {getUnitSquareTriangleMeshes()};
    meshes.addExpectedErrors("L2 error", {
                                             2.32455311e-01,  //------
                                             7.78569137e-02,  //  2.99
                                             1.97411586e-02,  //  3.94
                                             4.94737345e-03,  //  3.99
                                             1.23692989e-03,  //  4.00
                                             3.09147735e-04,  //  4.00
                                             7.72701909e-05,  //  4.00)
                                         });
    meshes.addExpectedErrors("Outflux", {
                                            3.39100853e-01,  //------
                                            3.20113283e-01,  //  1.06
                                            1.02324552e-01,  //  3.13
                                            2.70451883e-02,  //  3.78
                                            6.85129516e-03,  //  3.95
                                            1.71800404e-03,  //  3.99
                                            4.29764608e-04,  //  4.00
                                        });
    meshes.addExpectedErrors("Influx", {
                                           -3.45303837e-01,  //------
                                           -3.15202590e-01,  //  1.10
                                           -1.01889058e-01,  //  3.09
                                           -2.70226520e-02,  //  3.77
                                           -6.85008794e-03,  //  3.94
                                           -1.71793569e-03,  //  3.99
                                           -4.29760562e-04,  //  4.00
                                       });
    runConvergenceTest(meshes, ignoreFailures,
                       [&dgmax](std::string meshFile, std::size_t order) {
                           return solve(meshFile, order, dgmax,
                                        "planewave-solution-dgmax");
                       });

    ConvergenceTestSet meshes2 = {getUnitSquareTriangleMeshes(0, 6)};
    meshes2.addExpectedErrors("L2 error", {
                                              2.50221220e-01,  //------
                                              7.84254413e-02,  //  3.19
                                              1.97608759e-02,  //  3.97
                                              4.94782994e-03,  //  3.99
                                              1.23693191e-03,  //  4.00
                                              3.09147016e-04,  //  4.00
                                          });
    meshes2.addExpectedErrors("Outflux", {
                                             3.34694680e-01,  //------
                                             3.16171119e-01,  //  1.06
                                             1.01976385e-01,  //  3.10
                                             2.70322081e-02,  //  3.77
                                             6.85283654e-03,  //  3.94
                                             1.71869399e-03,  //  3.99
                                         });
    meshes2.addExpectedErrors("Influx", {
                                            -3.17129614e-01,  //------
                                            -3.12017624e-01,  //  1.02
                                            -1.01503813e-01,  //  3.07
                                            -2.70028573e-02,  //  3.76
                                            -6.85106208e-03,  //  3.94
                                            -1.71858594e-03,  //  3.99
                                        });

    DivDGMaxDiscretizationBase::Stab stab;
    stab.stab1 = 105;
    stab.stab2 = 0;
    stab.stab3 = 10;
    stab.setAllFluxeTypes(DivDGMaxDiscretizationBase::FluxType::IP);
    auto divdgmax = std::make_shared<DivDGMaxDiscretization<2>>(2, stab);

    runConvergenceTest(meshes2, ignoreFailures,
                       [&divdgmax](std::string meshFile, std::size_t order) {
                           return solve(meshFile, order, divdgmax,
                                        "planewave-solution-divdgmax");
                       });
}
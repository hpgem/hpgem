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
#include <ProblemTypes/Harmonic/ExactFieldHarmonicProblem.h>
#include <Utils/PlaneWave.h>

#include <petsc.h>

using namespace DGMax;

using Vec2 = LinearAlgebra::SmallVectorC<2>;
using Point = Geometry::PointPhysical<2>;

struct ProblemData {
   public:
    // Different from 1 to increase the likelihood of detecting a bug
    constexpr static const Material material = Material(2.0, 1.2);
    // Not on the dispersion curve to require a source term
    // On dispersion would be k.l2Norm()/sqrt(epsilon)
    constexpr static const double omega = 1.5;
    constexpr static const double phase = 0.5;

    static const LinearAlgebra::SmallVector<2> k;
    static const LinearAlgebra::SmallVectorC<2> E0;

    ProblemData()
        : infos(material),
          structureDescription(StructureDescription::fromFunction(
              [this](const Base::Element*) { return &this->infos; })),
          problem(omega, std::make_shared<DGMax::PlaneWave<2>>(k, E0, phase)) {
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
    ExactFieldHarmonicProblem<2> problem;
};

class Driver : public AbstractHarmonicSolverDriver<2> {
   public:
    Driver(ExactFieldHarmonicProblem<2>& problem)
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
            // Get one of the element infos
            auto point = face->referenceToPhysical(midPoint);
            auto materialTensor =
                ElementInfos::get(*face->getPtrElementLeft())
                    .getMaterialConstantCurl(point, problem_->omega());
            // Point flux at the middle of the face
            double expectedFlux =
                normal * problem_->getField()->localFlux(point, materialTensor,
                                                         problem_->omega());
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
    ExactFieldHarmonicProblem<2>* problem_;
    Output::VTKSpecificTimeWriter<2>* plotter_;

    std::vector<double> errors_;
};

const Material ProblemData::material;
const double ProblemData::phase;

// k_ is chosen such that:
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
    fileName << prefix << level;
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
                                             2.31544492e-01,  //------
                                             7.78755325e-02,  //  2.97
                                             1.97460447e-02,  //  3.94
                                             4.94781906e-03,  //  3.99
                                             1.23696273e-03,  //  4.00
                                             3.09149950e-04,  //  4.00
                                             7.72703345e-05,  //  4.00
                                         });
    meshes.addExpectedErrors("Outflux", {
                                            2.86300153e-01,  //------
                                            2.66147184e-01,  //  1.08
                                            8.53078798e-02,  //  3.12
                                            2.25630286e-02,  //  3.78
                                            5.71747856e-03,  //  3.95
                                            1.43387188e-03,  //  3.99
                                            3.58708541e-04,  //  4.00
                                        });
    meshes.addExpectedErrors("Influx", {
                                           -2.75684494e-01,  //------
                                           -2.61542159e-01,  //  1.05
                                           -8.48998638e-02,  //  3.08
                                           -2.25414966e-02,  //  3.77
                                           -5.71630571e-03,  //  3.94
                                           -1.43380470e-03,  //  3.99
                                           -3.58704496e-04,  //  4.00
                                       });
    runConvergenceTest(meshes, ignoreFailures,
                       [&dgmax](std::string meshFile, std::size_t order) {
                           return solve(meshFile, order, dgmax,
                                        "planewave-solution-dgmax");
                       });

    ConvergenceTestSet meshes2 = {getUnitSquareTriangleMeshes(0, 6)};
    meshes2.addExpectedErrors("L2 error", {
                                              2.87642577e-01,  //------
                                              7.96248160e-02,  //  3.61
                                              1.97976622e-02,  //  4.02
                                              4.94929023e-03,  //  4.00
                                              1.23700845e-03,  //  4.00
                                              3.09151533e-04,  //  4.00
                                          });
    meshes2.addExpectedErrors("Outflux", {
                                             2.52041706e-01,  //------
                                             2.57375069e-01,  //  0.98
                                             8.45657740e-02,  //  3.04
                                             2.25232154e-02,  //  3.75
                                             5.71665736e-03,  //  3.94
                                             1.43422355e-03,  //  3.99
                                         });
    meshes2.addExpectedErrors("Influx", {
                                            -1.69073609e-01,  //------
                                            -2.54790852e-01,  //  0.66
                                            -8.40758550e-02,  //  3.03
                                            -2.24894193e-02,  //  3.74
                                            -5.71456834e-03,  //  3.94
                                            -1.43409513e-03,  //  3.98
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
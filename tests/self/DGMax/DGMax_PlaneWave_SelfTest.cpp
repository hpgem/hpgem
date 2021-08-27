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
#include <Utils/PredefinedStructure.h>

#include <petsc.h>

using namespace DGMax;

using Vec2 = LinearAlgebra::SmallVector<2>;
using Point = Geometry::PointPhysical<2>;

class TestingProblem : public ExactHarmonicProblem<2> {
   public:
    TestingProblem()
        // k_ is choosen such that:
        //  1. Not too small, as that would only show the linear part
        //  2. Not too large, as multiple period will only converge on very
        //     finer meshes.
        //  3. Not in the cardinal directions or along the diagonals, as those
        //     are the normal directions in the mesh.
        : k_({1 * M_PI, 0.5 * M_PI}), E0_({0.5, -1}), bface(false){};

    double omega() const final { return 3.0; }

    double modulation(const Point& p) const {
        return std::cos(p.getCoordinates() * k_);
    }

    Vec2 sourceTerm(const Point& p) const final {
        double omega2 = omega();
        omega2 *= omega2;
        return E0_ * modulation(p) * (k_.l2NormSquared() - omega2);
    }
    Vec2 exactSolution(const Point& p) const final {
        return E0_ * modulation(p);
    }
    Vec2 exactSolutionCurl(const Point& p) const final {
        return -k_.crossProduct(E0_) * std::sin(p.getCoordinates() * k_);
    }

    BoundaryConditionType getBoundaryConditionType(
        const Base::Face& face) const final {

        bface.setFace(&face);
        bface.setPointReference(face.getReferenceGeometry()->getCenter());
        Vec2 normal = bface.getUnitNormalVector();
        double nx = std::abs(normal[0]);
        if (std::abs(nx - 1.0) < 1e-8) {
            return BoundaryConditionType::NEUMANN;
        } else {
            return BoundaryConditionType::NEUMANN;
        }
    }

   private:
    Vec2 k_;
    /// Direction of the field
    Vec2 E0_;
    /// Used for computing the boundary condition type
    mutable Base::PhysicalFace<2> bface;
};

double solve(std::string meshFile, std::size_t level) {
    Base::ConfigurationData config(1);

    PredefinedStructureDescription structure =
        PredefinedStructureDescription(PredefinedStructure::VACUUM, 2);

    auto mesh = DGMax::readMesh<2>(meshFile, &config, structure, 2);
    DGMaxHarmonic<2> solver(*mesh, 100, 2);

    TestingProblem problem;

    solver.solve(problem);
    std::stringstream fileName;
    fileName << "solution-" << level;
    Output::VTKSpecificTimeWriter<2> output(fileName.str(), mesh.get(), 0, 2);
    solver.writeVTK(output);
    auto errors = solver.computeError({DGMaxDiscretizationBase::L2}, problem);
    return errors[DGMaxDiscretizationBase::L2];
}

int main(int argc, char** argv) {

    Base::parse_options(argc, argv);
    initDGMaxLogging();

    bool ignoreFailures = false;

    // Default the solver if not specified to a direct LU solver
    std::map<std::string, std::string> defaultOptions = {
        {"-ksp_type", "preonly"}, {"-pc_type", "lu"}};
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
    ConvergenceTestSet meshes = {getUnitSquareTriangleMeshes(),
                                 {
                                     1.42329582e-01,  //------
                                     5.25619738e-02,  //  2.71
                                     1.34616005e-02,  //  3.90
                                     3.43157573e-03,  //  3.92
                                     8.66123561e-04,  //  3.96
                                     2.17527652e-04,  //  3.98
                                     5.45037524e-05,  //  3.99
                                 }};
    runConvergenceTest(meshes, ignoreFailures, &solve);
}
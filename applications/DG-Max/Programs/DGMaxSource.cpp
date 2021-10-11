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

#include <chrono>
#include <exception>
#include <Base/CommandLineOptions.h>
#include <Base/MeshFileInformation.h>
#include <Output/VTKSpecificTimeWriter.h>

#include <DGMaxLogger.h>
#include <DGMaxProgramUtils.h>
#include <Algorithms/DGMaxHarmonic.h>
#include <Algorithms/DivDGMaxHarmonic.h>
#include <ProblemTypes/Harmonic/SampleHarmonicProblems.h>

auto& meshFileName = Base::register_argument<std::string>(
    'm', "meshFile", "The hpgem meshfile to use", true);

auto& structureString = Base::register_argument<std::string>(
    '\0', "structure", "Structure to use", false, "0");

auto& order = Base::register_argument<std::size_t>(
    'p', "order", "Polynomial order of the solution", true);

auto& method = Base::register_argument<std::string>(
    '\0', "method", "Method to use, either DivDGMax (default) or DGMax", false,
    "DivDGMax");

template <std::size_t dim>
void runWithDimension();

int main(int argc, char** argv) {
    using namespace hpgem;
    using namespace DGMax;

    Base::parse_options(argc, argv);
    initDGMaxLogging();
    printArguments(argc, argv);

    auto start = std::chrono::high_resolution_clock::now();

    Base::MeshFileInformation meshInfo =
        Base::MeshFileInformation::readInformation(meshFileName.getValue());
    try {
        switch (meshInfo.dimension) {
            case 2:
                runWithDimension<2>();
                break;
            case 3:
                runWithDimension<3>();
                break;
            default:
                logger.assert_always(false, "Only dimension 2 & 3 supported");
        }
    } catch (const std::exception& e) {
        DGMaxLogger(ERROR, e.what());
        return 1;
    } catch (const char* message) {
        DGMaxLogger(ERROR, message);
        return 1;
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedTime = end - start;
    DGMaxLogger(INFO, "Runtime %s", elapsedTime.count());
    return 0;
}

template <std::size_t dim>
void writeMesh(std::string, const Base::MeshManipulator<dim>* mesh);

template <std::size_t dim>
class TestingProblem : public HarmonicProblem<dim> {

   public:
    TestingProblem() : pface(false){};

    double omega() const final { return 4.0e-2; }

    LinearAlgebra::SmallVector<dim> sourceTerm(
        const Geometry::PointPhysical<dim>& point) const final {
        LinearAlgebra::SmallVector<dim> result;

        result.set(0.0);
        if (point[1] > 0 && point[1] < 50) result[0] = 1.0;

        return result;
    }

    LinearAlgebra::SmallVector<dim> boundaryCondition(
        Base::PhysicalFace<dim>& face) const override {
        return {};
    }

    DGMax::BoundaryConditionType getBoundaryConditionType(
        const Base::Face& face) const final {
        pface.setFace(&face);
        pface.setPointReference(face.getReferenceGeometry()->getCenter());
        if (face.isInternal()) {
            return DGMax::BoundaryConditionType::INTERNAL;
        }

        LinearAlgebra::SmallVector<dim> normal = pface.getUnitNormalVector();

        double nx = std::abs(std::abs(normal[0]) - 1.0);
        if (nx < 1e-8) {
            return DGMax::BoundaryConditionType::DIRICHLET;
        } else {
            return DGMax::BoundaryConditionType::NEUMANN;
        }
    }

   private:
    mutable Base::PhysicalFace<dim> pface;
};

template <std::size_t dim>
void runWithDimension() {
    bool divdgmax;
    std::size_t unknowns;
    if (method.getValue() == "DivDGMax") {
        divdgmax = true;
        unknowns = 2;
    } else if (method.getValue() == "DGMax") {
        divdgmax = false;
        unknowns = 1;
    } else {
        logger(ERROR, "Unkown method %", method.getValue());
        return;
    }

    std::size_t numberOfElementMatrices = 2;

    Base::ConfigurationData configData(unknowns, 1);
    auto structure =
        DGMax::determineStructureDescription(structureString.getValue(), dim);

    auto mesh = DGMax::readMesh<dim>(meshFileName.getValue(), &configData,
                                     *structure, numberOfElementMatrices);
    logger(INFO, "Loaded mesh % with % local elements", meshFileName.getValue(),
           mesh->getNumberOfElements());
    writeMesh<dim>("mesh", mesh.get());

    std::unique_ptr<DGMax::AbstractHarmonicSolver<dim>> solver;

    if (divdgmax) {
        // Placeholder for more complicate fluxes
        DivDGMaxDiscretizationBase::Stab stab;
        stab.stab1 = 5;
        stab.stab2 = 0;
        stab.stab3 = 5;
        stab.setAllFluxeTypes(DivDGMaxDiscretization<dim>::FluxType::BREZZI);

        solver = std::make_unique<DivDGMaxHarmonic<dim>>(*mesh, stab,
                                                         order.getValue());
    } else {
        double stab = 100;
        solver =
            std::make_unique<DGMaxHarmonic<dim>>(*mesh, stab, order.getValue());
    }

    // Problem definition
    std::unique_ptr<HarmonicProblem<dim>> problem;
    problem = std::make_unique<SampleHarmonicProblems<dim>>(
        SampleHarmonicProblems<dim>::Problem::CONSTANT, 1);
    problem = std::make_unique<TestingProblem<dim>>();

    solver->solve(*problem);

    // Output
    Output::VTKSpecificTimeWriter<dim> output("harmonic", mesh.get(), 0,
                                              order.getValue());
    output.write(
        [](Base::Element* element, const Geometry::PointReference<dim>&,
           std::size_t) {
            const ElementInfos* elementInfos =
                dynamic_cast<ElementInfos*>(element->getUserData());
            logger.assert_debug(elementInfos != nullptr,
                                "Incorrect user data type");
            return elementInfos->epsilon_;
        },
        "epsilon");
    output.write(
        [&problem](Base::Element* element,
                   const Geometry::PointReference<dim>& p, std::size_t) {
            return problem->sourceTerm(element->referenceToPhysical(p));
        },
        "source");
    solver->writeVTK(output);
}

template <std::size_t DIM>
void writeMesh(std::string fileName, const Base::MeshManipulator<DIM>* mesh) {
    Output::VTKSpecificTimeWriter<DIM> writer(fileName, mesh);
    writer.write(
        [](Base::Element* element, const Geometry::PointReference<DIM>&,
           std::size_t) {
            const ElementInfos* elementInfos =
                dynamic_cast<ElementInfos*>(element->getUserData());
            logger.assert_debug(elementInfos != nullptr,
                                "Incorrect user data type");
            return elementInfos->epsilon_;
        },
        "epsilon");
}
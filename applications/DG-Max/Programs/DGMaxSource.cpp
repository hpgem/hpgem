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
#include <Algorithms/HarmonicSolver.h>
#include <Algorithms/DGMaxDiscretization.h>
#include <Algorithms/DivDGMaxDiscretization.h>


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
    TestingProblem(double omega) : omega_(omega){};

    double omega() const final { return omega_; }

    LinearAlgebra::SmallVectorC<dim> sourceTerm(
        const Geometry::PointPhysical<dim>& point) const final {
        LinearAlgebra::SmallVector<dim> result;

        result.set(0.0);
        // if (point[1] > 0 && point[1] < 50) result[0] = 1.0;

        return result;
    }

    LinearAlgebra::SmallVectorC<dim> boundaryCondition(
        Base::PhysicalFace<dim>& face) const override {

        auto bct = getBoundaryConditionType(*face.getFace());

        if (bct == DGMax::BoundaryConditionType::DIRICHLET) {
            return {};
        }
        LinearAlgebra::SmallVectorC<dim> normal = face.getUnitNormalVector();

        // E = E_0 exp(ik.x)
        // g_N = Curl E + i kappa sqrt(epsilon) E x n
        //     = i[(k x E_0) + kappa sqrt(epsilon) (E_0 x n)] exp(ik.x)

        using namespace std::complex_literals;
        LinearAlgebra::SmallVectorC<dim> E0{1.0, 0}, k{0, omega()};

        auto* userData = dynamic_cast<ElementInfos*>(
            face.getFace()->getPtrElementLeft()->getUserData());
        logger.assert_always(userData != nullptr, "Incorrect data");

        k[1] *= std::sqrt(userData->epsilon_);

        // (kxE_0) + kappa * sqrt(epsilon) (E_0 x n)
        LinearAlgebra::SmallVectorC<dim> result = k.crossProduct(E0);
        result +=
            omega() * std::sqrt(userData->epsilon_) * E0.crossProduct(normal);

        // i exp(ik.x)
        result *=
            1i * std::exp(1i * face.getPointPhysical().getCoordinates() * k);
        return result;
    }

    DGMax::BoundaryConditionType getBoundaryConditionType(
        const Base::Face& face) const final {

        if (face.isInternal()) {
            return DGMax::BoundaryConditionType::INTERNAL;
        }

        LinearAlgebra::SmallVector<dim> normal = face.getNormalVector(
            face.getReferenceGeometry()->getCenter().castDimension<dim - 1>());
        normal /= normal.l2Norm();

        // Check if the normal matches: {1,0} or {-1,0}
        double nx = std::abs(std::abs(normal[0]) - 1.0);
        if (nx < 1e-8) {
            // For faces at x=0,1
            return DGMax::BoundaryConditionType::DIRICHLET;
        } else {
            // For faces at y=0,1
            return DGMax::BoundaryConditionType::SILVER_MULLER;
        }
    }

   private:
    double omega_;
};

template <std::size_t dim>
class Driver : public DGMax::AbstractHarmonicSolverDriver<dim> {
   public:
    Driver(Base::MeshManipulator<dim>& mesh)
        : mesh_(&mesh), problem_(), nextCalled_(0){};

    bool stop() const override { return nextCalled_ == 5; }
    void nextProblem() override {
        nextCalled_++;
        problem_ = std::make_shared<TestingProblem<dim>>(0.01 * nextCalled_);
    }

    const HarmonicProblem<dim>& currentProblem() const override {
        return *problem_;
    }

    bool hasChanged(HarmonicProblemChanges change) const override {
        if (nextCalled_ == 1) {
            return true;
        } else {
            switch (change) {
                case HarmonicProblemChanges::OMEGA:
                case HarmonicProblemChanges::BOUNDARY_CONDITION_VALUE:
                    return true;
                case HarmonicProblemChanges::BOUNDARY_CONDITION_TYPE:
                case HarmonicProblemChanges::CURRENT_SOURCE:
                    return false;
                default:
                    return true;
            }
        }
    }

    void handleResult(DGMax::AbstractHarmonicResult<dim>& result) override {
        std::stringstream outputFile;
        outputFile << "harmonic-" << nextCalled_;
        Output::VTKSpecificTimeWriter<dim> output(outputFile.str(), mesh_, 0,
                                                  order.getValue());
        result.writeVTK(output);
        // Additional output
        output.write(
            [this](Base::Element* element,
                   const Geometry::PointReference<dim>& p, std::size_t) {
                return problem_->sourceTerm(element->referenceToPhysical(p))
                    .real();
            },
            "source");

        auto* eproblem =
            dynamic_cast<ExactHarmonicProblem<dim>*>(problem_.get());
        if (eproblem != nullptr) {
            double l2Error = result.computeL2Error(*eproblem);
            DGMaxLogger(INFO, "L2 error %", l2Error);

            output.write(
                [&eproblem](Base::Element* element,
                            const Geometry::PointReference<dim>& p,
                            std::size_t) {
                    return eproblem
                        ->exactSolution(element->referenceToPhysical(p))
                        .real();
                },
                "SolutionReal");
            output.write(
                [&eproblem](Base::Element* element,
                            const Geometry::PointReference<dim>& p,
                            std::size_t) {
                    return eproblem
                        ->exactSolution(element->referenceToPhysical(p))
                        .imag();
                },
                "SolutionImag");
        }
    }

   private:
    Base::MeshManipulator<dim>* mesh_;
    std::shared_ptr<HarmonicProblem<dim>> problem_;
    int nextCalled_;
};

template <std::size_t dim>
void runWithDimension() {
    std::shared_ptr<DGMax::AbstractDiscretization<dim>> discretization;

    if (method.getValue() == "DivDGMax") {
        // Placeholder for more complicate fluxes
        DivDGMaxDiscretizationBase::Stab stab;
        stab.stab1 = 105;
        stab.stab2 = 0;
        stab.stab3 = 500;
        stab.setAllFluxeTypes(DivDGMaxDiscretization<dim>::FluxType::IP);
        discretization = std::make_shared<DivDGMaxDiscretization<dim>>(
            order.getValue(), stab);
    } else if (method.getValue() == "DGMax") {
        double stab = 100;
        discretization =
            std::make_shared<DGMaxDiscretization<dim>>(order.getValue(), stab);
    } else {
        logger(ERROR, "Unkown method %", method.getValue());
        return;
    }

    Base::ConfigurationData configData(discretization->getNumberOfUnknowns(),
                                       1);
    auto structure =
        DGMax::determineStructureDescription(structureString.getValue(), dim);

    auto mesh =
        DGMax::readMesh<dim>(meshFileName.getValue(), &configData, *structure,
                             discretization->getNumberOfElementMatrices());
    logger(INFO, "Loaded mesh % with % local elements", meshFileName.getValue(),
           mesh->getNumberOfElements());
    writeMesh<dim>("mesh", mesh.get());

    DGMax::HarmonicSolver<dim> solver(discretization);
    Driver<dim> driver(*mesh);
    solver.solve(*mesh, driver);
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
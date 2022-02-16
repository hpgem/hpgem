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
#include <iomanip>
#include <Base/CommandLineOptions.h>
#include <Base/MeshFileInformation.h>
#include <Output/VTKSpecificTimeWriter.h>

#include <DGMaxLogger.h>
#include <DGMaxProgramUtils.h>
#include <Algorithms/HarmonicSolver.h>
#include <Algorithms/DGMaxDiscretization.h>
#include <Algorithms/DivDGMaxDiscretization.h>
#include <Utils/PlaneWave.h>
#include <ProblemTypes/Harmonic/ScatteringProblem.h>
#include <ProblemTypes/Harmonic/InterfaceReflectionField.h>

#include <PMLElementInfos.h>

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

/**
 * Facets of the mesh through which the flux needs to be computed
 */
class FluxFacets {

   public:
    FluxFacets(const Base::MeshManipulatorBase& mesh);

    template <std::size_t dim>
    std::map<std::string, LinearAlgebra::SmallVector<4>> computeFluxes(
        DGMax::AbstractHarmonicResult<dim>& result, double wavenumber,
        const DGMax::FieldPattern<dim>* background) const;

    // Public to allow writing the header
    /**
     * Face for which an energy flux needs to be computed
     */
    struct FluxFace {
        /**
         * The face
         */
        Base::Face* face;
        /**
         * The side from which it needs to be computed
         */
        Base::Side side;
    };
    std::map<std::string, std::vector<FluxFace>> facets;
};

template <std::size_t dim>
class TestingProblem : public HarmonicProblem<dim> {

   public:
    TestingProblem(double omega) : omega_(omega){};

    double omega() const final { return omega_; }

    LinearAlgebra::SmallVectorC<dim> sourceTerm(
        const Base::Element&,
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

        auto* userData = dynamic_cast<ElementInfos*>(
            face.getFace()->getPtrElementLeft()->getUserData());
        logger.assert_always(userData != nullptr, "Incorrect data");

        DGMax::PlaneWave<dim> incidentWave =
            DGMax::PlaneWave<dim>::onDispersion(omega_, {0.0, 1.0}, {1.0, 0.0},
                                                userData->getMaterial(), 0.0);

        const auto& p = face.getPointPhysical();

        using namespace std::complex_literals;
        // Curl E/mu + i omega Z (E x n)
        return incidentWave.fieldCurl(p) / userData->getPermeability() +
               1i * omega_ * userData->getImpedance() *
                   incidentWave.field(p).crossProduct(normal);
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
    static const constexpr std::size_t NUMBER_OF_PROBLEMS = 301;

    // For plotting in combination with scattered fields
    struct Fields {
        LinearAlgebra::SmallVectorC<dim> field;
        LinearAlgebra::SmallVectorC<dim> incident;
    };

   public:
    Driver(Base::MeshManipulator<dim>& mesh)
        : mesh_(&mesh), problem_(), nextCalled_(0), fluxFacets_(mesh) {

        if (Base::MPIContainer::Instance().getProcessorID() == 0) {
            outfile.open("harmonic.csv");
            logger.assert_always(outfile.good(), "Opening output file failed");
            outfile << "wavenumber";
            for (const auto& facet : fluxFacets_.facets) {
                outfile << "," << facet.first;
                if (true) {
                    outfile << "," << facet.first << "-incident";
                    outfile << "," << facet.first << "-extinction";
                    outfile << "," << facet.first << "-total";
                }
            }
            outfile << std::endl;
        }
    };

    bool stop() const override { return nextCalled_ == NUMBER_OF_PROBLEMS; }

    size_t getExpectedNumberOfProblems() const override {
        return NUMBER_OF_PROBLEMS;
    }

    void nextProblem() override {
        nextCalled_++;
        problem_ =
            std::make_shared<TestingProblem<dim>>(0.01 + 0.0001 * nextCalled_);
        //        double omega = 2 * M_PI * (0.001 + 0.001 * (nextCalled_ - 1));

        // Inverse cm wavenumbers
        const double wavenumberMin = 2000;
        const double waveNumberMax = 8000;
        const double dWavenumber = waveNumberMax - wavenumberMin;
        // Conversion factor (mesh length-unit) = N cm
        const double meshConversion = 1e-6;

        // Linear spacing in wavenumber
        double omega = 2.0 * M_PI * meshConversion *
                       (wavenumberMin + dWavenumber * (nextCalled_ - 1) /
                                            (NUMBER_OF_PROBLEMS + 1));

        LinearAlgebra::SmallVector<dim> k;
        LinearAlgebra::SmallVectorC<dim> E0;
        k[1] = omega;  // Assume vacuum
        E0[0] = 1.0;
        DGMax::Material outside(1.0, 1.0);
        DGMax::Material inside(12.1, 1.0);
        problem_ = std::make_shared<DGMax::ScatteringProblem<dim>>(
            omega, std::make_shared<DGMax::InterfaceReflectionField<dim>>(
                       omega, 0.0, outside, inside, 0.0));
        //        problem_ = std::make_shared<DGMax::ScatteringProblem<dim>>(
        //            omega, std::make_shared<DGMax::PlaneWave<dim>>(k, E0,
        //            0.0));
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
                case HarmonicProblemChanges::CURRENT_SOURCE:
                    return true;
                case HarmonicProblemChanges::BOUNDARY_CONDITION_TYPE:
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
                return problem_
                    ->sourceTerm(*element, element->referenceToPhysical(p))
                    .real();
            },
            "source");

        plotResult(result, output);

        const DGMax::FieldPattern<dim>* background = nullptr;
        auto scatteringProblem =
            std::dynamic_pointer_cast<DGMax::ScatteringProblem<dim>>(problem_);
        if (scatteringProblem) {
            background = scatteringProblem->incidentFieldPattern().get();
        }

        // Fluxes
        std::map<std::string, LinearAlgebra::SmallVector<4>> fluxes =
            fluxFacets_.computeFluxes(result, problem_->omega(), background);

        outfile << std::setprecision(16) << problem_->omega();
        LinearAlgebra::SmallVector<4> totalFlux = {};
        for (const auto& entry : fluxes) {
            std::size_t nfluxes = problem_->isScatterFieldProblem() ? 4 : 1;
            for (std::size_t i = 0; i < nfluxes; ++i) {
                outfile << "," << std::setprecision(16) << entry.second[i];
            }

            // TODO: Change to only boundary
            totalFlux += entry.second;
        }
        outfile << std::endl;
        DGMaxLogger(INFO, "Total net flux %", totalFlux);
    }

   private:
    void plotResult(DGMax::AbstractHarmonicResult<dim>& result,
                    Output::VTKSpecificTimeWriter<dim>& output);

    Base::MeshManipulator<dim>* mesh_;
    std::shared_ptr<HarmonicProblem<dim>> problem_;
    std::ofstream outfile;
    int nextCalled_;
    FluxFacets fluxFacets_;
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

    std::vector<PMLElementInfos<dim>> pmls;
    {
        double strength = 1e-5;
        LinearAlgebra::SmallVector<dim> offset, direction;
        offset[1] = -300.0;
        direction[1] = -1.0;

        pmls.emplace_back(PMLElementInfos<dim>(DGMax::Material(), offset,
                                               direction, strength / 10));
        offset[1] = 719.0;
        direction[1] = 1.0;
        DGMax::Material silicon(12.1, 1);
        pmls.emplace_back(
            PMLElementInfos<dim>(silicon, offset, direction, strength));
    }
    for (Base::Element* element : mesh->getElementsList()) {
        const std::string& zoneName = element->getZone().getName();
        if (zoneName == "4") {
            element->setUserData(&pmls[0]);
        } else if (zoneName == "5") {
            element->setUserData(&pmls[1]);
        }
    }

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
            return elementInfos->getPermittivity();
        },
        "epsilon");
}

FluxFacets::FluxFacets(const Base::MeshManipulatorBase& mesh) {
    // When using more than 1 MPI rank, we would need to aggregate data. But the
    // first processor may not know all the regions and possible boundaries.
    logger.assert_always(Base::MPIContainer::Instance().getNumProcessors() == 1,
                         "FluxFacets not suitable for more than 1 MPI rank");
    for (Base::Face* face : mesh.getFacesList()) {
        using Base::Side;
        if (!face->isOwnedByCurrentProcessor()) {
            continue;
        }
        const Base::Zone& leftZone = face->getPtrElementLeft()->getZone();
        FluxFace fface;
        fface.face = face;
        std::string facetName;

        if (face->isInternal()) {
            const Base::Zone& rightZone = face->getPtrElementRight()->getZone();
            if (leftZone.getZoneId() == rightZone.getZoneId()) {
                // Only interested in zone boundaries
                continue;
            }
            // Stable direction
            fface.side = leftZone.getZoneId() < rightZone.getZoneId()
                             ? Side::LEFT
                             : Side::RIGHT;
            // Name
            std::stringstream name;
            name << "flux-"
                 << (fface.side == Side::LEFT ? leftZone : rightZone).getName()
                 << "--"
                 << (fface.side == Side::LEFT ? rightZone : leftZone).getName();
            facetName = name.str();
        } else {
            // Outward fluxes are always included
            fface.side = Side::LEFT;
            std::stringstream name;
            name << "flux-" << leftZone.getName();
            facetName = name.str();
        }
        facets[facetName].push_back(fface);
    }
}

template <std::size_t dim>
std::map<std::string, LinearAlgebra::SmallVector<4>> FluxFacets::computeFluxes(
    DGMax::AbstractHarmonicResult<dim>& result, double wavenumber,
    const DGMax::FieldPattern<dim>* background) const {
    std::map<std::string, LinearAlgebra::SmallVector<4>> fluxes;
    for (const auto& facet : facets) {
        LinearAlgebra::SmallVector<4> flux = {};
        for (const auto& face : facet.second) {
            flux += result.computeEnergyFlux(*face.face, face.side, wavenumber,
                                             background);
        }
        fluxes[facet.first] += flux;
    }
    // TODO: MPI
    return fluxes;
}

template <std::size_t dim>
void Driver<dim>::plotResult(DGMax::AbstractHarmonicResult<dim>& result,
                             Output::VTKSpecificTimeWriter<dim>& output) {
    // Plotting of extra quantities for various types of problems

    auto scatteringProblem =
        std::dynamic_pointer_cast<DGMax::ScatteringProblem<dim>>(problem_);
    if (scatteringProblem) {
        using VecR = LinearAlgebra::SmallVector<dim>;
        std::map<std::string, std::function<double(Fields&)>> scalars;
        std::map<std::string, std::function<VecR(Fields&)>> vectors;

        scalars["Eincident-mag"] = [](Fields& fields) {
            return fields.incident.l2Norm();
        };
        vectors["Eincident-real"] = [](Fields& fields) {
            return fields.incident.real();
        };
        vectors["Eincident-imag"] = [](Fields& fields) {
            return fields.incident.imag();
        };

        scalars["Etotal-mag"] = [](Fields& fields) {
            return (fields.incident + fields.field).l2Norm();
        };

        vectors["Etotal-real"] = [](Fields& fields) {
            return fields.incident.real() + fields.field.real();
        };
        vectors["Etotal-imag"] = [](Fields& fields) {
            return fields.incident.imag() + fields.field.imag();
        };

        output.template writeMultiple<Fields>(
            [&result, &scatteringProblem](
                Base::Element* element, const Geometry::PointReference<dim>& p,
                std::size_t) {
                Fields fields;
                fields.field = result.computeField(element, p);
                auto pphys = element->template referenceToPhysical(p);
                fields.incident = scatteringProblem->incidentField(pphys);
                return fields;
            },
            scalars, vectors);
    }

    auto* exactProblem =
        dynamic_cast<ExactHarmonicProblem<dim>*>(problem_.get());
    if (exactProblem != nullptr) {
        double l2Error = result.computeL2Error(*exactProblem);
        DGMaxLogger(INFO, "L2 error %", l2Error);

        output.write(
            [&exactProblem](Base::Element* element,
                            const Geometry::PointReference<dim>& p,
                            std::size_t) {
                return exactProblem
                    ->exactSolution(element->referenceToPhysical(p))
                    .real();
            },
            "SolutionReal");
        output.write(
            [&exactProblem](Base::Element* element,
                            const Geometry::PointReference<dim>& p,
                            std::size_t) {
                return exactProblem
                    ->exactSolution(element->referenceToPhysical(p))
                    .imag();
            },
            "SolutionImag");
    }
}
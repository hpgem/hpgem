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
#include <Utils/FluxFacets.h>
#include <ProblemTypes/Harmonic/ScatteringProblem.h>
#include <ProblemTypes/Harmonic/InterfaceReflectionField.h>

#include <PMLElementInfos.h>
#include <Utils/PMLTransmission.h>

auto& meshFileName = Base::register_argument<std::string>(
    'm', "meshFile", "The hpgem meshfile to use", true);

auto& structureString = Base::register_argument<std::string>(
    '\0', "structure", "Structure to use", false, "0");

auto& order = Base::register_argument<std::size_t>(
    'p', "order", "Polynomial order of the solution", true);

auto& method = Base::register_argument<std::string>(
    '\0', "method", "Method to use, either DivDGMax (default) or DGMax", false,
    "DivDGMax");

auto& wavenumberscm = Base::register_argument<std::string>(
    '\0', "wavenumberscm",
    "Wavenumbers to compute at in cm^{-1}.\n"
    "Either listed as w1,w2, ..,wn or as N equidistant values @N,wmin,wmax",
    true);
auto& meshsize = Base::register_argument<double>(
    '\0', "scalenm",
    "Length scale of the mesh in nm (i.e. 1 mesh unit is x nanometers)");

auto& pmlFile = Base::register_argument<std::string>(
    '\0', "pmlfile", "File with description of the PMLs", false);
auto& pmlString = Base::register_argument<std::string>(
    '\0', "pmls",
    "Inline description of the PMLs. Same format as a PML file but using ';' "
    "instead of newlines",
    false);

auto& fieldIndices = Base::register_argument<std::string>(
    '\0', "fields",
    "Computational indices at which to output the fields, comma separated "
    "(default = all)",
    false);

auto& computePMLTransmission = Base::register_argument<bool>(
    '\0', "pmltransmission", "Compute L2 intregral at far end of PML", false);

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

/// Parses the wavenumbers passed to the program.
///
/// \return A list of wavenumbers as pair. First the computational value, second
/// the input value (e.g. in cm^{-1})
std::vector<std::pair<double, double>> parseWaveNumbers();

/// Apply PMLs from command line arguments to the mesh
///
/// \tparam dim The dimension of the mesh
/// \param mesh The mesh
/// \return The PMLElementInfos for the PMLs (and referenced to by the Elements
/// in the mesh).
template <std::size_t dim>
std::vector<std::shared_ptr<PMLElementInfos<dim>>> applyPMLs(
    Base::MeshManipulator<dim>& mesh);

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
    // For plotting in combination with scattered fields
    struct Fields {
        LinearAlgebra::SmallVectorC<dim> field;
        LinearAlgebra::SmallVectorC<dim> incident;
        LinearAlgebra::SmallVectorC<dim> totalCurl;
        DGMax::Material material;
    };

   public:
    Driver(Base::MeshManipulator<dim>& mesh)
        : mesh_(&mesh),
          wavenumbers_(parseWaveNumbers()),
          problem_(),
          nextCalled_(-1),
          fluxFacets_(mesh),
          pmlTransmission_(mesh) {

        if (Base::MPIContainer::Instance().getProcessorID() == 0) {
            outfile.open("harmonic.csv");
            logger.assert_always(outfile.good(), "Opening output file failed");
            outfile << "wavenumber-computational,wavenumber";
            for (const auto& facetName : fluxFacets_.getFacetNames()) {
                outfile << "," << facetName;
                if (true) {
                    outfile << "," << facetName << "-incident";
                    outfile << "," << facetName << "-extinction";
                    outfile << "," << facetName << "-total";
                }
            }
            if (computePMLTransmission.isUsed()) {
                for (const auto& facetName : pmlTransmission_.getFacetNames()) {
                    outfile << "," << facetName << "-pmltransmission";
                }
            }
            outfile << std::endl;
        }
        if (fieldIndices.isUsed()) {
            auto indexStrings =
                DGMax::stringSplit(fieldIndices.getValue(), ',');
            for (const auto& indexString : indexStrings) {
                std::size_t pos;
                int index = std::stoi(indexString, &pos);
                DGMaxLogger.assert_always(pos == indexString.size(),
                                          "Could not parse index string '%'",
                                          indexString);
                fieldPrintIndices_.push_back(index);
            }
        }
    };

    bool stop() const override {
        return nextCalled_ == wavenumbers_.size() - 1;
    }

    size_t getExpectedNumberOfProblems() const override {
        return wavenumbers_.size();
    }

    void nextProblem() override {
        nextCalled_++;
        //        problem_ =
        //            std::make_shared<TestingProblem<dim>>(0.01 + 0.0001 *
        //            nextCalled_);
        //        double omega = 2 * M_PI * (0.001 + 0.001 * (nextCalled_ - 1));

        double omega = wavenumbers_[nextCalled_].first;
        DGMaxLogger(INFO,
                    "Current wavenumber: (computational %, real % cm^{-1})",
                    omega, wavenumbers_[nextCalled_].second);

        LinearAlgebra::SmallVector<dim> k;
        LinearAlgebra::SmallVectorC<dim> E0;
        k[1] = omega;  // Assume vacuum
        E0[0] = 1.0;
        DGMax::Material outside(1.0, 1.0);
        DGMax::Material inside(12.1, 1.0);
        problem_ = std::make_shared<DGMax::ScatteringProblem<dim>>(
            omega, std::make_shared<DGMax::InterfaceReflectionField<dim>>(
                       omega, 0.0, outside, inside, 0.0));
    }

    const HarmonicProblem<dim>& currentProblem() const override {
        return *problem_;
    }

    bool hasChanged(HarmonicProblemChanges change) const override {
        if (nextCalled_ == 0) {
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

        bool inPrintIndices =
            std::find(fieldPrintIndices_.begin(), fieldPrintIndices_.end(),
                      nextCalled_) != fieldPrintIndices_.end();
        if (inPrintIndices || !fieldIndices.isUsed()) {
            Output::VTKSpecificTimeWriter<dim> output(outputFile.str(), mesh_,
                                                      0, order.getValue());
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
        }

        const DGMax::FieldPattern<dim>* background = nullptr;
        auto scatteringProblem =
            std::dynamic_pointer_cast<DGMax::ScatteringProblem<dim>>(problem_);
        if (scatteringProblem) {
            background = scatteringProblem->incidentFieldPattern().get();
        }

        // Fluxes
        std::vector<LinearAlgebra::SmallVector<4>> fluxes =
            fluxFacets_.computeFluxes(result, problem_->omega(), background);

        outfile << std::setprecision(16) << problem_->omega();
        outfile << "," << std::setprecision(16)
                << wavenumbers_[nextCalled_].second;
        LinearAlgebra::SmallVector<4> totalFlux = {};
        for (const auto& entry : fluxes) {
            std::size_t nfluxes = problem_->isScatterFieldProblem() ? 4 : 1;
            for (std::size_t i = 0; i < nfluxes; ++i) {
                outfile << "," << std::setprecision(16) << entry[i];
            }
            totalFlux += entry;
        }
        if (computePMLTransmission.isUsed()) {
            std::vector<double> transmission =
                pmlTransmission_.pmlTransmission(result);
            double maxTransmission = 0.0, totalTransmission = 0.0;
            for (const auto& trans : transmission) {
                outfile << "," << std::setprecision(16) << trans;
                maxTransmission = std::max(maxTransmission, trans);
                totalTransmission += trans;
            }
            // For quick inspection
            DGMaxLogger(INFO,
                        "Maximum PML-transmission %, sum PML-transmission %",
                        maxTransmission, totalTransmission);
        }
        outfile << std::endl;
        DGMaxLogger(INFO, "Total net flux %", totalFlux);
    }

   private:
    void plotResult(DGMax::AbstractHarmonicResult<dim>& result,
                    Output::VTKSpecificTimeWriter<dim>& output);

    Base::MeshManipulator<dim>* mesh_;
    std::vector<std::pair<double, double>> wavenumbers_;
    std::vector<int> fieldPrintIndices_;
    std::shared_ptr<HarmonicProblem<dim>> problem_;
    std::ofstream outfile;
    int nextCalled_;
    DGMax::FluxFacets fluxFacets_;
    DGMax::PMLTransmission<dim> pmlTransmission_;
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

    auto pmlElementInfos = applyPMLs(*mesh);

    logger(INFO, "Loaded mesh % with % local elements", meshFileName.getValue(),
           mesh->getNumberOfElements());
    DGMax::writeMesh<dim>("mesh", *mesh);

    DGMax::HarmonicSolver<dim> solver(discretization);
    for (auto& pml : pmlElementInfos) {
        solver.addDispersive(pml);
    }

    Driver<dim> driver(*mesh);
    solver.solve(*mesh, driver);
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
        vectors["Poynting"] = [&scatteringProblem](Fields& fields) {
            // Real part of the Poynting vector
            auto totalField = fields.incident + fields.field;
            // E x (Curl E)^*
            auto rawPoyntingVector =
                -1.0 * LinearAlgebra::leftDoubledCrossProduct(
                           totalField, fields.totalCurl.conj());
            // Compensate for that -i mu omega Curl E = H
            return rawPoyntingVector.imag() /
                   (scatteringProblem->omega() *
                    fields.material.getPermittivity());
        };

        output.template writeMultiple<Fields>(
            [&result, &scatteringProblem](
                Base::Element* element, const Geometry::PointReference<dim>& p,
                std::size_t) {
                Fields fields;
                fields.field = result.computeField(element, p);
                auto pphys = element->template referenceToPhysical(p);
                fields.incident = scatteringProblem->incidentField(pphys);
                fields.totalCurl = result.computeFieldCurl(element, p) +
                                   scatteringProblem->incidentFieldCurl(pphys);
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

double parseDouble(const std::string& input) {
    std::size_t len;
    try {
        double result = std::stod(input, &len);
        logger.assert_always(len == input.size(),
                             "Failed to parse \"%\" as a number", input);
        return result;
    } catch (std::invalid_argument&) {
        DGMaxLogger.fail("Failed to parse \"%\" as a number", input);
    }
}

std::vector<std::pair<double, double>> parseWaveNumbers() {
    const std::string& input = wavenumberscm.getValue();
    std::vector<std::string> parts = DGMax::stringSplit(input, ',');
    std::vector<double> tempResult;
    if (parts.empty()) {
        DGMaxLogger.fail("Empty set of wavenumbers");
    }
    // Check the first value for an '@' to indicate equidistant spacing
    auto p0 = parts[0];
    DGMaxLogger.assert_always(!p0.empty(), "No first value");
    if (p0[0] == '@') {
        // Equidistant spaced numbers
        logger.assert_always(parts.size() == 3,
                             "Expected exactly 3 parts for a wavelength range");
        std::size_t len;
        unsigned long npoints = std::stoul(p0.substr(1), &len);
        DGMaxLogger.assert_always(
            len == p0.size() - 1,
            "Failed to convert \"%\" to a number of points", p0.substr(1));
        double w1 = parseDouble(parts[1]);
        double w2 = parseDouble(parts[2]);
        DGMaxLogger.assert_always(
            npoints >= 1, "Need at least 1 point for equidistant spacing.");
        double dw = (w2 - w1) / (npoints - 1);
        tempResult.reserve(npoints);
        for (std::size_t i = 0; i < npoints; ++i) {
            tempResult.push_back(w1 + i * dw);
        }
    } else {
        tempResult.reserve(parts.size());
        for (const std::string& waveLength : parts) {
            tempResult.push_back(parseDouble(waveLength));
        }
    }
    // Convert to computational wavenumbers
    // *1e-7    wavenumber (per cm) -> wavenumber (per nm)
    // *scalenm wavenumber (per nm) -> wavenumber (per unit length)
    // *2 PI    wavenumber -> angular wavenumber
    double factor = 1e-7 * meshsize.getValue() * 2 * M_PI;
    std::vector<std::pair<double, double>> result;
    result.reserve(tempResult.size());
    for (const double& waveLength : tempResult) {
        result.emplace_back(factor * waveLength, waveLength);
    }
    return result;
}

template <std::size_t dim>
std::vector<std::shared_ptr<PMLElementInfos<dim>>> applyPMLs(
    Base::MeshManipulator<dim>& mesh) {
    std::vector<DGMax::PMLZoneDescription<dim>> pmlDescriptions;

    // Step 1: Parse the input
    if (pmlString.isUsed()) {
        auto lines = DGMax::stringSplit(pmlString.getValue(), ';');
        for (const std::string& line : lines) {
            pmlDescriptions.push_back(
                DGMax::parsePMLZoneDescription<dim>(line));
        }
    }
    if (pmlFile.isUsed()) {
        std::ifstream file(pmlFile.getValue());
        DGMaxLogger.assert_always(
            file.good(), "Something went wrong in opening the PML file \"%\"",
            pmlFile.getValue());
        for (std::string line; std::getline(file, line);) {
            if (line.empty()) {
                continue;
            }
            pmlDescriptions.push_back(
                DGMax::parsePMLZoneDescription<dim>(line));
        }
    }
    // Step 2: Apply the input
    return DGMax::applyPMLs(mesh, pmlDescriptions);
}
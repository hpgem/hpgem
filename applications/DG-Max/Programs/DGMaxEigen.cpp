
#include <exception>
#include "Base/CommandLineOptions.h"
#include "Base/MeshFileInformation.h"
#include "Output/VTKSpecificTimeWriter.h"

#include "DGMaxLogger.h"
#include "DGMaxProgramUtils.h"
#include "Algorithms/DivDGMaxEigenvalue.h"
#include "Algorithms/DGMaxEigenvalue.h"
#include "Utils/KSpacePath.h"
#include "Utils/StructureDescription.h"
#include "Utils/PredefinedStructure.h"

using namespace hpgem;

// File name of the mesh file, e.g. -m mesh.hpgem
auto& meshFile = Base::register_argument<std::string>(
    'm', "meshFile", "The hpgem meshfile to use", true);
// Polynomial order of the basis functions, e.g. -p 2
auto& order = Base::register_argument<std::size_t>(
    'p', "order", "Polynomial order of the solution", true);
// Number of eigenvalues to compute, e.g. -e 40
auto& numEigenvalues = Base::register_argument<std::size_t>(
    'e', "eigenvalues", "The number of eigenvalues to compute", false, 24);

auto& method = Base::register_argument<std::string>(
    '\0', "method",
    "The method to be used, either 'DGMAX', 'DGMAX_PROJECT' or 'DIVDGMAX' "
    "(default)",
    false, "DIVDGMAX");

// Compute a single point --point 1,0.5,0 or a path of points
// [steps@]0,0:1,0:1,1
// Which corresponds to the point pi, 0.5pi, 0 in k-space and a path from 0,0
// via pi,0 to pi,pi in kspace taking each time "steps" steps, for each line,
// (excluding first point, including last).
// Untested with anything else than 0.0 as the first point of a path
auto& pointMode = Base::register_argument<std::string>(
    '\0', "points",
    "Compute a single point in or a set of lines through k-space", false);

// Number of steps to take in the k-space walk along the cube edge
// e.g. --steps 20
auto& steps = Base::register_argument<std::size_t>(
    '\0', "steps", "Steps for the k-space walk", false, 10);

// Penalty parameter, e.g. -a b5b0b5
auto& pparams = Base::register_argument<std::string>(
    'a', "penalty",
    "Penalty parameters and fluxes to use, e.g. b5b0b5.0\n"
    "\t  The letters are for the flux types and can be either b (Brezzi) "
    "or i (Interior penalty). The numbers are for the three penalty values.",
    false);

// Dimension, e.g. -d 2
auto& d = Base::register_argument<std::size_t>(
    'd', "dimension", "(deprecated) The dimension of the problem", false);

// Either a number for the predefined structures or a filename for zone based
// structure. See DGMax::determineStructureDescription for the exact format.
auto& structure = Base::register_argument<std::string>(
    '\0', "structure", "Structure to use", false, "0");

//
auto& fieldDir = Base::register_argument<std::string>(
    '\0', "fields", "Existing directory to output the fields to.", false, "");

// A natural mesh would be one where the unit cell in the mesh has a lattice
// constant of '1' in mesh coordinates. However, it may be more convenient to
// create a mesh at a different length scale. This option may be used to
// indicate and compensate for this difference in scale. Specifically:
//  1. The k-vectors are multiplied by 1/lengthscale
//  2. The resulting frequencies are multiplied by lengthscale
auto& lengthScale = Base::register_argument<double>(
    '\0', "lengthscale", "Length scale of the mesh", false, 1.0);

template <std::size_t DIM>
void runWithDimension();
double parseDGMaxPenaltyParameter();

DivDGMaxDiscretizationBase::Stab parsePenaltyParmaters();
template <std::size_t DIM>
KSpacePath<DIM> parsePath();

int main(int argc, char** argv) {
    registerLogLevelCommandLineFlag();
    Base::parse_options(argc, argv);
    initDGMaxLogging();
    DGMax::printArguments(argc, argv);

    time_t start, end;
    time(&start);

    const Base::MeshFileInformation info =
        Base::MeshFileInformation::readInformation(meshFile.getValue());
    const std::size_t dimension = info.dimension;
    // Check legacy dimension argument
    if (d.isUsed()) {
        logger.assert_always(
            d.getValue() == dimension,
            "Explicit dimension specified that does not match file contents");
    }
    try {
        switch (dimension) {
            case 2:
                runWithDimension<2>();
                break;
            case 3:
                runWithDimension<3>();
                break;
            default:
                logger.assert_always(false, "Invalid dimension %", dimension);
        }
    } catch (const std::exception& e) {
        DGMaxLogger(ERROR, e.what());
        exit(2);
    } catch (const char* message) {
        DGMaxLogger(ERROR, message);
        exit(1);
    }
    time(&end);
    DGMaxLogger(INFO, "Runtime %s", end - start);
    return 0;
}

template <std::size_t DIM>
class DGMaxEigenDriver : public AbstractEigenvalueSolverDriver<DIM> {

   public:
    DGMaxEigenDriver(KSpacePath<DIM>& path,
                     std::size_t targetNumberOfEigenvalues)
        : targetNumberOfEigenvalues_(targetNumberOfEigenvalues),
          currentPoint_(0),
          path_(path),
          frequencyResults_(path.totalNumberOfSteps()) {
        // Only one processor should write to the file
        Base::MPIContainer::Instance().onlyOnOneProcessor({[&]() {
            outFile.open("frequencies.csv");
            logger.assert_always(!outFile.fail(), "Output file opening failed");
            writeHeader(outFile, ',');
        }});
    }

    ~DGMaxEigenDriver() override {
        if (outFile.is_open()) {
            outFile.close();
        }
    }

    bool stop() const final {
        return currentPoint_ >= path_.totalNumberOfSteps();
    }

    void nextKPoint() final { ++currentPoint_; }

    LinearAlgebra::SmallVector<DIM> getCurrentKPoint() const final {
        logger.assert_debug(currentPoint_ < path_.totalNumberOfSteps(),
                            "Too large k-point index");
        return path_.k(currentPoint_);
    }

    std::size_t getNumberOfKPoints() const final {
        return path_.totalNumberOfSteps();
    }

    std::size_t getTargetNumberOfEigenvalues() const final {
        return targetNumberOfEigenvalues_;
    }

    void handleResult(AbstractEigenvalueResult<DIM>& result) final {
        frequencyResults_[currentPoint_] = result.getFrequencies();
        Base::MPIContainer::Instance().onlyOnOneProcessor({[&]() {
            writeFrequencies(outFile, currentPoint_,
                             frequencyResults_[currentPoint_], ',');
        }});

        if (fieldDir.isUsed()) {
            DGMaxLogger(INFO, "Writing field paterns");
            for (std::size_t i = 0; i < frequencyResults_[currentPoint_].size();
                 ++i) {
                std::stringstream outFile;
                // Note: Assumes the field directory is present
                if (fieldDir.hasArgument()) {
                    outFile << fieldDir.getValue() << "/";
                }

                outFile << "eigenfield-" << currentPoint_ << "-field-" << i;
                // Note the zero for timelevel is unused.
                Output::VTKSpecificTimeWriter<DIM> writer(
                    outFile.str(), result.getMesh(), 0, order.getValue());
                result.writeField(i, writer);
            }
        }

        writeOverlapIntegrals(result);
        writeKDerivativeFiles(result);
    }

    void printFrequencies() {
        writeHeader(std::cout, '\t');
        for (std::size_t i = 0; i < frequencyResults_.size(); ++i) {
            writeFrequencies(std::cout, i, frequencyResults_[i], '\t');
        }
    }

   private:
    std::size_t targetNumberOfEigenvalues_;
    std::size_t currentPoint_;
    KSpacePath<DIM> path_;
    std::vector<std::vector<double>> frequencyResults_;
    std::ofstream outFile;

    void writeHeader(std::ostream& stream, char separator) {
        // Mostly matching MPB
        // clang-format off
        stream << "freqs:"
               << separator << "k index"
               << separator << "kx/2pi"
               << separator << "ky/2pi"
               << separator << "kz/2pi"
               << separator << "kmag/2pi";
        // clang-format on

        // Add headers for the number of expected bands. The actual number may
        // be higher.
        for (std::size_t i = 0; i < targetNumberOfEigenvalues_; ++i) {
            stream << separator << "band " << i;
        }
        stream << std::endl;
    }

    void writeFrequencies(std::ostream& stream, std::size_t point,
                          std::vector<double>& frequencies, char separator) {
        auto k = path_.k(point);

        // Undo rescaling
        k *= lengthScale.getValue();

        // MPB like prefix. Importantly, we don't have a reciprocal lattice
        // defined so we output the kx-kz instead of k1-k3.

        // clang-format off
        stream << "freqs:"
               << separator << point + 1 //MPB uses 1 indexing
               << separator << k[0] / (2*M_PI)
               << separator << k[1] / (2*M_PI)
               << separator << (DIM == 3 ?  k[2] / (2*M_PI) : 0)
               << separator << k.l2Norm() / (2*M_PI);
        // clang-format on

        for (double frequency : frequencies) {
            // By convention the frequency is reported in units of
            // 2pi c/a, with c the speed of light and 'a' the lattice constant
            // (or similarly predefined length). For the computation we assume
            // c=1, and assume the length scale a matches that of the mesh. Thus
            // the distance between x=0 and x=1 is assumed to be 'a'.
            stream << separator
                   << frequency / (2 * M_PI) * lengthScale.getValue();
        }
        stream << std::endl;
    }

    void writeOverlapIntegrals(AbstractEigenvalueResult<DIM>& result) const {
        if (currentPoint_ == 0) {
            // No previous point
            return;
        }
        DGMaxLogger(INFO, "Writing overlap integrals");

        // Compute overlap integrals -> Needs all processors
        LinearAlgebra::MiddleSizeMatrix overlapIntegrals =
            result.computeFieldOverlap();
        Base::MPIContainer::Instance().onlyOnOneProcessor({[&]() {
            std::stringstream fileName;
            fileName << "overlap-" << currentPoint_ << ".csv";
            std::ofstream overlapFile(fileName.str());

            // Write header
            overlapFile << "Mode";
            for (std::size_t j = 0; j < overlapIntegrals.getNumberOfColumns();
                 ++j) {
                std::stringstream header;
                overlapFile << ","
                            << "Prev" << j;
            }
            overlapFile << std::endl;
            for (std::size_t i = 0; i < overlapIntegrals.getNumberOfRows();
                 ++i) {
                overlapFile << "New" << i;
                for (std::size_t j = 0;
                     j < overlapIntegrals.getNumberOfColumns(); ++j) {
                    overlapFile << "," << std::real(overlapIntegrals(i, j))
                                << "+" << std::imag(overlapIntegrals(i, j))
                                << "i";
                }
                overlapFile << std::endl;
            }
        }});
    }

    void writeKDerivativeFiles(AbstractEigenvalueResult<DIM>& result) const {
        if (!result.supportsWaveVectorDerivatives()) {
            return;
        }
        DGMaxLogger(INFO, "Writing wavevector derivatives");
        std::array<LinearAlgebra::MiddleSizeMatrix, DIM> derivatives =
            result.computeWaveVectorDerivatives();
        std::size_t numEigenvectors = derivatives[0].getNumberOfColumns();
        Base::MPIContainer::Instance().onlyOnOneProcessor({[&]() {
            std::stringstream fileName;
            fileName << "kderivatives-" << currentPoint_ << ".csv";
            std::ofstream kderivFile(fileName.str());
            // Write header
            kderivFile << "";
            for (std::size_t i = 0; i < numEigenvectors; ++i) {
                for (std::size_t kdir = 0; kdir < DIM; ++kdir) {
                    kderivFile << ",mode" << i << "d" << kdir;
                }
            }
            kderivFile << std::endl;
            // Body
            for (std::size_t i = 0; i < numEigenvectors; ++i) {
                kderivFile << "mode" << i;
                for (std::size_t j = 0; j < numEigenvectors; ++j) {
                    for (std::size_t kdir = 0; kdir < DIM; ++kdir) {
                        kderivFile << "," << std::real(derivatives[kdir](i, j))
                                   << std::showpos
                                   << std::imag(derivatives[kdir](i, j))
                                   << std::noshowpos << "i";
                    }
                }
                kderivFile << std::endl;
            }
        }});
    }
};

template <std::size_t DIM>
void runWithDimension() {
    bool useDivDGMax = true;
    DGMaxEigenvalueBase::ProjectorUse useProjector = DGMaxEigenvalueBase::NONE;
    std::size_t numberOfElementMatrices = 2;
    std::size_t unknowns = 0;
    // 2 unknowns, 1 time level
    if (method.getValue() == "DGMAX") {
        useDivDGMax = false;
        unknowns = 1;
    } else if (method.getValue() == "DGMAX_PROJECT" ||
               method.getValue() == "DGMAX_PROJECT1") {
        useDivDGMax = false;
        useProjector = method.getValue() == "DGMAX_PROJECT"
                           ? DGMaxEigenvalueBase::ALL
                           : DGMaxEigenvalueBase::INITIAL;
        unknowns = 2;
        // 1 more is needed for the projector operator
        numberOfElementMatrices = 3;
    } else if (method.getValue() == "DIVDGMAX") {
        useDivDGMax = true;
        unknowns = 2;
    } else {
        logger(ERROR,
               "Invalid method {}, should be either DGMAX, DGMAX_PROJECT, "
               "DGMAX_PROJECT1 or DIVDGMAX",
               method.getValue());
        return;
    }

    Base::ConfigurationData configData(unknowns, 1);

    std::unique_ptr<DGMax::StructureDescription> structureDesc =
        DGMax::determineStructureDescription(structure.getValue(), DIM);

    auto mesh = DGMax::readMesh<DIM>(meshFile.getValue(), &configData,
                                     *structureDesc, numberOfElementMatrices);
    logger(INFO, "Loaded mesh % with % local elements", meshFile.getValue(),
           mesh->getNumberOfElements());
    DGMax::writeMesh<DIM>("mesh", *mesh);
    // TODO: Parameterize

    KSpacePath<DIM> path = parsePath<DIM>();
    DGMaxEigenDriver<DIM> driver(path, numEigenvalues.getValue());

    // Method dependent solving
    if (useDivDGMax) {
        DivDGMaxDiscretizationBase::Stab stab = parsePenaltyParmaters();
        DivDGMaxEigenvalue<DIM> solver(*mesh, order.getValue(), stab);
        solver.solve(driver);
    } else {
        const double stab = parseDGMaxPenaltyParameter();
        DGMaxEigenvalueBase::SolverConfig config;
        config.stab_ = stab;
        config.useHermitian_ = true;
        config.shiftFactor_ = 0;
        config.useProjector_ = useProjector;
        DGMaxEigenvalue<DIM> solver(*mesh, order.getValue(), config);
        solver.solve(driver);
    }

    Base::MPIContainer::Instance().onlyOnOneProcessor(
        {[&]() { driver.printFrequencies(); }});
}

/// Parse DIM comma separated numbers as the coordinates of a point.
/// \tparam DIM The dimension of the point
/// \param pointString The string containing the point coordinates
/// \param start The starting index in pointString
/// \param point The point (out)
/// \return The first index in pointString after the number.
template <std::size_t DIM>
std::size_t parsePoint(const std::string& pointString, std::size_t start,
                       LinearAlgebra::SmallVector<DIM>& point) {
    for (std::size_t i = 0; i < DIM; ++i) {
        if (start >= pointString.size()) {
            throw std::invalid_argument(
                "Not enough coordinates for a reciprocal point");
        }
        std::size_t len = 0;
        point[i] = std::stod(pointString.substr(start), &len);
        if (len == 0) {
            throw std::invalid_argument("No value parsed");
        }
        start += len;
        if (i < DIM - 1) {
            // Skip the comma, space, whatever that ended the point
            start++;
        }
    }
    return start;
}

template <std::size_t DIM>
KSpacePath<DIM> parsePath() {
    if (pointMode.isUsed()) {
        DGMax::PointPath<DIM> path =
            DGMax::parsePath<DIM>(pointMode.getValue());
        // Compensate for factor of pi in the reciprocal lattice
        for (std::size_t i = 0; i < path.points_.size(); ++i) {
            path.points_[i] *= M_PI / lengthScale.getValue();
        }
        // Default steps to 1.
        if (path.steps_ < 0) {
            path.steps_ = 1;
        }
        return KSpacePath<DIM>(path.points_, (std::size_t)path.steps_);
    }
    if (!steps.isUsed()) {
        logger(INFO, "Using default number of steps %", steps.getValue());
    }
    return KSpacePath<DIM>::cubePath(steps.getValue(), false);
}

double parseDGMaxPenaltyParameter() {
    if (pparams.isUsed()) {
        std::size_t idx;
        try {
            double value = std::stod(pparams.getValue(), &idx);
            if (idx != pparams.getValue().size()) {
                throw std::invalid_argument(
                    "Invalid stabilization parameter, should be a single "
                    "number for DGMAX");
            }
            return value;
        } catch (const std::invalid_argument&) {
            throw std::invalid_argument(
                "Invalid stabilization parameter, should be a single number "
                "for DGMAX");
        }
    } else {
        // Default
        return 100;
    }
}

DivDGMaxDiscretizationBase::Stab parsePenaltyParmaters() {
    if (pparams.isUsed()) {
        DivDGMaxDiscretizationBase::Stab stab;
        std::string input = pparams.getValue();
        std::vector<bool> useBrezzi;
        std::vector<double> values;

        std::size_t index = 0;
        bool error = false;
        // Only need three parameters
        while (input.size() > index && !error && useBrezzi.size() < 3) {
            // Parse the flux type
            char fluxType = input[index++];
            if (fluxType == 'b')
                useBrezzi.emplace_back(true);
            else if (fluxType == 'i')
                useBrezzi.emplace_back(false);
            else {
                DGMaxLogger(ERROR, "Unknown flux type %", fluxType);
                error = true;
                break;
            }
            // Find the numeric digits
            std::size_t startIndex = index;
            while (input.size() > index) {
                char c = input[index];
                if (std::isdigit(c) || c == '.' || c == '-' || c == 'e') {
                    index++;
                } else {
                    // Not a valid part of the number
                    break;
                }
            }
            // Parse the number (if present)
            if (index != startIndex) {
                // There is some numeric content
                double value =
                    std::stod(input.substr(startIndex, index - startIndex));
                values.emplace_back(value);
            } else if (index == input.size()) {
                DGMaxLogger(ERROR,
                            "Not enough input for the stabilization parameter");
                error = true;
            } else {
                DGMaxLogger(
                    ERROR,
                    "No number at position % of the stabilization parameters",
                    index);
                error = true;
            }
        }

        // Check the validity of the result
        if (!error && useBrezzi.size() < 3) {
            DGMaxLogger(ERROR,
                        "Not enough stabilization parameters only parsed %",
                        useBrezzi.size());
            error = true;
        }
        if (!error && index != input.size()) {
            DGMaxLogger(ERROR, "Unconsumed stabilization parameter input at %",
                        index);
            error = true;
        }
        if (!error) {
            DivDGMaxDiscretizationBase::Stab result;
            result.stab1 = values[0];
            result.stab2 = values[1];
            result.stab3 = values[2];

            using FLUX = DivDGMaxDiscretizationBase::FluxType;

            result.fluxType1 = useBrezzi[0] ? FLUX::BREZZI : FLUX::IP;
            result.fluxType2 = useBrezzi[1] ? FLUX::BREZZI : FLUX::IP;
            result.fluxType3 = useBrezzi[2] ? FLUX::BREZZI : FLUX::IP;
            logger(INFO, "Using fluxes and stabilization: %", result);
            return result;
        }
        throw std::invalid_argument("Invalid stabilization parameter");

    } else {
        // Default values
        DivDGMaxDiscretizationBase::Stab stab;
        stab.stab1 = 5;
        stab.stab2 = 0;
        stab.stab3 = 5;
        stab.setAllFluxeTypes(DivDGMaxDiscretizationBase::FluxType::BREZZI);
        return stab;
    }
}

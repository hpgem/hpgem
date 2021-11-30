
#include "DGMaxProgramUtils.h"

#include "Base/MeshManipulator.h"

#include "DGMaxLogger.h"
#include "Utils/PredefinedStructure.h"
#include "Utils/ZoneStructureDescription.h"

#ifdef HPGEM_USE_MPI
// In case of MPI it is usefull to know where each process is located
#include <unistd.h>

#include "CMakeDefinitions.h"
#include <Base/CommandLineHelpers.h>

using namespace hpgem;

#endif

namespace DGMax {
void printArguments(int argc, char** argv) {
#ifdef HPGEM_USE_MPI
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // Apperently a HOST_NAME_MAX macro is defined on linux, but requires
    // quite a bit of macro fiddling, with the standard value of 64. Instead
    // we define it here, with a +1 for the null termination of the string.
    const std::size_t HOSTNAMEMAX = 65;
    char hostname[HOSTNAMEMAX];
    int ierr = gethostname(hostname, HOSTNAMEMAX);
    if (ierr == -1) {
        std::cerr << "Hostname error" << std::endl;
        // Almost surely a too long hostname, make sure it is null-terminated
        hostname[HOSTNAMEMAX - 1] = '\0';
    }
    logAll([&]() { DGMaxLogger(INFO, "Proc %/% on %", rank, size, hostname); });
#endif
    if (!loggingSuppressed()) {
        std::stringstream stream;
        stream << "Program arguments: " << std::endl;
        for (int i = 0; i < argc; ++i) {
            if (i != 0) stream << " ";
            stream << argv[i];
        }
        std::string message = stream.str();
        DGMaxLogger(INFO, message);
        DGMaxLogger(INFO, "git ref: % (%%)", getGITRef(),
                    isGITClean() ? "" : "dirty ", getGITHash());
    }
}

template <std::size_t DIM>
std::unique_ptr<Base::MeshManipulator<DIM>> readMesh(
    std::string fileName, Base::ConfigurationData* configData,
    StructureDescription& structureDescription,
    std::size_t numberOfElementMatrices) {

    // One for the double curl, one for the impedance part
    const std::size_t faceMatrices = 2;

    auto mesh = std::unique_ptr<Base::MeshManipulator<DIM>>(
        new Base::MeshManipulator<DIM>(configData, numberOfElementMatrices, 3,
                                       faceMatrices, 1));
    mesh->readMesh(fileName);
    for (typename Base::MeshManipulator<DIM>::ElementIterator it =
             mesh->elementColBegin(Base::IteratorType::GLOBAL);
         it != mesh->elementColEnd(Base::IteratorType::GLOBAL); ++it) {
        (*it)->setUserData(structureDescription.createElementInfo(*it));
        (*it)->setNumberOfTimeIntegrationVectors(1);
    }
    return mesh;
}

// Explicit instantiation of the 2,3D versions.
template std::unique_ptr<Base::MeshManipulator<2>> readMesh(
    std::string, Base::ConfigurationData*,
    StructureDescription& structureDescription, std::size_t);
template std::unique_ptr<Base::MeshManipulator<3>> readMesh(
    std::string, Base::ConfigurationData*,
    StructureDescription& structureDescription, std::size_t);

template <std::size_t DIM>
PointPath<DIM> parsePath(const std::string& path) {
    LinearAlgebra::SmallVector<DIM> point;
    std::size_t index = 0;  // Current parsing index

    int steps = -1;
    // Check if the string starts with number@, to denote step count
    std::size_t atIndex = path.find_first_of('@');
    if (atIndex != std::string::npos) {
        std::size_t len = 0;
        steps = std::stoul(path, &len);
        index = atIndex + 1;  // Start parsing points after the @
        if (len != atIndex) {
            throw std::invalid_argument("Left between number of steps and '@'");
        }
    }

    // Parse a string of points
    std::vector<LinearAlgebra::SmallVector<DIM>> points;
    while (index < path.length()) {
        LinearAlgebra::SmallVector<DIM> point;
        index = Base::parsePoint(path, index, point);
        points.push_back(point);
        // Strip character if needed
        while (
            index < path.length() &&
            (std::isspace(path[index])  // Allow spaces, just as sto*-functions
             || path[index] == ':'))  // Allow colon separation for readability
        {
            index++;
        }
    }
    return PointPath<DIM>(points, steps);
}

// Explicit instantiation

template PointPath<1> parsePath(const std::string& path);
template PointPath<2> parsePath(const std::string& path);
template PointPath<3> parsePath(const std::string& path);

std::unique_ptr<ZoneInfoStructureDefinition> readZonedDescription(
    std::istream& file, char separator) {
    std::string line;
    std::vector<std::regex> zoneRegexes;
    std::vector<double> zoneEpsilons;
    std::size_t lineNumber = 0;
    while (!file.eof()) {
        std::getline(file, line, separator);
        // Discard empty lines
        if (line.empty()) {
            continue;
        }
        // Find the position of the comma, the separator between the regex and
        // values.
        std::size_t commaLoc = line.find_last_of(',');
        logger.assert_always(commaLoc != std::string::npos,
                             "No comma found on line %: line", lineNumber,
                             line);
        std::string regexStr = line.substr(0, commaLoc);
        std::regex regex(regexStr);

        // Extract the epsilon value
        std::string epsilonStr = line.substr(commaLoc + 1);
        logger.assert_always(!epsilonStr.empty(),
                             "No value after the comma on line %: %",
                             lineNumber, line);
        std::size_t idx;
        double epsilon = 0.0;
        try {
            epsilon = std::stod(epsilonStr, &idx);
        } catch (std::invalid_argument&) {
            logger.assert_always(false, "Invalid epsilon value on line %: %",
                                 lineNumber, line);
        }
        logger.assert_always(idx == epsilonStr.size(),
                             "Trailing data after epsilon on line %: %",
                             lineNumber, line);
        zoneRegexes.push_back(regex);
        zoneEpsilons.push_back(epsilon);

        lineNumber++;
    }
    return std::make_unique<ZoneInfoStructureDefinition>(zoneRegexes,
                                                         zoneEpsilons);
}

/**
 * Tries to parse a string which contains a integer in base 10.
 * @param input The input string to parse
 * @param out The parsed number
 * @return Whether successful, that is whether the whole string was consumed.
 */
bool tryParseLong(const std::string& input, long& out) {
    const char* cinput = input.c_str();
    char* end;
    out = std::strtol(cinput, &end, 10);
    // End points to the first entry after the numerical value, if parsing was
    // successful this means that the complete string was consumed and the first
    // character is '\0'.
    return !*end;
}

std::unique_ptr<StructureDescription> determineStructureDescription(
    const std::string& input, std::size_t dim) {

    // Test to see if the input is a number, by trying to parse it
    long value;
    if (tryParseLong(input, value)) {
        DGMaxLogger(INFO, "Using predefined structure %", value);
        return std::make_unique<PredefinedStructureDescription>(
            DGMax::structureFromInt(value), dim);
    } else if (input.find(',') != std::string::npos) {
        // Inline regexes
        std::istringstream iinput(input);
        return readZonedDescription(iinput, ';');
    }

    std::ifstream file;
    file.open(input);
    if (file.good()) {
        DGMaxLogger(INFO, "Reading structure defined in file %", input);
        return readZonedDescription(file, '\n');
    }

    logger.assert_always(
        false, "Could not determine structure information type from input '%'",
        input);
    return nullptr;
}

}  // namespace DGMax

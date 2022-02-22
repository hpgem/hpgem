
#include "DGMaxProgramUtils.h"

#include "Base/MeshManipulator.h"

#include "DGMaxLogger.h"
#include "Utils/PredefinedStructure.h"
#include "Utils/ZoneStructureDescription.h"

#ifdef HPGEM_USE_MPI
// In case of MPI it is usefull to know where each process is located
#include <unistd.h>

#include "CMakeDefinitions.h"

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
        try {
            point[i] = std::stod(pointString.substr(start), &len);
        } catch (const std::invalid_argument&) {
            // No parse, i.e. len == 0
            throw std::invalid_argument("Number point parsing failed at '" +
                                        pointString.substr(start) +
                                        "', expected a coordinate");
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
        index = parsePoint(path, index, point);
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
        DGMaxLogger(INFO, "Adding zone regexp '%' with material %", regexStr,
                    epsilon);
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

std::vector<std::string> stringSplit(const std::string& input, char separator) {
    if (input.empty()) {
        return {};
    }
    std::size_t pos = 0;
    std::size_t next_pos;
    std::vector<std::string> result;
    while ((next_pos = input.find_first_of(separator, pos)) !=
           std::string::npos) {
        result.push_back(input.substr(pos, next_pos - pos));
        pos = next_pos + 1;
    }
    // Remainder after the last separator
    result.push_back(input.substr(pos));
    return result;
}

template <std::size_t dim>
std::vector<
    std::pair<Geometry::PointPhysical<dim>, Geometry::PointPhysical<dim>>>
    computeZoneBoundingBoxes(const Base::MeshManipulator<dim>& mesh) {

    logger.assert_always(Base::MPIContainer::Instance().getNumProcessors() == 1,
                         "Not MPI Enabled");

    std::size_t zoneCount = mesh.getZones().size();
    // Two vectors with the minimum/maximum coordinates for each zone.
    // Sequential layout (e.g. x_0, y_0, x_1, y_1, ... with subscripts for
    // zoneIds). Using vectors allows easier MPI communication.
    //
    // Default value +- infinity, as that is always updated by
    // std::min/std::max
    std::vector<double> mins(dim * zoneCount,
                             std::numeric_limits<double>::infinity());
    std::vector<double> maxs(dim * zoneCount,
                             -std::numeric_limits<double>::infinity());

    // Update mins/maxs for each element based on the node coordinates.
    for (const Base::Element* element : mesh.getElementsList()) {
        if (!element->isOwnedByCurrentProcessor()) {
            continue;
        }
        const Geometry::PhysicalGeometryBase* pgeom =
            element->getPhysicalGeometry();
        std::size_t zoneId = element->getZone().getZoneId();
        std::size_t nodeCount = pgeom->getNumberOfNodes();
        std::size_t offset = dim * zoneId;
        for (std::size_t i = 0; i < nodeCount; ++i) {
            const Geometry::PointPhysical<dim> node =
                pgeom->getLocalNodeCoordinates(i);
            for (std::size_t j = 0; j < dim; ++j) {
                mins[offset + j] = std::min(mins[offset + j], node[j]);
                maxs[offset + j] = std::max(maxs[offset + j], node[j]);
            }
        }
    }

#ifdef HPGEM_USE_MPI
    // We only know about the bounds due to the elements in our own part of the
    // mesh. Communicate with the other processors to get the global bounds.
    MPI_Allreduce(MPI_IN_PLACE, mins.data(), dim * zoneCount, MPI_DOUBLE,
                  MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, maxs.data(), dim * zoneCount, MPI_DOUBLE,
                  MPI_MAX, MPI_COMM_WORLD);
#endif

    // Convert the results to output types
    using Bounds =
        std::pair<Geometry::PointPhysical<dim>, Geometry::PointPhysical<dim>>;
    std::vector<Bounds> result(zoneCount);
    for (std::size_t i = 0; i < zoneCount; ++i) {
        Bounds& bounds = result[i];
        std::size_t offset = dim * i;
        for (std::size_t j = 0; j < dim; ++j) {
            bounds.first[j] = mins[offset + j];
            bounds.second[j] = maxs[offset + j];
        }
    }
    return result;
}

template
std::vector<std::pair<Geometry::PointPhysical<2>, Geometry::PointPhysical<2>>>
    computeZoneBoundingBoxes(const Base::MeshManipulator<2>& mesh);
template
std::vector<std::pair<Geometry::PointPhysical<3>, Geometry::PointPhysical<3>>>
    computeZoneBoundingBoxes(const Base::MeshManipulator<3>& mesh);

template <std::size_t dim>
PMLZoneDescription<dim> parsePMLZoneDescription(const std::string& input) {
    PMLZoneDescription<dim> result;

    std::vector<std::string> parts = stringSplit(input, ',');
    DGMaxLogger.assert_always(
        parts.size() == 2 + dim,
        "Expected 2+dim fields to a PML description but got % in \"%\"",
        parts.size(), input);
    DGMaxLogger.assert_always(!parts[0].empty(),
                              "Empty zone name for PML description");
    result.zoneName_ = std::move(parts[0]);

    DGMaxLogger.assert_always(
        parts[1].size() == dim,
        "Expected exactly % direction characters but got \"\"", dim, parts[1]);

    for (std::size_t j = 0; j < dim; ++j) {
        switch (parts[1][j]) {
            case '0':
                result.direction_[j] = 0.0;
                break;
            case '+':
                result.direction_[j] = 1.0;
                break;
            case '-':
                result.direction_[j] = -1.0;
                break;
            default:
                DGMaxLogger.fail(
                    "Unknown direction character \"%\" for the %-th direction",
                    parts[1][j], j);
        }
        {

            std::size_t p;
            double attenuation = std::stod(parts[2 + j], &p);
            DGMaxLogger.assert_always(
                p == parts[2 + j].length(),
                "Left over information after parsing attenuation\"%\"",
                parts[2 + j]);
            result.attenuation_[j] = attenuation;
            if (result.direction_[j] != 0) {
                DGMaxLogger.assert_always(attenuation <= 1.0 && attenuation > 0,
                                          "Invalid attenuation value %",
                                          attenuation);
            }
        }
    }

    return result;
}

template
PMLZoneDescription<2> parsePMLZoneDescription(const std::string& input);
template
PMLZoneDescription<3> parsePMLZoneDescription(const std::string& input);

template <std::size_t dim>
std::vector<std::shared_ptr<PMLElementInfos<dim>>> applyPMLs(
    Base::MeshManipulator<dim>& mesh,
    const std::vector<PMLZoneDescription<dim>>& pmls) {
    if (pmls.empty()) {
        return {};
    }
    const int NO_PML_NEEDED = std::numeric_limits<int>::min();

    auto boundingBoxes = computeZoneBoundingBoxes(mesh);
    // Actual PMLElementInfos used by this MPI rank
    std::vector<std::shared_ptr<PMLElementInfos<dim>>> pmlinfos;
    // Index of the PMLElementInfos for each zone.
    //   - Positive values are indices in pmlinfos
    //   - Negative values are indices in pmls (offset by 1)
    //   - NO_PML_NEEDED is a signalling value that there is no PML in the
    //   region
    const auto& zones = mesh.getZones();
    std::vector<int> pmlIndices(zones.size(), NO_PML_NEEDED);
    for (int i = 0; i < pmls.size(); ++i) {
        for (std::size_t j = 0; j < zones.size(); ++j) {
            if (pmls[i].zoneName_ == zones[j]->getName()) {
                pmlIndices[j] = -i - 1;
                break;
            }
        }
    }

    for (Base::Element* element : mesh.getElementsList()) {
        std::size_t zoneId = element->getZone().getZoneId();
        int index = pmlIndices[zoneId];
        if (index == NO_PML_NEEDED) {
            continue;
        }
        if (index < 0) {
            // PML needs to be generated
            const PMLZoneDescription<dim>& description = pmls[-index - 1];

            auto box = boundingBoxes[zoneId];
            // Size of the bounding box
            auto pmlThickness =
                box.second.getCoordinates() - box.first.getCoordinates();
            const Material& material =
                ElementInfos::get(*element).getMaterial();

            LinearAlgebra::SmallVector<dim> scaling =
                PMLElementInfos<dim>::computeScaling(
                    material, description.direction_, pmlThickness,
                    description.attenuation_);

            LinearAlgebra::SmallVector<dim> offset;
            for (std::size_t i = 0; i < dim; ++i) {
                // Should be +-1 or 0
                int direction = static_cast<int>(description.direction_[i]);
                if (direction < 0) {
                    // e.g. PML for x < 0
                    offset[i] = box.second[i];
                } else if (direction > 0) {
                    offset[i] = box.first[i];
                } else {
                    offset[i] = 0.0;
                }
            }
            pmlinfos.emplace_back(std::make_shared<PMLElementInfos<dim>>(
                material, offset, description.direction_, scaling));
            pmlIndices[zoneId] = pmlinfos.size() - 1;
        }
        // Now that we know that a PML is present, set it
        element->setUserData(pmlinfos[pmlIndices[zoneId]].get());
    }
    return pmlinfos;
}

template
std::vector<std::shared_ptr<PMLElementInfos<2>>> applyPMLs(
    Base::MeshManipulator<2>& mesh,
    const std::vector<PMLZoneDescription<2>>& pmls);
template
std::vector<std::shared_ptr<PMLElementInfos<3>>> applyPMLs(
    Base::MeshManipulator<3>& mesh,
    const std::vector<PMLZoneDescription<3>>& pmls);
}  // namespace DGMax

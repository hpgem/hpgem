
#include "DGMaxProgramUtils.h"

#include "Base/MeshManipulator.h"

#include "DGMaxLogger.h"
#include "ElementInfos.h"

#ifdef HPGEM_USE_MPI
// In case of MPI it is usefull to know where each process is located
#include <unistd.h>

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
    }
}

template <std::size_t DIM>
std::unique_ptr<Base::MeshManipulator<DIM>> readMesh(
    std::string fileName, Base::ConfigurationData* configData,
    ElementInfos::EpsilonFunc<DIM> epsilon) {
    auto mesh = std::unique_ptr<Base::MeshManipulator<DIM>>(
        new Base::MeshManipulator<DIM>(configData, 2, 3, 1, 1));
    mesh->readMesh(fileName);
    for (typename Base::MeshManipulator<DIM>::ElementIterator it =
             mesh->elementColBegin(Base::IteratorType::GLOBAL);
         it != mesh->elementColEnd(Base::IteratorType::GLOBAL); ++it) {
        (*it)->setUserData(ElementInfos::createStructure<DIM>(**it, epsilon));
    }
    return mesh;
}

// Explicit instantiation of the 2,3D versions.
template std::unique_ptr<Base::MeshManipulator<2>> readMesh(
    std::string, Base::ConfigurationData*, ElementInfos::EpsilonFunc<2>);
template std::unique_ptr<Base::MeshManipulator<3>> readMesh(
    std::string, Base::ConfigurationData*, ElementInfos::EpsilonFunc<3>);

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
}  // namespace DGMax

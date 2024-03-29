
#ifndef HPGEM_APP_DGMAXPROGRAMUTILS_H
#define HPGEM_APP_DGMAXPROGRAMUTILS_H

#include <memory>
#include "Utils/StructureDescription.h"
#include "PMLElementInfos.h"

namespace hpgem {
namespace Base {
class ConfigurationData;

template <std::size_t DIM>
class MeshManipulator;
}  // namespace Base
}  // namespace hpgem
using namespace hpgem;
namespace DGMax {

void printArguments(int argc, char** argv);

template <std::size_t DIM>
std::unique_ptr<Base::MeshManipulator<DIM>> readMesh(
    std::string fileName, Base::ConfigurationData* configData,
    StructureDescription& structureDescription,
    std::size_t numberOfElementMatrices = 2);

/// A path in either real or reciprocal space
template <std::size_t DIM>
struct PointPath {
    /// The points along the path
    std::vector<LinearAlgebra::SmallVector<DIM>> points_;
    /// Possible number of steps specified to take to that point, default -1
    /// if not specified
    int steps_;

    PointPath(std::vector<LinearAlgebra::SmallVector<DIM>> points,
              int steps = -1)
        : points_(points), steps_(steps) {}
};

/// \brief Parse a path description
/// Parse a path, a path follows the syntax
///    path = (steps '@')? point (':' point)*
///    point = coord (',' coord){dim-1}, where dim is the dimension
///    coord = double (
///    steps = unsigned integer
/// Thus a 2D path may look like 20@1,0:1,1 meaning two points 1,0 and 1,1 with
/// 20 steps \tparam DIM The dimensionality of the path \param path The path
/// string to parse \return The points on the path
template <std::size_t DIM>
PointPath<DIM> parsePath(const std::string& path);

/// Heuristically determine the way to determine the structure description.
/// It accepts three forms:
///   - A single integer: the index into the predefined structures
///     example: "1"
///   - A filename: The file should contain lines of the form 'regex,epsilon',
///     where regex is a regular expression that matches some of the zones and
///     epsilon is the corresponding epsilon. Empty lines are ignored.
///     example: "zones.csv" with as content of the file:
///       Silicon*,12.1
///       Pore*,1
///   - An inline version of file based case. In this case the 'regex,epsilon'
///     pairs are separated by semicolons.
///     example (equivalent to the file one): "Silicon*,12.1;Pore*,1"
///
/// \param input The index or file name
/// \param dim The dimension of the mesh
/// \return A structure description based on the input
std::unique_ptr<StructureDescription> determineStructureDescription(
    const std::string& input, std::size_t dim);

std::vector<std::string> stringSplit(const std::string& input, char separator);

/// Computes the bounding box for each zone
///
/// Computes the axis-aligned bounding box for each region. For linear elements
/// this is exact. For curvilinear elements on the boundary a small part of the
/// element may fall outside. It is guaranteed that all reference points are
/// inside the bounding box.
///
/// \tparam dim The dimension of the mesh
/// \param mesh The mesh
/// \return For each zone a vector with minimum and maximum coordinates. Result
/// for regions without elements are undefined.
template <std::size_t dim>
std::vector<
    std::pair<Geometry::PointPhysical<dim>, Geometry::PointPhysical<dim>>>
    computeZoneBoundingBoxes(const Base::MeshManipulator<dim>& mesh);

template <std::size_t dim>
struct PMLZoneDescription {
    std::string zoneName_;
    LinearAlgebra::SmallVector<dim> direction_;
    // Attenuation from the PML (excluding far side boundary condition) based on
    // a plane wave in the i-th direction. This includes both the path from the
    // incident face to the far end, and the way back (far end -> interface)
    LinearAlgebra::SmallVector<dim> attenuation_;
};

/// Parses the description of a single PML zone.
///
/// Parses the description of pml zones, it expects a comma separated line with
/// the following fields:
///  - The zone name of the region to which to apply the PML
///  - The direction, specified as a string of length dim containing the
///    characters {+,0,-}. For example in dim=3 this could be "+0-".
///  - dim fields with the attenuation in x,y(,z) direction specified as double
///    in the range (0, 1].
///
/// Example in 2D: 'PML-YPlus,0+,1,1e-2', describes a PML for the zone with name
/// PML-YPlus. The PML attenuates in the y-direction increasing in +y direction,
/// the attenuation for a plane wave in the y-direction is 1e-2. An attenuation
/// of 1 is specified in x direction, but this does not matter as the PML does
/// not attenuate in that direction.
///
/// \tparam dim The dimension of the mesh
/// \param input The line with input information
/// \return The zonedescription
template <std::size_t dim>
PMLZoneDescription<dim> parsePMLZoneDescription(const std::string& input);

template <std::size_t dim>
std::vector<std::shared_ptr<PMLElementInfos<dim>>> applyPMLs(
    Base::MeshManipulator<dim>& mesh,
    const std::vector<PMLZoneDescription<dim>>& pmls);

/// Write a mesh (including epsilon) to file
/// \tparam dim The dimension of the mesh (instantiated for 2,3)
/// \param fileName The prefix for the vtk files
/// \param mesh The mesh to write
template <std::size_t dim>
void writeMesh(const std::string& fileName, Base::MeshManipulator<dim>& mesh);

}  // namespace DGMax

#endif  // HPGEM_APP_DGMAXPROGRAMUTILS_H

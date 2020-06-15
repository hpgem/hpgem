
#ifndef HPGEM_PROGRAMUTILS_H
#define HPGEM_PROGRAMUTILS_H

#include <memory>
#include <ElementInfos.h>

namespace Base
{
    class ConfigurationData;

    template<std::size_t DIM>
    class MeshManipulator;
}

namespace DGMax
{

    void printArguments(int argc, char** argv);

    template<std::size_t DIM>
    std::unique_ptr<Base::MeshManipulator<DIM>>
            readMesh(std::string fileName, Base::ConfigurationData *configData,
                    ElementInfos::EpsilonFunc<DIM> epsilonFunc,
                    std::size_t numberOfElementMatrices = 2);


    /// A path in either real or reciprocal space
    template<std::size_t DIM>
    struct PointPath
    {
        /// The points along the path
        std::vector<LinearAlgebra::SmallVector<DIM>> points_;
        /// Possible number of steps specified to take to that point, default -1
        /// if not specified
        int steps_;

        PointPath(std::vector<LinearAlgebra::SmallVector<DIM>> points, int steps = -1)
                : points_ (points), steps_(steps)
        {}
    };

    /// \brief Parse a path description
    /// Parse a path, a path follows the syntax
    ///    path = (steps '@')? point (':' point)*
    ///    point = coord (',' coord){dim-1}, where dim is the dimension
    ///    coord = double (
    ///    steps = unsigned integer
    /// Thus a 2D path may look like 20@1,0:1,1 meaning two points 1,0 and 1,1 with 20 steps
    /// \tparam DIM The dimensionality of the path
    /// \param path The path string to parse
    /// \return The points on the path
    template<std::size_t DIM>
    PointPath<DIM> parsePath(const std::string& path);
}

#endif //HPGEM_PROGRAMUTILS_H


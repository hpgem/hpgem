
#ifndef HPGEM_PROGRAMUTILS_H
#define HPGEM_PROGRAMUTILS_H

#include <memory>

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
            createCubeMesh(std::size_t subdivisions, Base::ConfigurationData *configData);

    template<std::size_t DIM>
    std::unique_ptr<Base::MeshManipulator<DIM>>
            readMesh(std::string fileName, Base::ConfigurationData *configData);
}

#endif //HPGEM_PROGRAMUTILS_H


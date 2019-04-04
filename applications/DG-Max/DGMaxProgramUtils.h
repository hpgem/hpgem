
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
            createCubeMesh(std::size_t subdivisions, Base::ConfigurationData *configData,
                    ElementInfos::EpsilonFunc<DIM> epsilonFunc);

    template<std::size_t DIM>
    std::unique_ptr<Base::MeshManipulator<DIM>>
            readMesh(std::string fileName, Base::ConfigurationData *configData,
                    ElementInfos::EpsilonFunc<DIM> epsilonFunc);
}

#endif //HPGEM_PROGRAMUTILS_H



#include "DGMaxProgramUtils.h"

#include "Base/MeshManipulator.h"

#include "DGMaxLogger.h"
#include "ElementInfos.h"

namespace DGMax
{
    void printArguments(int argc, char** argv)
    {
#ifdef HPGEM_USE_MPI
        int rank;
        int size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        logAll([&](){DGMaxLogger(INFO, "Proc %/%", rank, size);});
#endif
        if (!loggingSuppressed())
        {
            std::stringstream stream;
            stream << "Program arguments: " << std::endl;
            for(int i = 0; i < argc; ++i)
            {
                if (i != 0) stream << " ";
                stream << argv[i];
            }
            std::string message = stream.str();
            DGMaxLogger(INFO, message);
        }
    }

    template<std::size_t DIM>
    std::unique_ptr<Base::MeshManipulator<DIM>> createCubeMesh(std::size_t subdivisions, Base::ConfigurationData* configData)
    {
        Geometry::PointPhysical<DIM> bottomLeft, topRight;
        std::vector<std::size_t> numElementsOneD (DIM);
        std::vector<bool> periodic(DIM, true);
        // Configure each dimension of the unit cube/square
        for (std::size_t i = 0; i < DIM; ++i) {
            bottomLeft[i] = 0;
            topRight[i] = 1;
            numElementsOneD[i] = subdivisions;
        }

        auto mesh = std::unique_ptr<Base::MeshManipulator<DIM>>(
                new Base::MeshManipulator<DIM>(configData, 2, 3, 1, 1));
        mesh->createTriangularMesh(bottomLeft, topRight, numElementsOneD, periodic);

        for (typename Base::MeshManipulator<DIM>::ElementIterator it = mesh->elementColBegin(Base::IteratorType::GLOBAL);
             it != mesh->elementColEnd(Base::IteratorType::GLOBAL); ++it)
        {
            (*it)->setUserData(new ElementInfos(**it));
        }
        return mesh;
    }

    template<std::size_t DIM>
    std::unique_ptr<Base::MeshManipulator<DIM>> readMesh(std::string fileName, Base::ConfigurationData* configData)
    {
        auto mesh = std::unique_ptr<Base::MeshManipulator<DIM>>(
                new Base::MeshManipulator<DIM>(configData, 2, 3, 1, 1));
        mesh->readMesh(fileName);
        for (typename Base::MeshManipulator<DIM>::ElementIterator it = mesh->elementColBegin(Base::IteratorType::GLOBAL);
             it != mesh->elementColEnd(Base::IteratorType::GLOBAL); ++it)
        {
            (*it)->setUserData(new ElementInfos(**it));
        }
        return mesh;
    }


    // Explicit instantiation of the 2,3D versions.
    template
    std::unique_ptr<Base::MeshManipulator<2>> createCubeMesh(std::size_t subdivisions, Base::ConfigurationData* configData);
    template
    std::unique_ptr<Base::MeshManipulator<3>> createCubeMesh(std::size_t subdivisions, Base::ConfigurationData* configData);
    template
    std::unique_ptr<Base::MeshManipulator<2>> readMesh(std::string fileName, Base::ConfigurationData* configData);
    template
    std::unique_ptr<Base::MeshManipulator<3>> readMesh(std::string fileName, Base::ConfigurationData* configData);


}
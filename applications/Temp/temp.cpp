
#include "Base/MeshManipulator.h"
#include "Base/CommandLineOptions.h"
#include "Geometry/Mappings/RefinementMapsForCube.h"

auto& meshName = Base::register_argument<std::string>(
    'm', "meshName", "name of the Mesh file", true);

int main(int argc, char* argv[]) {
    Base::parse_options(argc, argv);

    auto* configData = new Base::ConfigurationData(1, 0);
    auto* mesh = new Base::MeshManipulator<2>(configData, 1, 0, 0, 0);
    auto* refinementMapping = Geometry::RefinementMapForSquare3::instance();

    mesh->readMesh(meshName.getValue());
    auto& temp = mesh->getElementsList(Base::IteratorType::GLOBAL);
    for (Base::Element* ptrElement : temp) {
        std::cout << *ptrElement << std::endl;
    }

    std::cout << "Refining" << std::endl;

    mesh->refine(refinementMapping);
    for (auto it = mesh->elementColBegin(Base::IteratorType::GLOBAL);
         it != mesh->elementColEnd(Base::IteratorType::GLOBAL); ++it) {
        if (it.getTreeEntry()->hasChild())
            continue;
        std::cout << *(*it) << std::endl;
    }

    delete configData;
    delete mesh;

    return 0;
}
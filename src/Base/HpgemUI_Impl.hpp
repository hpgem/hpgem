#include "Base/Element.hpp"
#include "Integration/ReturnTrait1.hpp"

namespace Base
{
    template<unsigned int DIM>
    class HpgemUI;

    template<unsigned int DIM>
    HpgemUI<DIM>::HpgemUI(const GlobalData* global, const ConfigurationData* config):
        meshes_(),
        globalData_(global),
        configData_(config)
    {
    }
    
    template<unsigned int DIM>
    HpgemUI<DIM>::~HpgemUI()
    {
            for(int i = 0; i < meshes_.size() ; ++i)
                delete meshes_[i];
    }
    
    template<unsigned int DIM>
    bool
    HpgemUI<DIM>::initialiseMeshMover(const MeshMoverBaseT* meshMoverBase, unsigned int meshID)
    {
        meshes_[meshID]->setMeshMover(meshMoverBase);
        return true;
    }
    
    
    template<unsigned int DIM>
    typename HpgemUI<DIM>::MeshId
    HpgemUI<DIM>::addMesh(const RectangularMeshDescriptorT& meshDscr, const MeshType& meshType)
    {
        unsigned int numOfMeshes=meshes_.size();
        MeshManipulator<DIM>* mesh = new MeshManipulator<DIM>(configData_);
        
        if (meshType== RECTANGULAR)
        {   
            mesh->createRectangularMesh(meshDscr.bottomLeft_, meshDscr.topLeft_, meshDscr.numElementsInDIM_);
            meshes_.push_back(mesh);
        }
        else
        {
            std::cerr << "Other types are yet to be implemented! " << std::endl;
        }
        cout<<"I just created a mesh!!!"<<endl;
        mesh->outputMesh(std::cout);
        return numOfMeshes;
    }

}

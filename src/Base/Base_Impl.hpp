#include "Base/Element.hpp"
#include "Integration/ReturnTrait1.hpp"

namespace Base
{
    template<unsigned int DIM>
    class Base;



    template<unsigned int DIM>
    bool
    Base<DIM>::initialiseMeshMover(MeshMoverBaseT* meshMoverBase, int meshID=0)
    {
        meshes_[meshID]->setMeshMover(meshMoverBase);
        return true;
    }
    
    
    template<unsigned int DIM>
    void Base<DIM>::addMesh(std::string type, PointPhysicalT BottomLeft, PointPhysicalT TopRight, std::vector<unsigned int> linearNoElements)
    {
        int numOfMeshes=meshes_.size();
        MeshManipulator<DIM>* mesh = new MeshManipulator<DIM>();
        
        if (type=="Rectangular")
        {   
            mesh->createRectangularMesh(BottomLeft, TopRight, linearNoElements);
            meshes_.push_back(mesh);
        }
        else
        {
            std::cerr << "Error in addmesh command : Unknown mesh type " << type << std::endl;
        }
        cout<<"I just created a mesh!!!"<<endl;
        mesh->outputMesh(std::cout);
    }

}

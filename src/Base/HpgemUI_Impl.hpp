#include "Base/Element.hpp"
#include "Integration/ReturnTrait1.hpp"

namespace Base
{
    template<unsigned int DIM>
    class HpgemUI;

    template<unsigned int DIM>
    HpgemUI<DIM>::HpgemUI(GlobalData* const global, const ConfigurationData* config):
        meshes_(),
        globalData_(global),
        configData_(config)
    {
    }
    
    template<unsigned int DIM>
    HpgemUI<DIM>::~HpgemUI()
    {
        cout << "is this called yet?"<<endl;
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
        MeshManipulator<DIM>* mesh = new MeshManipulator<DIM>(configData_,1,1,1,2);
        
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
            //mesh->outputMesh(std::cout);
        return numOfMeshes;
    }
    
    template<unsigned int DIM>
    typename HpgemUI<DIM>::ConstElementIterator
    HpgemUI<DIM>::elementColBegin(MeshId mId)const
    {
        meshes_[mId]->elementColBegin();
    }
    template<unsigned int DIM>
    typename HpgemUI<DIM>::ConstElementIterator
    HpgemUI<DIM>::elementColEnd(MeshId mId)const
    {
        meshes_[mId]->elementColEnd();
    }
    template<unsigned int DIM>
    typename HpgemUI<DIM>::ElementIterator
    HpgemUI<DIM>::elementColBegin(MeshId mId)
    {
        meshes_[mId]->elementColBegin();
    }
    template<unsigned int DIM>
    typename HpgemUI<DIM>::ElementIterator
    HpgemUI<DIM>::elementColEnd(MeshId mId)
    {
        meshes_[mId]->elementColEnd();
    }
    
    template<unsigned int DIM>
    typename HpgemUI<DIM>::ConstFaceIterator
    HpgemUI<DIM>::faceColBegin(MeshId mId)const
    {
        meshes_[mId]->faceColBegin();
    }
    template<unsigned int DIM>
    typename HpgemUI<DIM>::ConstFaceIterator
    HpgemUI<DIM>::faceColEnd(MeshId mId)const
    {
        meshes_[mId]->faceColEnd();
    }
    template<unsigned int DIM>
    typename HpgemUI<DIM>::FaceIterator
    HpgemUI<DIM>::faceColBegin(MeshId mId)
    {
        meshes_[mId]->faceColBegin();
    }
    template<unsigned int DIM>
    typename HpgemUI<DIM>::FaceIterator
    HpgemUI<DIM>::faceColEnd(MeshId mId)
    {
        meshes_[mId]->faceColEnd();
    }

}

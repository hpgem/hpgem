#include "HpgemUI.hpp"

#include "Base/Element.hpp"
#include "Integration/ReturnTrait1.hpp"

namespace Base
{

    HpgemUI::HpgemUI(GlobalData* const global, const ConfigurationData* config):
        meshes_(),
        globalData_(global),
        configData_(config)
    {
    }
    
    HpgemUI::~HpgemUI()
    {
        cout << "is this called yet?"<<endl;
            for(int i = 0; i < meshes_.size() ; ++i)
                delete meshes_[i];
    }
    
    bool
    HpgemUI::initialiseMeshMover(const MeshMoverBaseT* meshMoverBase, unsigned int meshID)
    {
        meshes_[meshID]->setMeshMover(meshMoverBase);
        return true;
    }
    
    
    typename HpgemUI::MeshId
    HpgemUI::addMesh(const RectangularMeshDescriptorT& meshDscr, const MeshType& meshType,int nrOfElementMatrixes, int nrOfElementVectors, int nrOfFaceMatrixes, int nrOfFaceVectors)
    {
        unsigned int numOfMeshes=meshes_.size();
        MeshManipulator* mesh = new MeshManipulator(configData_,meshDscr.boundaryConditions_[0],meshDscr.boundaryConditions_[1],meshDscr.boundaryConditions_[2],configData_->polynomialOrder_,0,nrOfElementMatrixes,nrOfElementVectors,nrOfFaceMatrixes,nrOfFaceVectors);//
        
        if (meshType== RECTANGULAR)
        {   
            mesh->createRectangularMesh(meshDscr.bottomLeft_, meshDscr.topRight_, meshDscr.numElementsInDIM_);
            meshes_.push_back(mesh);
        }
        else if (meshType==TRIANGULAR)
        {
        	mesh->createTriangularMesh(meshDscr.bottomLeft_,meshDscr.topRight_,meshDscr.numElementsInDIM_);
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
    
    typename HpgemUI::ConstElementIterator
    HpgemUI::elementColBegin(MeshId mId)const
    {
        return meshes_[mId]->elementColBegin();
    }

    typename HpgemUI::ConstElementIterator
    HpgemUI::elementColEnd(MeshId mId)const
    {
        return meshes_[mId]->elementColEnd();
    }

    typename HpgemUI::ElementIterator
    HpgemUI::elementColBegin(MeshId mId)
    {
        return meshes_[mId]->elementColBegin();
    }

    typename HpgemUI::ElementIterator
    HpgemUI::elementColEnd(MeshId mId)
    {
        return meshes_[mId]->elementColEnd();
    }
    
    typename HpgemUI::ConstFaceIterator
    HpgemUI::faceColBegin(MeshId mId)const
    {
        return meshes_[mId]->faceColBegin();
    }

    typename HpgemUI::ConstFaceIterator
    HpgemUI::faceColEnd(MeshId mId)const
    {
        return meshes_[mId]->faceColEnd();
    }

    typename HpgemUI::FaceIterator
    HpgemUI::faceColBegin(MeshId mId)
    {
        return meshes_[mId]->faceColBegin();
    }

    typename HpgemUI::FaceIterator
    HpgemUI::faceColEnd(MeshId mId)
    {
        return meshes_[mId]->faceColEnd();
    }

}

#include"Integration/ReturnTrait1.hpp"
#include"Base/Element.hpp"

namespace Base
{
    class HpgemUISimplified;
    
    HpgemUISimplified::HpgemUISimplified(unsigned int DIM):HpgemUI(new GlobalData,new ConfigurationData(DIM,1,11,1)){}
    
    bool
    HpgemUISimplified::solve()
    {
        initialise();
        checkInitialisation();
        for (int i=0; i < HpgemUI::meshes_.size();i++)
        {
            HpgemUI::meshes_[i]->move(); // just for testing
        }
        doAllElementIntegration();
        doAllFaceIntegration();
        
        return true;
    }
    template<class T>
    struct ABC
    {
        void aas(T& a){}
    };
    
    void
    HpgemUISimplified::doAllFaceIntegration(unsigned int meshID)
    {
        bool useCache   = false;
        
        LinearAlgebra::Matrix fMatrixData;
        LinearAlgebra::NumericalVector fVectorData;
        FaceIntegralT   faceIntegral(useCache);
        
        for (MeshManipulator::FaceIterator citFe = Base::HpgemUI::faceColBegin(); citFe != Base::HpgemUI::faceColEnd(); ++citFe)
        {
        	int n=(*citFe)->getPtrElementLeft()->getNrOfUnknows()*(*citFe)->getPtrElementLeft()->getNrOfBasisFunctions();
        	if((*citFe)->isInternal())
        		n+=(*citFe)->getPtrElementRight()->getNrOfUnknows()*(*citFe)->getPtrElementRight()->getNrOfBasisFunctions();
        	fMatrixData.resize(n,n);
            faceIntegral.integrate<LinearAlgebra::Matrix>((*citFe), this, fMatrixData);
            (*citFe)->setFaceMatrix(fMatrixData);
            faceIntegral.integrate<LinearAlgebra::NumericalVector>((*citFe), this, fVectorData);
            (*citFe)->setFaceVector(fVectorData);
        }
    }
    
    void
    HpgemUISimplified::doAllElementIntegration(unsigned int meshID)
    {
        unsigned int ndof = HpgemUI::configData_->numberOfBasisFunctions_*HpgemUI::configData_->numberOfUnknowns_;
        LinearAlgebra::Matrix  	eMatrixData(ndof, ndof);
        LinearAlgebra::NumericalVector eVectorData(ndof);
        
        bool isUseCache(false);
        Integration::ElementIntegral 	elIntegral(isUseCache);
        
        for (ElementIterator it=HpgemUI::meshes_[meshID]->elementColBegin(); it!= HpgemUI::meshes_[meshID]->elementColEnd(); ++it)
        {
            elIntegral.integrate<LinearAlgebra::Matrix>((*it), this, eMatrixData);
            (*it)->setElementMatrix(eMatrixData);
            elIntegral.integrate<LinearAlgebra::NumericalVector>((*it), this, eVectorData);
            (*it)->setElementVector(eVectorData);

            //cout << result;
            //cout<< "#####################################END of ELEMENT######"<<endl;
        }
    }
    
    bool
    HpgemUISimplified::checkInitialisation()
    {
        if (HpgemUI::meshes_.size()==0)
        {

            std::cerr << "Error no mesh created : You need to create at least one mesh to solve a problem" << std::endl;

        }
    }
}

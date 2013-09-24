#include"Integration/ReturnTrait1.hpp"
#include"Base/Element.hpp"

namespace Base
{
    template<unsigned int DIM>
    class HpgemUISimplified;
    
    
    template<unsigned int DIM>
    bool
    HpgemUISimplified<DIM>::solve()
    {
        initialise();
        checkInitialisation();
        for (int i=0; i < HpgemUI<DIM>::meshes_.size();i++)
        {
            HpgemUI<DIM>::meshes_[i]->move(); // just for testing
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
    
    template<unsigned int DIM>
    void
    HpgemUISimplified<DIM>::doAllFaceIntegration(unsigned int meshID)
    {
        bool useCache   = false;
        unsigned int nb = 1;//put something in here
        
        
        typedef void  (HpgemUISimplified::*FaceIntegrand)(const FaceT*, const PointPhysicalT& normal, const PointReferenceOnTheFaceT&, LinearAlgebra::Matrix&);
        
        LinearAlgebra::NumericalVector fData(nb);//do not know if this is the one you want
        FaceIntegrand faceInteg = &HpgemUISimplified<DIM>::faceIntegrand;
        FaceIntegralT   faceIntegral(useCache);
        
        for (ConstFaceIterator citFe = Base::HpgemUI<DIM>::faceColBegin(); citFe != Base::HpgemUI<DIM>::faceColEnd(); ++citFe)
        {
            
            
            faceIntegral.integrate((*citFe), faceInteg, fData, this);
            
        }
    }
    
    template<unsigned int DIM>
    void
    HpgemUISimplified<DIM>::doAllElementIntegration(unsigned int meshID)
    {
        unsigned int ndof = HpgemUI<DIM>::configData_->numberOfBasisFunctions_;
        LinearAlgebra::Matrix  	eData(ndof, ndof);
        typedef void  (HpgemUISimplified<DIM>::*Function)(const Element<DIM>& , const PointReferenceT&, LinearAlgebra::Matrix&);
        Function f = &HpgemUISimplified<DIM>::elementIntegrand;
        
        
        bool isUseCache(false);
        Integration::ElementIntegral<DIM> 	elIntegral(isUseCache);
        
        for (ElementIterator it=HpgemUI<DIM>::meshes_[meshID]->elementColBegin(); it!= HpgemUI<DIM>::meshes_[meshID]->elementColEnd(); ++it)
        {
            
            elIntegral.integrate((*it), f, eData, this);
            //cout << result;
            
            cout<< "#####################################END of ELEMENT######"<<endl;
        }
    }
    
    template<unsigned int DIM>
    bool
    HpgemUISimplified<DIM>::checkInitialisation()
    {
        if (HpgemUI<DIM>::meshes_.size()==0)
        {

            std::cerr << "Error no mesh created : You need to create at least one mesh to solve a problem" << std::endl;

        }
    }
}

#include"Integration/ReturnTrait1.hpp"
#include"Base/Element.hpp"

namespace Base
{
    template<unsigned int DIM>
    class BaseSimplified;
    
    template<unsigned int DIM>
    bool
    BaseSimplified<DIM>::solve()
    {
        initialise();
        checkInitialisation();
        for (int i=0; i<Base<DIM>::meshes_.size();i++)
        {
            Base<DIM>::meshes_[i]->move(); // just for testing
        }
        doAllElementIntegration();
        
        return true;
    }
    template<class T>
    struct ABC
    {
        void aas(T& a){}
    };
    
    template<unsigned int DIM>
    void
    BaseSimplified<DIM>::doAllElementIntegration(unsigned int meshID)
    {
        LinearAlgebra::Matrix  	matrix(1,1);
        typedef void  (BaseSimplified<DIM>::*Function)(const Element<DIM>& , const PointReferenceT&, LinearAlgebra::Matrix&);
        Function f = &BaseSimplified<DIM>::elementIntegrand;
        
     //void  (BaseSimplified<DIM>::*f)(const Element<DIM>& , const Geometry::PointReference<DIM>&, LinearAlgebra::NumericalVector&)= &BaseSimplified<DIM>::elementIntegrand;
     //typedef  void (Base<DIM>::*Function)(const Element<DIM>& , const Geometry::PointReference<DIM>&, LinearAlgebra::NumericalVector&);
     //Function f= &Base<DIM>::elementIntegrand;
        
     //   void  (*f)(const Element<DIM>& , const Geometry::PointReference<DIM>&, LinearAlgebra::NumericalVector&) = &ElementIntegrand;
     //typename Integration::ReturnTrait1<Function>::ReturnType result;
        
        bool isUseCache(false);
        Integration::ElementIntegral<DIM> 	elIntegral(isUseCache);
        
        for (ElementIterator it=Base<DIM>::meshes_[meshID]->elementColBegin(); it!= Base<DIM>::meshes_[meshID]->elementColEnd(); ++it)
        {
            
            elIntegral.integrate(*(*it), f, matrix, this);
            //cout << result;
            
            cout<< "#####################################END of ELEMENT######"<<endl;
        }
    }
    
    template<unsigned int DIM>
    bool
    BaseSimplified<DIM>::checkInitialisation()
    {
        if (Base<DIM>::meshes_.size()==0)
        {

            std::cerr << "Error no mesh created : You need to create at least one mesh to solve a problem" << std::endl;

        }
    }
}

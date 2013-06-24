#ifndef BaseSimplifiedHPP
#define BaseSimplifiedHPP

#include "Base/HpgemUI.hpp"

namespace Base
{
    template <unsigned int DIM>
    class HpgemUISimplified : public HpgemUI<DIM>
    {
    
    public:
        
        typedef Base::Element<DIM>                                      ElementT;
        typedef Base::Face<DIM>                                         FaceT;
        typedef Geometry::PointPhysical<DIM>                            PointPhysicalT;
        typedef Geometry::PointReference<DIM>                           PointReferenceT;
        typedef typename MeshManipulator<DIM>::ConstElementIterator     ConstElementIterator;
        typedef typename MeshManipulator<DIM>::ElementIterator          ElementIterator;
        

        /// You need the basis functions before creating the mesh, because the mesh manipulator
        /// needs to create the elements, and they need some basis functions variables.
    public:

        /// \brief Where the user creates a mesh
        bool virtual initialise()=0;
        
        /// \brief User-defined element integrand
        virtual void elementIntegrand(const ElementT& element, const PointReferenceT& p, LinearAlgebra::Matrix& ret)=0;
    
        /// \brief User-defined face integrand
        virtual void faceIntegrand(const PointPhysicalT& normal, 
                                   const PointReferenceT& p,  LinearAlgebra::Matrix& ret)=0;
        
        /// \brief User-defined initial conditions
        virtual void initialConditions(const PointPhysicalT& p)=0;

        /// \brief Integrates, and other things.
        bool solve();
        
        ///Preforms all the element integrations
        void doAllElementIntegration(unsigned int meshID=0);


    private:
        
        //This is a function that checks the users defined initisation is fine.
        bool  checkInitialisation();
        
        
    };
}
#include "HpgemUISimplified_Imp.hpp"

#endif

#ifndef BaseSimplifiedHPP
#define BaseSimplifiedHPP

#include "Base/HpgemUI.hpp"
#include "Integration/ElementIntegrandBase.hpp"
#include "Integration/FaceIntegrandBase.hpp"
//#include "Base/MeshManipulator.hpp"

namespace Base
{
    class HpgemUISimplified : public HpgemUI,Integration::ElementIntegrandBase<LinearAlgebra::Matrix>,Integration::FaceIntegrandBase<LinearAlgebra::Matrix>,
    						Integration::FaceIntegrandBase<LinearAlgebra::NumericalVector>,Integration::ElementIntegrandBase<LinearAlgebra::NumericalVector>
    {
    
    public:
        
        typedef typename MeshManipulator::ConstElementIterator     ConstElementIterator;
        typedef typename MeshManipulator::ElementIterator          ElementIterator;
        typedef typename MeshManipulator::ConstFaceIterator        ConstFaceIterator;
        
        typedef Base::Element                                     ElementT;
        typedef Base::Face                                         FaceT;
        typedef Geometry::PointPhysical                            PointPhysicalT;
        typedef Geometry::PointReference                           PointReferenceT;
        typedef Geometry::PointReference                         PointReferenceOnTheFaceT;
        typedef Integration::FaceIntegral                          FaceIntegralT;
        

        /// You need the basis functions before creating the mesh, because the mesh manipulator
        /// needs to create the elements, and they need some basis functions variables.
    public:

        HpgemUISimplified(unsigned int DIMension);

        /// \brief Where the user creates a mesh
        bool virtual initialise()=0;
        
        /// \brief User-defined element integrand
        virtual void elementIntegrand(const ElementT* element, const PointReferenceT& p, LinearAlgebra::Matrix& ret)=0;

        virtual void elementIntegrand(const ElementT* element, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret)=0;
    
        /// \brief User-defined face integrand
        virtual void faceIntegrand(const FaceT* face, const PointPhysicalT& normal,
                                   const PointReferenceT& p,  LinearAlgebra::Matrix& ret)=0;
        
        virtual void faceIntegrand(const FaceT* face, const PointPhysicalT& normal,
                				   const PointReferenceT& p,  LinearAlgebra::NumericalVector& ret)=0;

        /// \brief User-defined initial conditions
        virtual void initialConditions(const PointPhysicalT& p)=0;

        /// \brief Integrates, and other things.
        bool solve();
        
        ///Preforms all the element integrations
        void doAllElementIntegration(unsigned int meshID=0);
        void doAllFaceIntegration(unsigned int meshID=0);


    private:
        
        //This is a function that checks the users defined initisation is fine.
        bool  checkInitialisation();
        
        
    };
}
#include "HpgemUISimplified_Imp.hpp"

#endif

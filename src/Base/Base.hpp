#include "MeshMoverBase.hpp"
#include "MeshManipulator.hpp"
#include "Output/TecplotDiscontinuousSolutionWriter.hpp"

namespace Base
{
    template <unsigned int DIM>
    class Base
    {
        
    public:
        typedef MeshManipulator<DIM>    MeshManipulatorT;
        typedef MeshMoverBase<DIM>      MeshMoverBaseT;
        typedef Geometry::PointReference<DIM-1>     PointReferenceT;
        typedef Geometry::PointPhysical<DIM>        PointPhysicalT;
        //typedef BasisFunctions<DIM>     BasisFunctionT;
        


        /// \note You need the basis functions before creating the mesh, because the mesh manipulator
        /// needs to create the elements, and they need some basis functions variables.
    public:
        
        
        Base():
            mesh_()
        {std::cout << "What the fuck" <<std::endl;}

        virtual ~Base() {};

        /// \brief Gives the pointer of meshMoverBase class to mesh.
        bool initialiseMeshMover(MeshMoverBaseT* meshMoverBase);

        /// \brief Creates mesh.
        bool virtual initialise()=0;
    
        /// \brief User-defined element integrand
        virtual void elementIntegrand(const PointReferenceT& p, LinearAlgebra::NumericalVector& ret)=0;
    
        /// \brief User-defined face integrand
        virtual void faceIntegrand(const PointPhysicalT& normal, 
                                   const PointReferenceT& p,  LinearAlgebra::NumericalVector& ret)=0;

        /// \brief Integrates, and other things.
        bool solve();

        /// \brief Virtual function that should be overwritten by specific problem, specifies initial conditions.
        //virtual void initialCondition() const;

    protected:
        MeshManipulatorT mesh_;
            //MeshMoverBaseT* meshMoverBase_;
//            BasisFunctionT basisFunction_;
//            GlobalData globalData_;
//            ConfigurationData configurationData_;
        
    private:
        
        //This is a function that checks the users defined initisation is fine.
        bool  checkInitialisation();
    };
};
#include "Base_Impl.hpp"

#include "MeshMoverBase.hpp"
#include "MeshManipulator.hpp"
#include "Output/TecplotDiscontinuousSolutionWriter.hpp"

namespace Base
{
    template <unsigned int DIM>
    class Base
    {
        typedef MeshManipulator<DIM>    MeshManipulatorT;
        typedef MeshMoverBase<DIM>      MeshMoverBaseT;
        //typedef BasisFunctions<DIM>     BasisFunctionT;

        /// \note You need the basis functions before creating the mesh, because the mesh manipulator
        /// needs to create the elements, and they need some basis functions variables.
        public:

            Base():
                mesh_()
            {}

            virtual ~Base() {};

            /// \brief Gives the pointer of meshMoverBase class to mesh.
            bool initialiseMeshMover(MeshMoverBaseT* meshMoverBase);

            /// \brief Creates mesh.
            bool virtual initialiseMesh()=0;

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
    };
};
#include "Base_Impl.hpp"

#ifndef BASE_HPP
#define BASE_HPP

#include "Base/MeshMoverBase.hpp"
#include "Base/MeshManipulator.hpp"
#include "Base/GlobalNamespaceBase.hpp"
#include "Base/RectangularMeshDescriptor.hpp"
#include "Output/TecplotDiscontinuousSolutionWriter.hpp"
#include "Integration/ElementIntegral.hpp"

#include "vector"

namespace Base
{
    template <unsigned int DIM>
    class HpgemUI
    {
    public:
        typedef typename MeshManipulator<DIM>::ConstElementIterator     ConstElementIterator;
        typedef typename MeshManipulator<DIM>::ElementIterator          ElementIterator;
        
    public:
        typedef unsigned int                                            MeshId;
        typedef std::vector<unsigned int>                               VectorOfUIntegers;
        typedef RectangularMeshDescriptor<DIM>                          RectangularMeshDescriptorT;
        typedef MeshManipulator<DIM>                                    MeshManipulatorT;
        typedef std::vector<MeshManipulatorT* >                         VectorOfMeshManipulatorT;
        typedef std::string                                             String;
        typedef MeshMoverBase<DIM>                                      MeshMoverBaseT;
        typedef Geometry::PointPhysical<DIM>                            PointPhysicalT;
     
        
        //typedef BasisFunctions<DIM>     BasisFunctionT;
        

    public:
        
        
        HpgemUI(GlobalData* global, ConfigurationData* config);

        virtual ~HpgemUI() ;

        /// \brief Gives the pointer of meshMoverBase class to mesh.
        virtual bool initialiseMeshMover(const MeshMoverBaseT* meshMoverBase, unsigned int meshID);

            /// Creating a mesh with in-house remesher.
        virtual MeshId addMesh(const RectangularMeshDescriptorT& meshDescriptor, const MeshType& meshType = RECTANGULAR);
            /// Reading a mesh from a file, currently only Centaur is supported.
        virtual MeshId addMesh(const String& fileName){}


        /// \brief Virtual function that should be overwritten by specific problem, specifies initial conditions.
        //virtual void initialCondition() const;

    protected:
        VectorOfMeshManipulatorT                meshes_;
       
        const GlobalData*                       globalData_;
        const ConfigurationData*                configData_;
    };
};
#include "HpgemUI_Impl.hpp"
#endif

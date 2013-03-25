#ifndef BASE_HPP
#define BASE_HPP

#include "MeshMoverBase.hpp"
#include "MeshManipulator.hpp"
#include "Output/TecplotDiscontinuousSolutionWriter.hpp"
#include "Integration/ElementIntegral.hpp"

#include "vector"

namespace Base
{
    template <unsigned int DIM>
    class Base
    {
    public:
        typedef typename MeshManipulator<DIM>::ConstElementIterator     ConstElementIterator;
        typedef typename MeshManipulator<DIM>::ElementIterator          ElementIterator;
        
    public:
        
        typedef std::vector<MeshManipulator<DIM>* > VectorOfMeshManipulatorT;
        typedef MeshMoverBase<DIM>                  MeshMoverBaseT;
        typedef Geometry::PointPhysical<DIM>        PointPhysicalT;
     
        
        //typedef BasisFunctions<DIM>     BasisFunctionT;
        


        /// You need the basis functions before creating the mesh, because the mesh manipulator
        /// needs to create the elements, and they need some basis functions variables.
    public:
        
        
        Base():meshes_()
        {}

        virtual ~Base() 
        {
            for(int i = 0; i < meshes_.size() ; ++i)
                delete meshes_[i];
        }

        /// \brief Gives the pointer of meshMoverBase class to mesh.
        bool initialiseMeshMover(MeshMoverBaseT* meshMoverBase, int meshID);

        
        void addMesh(std::string type, PointPhysicalT bottomLeft, PointPhysicalT topRight, std::vector<unsigned int> LinearNoElements);

        /// \brief Virtual function that should be overwritten by specific problem, specifies initial conditions.
        //virtual void initialCondition() const;

    protected:
        /// \todo change this to a vector of meshes.
        VectorOfMeshManipulatorT                meshes_;
            //MeshMoverBaseT* meshMoverBase_;
//            BasisFunctionT basisFunction_;
//            GlobalData globalData_;
//            ConfigurationData configurationData_;
        
        
        
        
    };
};
#include "Base_Impl.hpp"
#endif

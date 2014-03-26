#ifndef BASE_HPP
#define BASE_HPP

#include "Base/MeshMoverBase.hpp"
#include "Base/MeshManipulator.hpp"
#include "Base/GlobalNamespaceBase.hpp"
#include "Base/RectangularMeshDescriptor.hpp"
#include "Output/TecplotDiscontinuousSolutionWriter.hpp"
#include "Integration/ElementIntegral.hpp"
#include "Integration/FaceIntegral.hpp"

#include "vector"

namespace Base
{
    class HpgemUI
    {
    public:
        typedef typename MeshManipulator::ConstElementIterator     ConstElementIterator;
        typedef typename MeshManipulator::ElementIterator          ElementIterator;
        typedef typename MeshManipulator::ConstFaceIterator        ConstFaceIterator;
        typedef typename MeshManipulator::FaceIterator             FaceIterator;
        
    public:
        typedef Base::Element                                      ElementT;
        typedef Base::Face                                         FaceT;
        typedef RectangularMeshDescriptor                          RectangularMeshDescriptorT;
        typedef MeshManipulator                                    MeshManipulatorT;
        typedef MeshMoverBase                                      MeshMoverBaseT;
        typedef Geometry::PointPhysical                            PointPhysicalT;
        typedef Geometry::PointReference                           PointReferenceT;

        typedef unsigned int                                            MeshId;
        typedef std::vector<unsigned int>                               VectorOfUIntegers;
       
        typedef std::vector<MeshManipulatorT* >                         VectorOfMeshManipulatorT;
        typedef std::string                                             String;
                
        //typedef BasisFunctions<DIM>     BasisFunctionT;
        

    public:
             
        HpgemUI(GlobalData* const global, const ConfigurationData* config);

        virtual ~HpgemUI() ;
        
        

        /// \brief Gives the pointer of meshMoverBase class to mesh.
        virtual bool initialiseMeshMover(const MeshMoverBaseT* meshMoverBase, unsigned int meshID);

            /// Creating a mesh with in-house remesher.
        MeshId addMesh(const RectangularMeshDescriptorT& meshDescriptor, const MeshType& meshType = RECTANGULAR, int nrOfElementMatrixes=0, int nrOfElementVectors=0,int nrOfFaceMatrixes=0, int nrOfFaceVectors=0);
            /// Reading a mesh from a file, currently only Centaur is supported.
        MeshId addMesh(const String& fileName){}
        
        unsigned int getNumberOfElements(MeshId id)const {return meshes_[id]->getNumberOfElements();}
        
        ConstElementIterator    elementColBegin(MeshId mId=0)const;
        ConstElementIterator    elementColEnd(MeshId mId=0)const;
        
        ElementIterator         elementColBegin(MeshId mId=0);
        ElementIterator         elementColEnd(MeshId mId=0);
        
        ConstFaceIterator       faceColBegin(MeshId mId=0)const;
        ConstFaceIterator       faceColEnd(MeshId mId=0)const;
        
        FaceIterator            faceColBegin(MeshId mId=0);
        FaceIterator            faceColEnd(MeshId mId=0);
        

        /// \brief Virtual function that should be overwritten by specific problem, specifies initial conditions.
        //virtual void initialCondition() const;

    protected:
        VectorOfMeshManipulatorT                meshes_;
       
        GlobalData* const                       globalData_;
        const ConfigurationData* const          configData_;
    };
};
#endif

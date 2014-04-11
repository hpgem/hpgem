#ifndef MESHMANIPULATOR_H_
#define MESHMANIPULATOR_H_


#include "Geometry/GlobalNamespaceGeometry.hpp"
#include "Geometry/ReferenceGeometry.hpp"
#include "Geometry/PhysicalGeometry.hpp"
#include "Geometry/PointPhysical.hpp"
#include "Geometry/RefinementGeometry.hpp"

#include "Base/BasisFunctionSet.hpp"
#include "Base/AssembleBasisFunctionSet.hpp"
#include "OrientedBasisFunctionSet.hpp"
#include "Base/Face.hpp"
#include "Base/MeshMoverBase.hpp"
#include "Base/LevelTree.hpp"
#include "Base/GlobalData.hpp"
#include "Base/Element.hpp"
#include "Base/ConfigurationData.hpp"

#include "Integration/QuadratureRules/AllGaussQuadratureRules.hpp"



#include <vector>
#include <list>
#include <fstream>

namespace Base
{
    struct HalfFaceDescription {
        std::vector<unsigned int> nodeList;
        unsigned int elementNum;
        unsigned int localFaceIndex;
    };
    
    
    class MeshManipulator //: public MeshRefiner <DIM>
    {
    public:
        
        typedef unsigned int                                PointIndexT;
        typedef Element                                ElementT;
        typedef Face                                   FaceT;
        typedef MeshMoverBase                          MeshMoverBaseT;
        typedef Geometry::PointPhysical                PointPhysicalT;
        typedef Base::BasisFunctionSet                 BasisFunctionSetT;
        
        typedef LevelTree<ElementT>                         ElementLevelTreeT;
        typedef LevelTree<FaceT>                            FaceLevelTreeT;
        
        typedef std::list<FaceT*>                           ListOfFacesT;
        typedef std::list<ElementT*>                        ListOfElementsT;
        typedef std::vector<ElementT* >                     VectorOfElementPtrT;
        typedef std::vector<PointPhysicalT >                VectorOfPhysicalPointsT;
        typedef std::vector<PointIndexT>                    VectorOfPointIndicesT;
        typedef std::vector<const BasisFunctionSetT*>       CollectionOfBasisFunctionSets;
        typedef std::vector<ElementLevelTreeT*>             VecOfElementLevelTreePtrT;
        typedef std::vector<FaceLevelTreeT*>                VecOfFaceLevelTreePtrT;

        
        typedef typename ElementLevelTreeT::iterator        ElementIteratorT;
        typedef typename FaceLevelTreeT::iterator           FaceIteratorT;
        
        typedef typename ListOfElementsT::const_iterator    ConstElementIterator;
        typedef typename ListOfElementsT::iterator          ElementIterator;
        
        typedef typename ListOfFacesT::const_iterator       ConstFaceIterator;
        typedef typename ListOfFacesT::iterator             FaceIterator;

        
    public:
            /// idRangeBegin is the begining of the range, from where the Element's ids should be assigned.
            /// In case of multiple meshes, one has to take care of empty intersection of those ranges!!!
        MeshManipulator(const ConfigurationData* configData, bool xPer=0, bool yPer=0, bool zPer=0, unsigned int orderOfFEM=1, unsigned int idRangeBegin=0,
        		int nrOfElementMatrixes=0, int nrOfElementVectors=0, int nrOfFaceMatrixes=0, int nrOfFaceVectors=0);

        MeshManipulator(const MeshManipulator& other);

        virtual ~MeshManipulator();

  

        void                            createDefaultBasisFunctions(unsigned int order);

        ElementT*                       addElement(const VectorOfPointIndicesT& globalNodeIndexes);

        bool                            addFace(ElementT* leftElementPtr, unsigned int leftElementLocalFaceNo, 
                                                ElementT* rightElementPtr, unsigned int rightElementLocalFaceNo,
                                                const Geometry::FaceType& faceType=Geometry::WALL_BC);
        
        void                            addEdge(std::vector< Element*> elements, std::vector<unsigned int> localEdgeNrs);

        unsigned int                    getNumberOfElements(unsigned int meshId=0) const {return elements_.size(); }
        unsigned int                    getNumberOfFaces(unsigned int meshId=0) const {return faces_.size();}
        unsigned int                    getNumberOfEdges(unsigned int meshId=0) const {return edges_.size();}
        unsigned int                    getNumberOfNodes()const {return points_.size();}

        /// *****************Iteration through the Elements*******************
        ConstElementIterator            elementColBegin()const{return elements_.begin();}

        ConstElementIterator            elementColEnd()const{return elements_.end();}

        ElementIterator                 elementColBegin(){return elements_.begin();}
        ElementIterator                 elementColEnd(){return elements_.end();}
        
        ConstFaceIterator               faceColBegin()const{return faces_.begin();}
        
        ConstFaceIterator               faceColEnd()const{return faces_.end();}
        
        FaceIterator                    faceColBegin(){return faces_.begin();}
        FaceIterator                    faceColEnd(){return faces_.end();}

       std::list< Edge*>::const_iterator edgeColBegin()const{return edges_.begin();}

       std::list< Edge*>::const_iterator edgeColEnd()const{return edges_.end();}

        std::list< Edge*>::iterator      edgeColBegin(){return edges_.begin();}
        std::list< Edge*>::iterator      edgeColEnd(){return edges_.end();}
        /// *****************Iteration through the Elements*******************

        void                            createRectangularMesh(const PointPhysicalT& BottomLeft, const PointPhysicalT& TopRight, const VectorOfPointIndicesT& LinearNoElements);

	/**
	 * Crates a mesh of simplices for the specified cube
	 * \param [in] BottomLeft the bottomleft corner of the cube
	 * \param [in] TopRight The topRight corner of the cube
	 * \param [in] LinearNoElements A vector detailing the amount of refinement you want per direction
	 * This routine generates the same mesh structure as createRectangularMesh, but then refines each of the cubes into
	 * (DIM-1)^2+1 tetrahedra
	 */
        void                            createTriangularMesh(PointPhysicalT BottomLeft, PointPhysicalT TopRight, const VectorOfPointIndicesT& LinearNoElements);
	
        void                            readCentaurMesh(const std::string& filename);

        void                            outputMesh(ostream& os)const;


        //! Set MeshMoverBase object pointer, for moving meshes if needed
        void                            setMeshMover(const MeshMoverBaseT* const meshMover);

        void                            move();

        // ******************THESE SHOULD BE DELETED LATER***********************//actually, these should be replaced by iterable editions of the levelTree -FB
        //! Get const list of elements
        const ListOfElementsT&          getElementsList() const {return elements_; }

        //! Get non-const list of elements
        ListOfElementsT&                getElementsList() { return elements_; }

        //! Get const list of faces
        const ListOfFacesT&             getFacesList() const { return faces_; }

        //! Get non-const list of faces
        ListOfFacesT&                   getFacesList() { return faces_; }

        const std::list<Edge*>&         getEdgesList() const {return edges_;}

        std::list<Edge*>&               getEdgesList() {return edges_;}
        // ************************************************************************

  //---------------------------------------------------------------------
        //! Get the number of mesh-tree.
        int                             getNumberOfMeshes() const;
        
        //! Create a new (empty) mesh-tree.
        void                            createNewMeshTree();
        
        //! Get the element container of a specific mesh-tree.
        ElementLevelTreeT*              ElCont(int meshTreeIdx) const;
        
        //! Get the face container of a specific mesh-tree.
        FaceLevelTreeT*                 FaCont(int meshTreeIdx) const;

        //! Some mesh generator: centaur / rectangular / triangle / tetrahedra / triangular-prism.
        void                            someMeshGenerator(int meshTreeIdx);
        
        //! Set active mesh-tree.
        void                            setActiveMeshTree(unsigned int meshTreeIdx);
        
        //! Get active mesh-tree index.
        int                             getActiveMeshTree() const;

        //! Reset active mesh-tree.
        void                            resetActiveMeshTree();
        
        //! Get maximum h-level of a specific mesh-tree.
        unsigned int                    getMaxLevel(int meshTreeIdx) const;

        //! Set active level of a specific mesh-tree.
        void                            setActiveLevel(unsigned int meshTreeIdx, int level);
        
        //! Get active level of a specific mesh-tree.
        int                             getActiveLevel(int meshTreeIdx) const;
        
        //! Reset active level of a specific mesh-tree.
        void                            resetActiveLevel(int meshTreeIdx);
        
        //! Duplicate mesh contents including all refined meshes.
        void                            duplicate(unsigned int fromMeshTreeIdx, unsigned int toMeshTreeIdx, unsigned int upToLevel);

        //! Refine a specific mesh-tree.
        void                            doRefinement(unsigned int meshTreeIdx, int refinementType);

        //! Changes the default set of basisFunctions for this mesh and all of its elements. Ignores any conforming basisFunctionset that nay be linked to faces/edges/...
        void							setDefaultBasisFunctionSet(BasisFunctionSetT* bFSet);

        //! Adds vertex based degrees of freedom to the set of basisfunctions for this mesh and all of its vertices. This routine will assume that the first basisfunctionset connects to the first vertex (local numbering) and so on
        void                            addVertexBasisFunctionSet(CollectionOfBasisFunctionSets& bFsets);///\TODO support for mixed meshes

        //! Adds face based degrees of freedom to the set of basisfunctions for this mesh and all of its faces. This routine will assume that all needed orientations are available in the collection of basisfunctionsets
        void                            addFaceBasisFunctionSet(std::vector<const OrientedBasisFunctionSet*>& bFsets);///\TODO support for mixed meshes

        //! Adds edge based degrees of freedom to the set of basisfunctions for this mesh and all of its edges. This routine will assume that all needed orientations are available in the collection of basisfunctionsets
        void                            addEdgeBasisFunctionSet(std::vector<const OrientedBasisFunctionSet*>& bFsets);///\TODO support for mixed meshes

        int                             dimension(){return configData_->dimension_;}

        const std::vector<PointPhysicalT>& getNodes(){return points_;}
  //---------------------------------------------------------------------
    private:
        
        //!Does the actual reading for 2D centaur meshes
        void                            readCentaurMesh2D(std::ifstream& centaurFile);
        
        //!Does the actual reading for 3D centaur meshes
        void                            readCentaurMesh3D(std::ifstream& centaurFile);
	
        void                            faceFactory();
        void                            edgeFactory();
	
	//someone thinks its a good idea to declare HalfFaceDescription in an implemetnation file
	void                            findElementNumber(std::list<int>& a, std::list<int>& b, std::list<int>& c,int aNumber, int bNumber, int cNumber, std::list<int>& notOnFace, HalfFaceDescription& face, std::vector<Element*>& vectorOfElements);
	
	//!An alternative to faceFactory that only iterates over internal faces, use the boundary face information in the centaur file to construct the boundary faces
	void                            constructInternalFaces(std::vector<std::list<int> >& listOfElementsForEachNode, std::vector<Element*>& vectorOfElements);
        
        void                            rectangularCreateFaces1D(VectorOfElementPtrT& tempElementVector, const VectorOfPointIndicesT& linearNoElements);
        
        void                            rectangularCreateFaces2D(VectorOfElementPtrT& tempElementVector, const VectorOfPointIndicesT& linearNoElements);
        
        void                            rectangularCreateFaces3D(VectorOfElementPtrT& tempElementVector, const VectorOfPointIndicesT& linearNoElements);
 
	//! Make faces for the 1D contrarotating mesh
        void                            triangularCreateFaces1D(VectorOfElementPtrT& tempElementVector, const VectorOfPointIndicesT& linearNoElements);
 
	//!Make faces for the 2D triangular mesh
        void                            triangularCreateFaces2D(VectorOfElementPtrT& tempElementVector, const VectorOfPointIndicesT& linearNoElements);
 
	//!Make faces for the 3D tetrahedral mesh
        void                            triangularCreateFaces3D(VectorOfElementPtrT& tempElementVector, const VectorOfPointIndicesT& linearNoElements);
        
        //! Do refinement on the elements.
        void                            doElementRefinement(unsigned int meshTreeIdx);
        
        //! Do refinement on the faces.
        void                            doFaceRefinement(unsigned int meshTreeIdx);
        
        //! Check whether the two elements may be connected by a face or not.
        void                            pairingCheck(const ElementIterator elL, unsigned int locFaceNrL,
                                                     const ElementIterator elR, unsigned int locFaceNrR,
                                                     int& pairingValue, bool& sizeOrder);
                          
        //! Check whether the two elements may be connected by a face or not in periodic face case.
        void                            periodicPairingCheck(const FaceIteratorT fa, 
                                                             const ElementIteratorT elL, unsigned int localFaceNrL,
                                                             const ElementIteratorT elR, unsigned int localFaceNrR,
                                                             int& pairingValue, bool& sizeOrder);
  //---------------------------------------------------------------------
    private:
        
        const ConfigurationData*        configData_;

        unsigned int                    elementcounter_;

        static unsigned int             faceCounter_;
        static unsigned int             edgeCounter_;
        //! List of all elements. TODO: this should be replaced by the mesh-tree structure
        ListOfElementsT                 elements_;
        
        //! List of all faces. TODO: this should be replaced by the mesh-tree structure
        ListOfFacesT                    faces_;
        
        //! List of all edges.
        std::list< Edge*>                edges_;

        //! Global vector of physical nodes.
        VectorOfPhysicalPointsT         points_;

        //! Periodicity in x-direction.
        bool                            periodicX_;
        
        //! Periodicity in y-direction.
        bool                            periodicY_;
        
        //! Periodicity in z-direction.
        bool                            periodicZ_;
        
        /// Pointer to MeshMoverBase, in order to move points in the mesh, when needed by user.
        const MeshMoverBaseT*           meshMover_;

        CollectionOfBasisFunctionSets   collBasisFSet_;
        
        //const BasisFunctionSetT*        defaultSetOfBasisFunctions_;
        //! Collection of additional basis function set, if p-refinement is applied
        
        //! Active mesh-tree.
        int                             activeMeshTree_;
        
        //! Number of mesh-tree.
        int                             numMeshTree_;
    
        //! Vector elements LevelTree.
        VecOfElementLevelTreePtrT       vecOfElementTree_;
        
        //! Vector faces LevelTree.
        VecOfFaceLevelTreePtrT          vecOfFaceTree_;

        int                             numberOfElementMatrixes_;
        int                             numberOfFaceMatrixes_;
        int                             numberOfElementVectors_;
        int                             numberOfFaceVectors_;
    };
    
    
}





#endif /* MESHMANIPULATOR_H_ */

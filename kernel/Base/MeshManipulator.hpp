/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, University of Twente
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef MESHMANIPULATOR_H_
#define MESHMANIPULATOR_H_


#include <vector>
#include <list>
#include <fstream>

//for enum support...
#include "Geometry/FaceGeometry.hpp"
#include "FaceCacheData.hpp"
#include "Mesh.hpp"

namespace Geometry{
	class PointPhysical;
	class PointReference;
}

namespace Base
{
	class BasisFunctionSet;
	class OrientedBasisFunctionSet;
	class Face;
	class MeshMoverBase;
	template<class V>
	class LevelTree;
	class Element;
	class ConfigurationData;
	class Edge;

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

        
        //typedef typename ElementLevelTreeT::iterator        ElementIteratorT;
        //typedef typename FaceLevelTreeT::iterator           FaceIteratorT;
        
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

        unsigned int                    getNumberOfElements(unsigned int meshId=0) const {return theMesh_.getNumberOfElements(meshId);}
        unsigned int                    getNumberOfFaces(unsigned int meshId=0) const {return theMesh_.getNumberOfFaces(meshId);}
        unsigned int                    getNumberOfEdges(unsigned int meshId=0) const {return theMesh_.getNumberOfEdges(meshId);}
        unsigned int                    getNumberOfNodes()const {return theMesh_.getNumberOfNodes();}

        /// *****************Iteration through the Elements*******************
        ConstElementIterator            elementColBegin()const{return theMesh_.elementColBegin();}
        ConstElementIterator            elementColEnd()const{return theMesh_.elementColEnd();}

        ElementIterator                 elementColBegin(){return theMesh_.elementColBegin();}
        ElementIterator                 elementColEnd(){return theMesh_.elementColEnd();}
        
        ConstFaceIterator               faceColBegin()const{return theMesh_.faceColBegin();}
        ConstFaceIterator               faceColEnd()const{return theMesh_.faceColEnd();}
        
        FaceIterator                    faceColBegin(){return theMesh_.faceColBegin();}
        FaceIterator                    faceColEnd(){return theMesh_.faceColEnd();}

       std::list< Edge*>::const_iterator edgeColBegin()const{return theMesh_.edgeColBegin();}
       std::list< Edge*>::const_iterator edgeColEnd()const{return theMesh_.edgeColEnd();}

        std::list< Edge*>::iterator      edgeColBegin(){return theMesh_.edgeColBegin();}
        std::list< Edge*>::iterator      edgeColEnd(){return theMesh_.edgeColEnd();}
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

        void                            outputMesh(std::ostream& os)const;


        //! Set MeshMoverBase object pointer, for moving meshes if needed
        void                            setMeshMover(const MeshMoverBaseT* const meshMover);

        void                            move();

        // ******************THESE SHOULD BE DELETED LATER***********************//actually, these should be replaced by iterable editions of the levelTree -FB
        //! Get const list of elements
        const ListOfElementsT&          getElementsList() const {return theMesh_.getElementsList(); }

        //! Get non-const list of elements
        ListOfElementsT&                getElementsList() { return theMesh_.getElementsList(); }

        //! Get const list of faces
        const ListOfFacesT&             getFacesList() const { return theMesh_.getFacesList(); }

        //! Get non-const list of faces
        ListOfFacesT&                   getFacesList() { return theMesh_.getFacesList(); }

        const std::list<Edge*>&         getEdgesList() const {return theMesh_.getEdgesList();}

        std::list<Edge*>&               getEdgesList() {return theMesh_.getEdgesList();}
        // ************************************************************************

        //! Changes the default set of basisFunctions for this mesh and all of its elements. Ignores any conforming basisFunctionset that nay be linked to faces/edges/...
        void							setDefaultBasisFunctionSet(BasisFunctionSetT* bFSet);

        //! Adds vertex based degrees of freedom to the set of basisfunctions for this mesh and all of its vertices. This routine will assume that the first basisfunctionset connects to the first vertex (local numbering) and so on
        void                            addVertexBasisFunctionSet(CollectionOfBasisFunctionSets& bFsets);///\TODO support for mixed meshes

        //! Adds face based degrees of freedom to the set of basisfunctions for this mesh and all of its faces. This routine will assume that all needed orientations are available in the collection of basisfunctionsets
        void                            addFaceBasisFunctionSet(std::vector<const OrientedBasisFunctionSet*>& bFsets);///\TODO support for mixed meshes

        //! Adds edge based degrees of freedom to the set of basisfunctions for this mesh and all of its edges. This routine will assume that all needed orientations are available in the collection of basisfunctionsets
        void                            addEdgeBasisFunctionSet(std::vector<const OrientedBasisFunctionSet*>& bFsets);///\TODO support for mixed meshes

		int dimension();

        const std::vector<PointPhysicalT>& getNodes()const{return theMesh_.getNodes();}

        //routines that deal with level trees
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
        //void                            periodicPairingCheck(const FaceIteratorT fa,
        //                                                     const ElementIteratorT elL, unsigned int localFaceNrL,
        //                                                     const ElementIteratorT elR, unsigned int localFaceNrR,
        //                                                     int& pairingValue, bool& sizeOrder);
  //---------------------------------------------------------------------
    private:
        
        Mesh                            theMesh_;
        
        const ConfigurationData*        configData_;
        //! Periodicity in x-direction.
        bool                            periodicX_;
        
        //! Periodicity in y-direction.
        bool                            periodicY_;
        
        //! Periodicity in z-direction.
        bool                            periodicZ_;
        
        /// Pointer to MeshMoverBase, in order to move points in the mesh, when needed by user.
        const MeshMoverBaseT*           meshMover_;

        //! Collection of additional basis function set, if p-refinement is applied
        CollectionOfBasisFunctionSets   collBasisFSet_;
        
        //const BasisFunctionSetT*        defaultSetOfBasisFunctions_;
        
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

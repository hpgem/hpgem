#ifndef MESHMANIPULATOR_H_
#define MESHMANIPULATOR_H_

#include "Base/LevelTree.hpp"
#include "Base/Element.hpp"
#include "Geometry/GlobalNamespaceGeometry.hpp"
#include "Geometry/ReferenceGeometry.hpp"
#include "Geometry/PhysicalGeometry.hpp"
#include "Geometry/PointPhysical.hpp"
#include "Geometry/ReferenceSquare.hpp"
#include "Geometry/ReferenceCube.hpp"
#include "Geometry/ReferenceTriangle.hpp"
#include "Geometry/ReferenceLine.hpp"
#include "Geometry/RefinementGeometry.hpp"
#include "Base/BasisFunctionSet.hpp"
#include "Base/AssembleBasisFunctionSet.hpp"
#include "Base/Face.hpp"
#include "Base/MeshMoverBase.hpp"



//#include "Face.hpp"
#include <vector>
#include <list>
#include <fstream>

namespace Base
{
   

    

    template <unsigned int DIM>
    class MeshManipulator //: public MeshRefiner <DIM>
    {
    public:
        typedef Element<DIM>                                ElementT;
        typedef std::list<ElementT*>                        ListOfElementsT;
        typedef Face<DIM>                                   FaceT;
        typedef Base::BasisFunctionSet<DIM>                 BasisFunctionSetT;
        typedef std::list<FaceT>                            ListOfFacesT;
        typedef unsigned int                                PointIndexT;
        typedef Geometry::PointPhysical<DIM>                PointPhysicalT;
        typedef std::vector<PointPhysicalT >                VectorOfPhysicalPointsT;
        typedef std::vector<PointIndexT>                    VectorOfPointIndexesT;
        typedef std::vector<BasisFunctionSetT*>             CollectionOfBasisFunctionSets;

        typedef LevelTree<ElementT>                         ElementLevelTreeT;
        typedef LevelTree<FaceT>                            FaceLevelTreeT;
        typedef typename ElementLevelTreeT::iterator        ElementIteratorT;
        typedef typename FaceLevelTreeT::iterator           FaceIteratorT;
        typedef std::vector<ElementLevelTreeT*>             VecOfElementLevelTreePtrT;
        typedef std::vector<FaceLevelTreeT*>                VecOfFaceLevelTreePtrT;


    public:
        MeshManipulator(bool xPer=0, bool yPer=0, bool zPer=0):
            periodicX_(xPer),
            periodicY_(yPer),
            periodicZ_(zPer),
            activeMeshTree_(0), 
            numMeshTree_(0)
        {
            unsigned int order =1;
            createBasisFunctions(order);
            createNewMeshTree();
        }
        
        void createBasisFunctions(unsigned int order);
        
        virtual ~MeshManipulator()
        {
            for (typename ListOfElementsT::iterator it=elements_.begin(); it!=elements_.end();++it)
            {
                ElementT* el = *it;
                delete el;
            }
            
            for (typename CollectionOfBasisFunctionSets::iterator bit=collBasisFSet_.begin(); bit!=collBasisFSet_.end();++bit)
            {
                BasisFunctionSetT* bf = *bit;
                delete bf;
            }

            // Kill all faces in all mesh-tree
            while (!vecOfFaceTree_.empty())
            {
              delete vecOfFaceTree_.back();
              vecOfFaceTree_.pop_back();
            }

            // Kill all elements in all mesh-tree
            while (!vecOfElementTree_.empty())
            {
              delete vecOfElementTree_.back();
              vecOfElementTree_.pop_back();
            }
        }
        
        Base::Element<DIM>* addElement(const std::vector<unsigned int>& globalNodeIndexes,
                        const Geometry::ReferenceGeometry<DIM>* const referenceGeometry)
            {
                unsigned int numOfUnknowns      = 10;
                unsigned int numOfTimeLevels    = 10;
                unsigned int counter            = 10;
                
                const BasisFunctionSetT* const bf= collBasisFSet_[0];
                
                
                ElementT*  myElement = new ElementT(globalNodeIndexes,bf, points_, numOfUnknowns, numOfTimeLevels, counter);
                
                
                elements_.push_back(myElement);
                
                return myElement;
            }
        
        bool addFace(Base::Element<DIM>* leftElementPtr, unsigned int leftElementLocalFaceNo, Base::Element<DIM>* rightElementPtr, unsigned int rightElementLocalFaceNo, const Geometry::FaceType& faceType=Geometry::INTERNAL)
        {
             std::cout << "-----------------------\n";
            
             std::cout << "{Left Element" << (*leftElementPtr) <<" , FaceIndex=" <<leftElementLocalFaceNo<<"}"<<endl;
            if (rightElementPtr!=NULL)
            {
                 std::cout << "{Right Element" << (*rightElementPtr) << " , FaceIndex=" <<rightElementLocalFaceNo<<"}"<<endl;
                
                Face<DIM> face(leftElementPtr, leftElementLocalFaceNo, rightElementPtr, rightElementLocalFaceNo);
                
                faces_.push_back(face);
            }
            else 
            {
                 std::cout << "This is a boundary face" << std::endl;
                
                
               Face<DIM> bFace(leftElementPtr, leftElementLocalFaceNo, faceType);
                
                faces_.push_back(bFace);
            }
            
             std::cout << "-----------------------\n";

            
            
            
        }
        //bool readMesh(FILE file);
        
        //This now works in only 2D.
        void createRectangularMesh(PointPhysicalT BottomLeft, PointPhysicalT TopRight, std::vector<unsigned int> LinearNoElements);
        
        void readCentaurMesh(const std::string filename);
             
        void outputMesh(ostream& os)const
        {
            for (int i=0;i<points_.size();i++){os<<"Node " <<i<<" "<<points_[i]<<endl;}
            
            
            int elementNum=0;
            
            for (typename ListOfElementsT::const_iterator cit=elements_.begin(); cit !=elements_.end(); ++cit)
                {
                os << "Element " <<elementNum <<" " <<*(*cit)<<endl;
                elementNum++;
                }
            
            int faceNum=0;
            for (typename ListOfFacesT::const_iterator cit=faces_.begin(); cit !=faces_.end(); ++cit)
            {
                /// \bug need at output routine for the face to test this.
               // os << "Face " <<faceNum <<" " <<(*cit)<<endl;
                faceNum++;
            }
            
            
        }

        //! Set MeshMoverBase object pointer, for moving meshes if needed
        void setMeshMover(MeshMoverBase<DIM>* meshMover)
        {
            meshMover_ = meshMover;
        }

        void move()
        {
            for (int i=0;i<points_.size();i++)
            {
                cout << "Node " << i << " " << points_[i]; cout << endl;
                meshMover_->movePoint(&points_[i]);
                cout << "Node " << i << " " << points_[i]; cout << endl;
            }
        }

        //! Get const list of elements
        const ListOfElementsT& getElementsList() const { return elements_; }

        //! Get non-const list of elements
        ListOfElementsT& getElementsList() { return elements_; }

        //! Get const list of faces
        const ListOfFacesT& getFacesList() const { return faces_; }

        //! Get non-const list of faces
        ListOfFacesT& getFacesList() { return faces_; }


  //---------------------------------------------------------------------
        //! Get the number of mesh-tree.
        int size() const;
        
        //! Create a new (empty) mesh-tree.
        void createNewMeshTree();
        
        //! Get the element container of a specific mesh-tree.
        ElementLevelTreeT* ElCont(int meshTreeIdx) const;
        
        //! Get the face container of a specific mesh-tree.
        FaceLevelTreeT* FaCont(int meshTreeIdx) const;

        //! Some mesh generator: centaur / rectangular / triangle / tetrahedra / triangular-prism.
        void someMeshGenerator(int meshTreeIdx);
        
        //! Set active mesh-tree.
        void setActiveMeshTree(unsigned int meshTreeIdx);
        
        //! Get active mesh-tree index.
        int getActiveMeshTree() const;

        //! Reset active mesh-tree.
        void resetActiveMeshTree();
        
        //! Get maximum h-level of a specific mesh-tree.
        unsigned int getMaxLevel(int meshTreeIdx) const;

        //! Set active level of a specific mesh-tree.
        void setActiveLevel(unsigned int meshTreeIdx, int level);
        
        //! Get active level of a specific mesh-tree.
        int getActiveLevel(int meshTreeIdx) const;
        
        //! Reset active level of a specific mesh-tree.
        void resetActiveLevel(int meshTreeIdx);
        
        //! Duplicate mesh contents including all refined meshes.
        void duplicate(unsigned int fromMeshTreeIdx, unsigned int toMeshTreeIdx, unsigned int upToLevel);

        //! Refine a specific mesh-tree.
        void doRefinement(unsigned int meshTreeIdx, int refinementType);
  //---------------------------------------------------------------------
    

    private:
        
        //Does the actually reading for 2D centaur meshes
        void readCentaurMesh2D(std::ifstream& centaurFile);
        
        void faceFactory();
        
        void rectangularCreateFaces1D(std::vector<Base::Element<DIM>* >& tempElementVector, const std::vector<unsigned int>& linearNoElements);
        void rectangularCreateFaces2D(std::vector<Base::Element<DIM>* >& tempElementVector, const std::vector<unsigned int>& linearNoElements);
        void rectangularCreateFaces3D(std::vector<Base::Element<DIM>* >& tempElementVector, const std::vector<unsigned int>& linearNoElements);
        
  //---------------------------------------------------------------------
        //! Do refinement on the elements.
        void doElementRefinement(unsigned int meshTreeIdx);
        
        //! Do refinement on the faces.
        void doFaceRefinement(unsigned int meshTreeIdx);
        
        //! Check whether the two elements may be connected by a face or not.
        void pairingCheck(const ElementIteratorT elL, unsigned int locFaceNrL,
                          const ElementIteratorT elR, unsigned int locFaceNrR,
                          int& pairingValue, bool& sizeOrder);
                          
        //! Check whether the two elements may be connected by a face or not in periodic face case.
        void periodicPairingCheck(const FaceIteratorT fa, 
                          const ElementIteratorT elL, unsigned int localFaceNrL,
                          const ElementIteratorT elR, unsigned int localFaceNrR,
                          int& pairingValue, bool& sizeOrder);
  //---------------------------------------------------------------------
        
        //! List of all elements. TODO: this should be replaced by the mesh-tree structure
        ListOfElementsT                 elements_;
        
        //! List of all faces. TODO: this should be replaced by the mesh-tree structure
        ListOfFacesT                    faces_;
        
        //! Global vector of physical nodes.
        VectorOfPhysicalPointsT         points_;

        //! Periodicity in x-direction.
        bool                            periodicX_;
        
        //! Periodicity in y-direction.
        bool                            periodicY_;
        
        //! Periodicity in z-direction.
        bool                            periodicZ_;
        
        /// Pointer to MeshMoverBase, in order to move points in the mesh, when needed by user.
        MeshMoverBase<DIM>*             meshMover_;
        
        //! Collection of basis function set.
        CollectionOfBasisFunctionSets   collBasisFSet_;
        
        //! Active mesh-tree.
        int                             activeMeshTree_;
        
        //! Number of mesh-tree.
        int                             numMeshTree_;
    
        //! Vector elements LevelTree.
        VecOfElementLevelTreePtrT       vecOfElementTree_;
        
        //! Vector faces LevelTree.
        VecOfFaceLevelTreePtrT          vecOfFaceTree_;
    };
     
#include "MeshManipulator_Impl.hpp"

}



#endif /* MESHMANIPULATOR_H_ */

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

#include <stdlib.h>
#include <ctime>
#include <iostream>
#include <vector>
#include <typeinfo>
#include "Base/LevelTree.hpp"
#include "Base/TreeIterator.hpp"
#include "Base/TreeEntry.hpp"

//=========================================================== Element class definition
template <unsigned int dim>
class Element
{
  public:
    Element(int id);
    Element(const Element& other);
    ~Element();
    int getId();

    friend std::ostream& operator<<(std::ostream& os, const Element<dim>& e)
    {
        os<< e.id_ << " ";
        return os;
    }

  private:
    Element();
    int id_;
};

template<unsigned int dim>
Element<dim>::Element(int id): id_(id)  
{ 
//     std::cout << "Element " << id_ << " is created.\n"; 
}

template<unsigned int dim>
Element<dim>::Element(const Element& other): id_(other.id_) 
{ 
//     std::cout << "Element " << id_ << " is copied.\n"; 
}

template<unsigned int dim>
Element<dim>::~Element() 
{ 
//     std::cout << "Element " << id_ << " is destroyed.\n"; 
}

template<unsigned int dim>
int Element<dim>::getId() { return id_; }


//=========================================================== Face class definition
template <unsigned int dim>
class Face
{
  public:
    Face(int id);
    Face(const Face& other);
    ~Face();
    int getId();

    friend std::ostream& operator<<(std::ostream& os, const Face<dim>& e)
    {
        os<< e.id_ << " ";
        return os;
    }

  private:
    Face();
    int id_;
};

template<unsigned int dim>
Face<dim>::Face(int id): id_(id)  
{ 
//     std::cout << "Face " << id_ << " is created.\n"; 
}

template<unsigned int dim>
Face<dim>::Face(const Face& other): id_(other.id_) 
{ 
//     std::cout << "Face " << id_ << " is copied.\n"; 
}

template<unsigned int dim>
Face<dim>::~Face() 
{ 
//     std::cout << "Face " << id_ << " is destroyed.\n"; 
}

template<unsigned int dim>
int Face<dim>::getId() { return id_; }


//========================================================== MeshManipulator class definition
template <unsigned int dim>
class MeshManipulator
{
  public:
    typedef Element<dim>                            ElementT;
    typedef Face<dim>                               FaceT;
    typedef Base::LevelTree<ElementT>               ElementLevelTreeT;
    typedef Base::LevelTree<FaceT>                  FaceLevelTreeT;
    typedef typename ElementLevelTreeT::iterator    ElementIteratorT;
    typedef typename FaceLevelTreeT::iterator       FaceIteratorT;
    typedef std::vector<ElementLevelTreeT*>         VecOfElementLevelTreePtrT;
    typedef std::vector<FaceLevelTreeT*>            VecOfFaceLevelTreePtrT;
    
    MeshManipulator();
    ~MeshManipulator();

    //! Get the number of mesh-tree.
    int size() const;
    
    //! Create a new (empty) mesh-tree.
    void createNewMeshTree();
    
    //! Get the element container of a specific mesh-tree.
    ElementLevelTreeT* ElCont(int meshTreeIdx = -1) const;
    
    //! Get the face container of a specific mesh-tree.
    FaceLevelTreeT* FaCont(int meshTreeIdx = -1) const;

    //! Some mesh generator: centaur / rectangular / triangle / tetrahedra / triangular-prism.
    void someMeshGenerator(int meshTreeIdx = -1);
    
    //! Set active mesh-tree.
    void setActiveMeshTree(unsigned int meshTreeIdx);
    
    //! Get active mesh-tree index.
    unsigned int getActiveMeshTree() const;

    //! Reset active mesh-tree.
    void resetActiveMeshTree();
    
    //! Get maximum h-level of a specific mesh-tree.
    unsigned int getMaxLevel(int meshTreeIdx = -1) const;

    //! Set active level of a specific mesh-tree.
    void setActiveLevel(unsigned int meshTreeIdx, int level);
    
    //! Get active level of a specific mesh-tree.
    int getActiveLevel(int meshTreeIdx = -1) const;
    
    //! Reset active level of a specific mesh-tree.
    void resetActiveLevel();
    
    //! Duplicate mesh contents including all refined meshes.
    void duplicate(unsigned int fromMeshIdx, unsigned int toMeshIdx, unsigned int upToLevel = 0);

    //! Refine a specific mesh-tree.
    void doRefinement(unsigned int meshTreeIdx, int refinementType = -1);
    
  private:
    //! Do refinement on the elements.
    void doElementRefinement(unsigned int meshTreeIdx);
    
    //! Do refinement on the faces.
    void doFaceRefinement(unsigned int meshTreeIdx);
    
    //! Check whether the two elements may be connected by a face or not.
    void pairingCheck(const ElementIteratorT elL, unsigned int locFaceNrL,
                      const ElementIteratorT elR, unsigned int locFaceNrR,
                      int& pairingValue, bool& sizeOrder);
                      
    //! Check whether the two elements may be connected by a face or not in periodic face case.
    bool periodicPairingCheck(const FaceIteratorT fa, 
                      const ElementIteratorT elL, unsigned int localFaceNrL,
                      const ElementIteratorT elR, unsigned int localFaceNrR,
                      int& pairingValue, bool& sizeOrder);


  private:
    //! Active mesh-tree.
    int                         activeMeshTree_;
    
    //! Number of mesh-tree.
    int                         numMeshTree_;
    
    //! Vector elements LevelTree.
    VecOfElementLevelTreePtrT   vecOfElementTree_;
    
    //! Vector faces LevelTree.
    VecOfFaceLevelTreePtrT      vecOfFaceTree_;
};


template<unsigned int dim>
MeshManipulator<dim>::MeshManipulator() : activeMeshTree_(0), numMeshTree_(0)
{
    this->createNewMeshTree();
}

template<unsigned int dim>
MeshManipulator<dim>::~MeshManipulator()
{
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

template<unsigned int dim>
int MeshManipulator<dim>::size() const
{
  return vecOfElementTree_.size();
}

template<unsigned int dim>
void MeshManipulator<dim>::createNewMeshTree()
{
  vecOfElementTree_.push_back(new ElementLevelTreeT);
  vecOfFaceTree_.push_back(new FaceLevelTreeT);
  ++numMeshTree_;
}

template<unsigned int dim>
typename MeshManipulator<dim>::ElementLevelTreeT* MeshManipulator<dim>::ElCont(int meshTreeIdx) const
{
  if ((meshTreeIdx >= 0) && (meshTreeIdx < numMeshTree_))
  {
    return vecOfElementTree_[meshTreeIdx];
  }
  else if ((meshTreeIdx == -1) && (activeMeshTree_ >= 0))
  {
    return vecOfElementTree_[activeMeshTree_];
  }
  
  return NULL;
}

template<unsigned int dim>
typename MeshManipulator<dim>::FaceLevelTreeT* MeshManipulator<dim>::FaCont(int meshTreeIdx) const
{
  if ((meshTreeIdx >= 0) && (meshTreeIdx < numMeshTree_))
  {
    return vecOfFaceTree_[meshTreeIdx];
  }
  else if ((meshTreeIdx == -1) && (activeMeshTree_ >= 0))
  {
    return vecOfFaceTree_[activeMeshTree_];
  }
  
  return NULL;
}

template<unsigned int dim>
void MeshManipulator<dim>::someMeshGenerator(int meshTreeIdx)
{
  int createOnIndex;
  if ((meshTreeIdx >= 0) && (meshTreeIdx < numMeshTree_))
  {
    createOnIndex = meshTreeIdx;
  }
  else if ((meshTreeIdx == -1) && (activeMeshTree_ >= 0))
  {
    createOnIndex = activeMeshTree_;
  }
  else
    throw "MeshManipulator::someMeshGenerator: invalid mesh-tree index";
    
  const int numberOfElement = 1 + (rand() % 10);
  const int startElementId = (createOnIndex+1)*1000;
  for (int id=startElementId; id<startElementId+numberOfElement; ++id)
  {
    Element<dim> el(id);
    vecOfElementTree_[createOnIndex]->addEntry(el);
  }

  const int numberOfFace = 1 + (rand() % 10);
  const int startFaceId = (createOnIndex+1)*1000;
  for (int id=startFaceId; id<startFaceId+numberOfFace; ++id)
  {
    Face<dim> fa(id);
    vecOfFaceTree_[createOnIndex]->addEntry(fa);
  }
}


//===================================================== Displaying container contents
template <typename T>
void display(T& container)
{
  std::cout << "There are " << container.size() << " entries in the container.\n";
  std::cout << "Content IDs: ";
  for (typename T::iterator it=container.begin(); it!=container.end(); ++it)
    std::cout << (*it).getId() << " ";
  
  std::cout << std::endl;
}


//========================================================= Here is the main function
int main()
{
    const unsigned int dim = 2;
    
    long seed = std::time(NULL);
    srand ( seed );
//     std::cout << "random number seed = " << seed << std::endl;

    typedef Base::LevelTree<Element<dim> > LevelTreeT;
    typedef std::vector<LevelTreeT*>  vecLevelTreePtrT;
    
    std::cout << "\nInitial test: adding people into a LevelTree and display them.\n";
    std::cout << "----------------------------------------------------------------\n";
    
    Element<dim> Cathy(-1);
    Element<dim> Denny(-2);
    Element<dim> Victor(-3);
    LevelTreeT treeElement;
    
    std::cout << "Before any people are added...";
    display<LevelTreeT>(treeElement);
    treeElement.addEntry(Cathy);
    treeElement.addEntry(Denny);
    treeElement.addEntry(Victor);
    std::cout << "After three people are added...";
    display<LevelTreeT>(treeElement);

    std::cout << "\n\n";
    std::cout << "Multiple LevelTree test: adding LevelTree(s), populate and display their contents.\n";
    std::cout << "----------------------------------------------------------------------------------\n";
    
    // number of LevelTrees
    const int numOfLevelTrees = 5;
    
    // create vector of LevelTrees
    vecLevelTreePtrT vecOfLevelTree;
    for (int i=0; i<numOfLevelTrees; ++i)
    {
      vecOfLevelTree.push_back(new LevelTreeT);
    }

    for (int i=0; i<vecOfLevelTree.size(); ++i)
    {
        // working with LevelTree-i
        LevelTreeT& myTree = *(vecOfLevelTree[i]);

        const int numberOfElement = 1 + (rand() % 10);
        const int startId = (i+1)*1000;
        for (int i=startId; i<startId+numberOfElement; ++i)
        {
          Element<dim> student(i);
          myTree.addEntry(student);
        }
        display<LevelTreeT>(myTree);
    }

    // Destroy the vector of LevelTree before quiting.
//     std::cout << "Deleting vector of LevelTree....\n";
    while (!vecOfLevelTree.empty())
    {
      delete vecOfLevelTree.back();
      vecOfLevelTree.pop_back();
    }
//     std::cout << "vector of LevelTree is deleted.\n";
    
    std::cout << "\n\n";
    std::cout << "Overall test: working with multiple mesh-tree.\n";
    std::cout << "----------------------------------------------------------------------------------\n";

    int meshTreeIdx;

    MeshManipulator<dim> myMesh;
    std::cout << "Number of mesh-tree: " << myMesh.size() << std::endl;
    myMesh.someMeshGenerator();

    meshTreeIdx = 0;
    std::cout << "Number of elements in mesh-tree-" << meshTreeIdx << " is " << myMesh.ElCont(meshTreeIdx)->size() << std::endl;
    for (MeshManipulator<dim>::ElementIteratorT it=myMesh.ElCont()->begin(); it!=myMesh.ElCont()->end(); ++it)
    {
      std::cout << (*it).getId() << " ";
    }
    std::cout << std::endl;

    myMesh.createNewMeshTree();
    meshTreeIdx = 1;
    myMesh.someMeshGenerator(meshTreeIdx);
    std::cout << "Number of elements in mesh-tree-" << meshTreeIdx << " is " << myMesh.ElCont(meshTreeIdx)->size() << std::endl;
    for (MeshManipulator<dim>::ElementIteratorT it=myMesh.ElCont(meshTreeIdx)->begin(); 
                      it!=myMesh.ElCont(meshTreeIdx)->end(); ++it)
    {
      std::cout << (*it).getId() << " ";
    }
    std::cout << std::endl;

    return 0;
}

#include <stdlib.h>
#include <ctime>
#include <iostream>
#include <vector>
#include <typeinfo>
#include "LevelTree.hpp"


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
    typedef typename ElementLevelTreeT::iterator    ElementIterator;
    typedef typename FaceLevelTreeT::iterator       FaceIterator;
    typedef std::vector<ElementLevelTreeT*>         VecOfElementLevelTreePtrT;
    typedef std::vector<FaceLevelTreeT*>            VecOfFaceLevelTreePtrT;
    
    MeshManipulator() 
    {
        this->createNewMeshTree();
    }
    
    ~MeshManipulator()
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

    int size() const
    {
      return vecOfElementTree_.size();
    }
    
    void createNewMeshTree()
    {
      vecOfElementTree_.push_back(new ElementLevelTreeT);
      vecOfFaceTree_.push_back(new FaceLevelTreeT);
    }

    ElementLevelTreeT* ElCont(unsigned int meshTreeIdx=0)
    {
      return vecOfElementTree_[meshTreeIdx];
    }

    FaceLevelTreeT* FaCont(unsigned int meshTreeIdx=0)
    {
      return vecOfFaceTree_[meshTreeIdx];
    }

    void someMeshGenerator(unsigned int meshTreeIdx = 0)
    {
      if (meshTreeIdx >= this->size())
        throw "MeshManipulator::someMeshGenerator: invalid mesh-tree index";
        
      const int numberOfElement = 1 + (rand() % 10);
      const int startElementId = (meshTreeIdx+1)*1000;
      for (int id=startElementId; id<startElementId+numberOfElement; ++id)
      {
        Element<dim> el(id);
        vecOfElementTree_[meshTreeIdx]->addEntry(el);
      }

      const int numberOfFace = 1 + (rand() % 10);
      const int startFaceId = (meshTreeIdx+1)*1000;
      for (int id=startFaceId; id<startFaceId+numberOfFace; ++id)
      {
        Face<dim> fa(id);
        vecOfFaceTree_[meshTreeIdx]->addEntry(fa);
      }
    }

  private:
    VecOfElementLevelTreePtrT   vecOfElementTree_;
    VecOfFaceLevelTreePtrT      vecOfFaceTree_;
};



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
    std::cout << "Number of elements in mesh-tree-" << meshTreeIdx << " is " << myMesh.size() << std::endl;
    for (MeshManipulator<dim>::ElementIterator it=myMesh.ElCont()->begin(); it!=myMesh.ElCont()->end(); ++it)
    {
      std::cout << (*it).getId() << " ";
    }
    std::cout << std::endl;

    myMesh.createNewMeshTree();
    meshTreeIdx = 1;
    myMesh.someMeshGenerator();
    std::cout << "Number of elements in mesh-tree-" << meshTreeIdx << " is " << myMesh.size() << std::endl;
    for (MeshManipulator<dim>::ElementIterator it=myMesh.ElCont(meshTreeIdx)->begin(); 
                      it!=myMesh.ElCont(meshTreeIdx)->end(); ++it)
    {
      std::cout << (*it).getId() << " ";
    }
    std::cout << std::endl;
    
    return 0;
}

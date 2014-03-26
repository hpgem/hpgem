#include <stdlib.h>
#include <ctime>
#include <iostream>
#include <vector>
#include <typeinfo>
#include "LevelTree.hpp"


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

    typedef Base::LevelTree<Base::Element > LevelTreeT;
    typedef std::vector<LevelTreeT*>  vecLevelTreePtrT;
    
    std::cout << "\nInitial test: adding people into a LevelTree and display them.\n";
    std::cout << "----------------------------------------------------------------\n";
    
    Base::Element Cathy(-1);
    Base::Element Denny(-2);
    Base::Element Victor(-3);
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
          Base::Element student(i);
          myTree.addEntry(student);
        }
        display<LevelTreeT>(myTree);
    }

    // Destroy the vector of LevelTree before quiting.
    while (!vecOfLevelTree.empty())
    {
      delete vecOfLevelTree.back();
      vecOfLevelTree.pop_back();
    }
    
    std::cout << "\n\n";
    std::cout << "Overall test: working with multiple mesh-tree.\n";
    std::cout << "----------------------------------------------------------------------------------\n";
    
    int meshTreeIdx;
    
    Base::MeshManipulator myMesh;
    std::cout << "Number of mesh-tree: " << myMesh.size() << std::endl;
    myMesh.someMeshGenerator(0);
    
    meshTreeIdx = 0;
    std::cout << "Number of elements in mesh-tree-" << meshTreeIdx << " is " << myMesh.size() << std::endl;
    for (Base::MeshManipulator::ElementIterator it=myMesh.ElCont()->begin(); it!=myMesh.ElCont()->end(); ++it)
    {
      std::cout << (*it).getId() << " ";
    }
    std::cout << std::endl;

    myMesh.createNewMeshTree();
    meshTreeIdx = 1;
    myMesh.someMeshGenerator(0);
    std::cout << "Number of elements in mesh-tree-" << meshTreeIdx << " is " << myMesh.size() << std::endl;
    for (Base::MeshManipulator::ElementIterator it=myMesh.ElCont(meshTreeIdx)->begin();
                      it!=myMesh.ElCont(meshTreeIdx)->end(); ++it)
    {
      std::cout << (*it).getId() << " ";
    }
    std::cout << std::endl;
    
    return 0;
}

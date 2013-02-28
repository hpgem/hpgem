#include "Base/MeshManipulator.hpp"
#include "Geometry/PointPhysical.hpp"
#include <vector>
#include <iostream>

using namespace std;


int main()
{
    //first 1D test
    {//
//        Geometry::PointPhysical<1> bottomLeft, topLeft;
//        std::vector<unsigned int> numElementsOneD(1);
//    
//        bottomLeft[0]=0;
//    
//        topLeft[0]=2;
//    
//        numElementsOneD[0]=2;
//    
//        Base::MeshManipulator<1> myOneDDemoMesh(1);
//    
//    
//    
//        myOneDDemoMesh.createRectangularMesh(bottomLeft,topLeft,numElementsOneD);
//    
//        myOneDDemoMesh.outputMesh(std::cout);
        
    }
    cout<<"asdasdas"<<endl;    //Now 2D test
    {
      //  Geometry::PointPhysical<2> bottomLeft, topLeft;
//        std::vector<unsigned int> numElementsTwoD(2);
//    
//        bottomLeft[0]=0;
//        bottomLeft[1]=0;
//    
//        topLeft[0]=2;
//        topLeft[1]=2;
//    
//        numElementsTwoD[0]=30;
//        numElementsTwoD[1]=20;
//    
//    
//        Base::MeshManipulator<2> myTwoDDemoMesh(1,1);
//    
//        myTwoDDemoMesh.createRectangularMesh(bottomLeft,topLeft,numElementsTwoD);
//    
//        myTwoDDemoMesh.outputMesh(std::cout);

    }
    
    //Now 3D test
    {
        Geometry::PointPhysical<3> bottomLeft, topRight;
        std::vector<unsigned int> numElementsThreeD(3);
  
        bottomLeft[0]=0;
        bottomLeft[1]=0;
        bottomLeft[2]=0;
        
        topRight[0]=2;
        topRight[1]=2;
        topRight[2]=2;
        
        numElementsThreeD[0]=3;
        numElementsThreeD[1]=3;
        numElementsThreeD[2]=1;

        Base::MeshManipulator<3> myThreeDDemoMesh(1,1,1);
        myThreeDDemoMesh.createRectangularMesh(bottomLeft,topRight,numElementsThreeD);
    
        myThreeDDemoMesh.outputMesh(std::cout);
    
    }

}

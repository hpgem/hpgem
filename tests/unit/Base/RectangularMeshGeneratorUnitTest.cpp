#include "Base/MeshManipulator.hpp"
#include "Geometry/PointPhysical.hpp"
#include <vector>
#include <iostream>

using namespace std;


int main()
{
    //first 1D test
//    {//
//        
//        Base::ConfigurationData config(1,1,1);
//        
//        config.numberOfUnknowns_       = 1;
//        config.numberOfTimeLevels_     = 1;
//        config.numberOfBasisFunctions_ = 1;
//
//        Geometry::PointPhysical<1> bottomLeft, topLeft;
//        std::vector<unsigned int> numElementsOneD(1);
//
//        bottomLeft[0]=0;
//
//        topLeft[0]=1;
//    
//        numElementsOneD[0]=10;
//    
//        Base::MeshManipulator<1> myOneDDemoMesh(&config, 1);
//    
//    
//    
//        myOneDDemoMesh.createRectangularMesh(bottomLeft,topLeft,numElementsOneD);
//    
//        myOneDDemoMesh.outputMesh(std::cout);
        
//    }
      //Now 2D test
    {
        
//        Base::ConfigurationData config(1,1,1);
//        
//        Geometry::PointPhysical<2> bottomLeft, topLeft;
//        std::vector<unsigned int> numElementsTwoD(2);
//    
//        bottomLeft[0]=0;
//        bottomLeft[1]=0;
//    
//        topLeft[0]=1;
//        topLeft[1]=1;
//    
//        numElementsTwoD[0]=2;
//        numElementsTwoD[1]=2;
//    
//    
//        Base::MeshManipulator<2> myTwoDDemoMesh(&config, 1,1);
//    
//        myTwoDDemoMesh.createRectangularMesh(bottomLeft,topLeft,numElementsTwoD);
//    
//        myTwoDDemoMesh.outputMesh(std::cout);

    }
    
    //Now 3D test
    {
        Base::ConfigurationData config(3,1,1,1);
            //
            //        config.numberOfUnknowns_       = 1;
            //        config.numberOfTimeLevels_     = 1;
            //        config.numberOfBasisFunctions_ = 1;
        
         Geometry::PointPhysical bottomLeft(3), topRight(3);
         std::vector<unsigned int> numElementsThreeD(3);
   
         bottomLeft[0]=0;
         bottomLeft[1]=0;
         bottomLeft[2]=0;
         
         topRight[0]=1;
         topRight[1]=1;
         topRight[2]=1;
         
         numElementsThreeD[0]=3;
         numElementsThreeD[1]=3;
         numElementsThreeD[2]=3;
 
         Base::MeshManipulator myThreeDDemoMesh(&config, 1,1,1);
         myThreeDDemoMesh.createRectangularMesh(bottomLeft,topRight,numElementsThreeD);
     
         myThreeDDemoMesh.outputMesh(std::cout);
//
    }

}

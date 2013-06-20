#include "Base/MeshManipulator.hpp"
#include "Geometry/PointPhysical.hpp"
#include <vector>
#include "CMakeDefinitions.hpp"
#include <sstream>
#include "Output/TecplotDiscontinuousSolutionWriter.hpp"

using namespace std;


int main()
{
    Base::ConfigurationData config;
    
    config.numberOfUnknowns_       = 1;
    config.numberOfTimeLevels_     = 1;
    config.numberOfBasisFunctions_ = 1;
    
    Base::MeshManipulator<2> myTwoDDemoMesh(&config, 1,1);
    
    std::stringstream filename;
    
    filename << getCMAKE_hpGEM_SOURCE_DIR() << "/tests/files/centaurQuadMinimum.hyb";

    myTwoDDemoMesh.readCentaurMesh(filename.str());
    
    myTwoDDemoMesh.outputMesh(std::cout);
    
    
    std::ofstream file2D;
    file2D.open ("SavedCentaurQuadMinimum.dat");
    
    int dimensionsToWrite[2] = {0,1};
    
    Output::TecplotDiscontinuousSolutionWriter<2> out(file2D,"QuadMinimum Test Mesh",dimensionsToWrite,"xy");
    out.write(&myTwoDDemoMesh,"holi",false);
    
    file2D.close();
    
    //Now do it again with a more complicated mesh
    
    Base::MeshManipulator<2> triQuadTwoDDemoMesh(&config,1,1);
    
    std::stringstream filename2;
    
    filename2 << getCMAKE_hpGEM_SOURCE_DIR() << "/tests/files/centaurMinimum.hyb";
    
    cout << "filename " << filename2.str() <<endl;
    
    triQuadTwoDDemoMesh.readCentaurMesh(filename2.str());
    
    triQuadTwoDDemoMesh.outputMesh(std::cout);
    
    file2D.open ("SavedCentaurMinimum.dat");
    
    Output::TecplotDiscontinuousSolutionWriter<2> out2(file2D,"Minimum Test Mesh",dimensionsToWrite,"xy");
    out2.write(&triQuadTwoDDemoMesh,"holi",false);
    
    file2D.close();
    

}

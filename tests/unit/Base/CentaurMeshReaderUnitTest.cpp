#include "Base/MeshManipulator.hpp"
#include "Geometry/PointPhysical.hpp"
#include "Base/Element.hpp"
#include "Geometry/PhysicalGeometry.hpp"

#include <vector>
#include "CMakeDefinitions.hpp"
#include <sstream>
#include "Output/TecplotDiscontinuousSolutionWriter.hpp"

using namespace std;

class Dummy
{
public:
    Dummy(){}
    void operator()(const Base::Element<2>& el, const Geometry::PointReference<2>& p, ostream& os)
    {
    }
};

int main()
{
    Base::ConfigurationData config(3,1,1);
    
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
    
    
    Output::TecplotDiscontinuousSolutionWriter<2> out(file2D,"QuadMinimum Test Mesh","01","xy");
    Dummy d;
    out.write(&myTwoDDemoMesh,"holi",false, d);
    
    file2D.close();
    
    //Now do it again with a more complicated mesh
    
    Base::MeshManipulator<2> triQuadTwoDDemoMesh(&config,1,1);
    
    std::stringstream filename2;
    
    filename2 << getCMAKE_hpGEM_SOURCE_DIR() << "/tests/files/centaurMinimum.hyb";
    
    cout << "filename " << filename2.str() <<endl;
    
    triQuadTwoDDemoMesh.readCentaurMesh(filename2.str());
    
    triQuadTwoDDemoMesh.outputMesh(std::cout);
    
    file2D.open ("SavedCentaurMinimum.dat");
    
    Output::TecplotDiscontinuousSolutionWriter<2> out2(file2D,"Minimum Test Mesh","01","x,y");
    out2.write(&triQuadTwoDDemoMesh,"holi",false,d);
    
    file2D.close();
    

}

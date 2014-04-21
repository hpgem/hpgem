/*
 * Problem.cpp
 *
 *  Created on: Feb 21, 2013
 *      Author: nicorivas
 */

#include "Base/HpgemUI.hpp"
#include "Base/RectangularMeshDescriptor.hpp"
#include "MeshMover.hpp"
#include "Base/GlobalData.hpp"
#include "Base/ConfigurationData.hpp"
#include "Output/TecplotSingleElementWriter.hpp"

#include "Output/TecplotDiscontinuousSolutionWriter.hpp"
#include "Output/TecplotPhysicalGeometryIterator.hpp"
using Base::RectangularMeshDescriptor;
using Base::ConfigurationData;
using Base::GlobalData;

const unsigned int DIM = 2;

//Note: the intended use of the prototype classes is to merge Dummy with MeshMoverExampleProblem
class Dummy:public Output::TecplotSingleElementWriter
{
public:
    Dummy(){}
    void writeToTecplotFile(const Base::Element* el, const Geometry::PointReference& p, ostream& os)
    {
    }
};

class MeshMoverExampleProblem : public Base::HpgemUI
{
    
public:
    MeshMoverExampleProblem(GlobalData* const global, const ConfigurationData* config):
        Base::HpgemUI(global, config)
    {
    }
    
    bool initialise()
    {
        
        RectangularMeshDescriptor rectangularMesh(2);
        
        rectangularMesh.bottomLeft_[0] = 0;
        rectangularMesh.bottomLeft_[1] = 0;
        rectangularMesh.topRight_[0] = 1;
        rectangularMesh.topRight_[1] = 1;
        rectangularMesh.numElementsInDIM_[0] = 8;
        rectangularMesh.numElementsInDIM_[1] = 8;
        
        Base::HpgemUI::MeshId id = addMesh(rectangularMesh);

        //Set up the move of the mesh;
        const MeshMover* meshMover= new MeshMover;
        initialiseMeshMover(meshMover, id);
        
        return true;
    }
    
    void output()
    {
        std::ofstream file2D;
        file2D.open ("out.dat");
        int dimensionsToWrite[2] = {0,1};
        Output::TecplotDiscontinuousSolutionWriter out(file2D,"RectangularMesh","01","xy");
        
        Dummy d;
        out.write(meshes_[0],"holi",false, &d);
    }
    
    void solve()
    {
     
        meshes_[0]->move();
        
    }
    
};

int main(int argc, char **argv)
{
   
    Base::GlobalData globalData;
    
    Base::ConfigurationData config(1,1,1);
    
    config.numberOfUnknowns_       = 1;
    config.numberOfTimeLevels_     = 1;
    config.numberOfBasisFunctions_ = 1;
    
    globalData.numberOfUnknowns_=10;
    globalData.numberOfTimeLevels_=1;

    MeshMoverExampleProblem problem(&globalData, &config);
    
    problem.initialise();
    
    problem.solve();
    
    problem.output();
}

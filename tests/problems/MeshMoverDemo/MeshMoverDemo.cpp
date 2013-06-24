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


#include "Output/TecplotDiscontinuousSolutionWriter.hpp"
#include "Output/TecplotPhysicalGeometryIterator.hpp"
using Base::RectangularMeshDescriptor;
using Base::ConfigurationData;
using Base::GlobalData;

const unsigned int DIM = 2;
class MeshMoverExampleProblem : public Base::HpgemUI<DIM>
{
    
public:
    MeshMoverExampleProblem(GlobalData* global, const ConfigurationData* config):
        HpgemUI(global, config)
    {
    }
    
    bool initialise()
    {
        
        RectangularMeshDescriptor<DIM> rectangularMesh;
        
        rectangularMesh.bottomLeft_[0] = 0;
        rectangularMesh.bottomLeft_[1] = 0;
        rectangularMesh.topLeft_[0] = 1;
        rectangularMesh.topLeft_[1] = 1;
        rectangularMesh.numElementsInDIM_[0] = 8;
        rectangularMesh.numElementsInDIM_[1] = 8;
        
        typename HpgemUI<DIM>::MeshId id = addMesh(rectangularMesh);

        //Set up the move of the mesh;
        const MeshMover<2>* meshMover= new MeshMover<2>;
        initialiseMeshMover(meshMover, id);
        
        return true;
    }
    
    void output()
    {
        std::ofstream file2D;
        file2D.open ("out.dat");
        int dimensionsToWrite[2] = {0,1};
        Output::TecplotDiscontinuousSolutionWriter<2> out(file2D,"RectangularMesh",dimensionsToWrite,"xy");
        out.write(meshes_[0],"holi",false);
    }
    
    void solve()
    {
     
        meshes_[0]->move();
        
    }
    
};

int main(int argc, char **argv)
{
   
    Base::GlobalData globalData;
    
    Base::ConfigurationData config;
    
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

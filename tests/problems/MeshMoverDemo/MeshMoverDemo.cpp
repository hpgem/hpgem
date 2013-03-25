/*
 * Problem.cpp
 *
 *  Created on: Feb 21, 2013
 *      Author: nicorivas
 */

#include "Base/Base.hpp"
#include "MeshMover.hpp"
#include "Output/TecplotDiscontinuousSolutionWriter.hpp"
#include "Output/TecplotPhysicalGeometryIterator.hpp"

class MeshMoverExampleProblem : public Base::Base<2>
{
    
public:

    
    bool initialise()
    {
        Geometry::PointPhysical<2> bottomLeft, topRight;
        std::vector<unsigned int> numElementsOneD(2);
        bottomLeft[0] = 0;
        bottomLeft[1] = 0;
        topRight[0] = 1;
        topRight[1] = 1;
        numElementsOneD[0] = 8;
        numElementsOneD[1] = 8;
        addMesh("Rectangular",bottomLeft, topRight, numElementsOneD);
        
        //Set up the move of the mesh;
        Base::MeshMover<2>* meshMover= new Base::MeshMover<2>;
        initialiseMeshMover(meshMover);
        
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
    MeshMoverExampleProblem problem;

    problem.initialise();
    
    problem.solve();
    
    problem.output();
}

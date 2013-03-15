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
    
    bool initialiseMesh()
    {
        Geometry::PointPhysical<2> bottomLeft, topLeft;
        std::vector<unsigned int> numElementsOneD(2);
        bottomLeft[0] = 0;
        bottomLeft[1] = 0;
        topLeft[0] = 1;
        topLeft[1] = 1;
        numElementsOneD[0] = 8;
        numElementsOneD[1] = 8;
        mesh_.createRectangularMesh(bottomLeft, topLeft, numElementsOneD);
        return true;
    }
    
  void output()
    {
        std::ofstream file2D;
        file2D.open ("out.dat");
        int dimensionsToWrite[2] = {0,1};
        Output::TecplotDiscontinuousSolutionWriter<2> out(file2D,"RectangularMesh",dimensionsToWrite,"xy");
        out.write(mesh_,"holi",false);
    }
    
};

int main(int argc, char **argv)
{
    MeshMoverExampleProblem problem;

    problem.initialiseMesh();

    Base::MeshMover<2> meshMover;

    problem.initialiseMeshMover(&meshMover);


    problem.solve();
    
    problem.output();
}

/*
 * Problem.cpp
 *
 *  Created on: March 15, 2013
 *      Author: the ghost of nicorivas
 */

#include "Base/BaseSimplified.hpp"
#include "Output/TecplotDiscontinuousSolutionWriter.hpp"
#include "Output/TecplotPhysicalGeometryIterator.hpp"
#include "LinearAlgebra/NumericalVector.hpp"

class SimpleDemoProblem : public Base::BaseSimplified<2>
{
    
public:
    
    bool initialise()
    {
        Geometry::PointPhysical<2> bottomLeft, topLeft;
        std::vector<unsigned int> numElementsOneD(2);
        bottomLeft[0] = 0;
        bottomLeft[1] = 0;
        topLeft[0] = 1;
        topLeft[1] = 1;
        numElementsOneD[0] = 8;
        numElementsOneD[1] = 8;
        addMesh("Rectangular",bottomLeft, topLeft, numElementsOneD);
        return true;
    }

    void elementIntegrand(const ElementT& element, const PointReferenceT& p, LinearAlgebra::Matrix& ret)
    {
        
//         double phi=element.getData(0,0)
//         ret[0]=
        
        
    }
    
    void faceIntegrand(const PointPhysicalT& normal, 
                       const PointReferenceT& p,  LinearAlgebra::NumericalVector& ret)
    {
    
        
    
    }
    
    void initialConditions(const PointPhysicalT& p){}
    
    void output()
    {
        std::ofstream file2D;
        file2D.open ("out.dat");
        int dimensionsToWrite[2] = {0,1};
        Output::TecplotDiscontinuousSolutionWriter<2> out(file2D,"RectangularMesh",dimensionsToWrite,"xy");
        out.write(meshes_[0],"holi",false);
    }
    
};

int main(int argc, char **argv)
{
//     SimpleDemoProblem problem;
//     
//     problem.solve();
//     
//     problem.output();
}



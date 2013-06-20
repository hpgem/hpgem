/*
 * Problem.cpp
 *
 *  Created on: March 15, 2013
 *      Author: the ghost of nicorivas
 */

#include "Base/HpgemUISimplified.hpp"
#include "Base/RectangularMeshDescriptor.hpp"
#include "Output/TecplotDiscontinuousSolutionWriter.hpp"
#include "Output/TecplotPhysicalGeometryIterator.hpp"
#include "LinearAlgebra/NumericalVector.hpp"
using Base::RectangularMeshDescriptor;
using Base::HpgemUISimplified;

const unsigned int DIM = 2;
class SimpleDemoProblem : public HpgemUISimplified<DIM>
{
    
public:
    
    bool initialise()
    {
        RectangularMeshDescriptor<DIM> rectangularMesh;
        
        rectangularMesh.bottomLeft_[0] = 0;
        rectangularMesh.bottomLeft_[1] = 0;
        rectangularMesh.topLeft_[0] = 1;
        rectangularMesh.topLeft_[1] = 1;
        rectangularMesh.numElementsInDIM_[0] = 8;
        rectangularMesh.numElementsInDIM_[1] = 8;
        
        addMesh(rectangularMesh);
        
        return true;
    }

    void elementIntegrand(const ElementT& element, const PointReferenceT& p, LinearAlgebra::Matrix& ret)
    {
        
     //   double phi=element.getSolution(0);
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



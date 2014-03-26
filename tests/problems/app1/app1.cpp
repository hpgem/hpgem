/*
 * Problem.cpp
 *
 *  Created on: March 15, 2013
 *      Author: the ghost of nicorivas
 *
 * This a simple 1D advection example.
 */

#include "Base/HpgemUISimplified.hpp"
#include "Base/RectangularMeshDescriptor.hpp"
#include "Output/TecplotDiscontinuousSolutionWriter.hpp"
#include "Output/TecplotPhysicalGeometryIterator.hpp"
#include "LinearAlgebra/NumericalVector.hpp"
#include "Base/PhysGradientOfBasisFunction.hpp"
#include "Base/Norm2.hpp"
#include "Output/TecplotSingleElementWriter.hpp"
using Base::RectangularMeshDescriptor;
using Base::HpgemUISimplified;

const unsigned int DIM = 2;

//Note: the intended use of the prototype classes is to merge Dummy with SimpleDemoProblem
class Dummy: public Output::TecplotSingleElementWriter
{
public:
    Dummy(){}
    void writeToTecplotFile(const Base::Element* el, const Geometry::PointReference& p, ostream& os)
    {
    }
};

class SimpleDemoProblem : public HpgemUISimplified
{
    
public:
    
    bool initialise()
    {
        RectangularMeshDescriptor rectangularMesh(2);
        
        rectangularMesh.bottomLeft_[0] = 0;
        rectangularMesh.bottomLeft_[1] = 0;
        rectangularMesh.topRight_[0] = 1;
        rectangularMesh.topRight_[1] = 1;
        rectangularMesh.numElementsInDIM_[0] = 8;
        rectangularMesh.numElementsInDIM_[1] = 8;
        
        addMesh(rectangularMesh);
        
        return true;
    }

    void elementIntegrand(const ElementT* element, const PointReferenceT& p, LinearAlgebra::Matrix& ret)
    {
        
        unsigned int numberOfDegreesOfFreedom=element->getNrOfBasisFunctions();
       
        LinearAlgebra::NumericalVector sol;
        element->getSolution(0,p,sol);
        
        //This is the grad of the basic function.
        NumericalVector grads(2);
        
        for (unsigned int i=0; i < numberOfDegreesOfFreedom; ++i)
        {
            element->basisFunctionDeriv(i,p,grads);
            
            ret(i,0) = sol(i) * grads[i];
            
        }
    
        
    }
    
    void faceIntegrand(const FaceT* face, const PointPhysicalT& normal, 
                       const PointReferenceT& p,  LinearAlgebra::NumericalVector& ret)
    {
        
        if (face->isInternal())
        {
         
            const double magn                     = Utilities::norm2(normal);
            unsigned int numberOfDegreesOfFreedom = face->getPtrElementLeft()->getNrOfBasisFunctions();
            
        }
        else
        {
            //here you have to implement the boundary conditions
        }

    
    }
    
    double initialConditions(const PointPhysicalT& p){return 0;}
    
    void output()
    {
        std::ofstream file2D;
        file2D.open ("out.dat");
        int dimensionsToWrite[2] = {0,1};
        Output::TecplotDiscontinuousSolutionWriter out(file2D,"RectangularMesh","01","xy");
        Dummy d;
        out.write(meshes_[0],"holi",false, &d);
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



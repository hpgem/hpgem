/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, University of Twente
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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



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

//this file should contain all relevant information about how the integrands look like and what problem is solved

#define _USE_MATH_DEFINES
#include <cstdlib>
#include "Output/TecplotDiscontinuousSolutionWriter.h"
#include "Output/TecplotPhysicalGeometryIterator.h"
#include "LinearAlgebra/MiddleSizeVector.h"
#include "LinearAlgebra/SmallVector.h"
#include "BasisFunctionCollection_Curl.h"
#include <iostream>
#include "Base/L2Norm.h"
#include "BaseExtended.h"
#include "math.h"
#include <ctime>
#include "ElementInfos.h"
#include "Base/HpgemAPISimplified.h"
#include "Integration/FaceIntegral.h"
#include "Integration/FaceIntegral_Impl.h"
#include "Integration/ElementIntegral.h"
#include "Integration/ElementIntegral_Impl.h"
#include "Output/TecplotSingleElementWriter.h"
#include "Geometry/PointPhysical.h"
#include "Base/FaceCacheData.h"
#include "Base/ElementCacheData.h"
#include "Base/ConfigurationData.h"
#include "Base/HCurlConformingTransformation.h"
#include "Integration/ElementIntegrandBase.h"
#include "Integration/FaceIntegrandBase.h"
#include "Base/Face.h"
#include "Base/Element.h"
#include "LinearAlgebra/MiddleSizeMatrix.h"
#include "Utilities/GlobalMatrix.h"
#include "Utilities/GlobalVector.h"
#include "petscksp.h"
#include "matrixAssembly.h"

using BaseMeshManipulatorT = Base::MeshManipulator<DIM>;
using basisFunctionT = Base::threeDBasisFunction;
using FaceIterator = Base::MeshManipulator<DIM>::FaceIterator;

/**
 * This class should provide problem specific information about the maxwell equations.
 */
class DGMax : public hpGemUIExtentions
{
private:
    
    using PointPhysicalT = Geometry::PointPhysical<DIM>;
    using ElementT = Base::Element;
    using FaceT = Base::Face;

    using FaceIterator = Base::MeshManipulatorBase::FaceIterator;

public:
    
    DGMax(MaxwellData* globalConfig, Base::ConfigurationData* elementConfig, MatrixAssembly* fill)
            : hpGemUIExtentions(globalConfig, elementConfig, fill)
    {
    }
    
    /**
     * set up the mesh and complete initialisation of the global data and the configuration data
     */
    
    bool initialise()
    {
        int n = getData()->NumberOfIntervals_;
        Geometry::PointPhysical<DIM> bottomLeft, topRight;
        std::vector<std::size_t> numElementsOneD(3);
        bottomLeft[0] = 0;
        bottomLeft[1] = 0;
        bottomLeft[2] = 0;
        topRight[0] = 1;
        topRight[1] = 1;
        topRight[2] = 1;
        numElementsOneD[0] = n;
        numElementsOneD[1] = n;
        numElementsOneD[2] = n;
        
        // this is the old way of creating the mesh.
        //BaseMeshManipulatorT* mesh = new MyMeshManipulator(getConfigData(), Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, getData()->PolynomialOrder_, 0, 2, 3, 1, 1, true);
        
        // this is the way the mesh is created in hpGEM.
        BaseMeshManipulatorT* mesh = new Base::MeshManipulator<3>(getConfigData(), Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, Base::BoundaryType::SOLID_WALL, getData()->PolynomialOrder_, 0, 2, 3, 1, 1);
        mesh->createTriangularMesh(bottomLeft, topRight, numElementsOneD);
        mesh->useNedelecDGBasisFunctions();
        //mesh->useAinsworthCoyleDGBasisFunctions();
        const_cast<MaxwellData*>(getData())->numberOfUnknowns_ = (*mesh->elementColBegin())->getNrOfBasisFunctions();
        setConfigData();
        //mesh->readCentaurMesh("Cylinder3.hyb");
        //mesh->readCentaurMesh("input_basic2.hyb");
        //mesh->readCentaurMesh("Cube_final.hyb");
        
        addMesh(mesh);
        setNumberOfTimeIntegrationVectorsGlobally(21);
        const_cast<MaxwellData*>(getData())->numberOfUnknowns_ *= mesh->getElementsList().size();
        for (Base::MeshManipulator<DIM>::ElementIterator it = mesh->elementColBegin(); it != mesh->elementColEnd(); ++it)
        {
            (*it)->setUserData(new ElementInfos(**it));
        }
        return true;
    }
};

/**
 * Computes element contributions to the mass matrix i.e. phi_i * phi_j
 * returns the contibutions at this gauss point to the entire element matrix in one go
 */

void MatrixAssemblyIP::CompleteElementIntegrationIP(hpGemUIExtentions* matrixContainer)
{
    LinearAlgebra::MiddleSizeMatrix matrix1(1, 1), matrix2(1, 1);
    LinearAlgebra::MiddleSizeVector vector1(1), vector2(1), vector3(1);
    Integration::ElementIntegral<DIM> elIntegral(false);
    
    elIntegral.setTransformation(std::shared_ptr<Base::CoordinateTransformation<DIM>> (new Base::HCurlConformingTransformation<DIM>()));
    for (hpGemUIExtentions::ElementIterator it = matrixContainer->elementColBegin(); it != matrixContainer->elementColEnd(); ++it)
    {
        matrix1.resize((*it)->getNrOfBasisFunctions(), (*it)->getNrOfBasisFunctions());
        matrix1 = elIntegral.integrate<LinearAlgebra::MiddleSizeMatrix>((*it), &(matrixContainer->elementMassIntegrand));
        if(matrixContainer->MHasToBeInverted_)
        {
            matrix1 = matrix1.inverse();
        }
        (*it)->setElementMatrix(matrix1, 0);
        
        matrix2.resize((*it)->getNrOfBasisFunctions(), (*it)->getNrOfBasisFunctions());
        matrix2 = elIntegral.integrate<LinearAlgebra::MiddleSizeMatrix>((*it), &(matrixContainer->elementStiffnessIntegrand));
        (*it)->setElementMatrix(matrix2, 1);
        
        vector1.resize((*it)->getNrOfBasisFunctions());
        vector1 = elIntegral.integrate<LinearAlgebra::MiddleSizeVector>((*it), &(matrixContainer->initialConditionsIntegrand));
        (*it)->setElementVector(vector1, 0);
        
        vector2.resize((*it)->getNrOfBasisFunctions());
        vector2 = elIntegral.integrate<LinearAlgebra::MiddleSizeVector>((*it), &(matrixContainer->initialConditionsDerivIntegrand));
        (*it)->setElementVector(vector2, 1);
        
        vector3.resize((*it)->getNrOfBasisFunctions());
        vector3 = elIntegral.integrate<LinearAlgebra::MiddleSizeVector>((*it), &(matrixContainer->elementSpaceIntegrand));
        (*it)->setElementVector(vector3, 2);
        
    }
   
}

void MatrixAssemblyIP::CompleteFaceIntegrationIP(hpGemUIExtentions* matrixContainer)
{
    LinearAlgebra::MiddleSizeMatrix matrix(1, 1), matrix1(1, 1), matrix2(1, 1);
    LinearAlgebra::MiddleSizeVector vector1(1);
    Integration::FaceIntegral<DIM> faIntegral(false);
    
    faIntegral.setTransformation(std::shared_ptr<Base::CoordinateTransformation<DIM>> (new Base::HCurlConformingTransformation<DIM>()));
    for (hpGemUIExtentions::FaceIterator it = matrixContainer->faceColBegin(); it != matrixContainer->faceColEnd(); ++it)
    {
        
        
        if ((*it)->isInternal())
        {
            matrix1.resize((*it)->getPtrElementLeft()->getNrOfBasisFunctions() + (*it)->getPtrElementRight()->getNrOfBasisFunctions(), (*it)->getPtrElementLeft()->getNrOfBasisFunctions() + (*it)->getPtrElementRight()->getNrOfBasisFunctions());
           
        
            matrix2.resize((*it)->getPtrElementLeft()->getNrOfBasisFunctions() + (*it)->getPtrElementRight()->getNrOfBasisFunctions(), (*it)->getPtrElementLeft()->getNrOfBasisFunctions() + (*it)->getPtrElementRight()->getNrOfBasisFunctions());
            
            matrix.resize((*it)->getPtrElementLeft()->getNrOfBasisFunctions() + (*it)->getPtrElementRight()->getNrOfBasisFunctions(), (*it)->getPtrElementLeft()->getNrOfBasisFunctions() + (*it)->getPtrElementRight()->getNrOfBasisFunctions());
            
            vector1.resize((*it)->getPtrElementLeft()->getNrOfBasisFunctions() + (*it)->getPtrElementRight()->getNrOfBasisFunctions());
           
        }
        
        else
        {
            matrix1.resize((*it)->getPtrElementLeft()->getNrOfBasisFunctions(), (*it)->getPtrElementLeft()->getNrOfBasisFunctions());
            
            matrix2.resize((*it)->getPtrElementLeft()->getNrOfBasisFunctions(), (*it)->getPtrElementLeft()->getNrOfBasisFunctions());
            
            matrix.resize((*it)->getPtrElementLeft()->getNrOfBasisFunctions(), (*it)->getPtrElementLeft()->getNrOfBasisFunctions());
            
            vector1.resize((*it)->getPtrElementLeft()->getNrOfBasisFunctions());
            

        }
        matrix1 = faIntegral.integrate<LinearAlgebra::MiddleSizeMatrix>((*it), &(matrixContainer->faceStiffnessIntegrand));
        matrix2 = faIntegral.integrate<LinearAlgebra::MiddleSizeMatrix>((*it), &(matrixContainer->faceStiffnessIntegrandIP));
        matrix = matrix1 + matrix2;
        (*it)->setFaceMatrix(matrix, 0);
        vector1 = faIntegral.integrate<LinearAlgebra::MiddleSizeVector>((*it), &(matrixContainer->faceSpaceIntegrandIP));
        (*it)->setFaceVector(vector1, 0);
    }
}

void MatrixAssemblyIP::fillMatrices(hpGemUIExtentions* matrixContainer)
{
    
    CompleteElementIntegrationIP(matrixContainer);
    std::cout << "CompleteElementIntegrationIP done" << std::endl;
    CompleteFaceIntegrationIP(matrixContainer);
    std::cout << "CompleteFaceIntegrationIP done" << std::endl;
    std::cout << "fillMatricesIP done" << std::endl;
    
}

auto& numElements = Base::register_argument<std::size_t>('n', "numElems", "number of elements per dimension", true);
auto& p = Base::register_argument<std::size_t>('p', "order", "polynomial order of the solution", true);

/**
 * main routine. Suggested usage: make a problem, solve it, tell the user your results.
 */
int main(int argc, char** argv)
{
    Base::parse_options(argc, argv);
    std::cout<<"This is the parallel version"<<std::endl;
    
    //set up timings
    time_t start, end, initialised, solved;
    time(&start);
    logger(INFO, "using % elements", std::pow(numElements.getValue(), 3) * 5);
    logger(INFO, "using polynomial order: %", p.getValue());
    //set up problem and decide flux type 
    DGMax problem(new MaxwellData(numElements.getValue(), p.getValue()), new Base::ConfigurationData(3, 1, p.getValue(), 1), new MatrixAssemblyIP);
    try
    {
        problem.initialise();
        time(&initialised);
        //choose what problem to solve
        //problem.solveEigenvalues();
        problem.solveHarmonic();
        //problem.solveTimeDependent(false,true);
        std::cout << "solved for Harmonic" << std::endl;
        time(&solved);
        char filename[] = "output.dat";
        problem.makeOutput(filename); //Issue with parallelisation occurs in the makeOutput function
        time(&end);
        
        //display timing data
        std::cout << "Initialisation took " << difftime(initialised, start) << " seconds." << std::endl;
        std::cout << "Solving the problem took " << difftime(solved, initialised) << " seconds." << std::endl;
        std::cout << "The rest took " << difftime(end, solved) << " seconds." << std::endl;
        
    }
    catch (const char* message)
    {
        std::cout << message;
    }
    return 0;
}

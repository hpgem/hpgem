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
#include "LinearAlgebra/NumericalVector.h"
#include "BasisFunctionCollection_Curl.h"
#include <iostream>
#include "Base/L2Norm.h"
#include "BaseExtended.h"
#include "math.h"
#include <ctime>
#include "ElementInfos.h"
#include "Base/HpgemUI.h"
#include "Base/HpgemUISimplified.h"
#include "Integration/ElementIntegrandBase.h"
#include "Integration/FaceIntegrandBase.h"
#include "Integration/FaceIntegral.h"
#include "Integration/ElementIntegral.h"
#include "Output/TecplotSingleElementWriter.h"
#include "Geometry/PointPhysical.h"
#include "Base/FaceCacheData.h"
#include "Base/ElementCacheData.h"
#include "Base/ConfigurationData.h"
#include "Base/ShortTermStorageElementHcurl.h"
#include "Base/ShortTermStorageFaceHcurl.h"
#include "Integration/ElementIntegrandBase.h"
#include "Integration/FaceIntegrandBase.h"
#include "Base/Face.h"
#include "Base/Element.h"
#include "LinearAlgebra/NumericalVector.h"
#include "Utilities/GlobalMatrix.h"
#include "Utilities/GlobalVector.h"
#include "petscksp.h"
#include "matrixAssembly.h"

using BaseMeshManipulatorT = Base::MeshManipulator;
using basisFunctionT = Base::threeDBasisFunction;
//typedef std::list<Base::Face<3> >::iterator FaceIteratorT;

using FaceIterator = Base::MeshManipulator::FaceIterator;

/**
 * This class should provide problem specific information about the maxwell equations.
 */
class DGMax : public hpGemUIExtentions
{
private:
    
    using PointPhysicalT = Geometry::PointPhysical;
    using ElementT = Base::Element;
    using FaceT = Base::Face;

    using FaceIterator = Base::MeshManipulator::FaceIterator;

public:
    
    DGMax(int argc, char** argv, MaxwellData* globalConfig, Base::ConfigurationData* elementConfig, MatrixAssembly* fill)
            : hpGemUIExtentions(argc, argv, globalConfig, elementConfig, fill)
    {
    }
    
    /**
     * set up the mesh and complete initialisation of the global data and the configuration data
     */
    bool initialise()
    {
        int n = getData()->NumberOfIntervals_;
        Geometry::PointPhysical bottomLeft(3), topRight(3);
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
        
        BaseMeshManipulatorT* mesh = new MyMeshManipulator(getConfigData(), false, false, false, getData()->PolynomialOrder_, 0, 2, 3, 1, 1);
        mesh->createTriangularMesh(bottomLeft, topRight, numElementsOneD);
        
        const_cast<MaxwellData*>(getData())->numberOfUnknowns_ = (*mesh->elementColBegin())->getNrOfBasisFunctions();
        
        setConfigData();
        //mesh->readCentaurMesh("Cylinder3.hyb");
        //mesh->readCentaurMesh("input_basic2.hyb");
        //mesh->readCentaurMesh("Cube_final.hyb");
        
        addMesh(mesh);
        const_cast<MaxwellData*>(getData())->numberOfUnknowns_ *= mesh->getElementsList().size();
        
        //addMesh("Triangular",bottomLeft, topRight, numElementsOneD);
        
        //actually store the full number of unknowns
        
        for (Base::MeshManipulator::ElementIterator it = mesh->elementColBegin(); it != mesh->elementColEnd(); ++it)
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
    
    LinearAlgebra::Matrix matrix1(1, 1), matrix2(1, 1);
    LinearAlgebra::MiddleSizeVector vector1(1), vector2(1), vector3(1);
    
    Base::ShortTermStorageElementBase* localElement_;
    localElement_ = new Base::ShortTermStorageElementHcurl(3); // creating object for H curl transformation
    Integration::ElementIntegral elIntegral(false);
    elIntegral.setStorageWrapper(localElement_);
    //std::cout<<"ElementHcurl"<<std::endl;
    
    for (hpGemUIExtentions::ElementIterator it = matrixContainer->elementColBegin(); it != matrixContainer->elementColEnd(); ++it)
    {
        
        matrix1.resize((*it)->getNrOfBasisFunctions(), (*it)->getNrOfBasisFunctions());
        //std::cout<<"Matrix Resize done"<<std::endl;
        elIntegral.integrate<LinearAlgebra::Matrix>((*it), &(matrixContainer->elementMassIntegrand), matrix1);
        //std::cout<<"Integral Called"<<std::endl;
        (*it)->setElementMatrix(matrix1, 0);
        
        //std::cout<<"FirstIntegrationDone"<<std::endl;
        
        matrix2.resize((*it)->getNrOfBasisFunctions(), (*it)->getNrOfBasisFunctions());
        elIntegral.integrate<LinearAlgebra::Matrix>((*it), &(matrixContainer->elementStiffnessIntegrand), matrix2);
        (*it)->setElementMatrix(matrix2, 1);
        
        vector1.resize((*it)->getNrOfBasisFunctions());
        elIntegral.integrate<LinearAlgebra::MiddleSizeVector>((*it), &(matrixContainer->initialConditionsIntegrand), vector1);
        (*it)->setElementVector(vector1, 0);
        
        vector2.resize((*it)->getNrOfBasisFunctions());
        elIntegral.integrate<LinearAlgebra::MiddleSizeVector>((*it), &(matrixContainer->initialConditionsDerivIntegrand), vector2);
        (*it)->setElementVector(vector2, 1);
        
        vector3.resize((*it)->getNrOfBasisFunctions());
        elIntegral.integrate<LinearAlgebra::MiddleSizeVector>((*it), &(matrixContainer->elementSpaceIntegrand), vector3);
        (*it)->setElementVector(vector3, 2);
        
    }
}

void MatrixAssemblyIP::CompleteFaceIntegrationIP(hpGemUIExtentions* matrixContainer)
{
    LinearAlgebra::Matrix matrix(1, 1), matrix1(1, 1), matrix2(1, 1);
    LinearAlgebra::MiddleSizeVector vector0(1), vector1(1);
    //std::cout<<"Complete Face Integration IP started"<<std::endl;
    
    Base::ShortTermStorageFaceBase* localFace_;
    localFace_ = new Base::ShortTermStorageFaceHcurl(3); // creating object for H curl transformation
    Integration::FaceIntegral faIntegral(false);
    faIntegral.setStorageWrapper(localFace_);
    
    //std::cout<<"FaceHcurl"<<std::endl;
    
    for (hpGemUIExtentions::FaceIterator it = matrixContainer->faceColBegin(); it != matrixContainer->faceColEnd(); ++it)
    {
        //std::cout<<"Iteration over face started"<<std::endl;
        if ((*it)->isInternal())
        {
            //std::cout<<"Internal face"<<std::endl;
            //matrix1.resize((*it)->getPtrElementLeft()->getNrOfBasisFunctions() + (*it)->getPtrElementRight()->getNrOfBasisFunctions(), (*it)->getPtrElementLeft()->getNrOfBasisFunctions() + (*it)->getPtrElementRight()->getNrOfBasisFunctions());
            matrix1.resize((*it)->getNrOfBasisFunctions(), (*it)->getNrOfBasisFunctions());
            
            faIntegral.integrate<LinearAlgebra::Matrix>((*it), &(matrixContainer->faceStiffnessIntegrand), matrix1);
            
            //matrix2.resize((*it)->getPtrElementLeft()->getNrOfBasisFunctions() + (*it)->getPtrElementRight()->getNrOfBasisFunctions(), (*it)->getPtrElementLeft()->getNrOfBasisFunctions() + (*it)->getPtrElementRight()->getNrOfBasisFunctions());
            matrix1.resize((*it)->getNrOfBasisFunctions(), (*it)->getNrOfBasisFunctions());
            
            faIntegral.integrate<LinearAlgebra::Matrix>((*it), &(matrixContainer->faceStiffnessIntegrandIP), matrix2);
            
            matrix = matrix1 + matrix2;
            
            (*it)->setFaceMatrix(matrix, 0);
            
            vector0.resize((*it)->getNrOfBasisFunctions());
            for (int i = 0; i < vector0.size(); i++)
                vector0(i) = 0;
            
            //(*it)->setFaceVector(vector0, 0);
            (*it)->setFaceVector(vector0, 0);
            
        }
        
        else
        {
            //std::cout<<"Boundary face"<<std::endl;
            matrix1.resize((*it)->getPtrElementLeft()->getNrOfBasisFunctions(), (*it)->getPtrElementLeft()->getNrOfBasisFunctions());
            
            faIntegral.integrate<LinearAlgebra::Matrix>((*it), &(matrixContainer->faceStiffnessIntegrand), matrix1);
            
            matrix2.resize((*it)->getPtrElementLeft()->getNrOfBasisFunctions(), (*it)->getPtrElementLeft()->getNrOfBasisFunctions());
            
            faIntegral.integrate<LinearAlgebra::Matrix>((*it), &(matrixContainer->faceStiffnessIntegrandIP), matrix2);
            
            matrix = matrix1 + matrix2;
            
            (*it)->setFaceMatrix(matrix, 0);
            
            // std::cout<<"Matrix1 has be set"<<std::endl;
            vector1.resize((*it)->getPtrElementLeft()->getNrOfBasisFunctions());
            
            faIntegral.integrate<LinearAlgebra::MiddleSizeVector>((*it), &(matrixContainer->faceSpaceIntegrandIP), vector1);
            
            (*it)->setFaceVector(vector1, 0);
            
            //std::cout<<"Vector set"<<std::endl;
            
        }
        
    }
    //std::cout<<"Iteration over face ended"<<std::endl;
}

void MatrixAssemblyIP::fillMatrices(hpGemUIExtentions* matrixContainer)
{
    
    CompleteElementIntegrationIP(matrixContainer);
    CompleteFaceIntegrationIP(matrixContainer);
    std::cout << "fillMatricesIP done" << std::endl;
    
}

/**
 * main routine. Suggested usage: make a problem, solve it, tell the user your results.
 */
int main(int argc, char** argv)
{
    
    //set up timings
    time_t start, end, initialised, solved;
    time(&start);
    
    //read input data
    int elements = 1;
    if (argc > 1)
    {
        elements = std::atoi(argv[1]);
        std::cout << "using " << elements * elements * elements * 5 << " elements" << std::endl;
    }
    int order = 1;
    if (argc > 2)
    {
        order = std::atoi(argv[2]);
        std::cout << "using polynomial order: " << order << std::endl;
    }
    else
    {
        std::cout << "usage:./Maxwell.out <elements> <order> [<petsc-args>]";
        exit(1);
    }
    //set up problem and decide flux type
    DGMax problem(argc - 2, &argv[2], new MaxwellData(elements, order), new Base::ConfigurationData(3, 1, order, 1), new MatrixAssemblyIP);
    //Base::ConfigurationData(3,1,order,61) for solveEigenvalues()
    try
    {
        problem.initialise();
        time(&initialised);
        
        //choose what problem to solve
        //problem.solveEigenvalues();
        problem.solveHarmonic();
        std::cout << "solved for Harmonic" << std::endl;
        time(&solved);
        char filename[] = "output.dat";
        problem.makeOutput(filename);
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

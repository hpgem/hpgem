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

#include "HpgemAPILinearSteadyState.h"

#include "Base/CommandLineOptions.h"
#include "Base/ConfigurationData.h"
#include "Base/Element.h"
#include "Base/Face.h"
#include "Base/MpiContainer.h"
#include "Base/RectangularMeshDescriptor.h"
#include "Base/TimeIntegration/AllTimeIntegrators.h"
#include "Geometry/PointReference.h"
#include "Integration/ElementIntegral.h"
#include "Integration/FaceIntegral.h"
#include "Output/TecplotDiscontinuousSolutionWriter.h"
#include "Output/TecplotSingleElementWriter.h"
#include "Utilities/BasisFunctions1DH1ConformingLine.h"
#include "Utilities/BasisFunctions2DH1ConformingSquare.h"
#include "Utilities/BasisFunctions2DH1ConformingTriangle.h"
#include "Utilities/BasisFunctions3DH1ConformingCube.h"
#include "Utilities/BasisFunctions3DH1ConformingTetrahedron.h"
#include "Utilities/GlobalMatrix.h"
#include "Utilities/GlobalVector.h"

#if defined(HPGEM_USE_PETSC) || defined(HPGEM_USE_COMPLEX_PETSC)
    #include "petscksp.h"
#endif

#include "Logger.h"

namespace Base
{
    /// \param[in] dimension Dimension of the domain
    /// \param[in] numOfVariables Number of variables in the PDE
    /// \param[in] polynomialOrder Polynomial order of the basis functions
    /// \param[in] useSourceTerm Boolean to indicate if there is a source term.
    /// \param[in] useSourceTerm Boolean to indicate if there is a source term at the domain boundary.
    HpgemAPILinearSteadyState::HpgemAPILinearSteadyState
    (
     const std::size_t dimension,
     const std::size_t numOfVariables,
     const std::size_t polynomialOrder,
     const bool useSourceTerm,
     const bool useSourceTermAtBoundary
     ) :
    HpgemAPILinear(dimension, numOfVariables, polynomialOrder, Base::AllTimeIntegrators::Instance().getRule(1, 1), 1, useSourceTerm, useSourceTermAtBoundary),
    sourceElementVectorID_(0),
    sourceFaceVectorID_(0)
    {
    }
    
    void HpgemAPILinearSteadyState::createMesh(const std::size_t numOfElementsPerDirection, const Base::MeshType meshType)
    {
        const Base::RectangularMeshDescriptor description = createMeshDescription(numOfElementsPerDirection);
        
        // Set the number of Element/Face Matrices/Vectors.
        std::size_t numOfElementMatrices = 2;   // Mass matrix and stiffness matrix
        std::size_t numOfElementVectors = 1;    // Source term vector
        std::size_t numOfFaceMatrices = 1;      // Stiffness matrix
        std::size_t numOfFaceVectors = 1;       // Source term vector at boundary

        // Create mesh and set basis functions.
        addMesh(description, meshType, numOfElementMatrices, numOfElementVectors, numOfFaceMatrices, numOfFaceVectors);
        meshes_[0]->useDefaultDGBasisFunctions();
        
        std::size_t nElements = meshes_[0]->getNumberOfElements();
        logger(VERBOSE, "Total number of elements: %", nElements);
    }
    
    void HpgemAPILinearSteadyState::createSourceTerms()
    {
        logger(INFO, "- Creating source term vectors for the elements.");
        for (Base::Element *ptrElement : meshes_[0]->getElementsList())
        {
            LinearAlgebra::NumericalVector sourceTerm(integrateSourceTermAtElement(ptrElement, 0, 0));
            ptrElement->setElementVector(sourceTerm, sourceElementVectorID_);
        }
        
        logger(INFO, "- Creating source term vectors for the boundary faces.");
        for (Base::Face *ptrFace : meshes_[0]->getFacesList())
        {
            LinearAlgebra::NumericalVector sourceTerm(integrateSourceTermAtFace(ptrFace, 0, 0));
            ptrFace->setFaceVector(sourceTerm, sourceFaceVectorID_);
        }
    }
    
    void HpgemAPILinearSteadyState::tasksBeforeSolving()
    {
        logger(INFO, "Computing stiffness matrices.");
        createStiffnessMatrices();
        logger(INFO, "Computing source terms.");
        createSourceTerms();
    }
    
    /// \details Solve the linear problem \f$Ax=b\f$, where \f$A=-S\f$ with \f$S\f$ being the stiffness matrix. At the moment it is not possible to solve the steady-state problem using an intial solution.
    void HpgemAPILinearSteadyState::solveSteadyStateWithPetsc()
    {
#if defined(HPGEM_USE_PETSC) || defined(HPGEM_USE_COMPLEX_PETSC)
        // Create output files for Paraview.
        std::string outputFileNameVTK = outputFileName_;
        
        registerVTKWriteFunctions();
        Output::VTKTimeDependentWriter VTKWriter(outputFileNameVTK, meshes_[0]);
        
        // Create output files for Tecplot.
#ifdef HPGEM_USE_MPI
        std::string outputFileName = outputFileName_ + "." + std::to_string(Base::MPIContainer::Instance().getProcessorID());
#else
        std::string outputFileName = outputFileName_;
#endif
        std::string outputFileNameTecplot = outputFileName + ".dat";
        std::string dimensionsToWrite = "";
        for(std::size_t i = 0; i < configData_->dimension_; i++)
        {
            dimensionsToWrite = dimensionsToWrite + std::to_string(i);
        }
        
        std::string variableString = variableNames_[0];
        for(std::size_t iV = 1; iV < variableNames_.size(); iV++)
        {
            variableString = variableString + "," + variableNames_[iV];
        }
        
        std::ofstream outputFile(outputFileNameTecplot);
        Output::TecplotDiscontinuousSolutionWriter tecplotWriter(outputFile, internalFileTitle_, dimensionsToWrite, variableString);
        
        // Create and Store things before solving the problem.
        tasksBeforeSolving();
        
        // Solve the linear problem
        //Assemble the matrix A of the system Ax = b.
        Utilities::GlobalPetscMatrix A(HpgemAPIBase::meshes_[0], stiffnessElementMatrixID_, stiffnessFaceMatrixID_);
        MatScale(A,-1);
        //Declare the vectors x and b of the system Ax = b.
        Utilities::GlobalPetscVector b(HpgemAPIBase::meshes_[0], sourceElementVectorID_, sourceFaceVectorID_), x(HpgemAPIBase::meshes_[0]);
        
        //Assemble the vector b. This is needed because Petsc assumes you don't know
        //yet whether a vector is a variable or right-hand side the moment it is
        //declared.
        b.assemble();
        
        //Make the Krylov supspace method
        KSP ksp;
        KSPCreate(PETSC_COMM_WORLD, &ksp);
        //Tell ksp that it will solve the system Ax = b.
        KSPSetOperators(ksp, A, A);
        KSPSetFromOptions(ksp);
        KSPSolve(ksp, b, x);
        //Do PETSc magic, including solving.
        KSPConvergedReason converge;
        KSPGetConvergedReason(ksp, &converge);
        int iterations;
        KSPGetIterationNumber(ksp, &iterations);
        logger(INFO, "KSP solver ended because of % in % iterations.", KSPConvergedReasons[converge], iterations);
        
        x.writeTimeLevelData(solutionTimeLevel_);
        
        tecplotWriter.write(meshes_[0], solutionTitle_, false, this, 0);
        VTKWrite(VTKWriter, 0, solutionTimeLevel_);
        return;
#endif
        logger(ERROR, "Petsc is needed to solve the steady state problem using this function (solveSteadyStateWithPetsc).");
    }
}




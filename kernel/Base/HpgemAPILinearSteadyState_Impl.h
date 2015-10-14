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
    /// \param[in] numberOfVariables Number of variables in the PDE
    /// \param[in] polynomialOrder Polynomial order of the basis functions
    /// \param[in] useSourceTerm Boolean to indicate if there is a source term.
    /// \param[in] useSourceTermAtBoundary Boolean to indicate if there is a source term at the domain boundary.
    template<std::size_t DIM>
    HpgemAPILinearSteadyState<DIM>::HpgemAPILinearSteadyState
    (
     const std::size_t numberOfVariables,
     const std::size_t polynomialOrder,
     const bool useSourceTerm,
     const bool useSourceTermAtBoundary
     ) :
    HpgemAPILinear<DIM>(numberOfVariables, polynomialOrder, 1, 0, useSourceTerm, useSourceTermAtBoundary),
    sourceElementVectorID_(0),
    sourceFaceVectorID_(0)
    {
    }

    template<std::size_t DIM>
    void HpgemAPILinearSteadyState<DIM>::createMesh(const std::size_t numberOfElementsPerDirection, const Base::MeshType meshType)
    {
        const Base::RectangularMeshDescriptor<DIM> description = this->createMeshDescription(numberOfElementsPerDirection);
        
        // Set the number of Element/Face Matrices/Vectors.
        std::size_t numberOfElementMatrices = 2;   // Mass matrix and stiffness matrix
        std::size_t numberOfElementVectors = 1;    // Source term vector
        std::size_t numberOfFaceMatrices = 1;      // Stiffness matrix
        std::size_t numberOfFaceVectors = 1;       // Source term vector at boundary

        // Create mesh and set basis functions.
        this->addMesh(description, meshType, numberOfElementMatrices, numberOfElementVectors, numberOfFaceMatrices, numberOfFaceVectors);
        this->meshes_[0]->useDefaultDGBasisFunctions();
        
        // Set the number of time integration vectors according to the size of the Butcher tableau.
        this->setNumberOfTimeIntegrationVectorsGlobally(this->globalNumberOfTimeIntegrationVectors_);
        
        // Plot info about the mesh
        std::size_t nElements = this->meshes_[0]->getNumberOfElements();
        logger(VERBOSE, "Total number of elements: %", nElements);
    }

    template<std::size_t DIM>
    void HpgemAPILinearSteadyState<DIM>::createSourceTerms()
    {
        logger(INFO, "- Creating source term vectors for the elements.");
        for (Base::Element *ptrElement : this->meshes_[0]->getElementsList())
        {
            LinearAlgebra::MiddleSizeVector sourceTerm(this->integrateSourceTermAtElement(ptrElement, 0, 0));
            ptrElement->setElementVector(sourceTerm, sourceElementVectorID_);
        }
        
        logger(INFO, "- Creating source term vectors for the boundary faces.");
        for (Base::Face *ptrFace : this->meshes_[0]->getFacesList())
        {
            LinearAlgebra::MiddleSizeVector sourceTerm(this->integrateSourceTermAtFace(ptrFace, 0, 0));
            ptrFace->setFaceVector(sourceTerm, sourceFaceVectorID_);
        }
    }

    template<std::size_t DIM>
    void HpgemAPILinearSteadyState<DIM>::tasksBeforeSolving()
    {
        logger(INFO, "Computing stiffness matrices.");
        this->createStiffnessMatrices();
        logger(INFO, "Computing source terms.");
        createSourceTerms();
    }
    
    /// \details Solve the linear problem \f$Ax=b\f$, where \f$A=-S\f$ with \f$S\f$ being the stiffness matrix. At the moment it is not possible to solve the steady-state problem using an intial solution.
    template<std::size_t DIM>
    void HpgemAPILinearSteadyState<DIM>::solveSteadyStateWithPetsc(bool doComputeError)
    {
#if defined(HPGEM_USE_PETSC) || defined(HPGEM_USE_COMPLEX_PETSC)
        // Create output files for Paraview.
        std::string outputFileNameVTK = this->outputFileName_;
        
        this->registerVTKWriteFunctions();
        Output::VTKTimeDependentWriter<DIM> VTKWriter(outputFileNameVTK, this->meshes_[0]);
        
        // Create output files for Tecplot.
#ifdef HPGEM_USE_MPI
        std::string outputFileName = this->outputFileName_ + "." + std::to_string(Base::MPIContainer::Instance().getProcessorID());
#else
        std::string outputFileName = this->outputFileName_;
#endif
        std::string outputFileNameTecplot = outputFileName + ".dat";
        std::string dimensionsToWrite = "";
        for(std::size_t i = 0; i < this->configData_->dimension_; i++)
        {
            dimensionsToWrite = dimensionsToWrite + std::to_string(i);
        }
        
        std::string variableString = this->variableNames_[0];
        for(std::size_t iV = 1; iV < this->variableNames_.size(); iV++)
        {
            variableString = variableString + "," + this->variableNames_[iV];
        }
        
        std::ofstream outputFile(outputFileNameTecplot);
        Output::TecplotDiscontinuousSolutionWriter<DIM> tecplotWriter(outputFile, this->internalFileTitle_, dimensionsToWrite, variableString);
        
        // Create and Store things before solving the problem.
        tasksBeforeSolving();
        
        // Solve the linear problem
        //Assemble the matrix A of the system Ax = b.
        Utilities::GlobalPetscMatrix A(HpgemAPIBase<DIM>::meshes_[0], this->stiffnessElementMatrixID_, this->stiffnessFaceMatrixID_);
        MatScale(A,-1);
        //Declare the vectors x and b of the system Ax = b.
        Utilities::GlobalPetscVector b(HpgemAPIBase<DIM>::meshes_[0], sourceElementVectorID_, sourceFaceVectorID_), x(HpgemAPIBase<DIM>::meshes_[0]);

        //Assemble the vector b. This is needed because Petsc assumes you don't know
        //yet whether a vector is a variable or right-hand side the moment it is
        //declared.
        b.assemble();
        
        //Make the Krylov supspace method
        KSP ksp;
        KSPCreate(PETSC_COMM_WORLD, &ksp);
        KSPSetTolerances(ksp, 1e-12, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
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
        
        x.writeTimeIntegrationVector(this->solutionVectorId_);
        
        tecplotWriter.write(this->meshes_[0], this->solutionTitle_, false, this, 0);
        this->VTKWrite(VTKWriter, 0, this->solutionVectorId_);
        
        // Compute the energy norm of the error
        if(doComputeError)
        {
            LinearAlgebra::MiddleSizeVector::type totalError = this->computeTotalError(this->solutionVectorId_, 0);
            logger(INFO, "Total error: %.", totalError);
            LinearAlgebra::MiddleSizeVector maxError = this->computeMaxError(this->solutionVectorId_, 0);
            logger.assert(maxError.size() == this->configData_->numberOfUnknowns_, "Size of maxError (%) not equal to the number of variables (%)", maxError.size(), this->configData_->numberOfUnknowns_);
            for(std::size_t iV = 0; iV < this->configData_->numberOfUnknowns_; iV ++)
            {
                logger(INFO, "Maximum error %: %", this->variableNames_[iV], maxError(iV));
            }
        }
        
        return;
#endif
        logger(ERROR, "Petsc is needed to solve the steady state problem using this function (solveSteadyStateWithPetsc). Please put if(hpGEM_USE_PETSC) in the CMakeLists.txt of your application to make this clearer to other users");
    }
}




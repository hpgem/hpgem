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

#include "HpgemAPINonLinearSteadyState.h"
#include "Utilities/GlobalVector.h"

#if defined(HPGEM_USE_SUNDIALS)
	#include "kinsol/kinsol.h"
	#include "kinsol/kinsol_spgmr.h"
	#include "nvector/nvector_serial.h"
#endif

#include "Logger.h"

namespace Base
{
	//todo:write documentation
	//todo: Rewrite in such a way it is readable
	//todo: Butcher table is not required, remove this from the constructor
	//todo: Add a bool that can switch if you want to store the intermediate solutions or only the last solution
	//todo: Initial condition is optional; i.e. if not specified everything is 0
	//todo: Possibly add restart from data file to create a p-multigrid solution strategy


	//note: compute error is always on, compute both faces is always on
	//NOTE: Temporary just 10 coefficients for testing
	template<std::size_t DIM>
	HpgemAPINonLinearSteadyState<DIM>::HpgemAPINonLinearSteadyState
	(
			const std::size_t numberOfVariables,
			const std::size_t polynomialOrder,
			const TimeIntegration::ButcherTableau * const ptrButcherTableau,
			const bool computeBothFaces
	) :
	HpgemAPISimplified<DIM>(numberOfVariables, polynomialOrder, ptrButcherTableau, 1, computeBothFaces),
    sourceElementVectorID_(0),
    sourceFaceVectorID_(0)
	{
	}

#if defined(HPGEM_USE_SUNDIALS)
	//static int KINSol function required function to compute the right hand side.
	//The user_data is recast into the API class and a corresponding memeber function is then called to compute the RHS
	template<std::size_t DIM>
	int HpgemAPINonLinearSteadyState<DIM>::func(N_Vector u, N_Vector fval, void *user_data)
	{
		return static_cast<HpgemAPINonLinearSteadyState<DIM>*>(user_data)->computeRHS(u, fval);
	}

	//Member function computing the rhs
	//todo: Add the option to write intermediate solutions to VTK output
	template<std::size_t DIM>
	int HpgemAPINonLinearSteadyState<DIM>::computeRHS(N_Vector u, N_Vector fval)
	{

		std::cout << "YOLO" << std::endl;

		//Pass the solution vector u to the GlobalSundialsVector
		globalVector_->setVector(u);
		//The GlobalSundialsVector will then put the data in the hpGEM data structure
		globalVector_->writeTimeIntegrationVector(this->solutionVectorId_);

		//Compute RHS with the new data
		//todo: mention that computeRightHandSide is very confusing on naming
		this->computeRightHandSide(0, 1, 0);

		//Compute intermediate solutions
		if (doOutputIntermediateSolutions_)
		{
			logger(INFO,"Writing intermediate solution:");
			tecplotWriter_->write(this->meshes_[0], this->solutionTitle_, false, this, 0);
			this->VTKWrite(*(this->VTKWriter_), 0, this->solutionVectorId_);
		}

		//Set the correct N_Vector pointer
		globalVector_->setVector(fval);
		//Pass the new solution back to KINsol
		globalVector_->constructFromTimeIntegrationVector(0);

		return 0;
	}
#endif

	template<std::size_t DIM>
	bool HpgemAPINonLinearSteadyState<DIM>::solve(bool doComputeInitialCondition, bool doComputeError)
	{
#if defined(HPGEM_USE_SUNDIALS)
		std::cout << "Hello!" << std::endl;
		int flag;
		int maxl;
		int globalStrategy = KIN_NONE; // For now: nothing special, in future maybe linesearch

        // Create output files for Paraview.
        std::string outputFileNameVTK = this->outputFileName_;

        this->registerVTKWriteFunctions();
        Output::VTKTimeDependentWriter<DIM> VTKWriter(outputFileNameVTK, this->meshes_[0]);
        VTKWriter_ = &VTKWriter;

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
        tecplotWriter_ = &tecplotWriter;

		//Create KINsol memory block
		void *kmem = KINCreate();
		logger.assert_always(kmem != nullptr,"Failed creating KINmem.");

		//Initialise Global Solution Vector and template N_Vector u
		Utilities::GlobalSundialsVector globalVector(HpgemAPIBase<DIM>::meshes_[0], sourceElementVectorID_, sourceFaceVectorID_);
		globalVector_ = &globalVector;
		N_Vector u = globalVector;
		int numberOfDOF = NV_LENGTH_S(u);

		//Create the scale vector. For now this is just all set to 1
		N_Vector scale = N_VNew_Serial(numberOfDOF);
		double *ptrScale = NV_DATA_S(scale);
		for (std::size_t i = 0; i < numberOfDOF; i++)
		{
			ptrScale[i] = 1.0;
		}

		//If an initial condition is supplied, compute it and initialise u with this vector, else the vector is initialised as 0
		if (doComputeInitialCondition == true)
		{
			//Let hpGEM compute the initial coefficients
			this->setInitialSolution(this->solutionVectorId_, 0, 0);

			//obtain globalvector for initial conditions
			globalVector_->setVector(u);
			globalVector_->constructFromTimeIntegrationVector(0);
		}

		///DEBUG
		//Write final solution to VTK
		tecplotWriter.write(this->meshes_[0], this->solutionTitle_, false, this, 0);
		this->VTKWrite(VTKWriter, 0, this->solutionVectorId_);
		///DEBUG

		//Give userdata to KINSol, as a void pointer
		void * userdata = static_cast<void*>(this);
		flag = KINSetUserData(kmem, userdata);
		logger.assert_always(flag >= 0, "User data not correctly passed to KINSol.");

		//Initialise kmem
		flag = KINInit(kmem, func, u);
		logger.assert_always(flag >= 0, "Initialisation failed with flag %.", flag);

		//Set type of solver
		maxl = 15;
		flag = KINSpgmr(kmem, maxl);
		logger.assert(flag >= 0, "Failed to attach GMRS solver.");

		//Set parameters

		//Call the KINsol function
		flag = KINSol(kmem, u, globalStrategy, scale, scale);
		logger.assert(flag >= 0, "Failed to solve the problem with flag %.",flag);
		logger(WARN,"Warning: An event happened during the problem solving, with flag %.", flag);

		//Write final solution to VTK
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

		//Free memory
		KINFree(&kmem);
		N_VDestroy_Serial(scale);

		return false;
#endif
		logger(ERROR, "Sundials is required to solve the non linear steady state problem using this function.");
	}
}

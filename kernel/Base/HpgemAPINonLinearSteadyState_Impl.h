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
	template<std::size_t DIM>
	int HpgemAPINonLinearSteadyState<DIM>::func(N_Vector u, N_Vector fval, void *user_data)
	{
		return static_cast<HpgemAPINonLinearSteadyState<DIM>*>(user_data)->computeRHS(u, fval);
	}


	//Member function computing the rhs

	template<std::size_t DIM>
	int HpgemAPINonLinearSteadyState<DIM>::computeRHS(N_Vector u, N_Vector fval)
	{

		std::cout << "YOLO" << std::endl;

		//todo: write this functionality
		//Utilities::GlobalSundialsVector solution = u;

		//Compute the right hand side with the new solution
		//computeRightHandSide(inputVectorIds, coefficientsInputVectors, resultVectorId, 0);

		//Pass the new solution back to KINsol
		//fval =

		return 0;
	}
#endif

	template<std::size_t DIM>
	bool HpgemAPINonLinearSteadyState<DIM>::solve()
	{
#if defined(HPGEM_USE_SUNDIALS)
		std::cout << "Hello!" << std::endl;
		int flag;
		int maxl;
		int globalStrategy = KIN_NONE; // For now: nothing special, in future maybe linesearch

		//Create KINsol memory block
		void *kmem = KINCreate();
		logger.assert_always(kmem == nullptr,"Failed creating KINmem.");

		//Initialise KINsol
		Utilities::GlobalSundialsVector uGV(HpgemAPIBase<DIM>::meshes_[0], sourceElementVectorID_, sourceFaceVectorID_);
		Utilities::GlobalSundialsVector scaleGV(HpgemAPIBase<DIM>::meshes_[0]);
		//If an initial condition is supplied, compute it and initialise u with this vector, else the vector is initialised as 0
		if (computeInitialCondition_ == true)
		{
			//Let hpGEM compute the initial coefficients
			this->setInitialSolution(this->solutionVectorId_, 0, 0);

			//obtain globalvector for initial conditions
		}

		//compute scale vector: currently set everything to 1
		N_Vector scale = scaleGV;
		N_Vector u = uGV;
		int length = NV_LENGTH_S(scale);
		double *ptrScale = NV_DATA_S(scale);
		for (std::size_t i = 0; i < length; i++)
		{
			ptrScale[i] = 1.0;
		}

		//Set the pointer of the class
		void * userdata = static_cast<void*>(this);
		flag = KINSetUserData(kmem, userdata);
		logger.assert_always(flag < 0, "User data not correctly passed to KINSol.");

		//Initialise kmem
		flag = KINInit(kmem, func, u);
		logger.assert_always(flag < 0, "Initialisation failed with flag %.", flag);

		//Set type of solver
		maxl = 15;
		flag = KINSpgmr(kmem, maxl);
		logger.assert_always(flag < 0, "Failed to attach GMRS solver.");

		//Set parameters

		//Call the KINsol function
		flag = KINSol(kmem, u, globalStrategy, scale, scale);
		if (flag < 0)
		{
			logger(ERROR,"Failed to solve the problem with flag %.",flag);
		}
		else if (flag > 0)
		{
			logger(WARN,"An error occurred during the problem solving, with flag %.", flag);
		}


		//Free memory
		KINFree(&kmem);

		return false;
#endif
		logger(ERROR, "Sundials is required to solve the non linear steady state problem using this function.");
	}
}

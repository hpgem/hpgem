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

#ifndef BaseNonLinearSteadyStateH
#define BaseNonLinearSteadyStateH

#include "Base/HpgemAPISimplified.h"
#include "Base/TimeIntegration/AllTimeIntegrators.h"
#include "Utilities/GlobalVector.h"

#if defined(HPGEM_USE_SUNDIALS)
	#include "nvector/nvector_serial.h"
#endif


namespace Base
{
	//todo: Write documentation

	template<std::size_t DIM>
	class HpgemAPINonLinearSteadyState : public HpgemAPISimplified<DIM>
	{
	public:

		//Constructor
		HpgemAPINonLinearSteadyState
		(
				const std::size_t numberOfUnknowns,
				const std::size_t polynomialOrder,
				const TimeIntegration::ButcherTableau * const ptrButcherTableau,
				const bool computeBothFaces
		);

#if defined(HPGEM_USE_SUNDIALS)
		//todo: Write a function that computes the RHS, creates the vector and can be used by KINsol
		static int func(N_Vector cc, N_Vector fval, void *user_data);

		//Member function computing the rhs
		int computeRHS(N_Vector u, N_Vector fval);
#endif

		//todo: Write a function that computes the Jacobian, creates a sparse matrix and can be used by KINsol

		//todo: Write a function that solves a steady state solution of a non linear PDE using KINsol
		virtual bool solve(bool doComputeInitialCondition, bool doComputeError);

    protected:
        /// Index to indicate where the vectors for the source terms for the elements are stored.
        const std::size_t sourceElementVectorID_;

        /// Index to indicate where the vectors for the source terms for the faces are stored.
        const std::size_t sourceFaceVectorID_;

    private:

        //Flag to output intermediate solutions
        bool doOutputIntermediateSolutions_ = true;

#if defined(HPGEM_USE_SUNDIALS)
        //Global Vector and global matrix
        Utilities::GlobalSundialsVector *globalVector_;
        //Utilities::GlobalSundialsMatrix jacobianMatrix_;
#endif

        //tecplotwriter
        Output::TecplotDiscontinuousSolutionWriter<DIM> *tecplotWriter_;
        Output::VTKTimeDependentWriter<DIM> *VTKWriter_;

	};

/*    std::ostream& operator<<(std::ostream &os, const N_Vector &A)
    {
    	std::size_t length = NV_LENGTH_S(A);
    	double *data = NV_DATA_S(A);
            os << '[';
        for (std::size_t i = 0; i < length; i++)
        {
            if(i != 0)
            {
                os << ", ";
            }
            os << data[i];
        }
        os << ']';
        return os;
    }*/

}

#include "HpgemAPINonLinearSteadyState_Impl.h"

#endif

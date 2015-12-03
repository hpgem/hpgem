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
#include "Base/FaceMatrix.h"
#include "Utilities/GlobalVector.h"

#if defined(HPGEM_USE_SUNDIALS)
	#include "nvector/nvector_serial.h"
#endif


namespace Base
{
	/// \brief Interface for solving steady-state solutions of non-linear PDE's. At the moment this class can only solve steady state problems using Sundials.
	/** The API is well suited to solve steady-state problems of non-Linear PDE's. The API solves problems of the form \f[ F(u) = 0 \f]. Where the function F(u) is
	 * the right hand side of the system of equations that has to be solved and u are the solution coefficients.
	 * To use the API, only the functions are required to compute the RHS must be user supplied. Jacobian information is not required, since internally a finite difference approximation
	 * will be constructed. In future however Jacobian information will be added to the API
	 */

	/** \details To solve some steady-state non-linear PDE with this class you should at least do the following:
	 * \li Create your own class that inherits this class.
	 * \li Implement the function 'createMeshDescription' to create a mesh description (e.g. domain, number of elements, etc.).
	 * \li Implement the function 'computeRightHandSideAtFace' to compute the face contribution to the rhs vector \f[ F(u) \f].
	 * \li Implement the function 'computeRightHandSideAtElement' to get the element contribution to the rhs vector \f[ F(u) \f].
	 */
/** \details To solve the PDE do the following in the main routine:
 * \li Create an object of your own class, that inherits from this class and has the necessary functions implemented (see list above).
 * \li Call the function 'CreateMesh' to create the mesh.
 * \li Call the function 'setOutputNames' to set the names for the output files.
 * \li Call the function 'solve'.
 */
/** \details An example using this interface is still under construction.
 */

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

        HpgemAPINonLinearSteadyState(const HpgemAPINonLinearSteadyState &other) = delete;

        /// \brief Create the mesh.
        virtual void createMesh(const std::size_t numberOfElementsPerDirection, const Base::MeshType meshType) override;

#if defined(HPGEM_USE_SUNDIALS)
		/// \brief This function is called by the Sundials library and then calls the function int computeRHS(N_Vector u, N_Vector fval);
		static int func(N_Vector cc, N_Vector fval, void *user_data);

		//todo: write this function
		static int jtimes(N_Vector v, N_Vector Jv, N_Vector u, booleantype *new_u, void *user_data);

		/// \brief This function computes the right hand side of the system of equations. N_Vector u contains the supplied solution coefficients
		/// N_Vector fval will contain the computed rhs
		int computeRHS(N_Vector u, N_Vector fval);

		//todo: write this function
		int computeJacTimesVector(N_Vector v, N_Vector Jv, N_Vector u, booleantype *new_u);
#endif

		/// \brief This function computes the local and non-local Jacobian contributions from face i of the face integrals
		virtual Base::FaceMatrix computeJacobianAtFace(Base::Face *ptrFace, LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeft, LinearAlgebra::MiddleSizeVector &solutionCoefficientsRight, const double time)
		{
            logger(ERROR, "No function computeJacobianFaceAtElement() implemented to compute the local and non local contributions from face integrals.");
            Base::FaceMatrix jacobianAtFace;
            return jacobianAtFace;
		}

		/// \brief This function computes the local Jacobian contributions from the element integrals
		virtual LinearAlgebra::MiddleSizeMatrix computeJacobianAtElement(Base::Element *ptrElement, const LinearAlgebra::MiddleSizeVector &solutionCoefficients, const double time)
		{
            logger(ERROR, "No function computeJacobianAtElement() implemented to compute the local contributions of the Jacobian matrix from element integrals.");
            LinearAlgebra::MiddleSizeMatrix jacobianAtElement(1,1);
            return jacobianAtElement;
		}


		/// \brief Solve the steady-state problem using Sundials
		virtual bool solve(bool doComputeInitialCondition, bool doComputeError, bool doUseJacobian);

    protected:
        /// Index to indicate where the vectors for the source terms for the elements are stored.
        const std::size_t sourceElementVectorID_;

        /// Index to indicate where the vectors for the source terms for the faces are stored.
        const std::size_t sourceFaceVectorID_;

    private:

        /// \brief for Debugging purposes or curiousity. If the flag is true it will output intermediate solutions.
        bool doOutputIntermediateSolutions_ = true;

        /// \brief Location of the computed RHS
        const std::size_t solutionRHS_ = 1;
        /// \brief Location of the jacobian times vector v in the solution
        const std::size_t jTimesvecID_ = 2;
        /// \brief Location of the vector V which will be multiplied by the jacobian matrices
        const std::size_t VectorVID_ = 3;

        /// \brief Location of the Jacobian Element Matrix in the elementMatrix vector
        const std::size_t jacobianElementMatrixID_ = 0;
        /// \brief Location of the faceMatrix
        const std::size_t jacobianFaceMatrixID_ = 0;

#if defined(HPGEM_USE_SUNDIALS)
        /// GlobalVector that performs operations on a given N_Vector. Either writes solutionCoefficients to the hpGEM data structure,
        /// or obtains the RHS from the hpGEM data structure
        Utilities::GlobalSundialsVector *globalVector_;

#endif

        /// tecplotwriter
        Output::TecplotDiscontinuousSolutionWriter<DIM> *tecplotWriter_;
        Output::VTKTimeDependentWriter<DIM> *VTKWriter_;
        double step_ = 0.0;

	};
}

#include "HpgemAPINonLinearSteadyState_Impl.h"

#endif
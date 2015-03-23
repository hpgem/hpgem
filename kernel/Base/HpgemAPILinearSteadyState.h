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

#ifndef BaseLinearSteadyStateH
#define BaseLinearSteadyStateH

#include "Base/ConfigurationData.h"
#include "Base/Element.h"
#include "Base/Face.h"
#include "Base/HpgemAPILinear.h"
#include "Base/RectangularMeshDescriptor.h"
#include "Base/TimeIntegration/AllTimeIntegrators.h"
#include "Integration/ElementIntegrandBase.h"
#include "Integration/FaceIntegrandBase.h"
#include "Output/TecplotSingleElementWriter.h"
#include "Output/VTKTimeDependentWriter.h"
#include <functional>

namespace Base
{
    /// \brief Interface for solving steady-state solutions of linear PDE's.
    /** At the moment this class is well-suited for steady-state problems corresponding to PDE's of the form \f[ l(\partial_t^k \vec{u}) = f + g(\vec{u}) \f], where \f$ \vec{u} \f$ is some vector function, \f$ l(\partial_t^k \vec{u})\f$ is some linear function, applied on the k-th order time-derivative of \f$ u \f$, and \f$ f(t) + g(\vec{u}) \f$ is the right-hand side, which is some linear function of \f$ \vec{u} \f$ (that can depend on arbitrary order spatial derivatives of \f$\vec{u}\f$) plus some source term \f$ f(t) \f$. The steady state problem has the form \f[ -g(\vec{u}) = f \f]. The resulting set of ODE's will have the form \f[ -Su = f \f], where \f$S\f$ is the stiffness matrix and \f$f\f$ is the source term. At the moment this class can only solve steady state problems using Petsc.
     */
    /** \details To solve some linear time depent PDE with this class you should at least do the following:
     * \li Create your own class that inherits this class.
     * \li Implement the function 'createMeshDescription' to create a mesh description (e.g. domain, number of elements, etc.).
     * \li Implement the function 'getSourceTerm' to define the source term (e.g. external force) if there is one.
     * \li Implement the function 'getSourceTermAtBoundary' to get the source term (e.g. because of a boundary force / boundary condition) if there is one.
     * \li Implement the function 'integrateSourceTermAtFace' for integrating the source term at a boundary face.
     * \li Implement the functions 'computeStiffnessMatrixAtElement' and 'computeStiffnessMatrixAtFace' for computing the stiffness matrix at an element or face. One can also choose to implement the functions 'computeIntegrandStiffnessMatrixAtElement' and 'computeIntegrandStiffnessMatrixAtFace' for computing the integrands at the face and element. The integration will be done by an automatic routine.
     */
    /** \details To solve the PDE do the following in the main routine:
     * \li Create an object of your own class, that inherits from this class and has the necessary functions implemented (see list above).
     * \li Call the function 'CreateMesh' to create the mesh.
     * \li Call the function 'setOutputNames' to set the names for the output files.
     * \li Call the function 'solveSteadyStateWithPetsc'.
     */
    /** \details Some other thinsgs you can do:
     * \li Implement the function 'getExactSolution' if you know the analytic solution and want to compute the error.
     * \li Implement the function 'integrateInitialSolutionAtElement' for integrating the initial solution at the element (by default this function computes the standard L2 inner product).
     * \li Implement the function 'computeMassMatrixAtElement' if you want to compute the mass matrix (by default a mass matrix is computed based on the L2 inner product).
     * \li Implement the function 'integrateSourceTermAtElement' to compute the source term at an element (if there is a source term, by default the L2 inner product is computed).
     * \li Implement the function 'integrateErrorAtElement' to compute the square of some user-defined norm of the error at an element (by default the L2-norm is computed).
     * \li Override the function 'writeToTecplotFile' to determine what data to write to the output file.
     */
    /** \details For an example of using this interface see the application 'TutorialPoisson'.
     */
    
    class HpgemAPILinearSteadyState : public HpgemAPILinear
    {
    public:
        // Constructor
        HpgemAPILinearSteadyState
        (
         const std::size_t dimension,
         const std::size_t numberOfUnknowns,
         const std::size_t polynomialOrder,
         const bool useSourceTerm = false,
         const bool useSourceTermAtBoundary = false
         );
        
        /// \brief Create the mesh.
        virtual void createMesh(const std::size_t numOfElementsPerDirection, const Base::MeshType meshType) override;
        
        /// \brief Compute the source term at a given physical point.
        virtual LinearAlgebra::NumericalVector getSourceTerm(const PointPhysicalT &pPhys)
        {
            logger(ERROR, "No source term implemented.");
            LinearAlgebra::NumericalVector sourceTerm;
            return sourceTerm;
        }
        
        /// \brief Compute the source term at a given physical point.
        virtual LinearAlgebra::NumericalVector getSourceTerm(const PointPhysicalT &pPhys, const double &time, const std::size_t orderTimeDerivative) override
        {
            return getSourceTerm(pPhys);
        }
        
        /// \brief Get the source term at the boundary at a given physical point.
        /// \details The source term at the boundary can be a result of certain boundary conditions (e.g. Neumann boundary conditions).
        virtual LinearAlgebra::NumericalVector getSourceTermAtBoundary(const PointPhysicalT &pPhys)
        {
            
            logger(ERROR, "No source term at the boundary implemented.");
            LinearAlgebra::NumericalVector sourceTerm;
            return sourceTerm;
        }
        
        /// \brief Get the source term at the boundary at a given physical point.
        /// \details The source term at the boundary can be a result of certain boundary conditions (e.g. Neumann boundary conditions).
        virtual LinearAlgebra::NumericalVector getSourceTermAtBoundary(const PointPhysicalT &pPhys, const double &time, const std::size_t orderTimeDerivative) override
        {
            
            return getSourceTermAtBoundary(pPhys);
        }
        
        /// \brief Create and store the source terms.
        virtual void createSourceTerms();
        
        /// \brief Create and store things before solving the problem.
        virtual void tasksBeforeSolving() override;
        
        /// \brief Solve the steady-state problem using Petsc.
        virtual void solveSteadyStateWithPetsc();
        
    protected:
        /// Index to indicate where the vectors for the source terms for the elements are stored.
        const std::size_t sourceElementVectorID_;
        
        /// Index to indicate where the vectors for the source terms for the faces are stored.
        const std::size_t sourceFaceVectorID_;
    };
}

#endif




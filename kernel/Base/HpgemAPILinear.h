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

#ifndef BaseLinearH
#define BaseLinearH

#include "Base/ConfigurationData.h"
#include "Base/Element.h"
#include "Base/Face.h"
#include "Base/HpgemAPISimplified.h"
#include "Base/RectangularMeshDescriptor.h"
#include "Base/TimeIntegration/AllTimeIntegrators.h"
#include "Integration/ElementIntegrandBase.h"
#include "Integration/FaceIntegrandBase.h"
#include "Output/TecplotSingleElementWriter.h"
#include "Output/VTKTimeDependentWriter.h"
#include <functional>

namespace Base
{
    /// \brief Simplified Interface for solving linear PDE's.
    /** At the moment this class is well-suited for problems of the form \f[ l(\partial_t^k \vec{u},t) = f(\vec{u},t) \f], where \f$ \vec{u} \f$ is some vector function, \f$ l(\partial_t^k \vec{u},t)\f$ is some linear function, applied on the k-th order time-derivative of \f$ u \f$, and \f$ f(\vec{u},t) \f$ is the right-hand side, which is some linear function of \f$ \vec{u} \f$ (that can depend on arbitrary order spatial derivatives of \f$\vec{u}\f$) plus some source term. The resulting set of ODE's will have the form \f[ M\partial_t^ku = Au + f(t)\f], where \f$A\f$ is the stiffness matrix and \f$f(t)\f$ is the source term. If you do not want to store all element- and face matrices for the mass matrix and stiffness matrix, it is advised to use the superclass HpgemAPISimplified instead.
     */
    /** \details To solve some linear time depent PDE with this class you should at least do the following:
     * \li Create your own class that inherits this class.
     * \li Implement the function 'createMeshDescription' to create a mesh description (e.g. domain, number of elements, etc.).
     * \li Implement the function 'getInitialConditions' to define the initial condition(s) of your problem.
     * \li Implement the function 'getSourceTerm' to define the source term (e.g. external force) if there is one.
     * \li Implement the functions 'computeStiffnessMatrixAtElement' and 'computeStiffnessMatrixAtFace' for computing the stiffness matrix at an element or face. One can also choose to implement the functions 'computeIntegrandStiffnessMatrixAtElement' and 'computeIntegrandStiffnessMatrixAtFace' for computing the integrands at the face and element. The integration will be done by an automatic routine. 
     * \li Implement the function 'integrateSourceTermAtElement' to compute the source term at an element (if there is a source term).
     */
    /** \details To solve the PDE do the following in the main routine:
     * \li Create an object of your own class, that inherits from this class and has the necessary functions implemented (see list above).
     * \li Call the function 'CreateMesh' to create the mesh.
     * \li Call the function 'setOutputNames' to set the names for the output files.
     * \li Call the function 'solve'.
     */
    /** \details Some other thinsgs you can do:
     * \li Implement the function 'getExactSolution' if you know the analytic solution and want to compute the error.
     * \li Implement the function 'integrateInitialSolutionAtElement' for integrating the initial solution at the element (by default this function computes the standard L2 inner product).
     * \li Implement the function 'computeMassMatrixAtElement' if you want to compute the mass matrix (by default a mass matrix is computed based on the L2 norm).
     * \li Implement the function 'integrateErrorAtElement' to compute the square of some user-defined norm of the error at an element (by default the L2-norm is computed).
     * \li Override the function 'writeToTecplotFile' to determine what data to write to the output file.
     * \li Override the function 'showProgress' to determine how you want to show the progress of the time integration routine.
     * \li Override the function 'solve' when using another time integration routine than a Runge-Kutta integration method.
     */
    /** \details For an example of using this interface see the application 'ExampleMultipleVariableProblem'.
     */
    
    class HpgemAPILinear : public HpgemAPISimplified
    {
    public:
        // Constructor
        HpgemAPILinear
        (
         const std::size_t dimension,
         const std::size_t numberOfUnknowns,
         const std::size_t polynomialOrder,
         const Base::ButcherTableau * const ptrButcherTableau = Base::AllTimeIntegrators::Instance().getRule(4, 4),
         const std::size_t numOfTimeLevels = 1,
         const bool useSourceTerm = false
         );
        
        /// \brief Create the mesh.
        virtual void createMesh(const std::size_t numOfElementsPerDirection, const Base::MeshType meshType);
        
        /// \brief Compute the source term at a given physical point.
        virtual double getSourceTerm(const PointPhysicalT &pPhys, const double &time, const std::size_t orderTimeDerivative)
        {
            logger(ERROR, "No source term implemented.");
            return 0.0;
        }
        
        /// \brief Compute and store the mass matrices.
        virtual void createMassMatrices();
        
        /// \brief Solve the mass matrix equations.
        virtual void solveMassMatrixEquations(const std::size_t timeLevel) override;
        
        /// \brief Compute the integrand for the stiffness matrix.
        virtual LinearAlgebra::Matrix computeIntegrandStiffnessMatrixAtElement(const Base::Element *ptrElement, const Geometry::PointReference &pRef)
        {
            logger(ERROR, "No function for computing the integrand for the stiffness matrix at an element implemented.");
            LinearAlgebra::Matrix integrandStiffnessMatrix;
            return integrandStiffnessMatrix;
        }
        
        /// \brief Compute the stiffness matrix corresponding to an element.
        virtual LinearAlgebra::Matrix computeStiffnessMatrixAtElement(Base::Element *ptrElement);
        
        /// \brief Compute the integrand for the stiffness matrix.
        virtual Base::FaceMatrix computeIntegrandStiffnessMatrixAtFace(const Base::Face *ptrFace, const LinearAlgebra::NumericalVector &normal, const Geometry::PointReference &pRef)
        {
            logger(ERROR, "No function for computing the integrand for the stiffness matrix at a face implemented.");
            Base::FaceMatrix integrandStiffnessMatrix;
            return integrandStiffnessMatrix;
        }
        
        /// \brief Compute the stiffness matrix corresponding to a face.
        virtual Base::FaceMatrix computeStiffnessMatrixAtFace(Base::Face *ptrFace);
        
        /// \brief Compute and store stiffness matrices for computing the right hand side.
        virtual void createStiffnessMatrices();
        
        /// \brief Integrate the source term at a single element.
        virtual LinearAlgebra::NumericalVector integrateSourceTermAtElement(Base::Element * ptrElement, const double startTime, const std::size_t orderTimeDerivative)
        {
            logger(ERROR, "No function for computing the integral for the source term at an element implemented.");
            LinearAlgebra::NumericalVector integralSourceTerm;
            return integralSourceTerm;
        }
        
        /// \brief Multiply the stiffness matrices with the solution at time level 'timeLevelIn' and store the result at time level 'timeLevelResult'.
        virtual void multiplyStiffnessMatrices(const std::size_t timeLevelIn, const std::size_t timeLevelResult);
        
        /// \brief Multiply the stiffness matrices with the linear combination of solutions at time level 'timeLevelIn' with coefficients given in coefficientsTimeLevels. Store the result at time level 'timeLevelResult'.
        virtual void multiplyStiffnessMatrices(const std::vector<std::size_t> timeLevelsIn, const std::vector<double> coefficientsTimeLevels, const std::size_t timeLevelResult);
        
        /// \brief Add the source term to the solution at time level 'timeLevelResult'.
        virtual void addSourceTerm(const std::size_t timeLevelResult, const double time, const std::size_t orderTimeDerivative);
        
        /// \brief Compute the right hand side for the solution at time level 'timeLevelIn' and store the result at time level 'timeLevelResult'.
        virtual void computeRightHandSide(const std::size_t timeLevelIn, const std::size_t timeLevelResult, const double time) override;
        
        /// \brief Compute the right hand side for the linear combination of solutions at time level 'timeLevelIn' with coefficients given in coefficientsTimeLevels. Store the result at time level 'timeLevelResult'.
        virtual void computeRightHandSide(const std::vector<std::size_t> timeLevelsIn, const std::vector<double> coefficientsTimeLevels, const std::size_t timeLevelResult, const double time) override;
        
        /// \brief Create and Store things before solving the problem.
        virtual void tasksBeforeSolving() override;
        
    protected:
        /// Boolean to indicate if there is a source term.
        const bool useSourceTerm_;
        
        /// Index to indicate where the mass matrix is stored
        const std::size_t massMatrixID_;
        
        /// Index to indicate where the stiffness matrix for the elements is stored
        const std::size_t stiffnessElementMatrixID_;
        
        /// Index to indicate where the stiffness matrix for the elements is stored
        const std::size_t stiffnessFaceMatrixID_;
    };
}

#endif



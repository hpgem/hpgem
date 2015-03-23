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

#ifndef BaseSimplifiedHPP
#define BaseSimplifiedHPP

#include "Base/ConfigurationData.h"
#include "Base/Element.h"
#include "Base/Face.h"
#include "Base/HpgemAPIBase.h"
#include "Base/RectangularMeshDescriptor.h"
#include "Base/TimeIntegration/AllTimeIntegrators.h"
#include "Integration/ElementIntegrandBase.h"
#include "Integration/FaceIntegrandBase.h"
#include "Integration/ElementIntegral.h"
#include "Integration/FaceIntegral.h"
#include "Output/TecplotSingleElementWriter.h"
#include "Output/VTKTimeDependentWriter.h"
#include <functional>

namespace Integration
{
    class FaceIntegral;
}

namespace Base
{
    extern CommandLineOption<double>& startTime;
    extern CommandLineOption<double>& endTime;
    extern CommandLineOption<double>& dt;
    extern CommandLineOption<std::size_t>& numberOfSnapshots;
    
    /// \brief Simplified Interface for solving PDE's.
    /** This class is well-suited for problems of the form \f[ l(\partial_t^k \vec{u}) = f(\vec{u},t) \f], where \f$ \vec{u} \f$ is some vector function, \f$ l(\partial_t^k \vec{u})\f$ is some linear function, applied on the k-th order time-derivative of \f$ u \f$, and \f$ f(\vec{u},t) \f$ is some function of \f$ \vec{u} \f$ that can depend on arbitrary order spatial derivatives of \f$\vec{u}\f$. This last term will be referred to as the right-hand side. The resulting set of ODE's will have the form \f[ M\partial_t^ku = f(u,t)\f], where \f$ M\f$ is the mass matrix, \f$u\f$ is the numerical solution vector and \f$f(u)\f$ is the right-hand side.
     */
    /** \details To solve some linear time depent PDE with this class you should at least do the following:
     * \li Create your own class that inherits this class.
     * \li Implement the function 'createMeshDescription' to create a mesh description (e.g. domain, number of elements, etc.).
     * \li Implement the function 'getInitialConditions' to define the initial condition(s) of your problem.
     * \li Implement the functions 'computeRightHandSideAtElement' and 'computeRightHandSideAtFace' to compute the right-hand side corresponding to an element or face.
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
     * \li Override the function 'solveMassMatrixEquationsAtElement' if you want to solve the mass matrix equtions without computing the mass matrix first.
     * \li Implement the function 'integrateErrorAtElement' to compute the square of some user-defined norm of the error at an element (by default the L2-norm is computed).
     * \li Override the function 'writeToTecplotFile' to determine what data to write to the output file.
     * \li Override the function 'showProgress' to determine how you want to show the progress of the time integration routine.
     * \li Override the function 'solve' when using another time integration routine than a Runge-Kutta integration method.
     */
    /** \details For an example of using this interface see the application 'ExampleMultipleVariableProblem'.
     */
    
    class HpgemAPISimplified : public HpgemAPIBase, public Output::TecplotSingleElementWriter
    {
    public:
        // Constructor
        HpgemAPISimplified
        (
         const std::size_t dimension,
         const std::size_t numberOfUnknowns,
         const std::size_t polynomialOrder,
         const Base::ButcherTableau * const ptrButcherTableau = Base::AllTimeIntegrators::Instance().getRule(4, 4),
         const std::size_t numOfTimeLevels = 1
         );

        /// \brief Create a mesh description
        virtual Base::RectangularMeshDescriptor createMeshDescription(const std::size_t numOfElementPerDirection)
        {
            logger(ERROR, "No routine for creating the domain implemented.");
            Base::RectangularMeshDescriptor description(configData_->dimension_);
            return description;
        }
        
        /// \brief Create the mesh.
        virtual void createMesh(const std::size_t numOfElementsPerDirection, const Base::MeshType meshType);
        
        /// \brief Compute the real solution at a given point in space and time.
        virtual LinearAlgebra::NumericalVector getExactSolution(const PointPhysicalT &pPhys, const double &time, const std::size_t orderTimeDerivative)
        {
            logger(ERROR, "No real solution implemented.");
            LinearAlgebra::NumericalVector realSolution(configData_->numberOfUnknowns_);
            return realSolution;
        }
        
        /// \brief Compute the initial solution at a given point in space and time.
        virtual LinearAlgebra::NumericalVector getInitialSolution(const PointPhysicalT &pPhys, const double &startTime, const std::size_t orderTimeDerivative)
        {
            logger(ERROR, "No initial solution implemented.");
            LinearAlgebra::NumericalVector realSolution(configData_->numberOfUnknowns_);
            return realSolution;
        }
        
        /// \brief Compute the mass matrix for a single element.
        virtual LinearAlgebra::Matrix computeMassMatrixAtElement(Base::Element *ptrElement);
        
        /// \brief Solve the mass matrix equations for a single element.
        /// \details Solve the equation \f$ Mu = r \f$ for \f$ u \f$ for a single element, where \f$ r \f$ is the right-hand sid and \f$ M \f$ is the mass matrix. The input is the right hand side here called 'solutionCoefficients' and the result is returned in this same vector.
        virtual void solveMassMatrixEquationsAtElement(Base::Element *ptrElement, LinearAlgebra::NumericalVector &solutionCoefficients);

        /// \brief Solve the mass matrix equations.
        virtual void solveMassMatrixEquations(const std::size_t timeLevel);

        /// \brief Integrate the initial solution at a single element.
        virtual LinearAlgebra::NumericalVector integrateInitialSolutionAtElement( Base::Element * ptrElement, const double startTime, const std::size_t orderTimeDerivative);
        
        /// \brief Integrate the initial solution.
        void integrateInitialSolution(const std::size_t timeLevelResult, const double startTime, const std::size_t orderTimeDerivative);

        /// \brief Integrate the square of some norm of the error on a single element.
        virtual LinearAlgebra::NumericalVector integrateErrorAtElement(Base::Element *ptrElement, LinearAlgebra::NumericalVector &solutionCoefficients, const double time);
        
        /// \brief Compute the (weighted) L2-norm of the error.
        double computeTotalError(const std::size_t solutionTimeLevel, const double time);
        
        /// \brief Compute the L-infinity norm (essential supremum) of the error at an element.
        virtual LinearAlgebra::NumericalVector computeMaxErrorAtElement(Base::Element *ptrElement, LinearAlgebra::NumericalVector &solutionCoefficients, const double time);
        
        /// \brief Compute the L-infinity norm (essential supremum) of the error.
        LinearAlgebra::NumericalVector computeMaxError(const std::size_t solutionTimeLevel, const double time);
        
        /// \brief Compute the right-hand side corresponding to an element
        virtual LinearAlgebra::NumericalVector computeRightHandSideAtElement
        (
         Base::Element *ptrElement,
         LinearAlgebra::NumericalVector &solutionCoefficients,
         const double time
        )
        {
            logger(ERROR, "No function for computing the right-hand side at an element implemented.");
            LinearAlgebra::NumericalVector rightHandSideElement;
            return rightHandSideElement;
        }
        
        /// \brief Compute the right-hand side corresponding to a face
        virtual LinearAlgebra::NumericalVector computeRightHandSideAtFace
        (
         Base::Face *ptrFace,
         const Base::Side side,
         LinearAlgebra::NumericalVector &solutionCoefficientsLeft,
         LinearAlgebra::NumericalVector &solutionCoefficientsRight,
         const double time
         )
        {
            logger(ERROR, "No function for computing the right-hand side at a face implemented.");
            LinearAlgebra::NumericalVector rightHandSideFace;
            return rightHandSideFace;
        }
        
        /// \brief Compute the right hand side for the solution at time level 'timeLevelIn' and store the result at time level 'timeLevelResult'. Make sure timeLevelIn is different from timeLevelResult.
        virtual void computeRightHandSide(const std::size_t timeLevelIn, const std::size_t timeLevelResult, const double time);

        /// \brief Get a linear combination of solutions at time level 'timeLevelIn' with coefficients given in coefficientsTimeLevels.
        LinearAlgebra::NumericalVector getSolutionCoefficients(const Base::Element *ptrElement, const std::vector<std::size_t> timeLevelsIn, const std::vector<double> coefficientsTimeLevels);

        /// \brief Compute the right hand side for the linear combination of solutions at time level 'timeLevelIn' with coefficients given in coefficientsTimeLevels. Store the result at time level 'timeLevelResult'.
        virtual void computeRightHandSide(const std::vector<std::size_t> timeLevelsIn, const std::vector<double> coefficientsTimeLevels, const std::size_t timeLevelResult, const double time);

        /// \brief Synchronize between the different submeshes (when using MPI)
        void synchronize(const std::size_t timeLevel);

        /// \brief Scale the solution coefficients of a given time level.
        void scaleTimeLevel(const std::size_t timeLevel, const double scale);

        /// \brief scale and add solution coefficients of a certain time level to the coefficients of another time level.
        void scaleAndAddTimeLevel(const std::size_t timeLevelToChange, const std::size_t timeLevelToAdd, const double scale);

        /// \brief Set the initial numerical solution (w at t=0).
        void setInitialSolution(const std::size_t solutionTimeLevel, const double startTime, const std::size_t orderTimeDerivative);

        /// \brief Compute the time derivative for a given time level.
        virtual void computeTimeDerivative(const std::size_t timeLevelIn, const std::size_t timeLevelResult, const double time);

        /// \brief Compute the time derivative for a given linear combination of solutions at different time levels.
        virtual void computeTimeDerivative(const std::vector<std::size_t> timeLevelsIn, const std::vector<double> coefficientsTimeLevels, const std::size_t timeLevelResult, const double time);
        
        /// \brief Compute one time step, using a Runge-Kutta scheme.
        virtual void computeOneTimeStep(double &time, const double dt);
        
        /// \brief Set output names.
        void setOutputNames(std::string outputFileName, std::string internalFileTitle, std::string solutionTitle, std::vector<std::string> variableNames);
        
        /// \brief Write output to a tecplot file.
        virtual void writeToTecplotFile(const ElementT *ptrElement, const PointReferenceT &pRef, std::ostream &out) override;
        
        void VTKWrite(Output::VTKTimeDependentWriter& out, double t, std::size_t timeLevel)
        {
            //you would say this could be done more efficiently, but p.first has different types each time
            for (auto p : VTKDoubleWrite_)
            {
                out.write(p.first, p.second, t, timeLevel);
            }
            for (auto p : VTKVectorWrite_)
            {
                out.write(p.first, p.second, t, timeLevel);
            }
            for (auto p : VTKMatrixWrite_)
            {
                out.write(p.first, p.second, t, timeLevel);
            }
        }
        
        virtual void registerVTKWriteFunctions();
        
        /// \brief Show the progress of the time integration.
        virtual void showProgress(const double time, const std::size_t timeStepID)
        {
            if (timeStepID % 10 == 0)
            {
                logger(VERBOSE, "% time steps computed.", timeStepID);
            }
        }
        
        /// \brief Create and Store things before solving the problem.
        virtual void tasksBeforeSolving()
        {
        }
        
        /// \brief Check things before solving (e.g. check if a mesh is created.)
        virtual bool checkBeforeSolving();
        
        /// \brief Solve the PDE, using a Runge-Kutta scheme.
        virtual bool solve(const double startTime, const double endTime, double dt, const std::size_t numOfOutputFrames, bool doComputeError);
        
    protected:
        /// Butcher tableau for time integration. The integration method is assumed to be explicit.
        const Base::ButcherTableau * const ptrButcherTableau_;
        
        /// Index to indicate where the coefficients for the solution are stored.
        std::size_t solutionTimeLevel_;
        
        /// Indices to indicate where the intermediate results are stored.
        std::vector<std::size_t> intermediateTimeLevels_;
        
        /// Name of the complete output file (including extensions like .dat).
        std::string outputFileName_;
        
        /// Title of the file as used by Tecplot internally.
        std::string internalFileTitle_;
        
        /// Title of the solution
        std::string solutionTitle_;
        
        /// Vector of variable names.
        std::vector<std::string> variableNames_;
        
        /// Integrator for the elements
        Integration::ElementIntegral elementIntegrator_;
        
        /// Integrator for the faces
        Integration::FaceIntegral faceIntegrator_;
        
    private:
        /// \brief Define how the solution should be written in the VTK files.
        /// \details For an example of using this function, see for example the application 'TutorialAdvection' to find out how to use this function.
        void registerVTKWriteFunction(std::function<double(Base::Element*, const Geometry::PointReference&, std::size_t)> function, std::string name)
        {
            VTKDoubleWrite_.push_back( {function, name});
        }
        
        void registerVTKWriteFunction(std::function<LinearAlgebra::NumericalVector(Base::Element*, const Geometry::PointReference&, std::size_t)> function, std::string name)
        {
            VTKVectorWrite_.push_back( {function, name});
        }
        
        void registerVTKWriteFunction(std::function<LinearAlgebra::Matrix(Base::Element*, const Geometry::PointReference&, std::size_t)> function, std::string name)
        {
            VTKMatrixWrite_.push_back( {function, name});
        }
        
        std::vector<std::pair<std::function<double(Base::Element*, const Geometry::PointReference&, std::size_t)>, std::string> > VTKDoubleWrite_;
        std::vector<std::pair<std::function<LinearAlgebra::NumericalVector(Base::Element*, const Geometry::PointReference&, std::size_t)>, std::string> > VTKVectorWrite_;
        std::vector<std::pair<std::function<LinearAlgebra::Matrix(Base::Element*, const Geometry::PointReference&, std::size_t)>, std::string> > VTKMatrixWrite_;
        
    };
}

#endif


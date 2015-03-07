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
#include "Base/HpgemUI.h"
#include "Integration/ElementIntegrandBase.h"
#include "Integration/FaceIntegrandBase.h"
#include "Output/TecplotSingleElementWriter.h"
#include "Output/VTKTimeDependentWriter.h"
#include <functional>

namespace Integration
{
    class FaceIntegral;
}

namespace Base
{
    /// \brief Simplified Interface for solving PDE's.
    /** \details To solve some linear time depent PDE you should do the following:
     * \li Create your own class that inherits this class.
     * \li Implement the function 'initialise' for creating the mesh.
     * \li Implement the function 'initialConditions' to define the initial conditions of your problem.
     * \li Implement the functions 'elementIntegrand' and 'faceIntegrand' for defining the integrands for element matrices and vectors and face matrices and vectors.
     * \li Implement the functions 'computeRhsLocal' and 'computeRhsFaces' for computing the right-hand-side corresponding to your time integration method.
     * \li Implement the function 'beforeTimeIntegration' when multiple element/face matrices/vectors are required.
     * \li Override the function 'solve' when using another time integration routine than forward Euler.
     * \li Implement the function 'writeToTecplotFile' to determine what data to write to the output file.
     */
    /** \details To solve the PDE do the following in the main routine:
     * \li Create an object of your own class.
     * \li Define how the solution should be written in the VTK files using the function 'registerVTKWriteFunction'.
     * \li Call the function 'solve'.
     */
    /** \details For an example of using this interface see the application 'TutorialAdvection'.
     */
    class HpgemAPISimplified : public HpgemUI, Output::TecplotSingleElementWriter
    {
    public:
        // Constructor
        HpgemAPISimplified(const std::size_t dimension, const std::size_t numberOfUnknowns, const std::size_t polynomialOrder, const std::size_t numberOfTimeLevels, const bool useMatrixStorage = false);

        /// \brief Create the mesh.
        virtual void createMesh(std::size_t numOfElementsPerDirection, Base::MeshType meshType = Base::MeshType::RECTANGULAR)
        {
            logger(ERROR, "No routine for creating the mesh implemented.");
        }
        
        /// \brief Compute the source term at a given physical point.
        virtual double getSourceTerm(const double &time, const PointPhysicalT &pPhys)
        {
            logger(ERROR, "No source term implemented.");
            return 0.0;
        }
        
        /// \brief Compute the real solution at a given point in space and time.
        virtual LinearAlgebra::NumericalVector getRealSolution(const double &time, const PointPhysicalT &pPhys)
        {
            logger(ERROR, "No real solution implemented.");
            LinearAlgebra::NumericalVector realSolution(configData_->numberOfUnknowns_);
            return realSolution;
        }
        
        /// \brief Compute the real solution at a given point in space and time.
        virtual LinearAlgebra::NumericalVector getInitialSolution(const double &startTime, const PointPhysicalT &pPhys)
        {
            logger(ERROR, "No initial solution implemented.");
            LinearAlgebra::NumericalVector realSolution(configData_->numberOfUnknowns_);
            return realSolution;
        }
        
        /// \brief Compute the mass matrix for a single element.
        virtual LinearAlgebra::Matrix computeMassMatrixAtElement(const Base::Element *ptrElement, const double time)
        {
            logger(ERROR, "No function for computing the mass matrix at an element implemented.");
            LinearAlgebra::Matrix massMatrix;
            return massMatrix;
        }
        
        /// \brief Compute and store the mass matrices.
        void createMassMatrices(const double time);

        /// \brief Solve the mass matrix equations.
        void solveMassMatrixEquations(const std::size_t timeLevel, const double time);

        /// \brief Integrate the initial solution at a single element.
        virtual LinearAlgebra::NumericalVector integrateInitialSolutionAtElement(const Base::Element * ptrElement, const double startTime, const std::size_t orderTimeDerivative)
        {
            logger(ERROR, "No function for computing the integral for an initial solution at an element implemented.");
            LinearAlgebra::NumericalVector integralInitialSolution;
            return integralInitialSolution;
        }
        
        /// \brief Integrate the initial solution.
        void integrateInitialSolution(const std::size_t timeLevelResult, const double startTime, const std::size_t orderTimeDerivative);

        /// \brief Compute the stiffness matrix corresponding to an element.
        virtual LinearAlgebra::Matrix computeStiffnessMatrixAtElement(const Base::Element *ptrElement, const double time)
        {
            logger(ERROR, "No function for computing the stiffness matrix at an element implemented.");
            LinearAlgebra::Matrix massMatrix;
            return massMatrix;
        }
        
        /// \brief Compute the stiffness matrix corresponding to a face.
        virtual Base::FaceMatrix computeStiffnessMatrixAtFace(const Base::Face *ptrFace, const double time)
        {
            logger(ERROR, "No function for computing the stiffness matrix at a face implemented.");
            Base::FaceMatrix stiffnessMatrix;
            return stiffnessMatrix;
        }
        
        /// \brief Compute and store stiffness matrices for computing the right hand side.
        void createStiffnessMatrices(const double time);

        /// \brief Integrate the energy of the error on a single element.
        /// \details The energy of the error is defined here as \f$ e^tMe\f$, where \f$ M\f$ is the mass matrix and \f$ e\f$ is the error. In most cases this is equivalent to the \f$ L^2 \f$ norm of the error.
        virtual LinearAlgebra::NumericalVector integrateEnergyErrorAtElement(const Base::Element *ptrElement, const double time, LinearAlgebra::NumericalVector &solutionCoefficients)
        {
            logger(ERROR, "No function for computing the error at an element implemented.");
            LinearAlgebra::NumericalVector errorAtElement(1);
            return errorAtElement;
        }
        
        /// \brief Compute the energy-norm of the error.
        double computeEnergyNormError(const std::size_t solutionTimeLevel, const double time);

        /// \brief Compute the right-hand side corresponding to an element
        virtual LinearAlgebra::NumericalVector computeRightHandSideAtElement(const Base::Element *ptrElement, const double time, LinearAlgebra::NumericalVector &solutionCoefficients)
        {
            logger(ERROR, "No function for computing the right-hand side at an element implemented.");
            LinearAlgebra::NumericalVector rightHandSideElement;
            return rightHandSideElement;
        }
        
        /// \brief Compute the right-hand side corresponding to a face
        virtual LinearAlgebra::NumericalVector computeRightHandSideAtFace(const Base::Face *ptrFace, const double time, const Base::Side side, LinearAlgebra::NumericalVector &solutionCoefficientsLeft, LinearAlgebra::NumericalVector &solutionCoefficientsRight)
        {
            logger(ERROR, "No function for computing the right-hand side at a face implemented.");
            LinearAlgebra::NumericalVector rightHandSideFace;
            return rightHandSideFace;
        }
        
        /// \brief Compute the right hand side for the solution at time level 'timeLevelIn' and store the result at time level 'timeLevelResult'. Make sure timeLevelIn is different from timeLevelResult.
        void computeRightHandSide(const std::size_t timeLevelIn, const std::size_t timeLevelResult, const double time);

        /// \brief Get a linear combination of solutions at time level 'timeLevelIn' with coefficients given in coefficientsTimeLevels.
        LinearAlgebra::NumericalVector getSolutionCoefficients(const Base::Element *ptrElement, const std::vector<std::size_t> timeLevelsIn, const std::vector<double> coefficientsTimeLevels);

        /// \brief Compute the right hand side for the linear combination of solutions at time level 'timeLevelIn' with coefficients given in coefficientsTimeLevels. Store the result at time level 'timeLevelResult'.
        void computeRightHandSide(const std::vector<std::size_t> timeLevelsIn, const std::vector<double> coefficientsTimeLevels, const std::size_t timeLevelResult, const double time);

        /// \brief Synchronize between the different submeshes (when using MPI)
        void synchronize(const std::size_t timeLevel);

        /// \brief Scale the solution coefficients of a given time level.
        void scaleTimeLevel(const std::size_t timeLevel, const double scale);

        /// \brief scale and add solution coefficients of a certain time level to the coefficients of another time level.
        void scaleAndAddTimeLevel(const std::size_t timeLevelToChange, const std::size_t timeLevelToAdd, const double scale);

        /// \brief Set the initial numerical solution (w at t=0).
        void setInitialSolution(const std::size_t solutionTimeLevel, const double startTime);

        /// \brief Compute the time derivative for a given time level.
        void computeTimeDerivative(const std::size_t timeLevelIn, const std::size_t timeLevelResult, const double time);

        /// \brief Compute the time derivative for a given linear combination of solutions at different time levels.
        void computeTimeDerivative(const std::vector<std::size_t> timeLevelsIn, const std::vector<double> coefficientsTimeLevels, const std::size_t timeLevelResult, const double time);

        /// \brief Write output to a tecplot file.
        virtual void writeToTecplotFile(const ElementT *ptrElement, const PointReferenceT &pRef, std::ostream &out)
        {
            logger(ERROR, "No function for writing to a tecplot file implemented.");
        }
        
        // \brief Solve the PDE.
        virtual bool solve(const double startTime, const double endTime, double dt, const std::size_t numOfOutputFrames = 1)
        {
            logger(ERROR, "No function for solving the problem implemented.");
            return false;
        }
        
    protected:
        /// Initial time. The PDE is solved over the time interval [startTime,endTime].
        double startTime_;

        /// Final time. The PDE is solved over the time interval [startTime,endTime].
        double endTime_;

        /// Boolean to indicate if matrices (e.g. mass matrix, stiffness matrix etc) should be stored.
        const bool useMatrixStorage_;

        /// Index to indicate where the mass matrix is stored
        const std::size_t massMatrixID_;

        /// Index to indicate where the stiffness matrix for the elements is stored
        const std::size_t stiffnessElementMatrixID_;

        /// Index to indicate where the stiffness matrix for the elements is stored
        const std::size_t stiffnessFaceMatrixID_;
        
    };
}

#endif


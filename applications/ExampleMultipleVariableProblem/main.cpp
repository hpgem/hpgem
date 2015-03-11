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

#include <fstream>

#include "Base/CommandLineOptions.h"
#include "Base/ConfigurationData.h"
#include "Base/Element.h"
#include "Base/Face.h"
#include "Base/HpgemAPISimplified.h"
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

#include "Logger.h"

/// \brief Class to demonstrate how a system of multiple PDE's can be solved, using the interface HpgemAPISimplified.
/** \details
 Currently only periodic boundary conditions are used in this example and we do not use a source term.
 
 This class can solve the scalar wave equation in 2D or 3D written as a first order system of PDE's. The original scalar wave equation is given by \f[ \partial_t^2 u = \nabla \cdot c \nabla u \f], where \f$u \f$ is the scalar variable and \f$c \f$ a material parameter corresponding to the velocity with which waves can propagate.
 
 For a first order scheme we define the scalar function \f$ v := \partial_t u \f$ and the vector function \f$ s := c \nabla u \f$. We can then obtain the equations \f[ \partial_t v = \nabla \cdot s \f] and \f[ c^{-1} \partial_t s = \nabla v \f]
 
 We define a new vector function \f$w = [w_0, w_1, w_2, ..] = [v, s_0, s_1, ..]\f$. We can then rewrite the system of PDE's as follows: \f[ \partial_t w_0 = \partial_i w_{i+1} \f] summing over i = 0 .. (DIM-1) and \f[ c^{-1} \partial_t w_{i+1} = \partial_i w_0 \f] for i = 0 .. (DIM-1).
 
 A boolean useMatrixStorage can be set to true if you want to store all mass and stiffness matrices. This only works for linear problems, but might reduce computational costs.
 
 This class consists of the following parts:
 \li A constructor to set the dimension, number of elements, polynomial order, butcher tableau, and boolean for storing matrices.
 \li The functions 'createMesh' and 'setMaterialParameter' are used to create the mesh and set the material parameters.
 \li The functions 'getCInv', 'getSourceTerm' and 'getRealSolution' return the material parameters, source term and analytic solution.
 \li The functions 'integrand...OnRefElement' and 'integrand...OnRefElement' compute the integrand for ... for the reference element/face. These functions are necessary to compute the mass matrix, stiffness matrix, initial solution and numerical error. When not storing matrices we compute the integrand for the right-hand side.
 \li The function 'solve' solves the PDE over the time interval [startTime, endTime].
 \li The function 'writeToTecplotFile' is used to determine which data should be written to the output file.
 */
/** \details To solve the problem, the following things are done in the main routine:
 \li An object of the class ExampleMultipleVariableProblem is created.
 \li The mesh is created with 'createMesh'.
 \li The material parameters are set using 'setMaterialParameters'.
 \li The names for the output files are set using 'setOutputNames'.
 \li The function 'solve' is then used to solve the PDE.
 */
class ExampleMultipleVariableProblem : public Base::HpgemAPISimplified
{
public:
    /// \param[in] dimension Dimension of the domain
    /// \param[in] numberOfVariables Number of variables in the PDE
    /// \param[in] polynomialOrder Polynomial order of the basis functions
    /// \param[in] useMatrixStorage Boolean to indicate if element and face matrices for the PDE should be stored
    /// \param[in] ptrButcherTableau Pointer to a Butcher Tableau used to do the time integration with a Runge-Kutta scheme. By default this is a RK4 scheme.
    ExampleMultipleVariableProblem
    (
     const std::size_t dimension,
     const std::size_t numOfVariables,
     const std::size_t polynomialOrder,
     const Base::ButcherTableau * const ptrButcherTableau,
     const bool useMatrixStorage
     ) :
    HpgemAPISimplified(dimension, numOfVariables, polynomialOrder, ptrButcherTableau, ptrButcherTableau->getNumStages() + 1, useMatrixStorage, false),
    DIM_(dimension),
    numOfVariables_(numOfVariables),
    elementIntegrator_(),
    faceIntegrator_(),
    cInv_(1.0)
    {
        solutionTimeLevel_ = 0;
        for (std::size_t i = 1; i < configData_->numberOfTimeLevels_; i++)
        {
            intermediateTimeLevels_.push_back(i);
        }
    }
    
    /// \brief Create a domain
    Base::RectangularMeshDescriptor createMeshDescription(const std::size_t numOfElementPerDirection) override
    {
        // Create the domain. In this case the domain is the square [0,1]^DIM and periodic.
        Base::RectangularMeshDescriptor description(DIM_);
        for (int i = 0; i < DIM_; ++i)
        {
            description.bottomLeft_[i] = 0;
            description.topRight_[i] = 1;
            description.numElementsInDIM_[i] = numOfElementPerDirection;
            description.boundaryConditions_[i] = Base::RectangularMeshDescriptor::PERIODIC;
        }
        
        return description;
    }
    
    /// \brief Set the material parameter.
    /// \param[in] c Material parameter corresponding to the speed with which waves can propagate.
    void setMaterialParameter(const double c)
    {
        cInv_ = 1.0 / c;
    }
    
    /// \brief Get the material parameter c^{-1} at a given physical point.
    double getCInv(const Geometry::PointPhysical &pPhys)
    {
        return cInv_;
    }
    
    /// \brief Compute the source term at a given physical point.
    double getSourceTerm(const PointPhysicalT &pPhys, const double &time, const std::size_t orderTimeDerivative = 0) override
    {
        return 0.0;
    }
    
    /// \brief Compute the real solution at a given point in space and time.
    LinearAlgebra::NumericalVector getRealSolution(const PointPhysicalT &pPhys, const double &time, const std::size_t orderTimeDerivative = 0) override
    {
        LinearAlgebra::NumericalVector realSolution(numOfVariables_);
        double c = std::sqrt(1.0 / cInv_); // Wave velocity.
                
        double xWave = 0; // The inner product of the physical point and some direction vector. In this case the direction (1,..,1).
        for (std::size_t iD = 0; iD < DIM_; iD++) // Index for the dimension.
        {
            xWave += pPhys[iD];
        }
        
        for (std::size_t iV = 0; iV < numOfVariables_; iV++) // iV is the index for the variable.
        {
            realSolution(iV) = 2 * M_PI * std::cos(2 * M_PI * (xWave - std::sqrt(DIM_) * c * time));
            if (iV == 0)
            {
                realSolution(iV) *= -std::sqrt(DIM_) * c;
            }
            else
            {
                realSolution(iV) /= cInv_;
            }
        }
        
        return realSolution;
    }
    
    /// \brief Compute the initial solution at a given point in space and time.
    LinearAlgebra::NumericalVector getInitialSolution(const PointPhysicalT &pPhys, const double &startTime, const std::size_t orderTimeDerivative = 0) override
    {
        return getRealSolution(pPhys, startTime, orderTimeDerivative);
    }
    
    /// \brief Compute the integrand for the mass matrix for the reference element.
    /// \details The integrand for the reference element is the same as the physical element, but scaled with the reference-to-physical element scale, which is the determinant of the jacobian of the reference-to-physical element mapping.
    LinearAlgebra::Matrix integrandMassMatrixOnRefElement(const Base::Element *ptrElement, const Geometry::PointReference &pRef)
    {
        std::size_t numOfBasisFunctions = ptrElement->getNrOfBasisFunctions();
        LinearAlgebra::Matrix integrand(numOfVariables_ * numOfBasisFunctions, numOfVariables_ * numOfBasisFunctions);
        Geometry::PointPhysical pPhys = ptrElement->referenceToPhysical(pRef);
        
        std::size_t iVB, jVB; // indices for both variable and basis function.
        for (std::size_t iV = 0; iV < numOfVariables_; iV++)
        {
            for (std::size_t iB = 0; iB < numOfBasisFunctions; iB++)
            {
                for (std::size_t jB = 0; jB < numOfBasisFunctions; jB++)
                {
                    iVB = ptrElement->convertToSingleIndex(iB, iV);
                    jVB = ptrElement->convertToSingleIndex(jB, iV);
                    integrand(iVB, jVB) = ptrElement->basisFunction(iB, pRef) * ptrElement->basisFunction(jB, pRef);
                    if (iV > 0)
                    {
                        integrand(iVB, jVB) *= getCInv(pPhys);
                    }
                }
            }
        }
        
        // Scale with the reference-to-physical element ratio.
        Geometry::Jacobian jac = ptrElement->calcJacobian(pRef);
        integrand *= jac.determinant();
        
        return integrand;
    }
    
    /// \brief Compute the integrand for the reference element for obtaining the initial solution.
    /// \details The integrand for the initial solution is the exact solution at time 0 multiplied by a test function. The integrand is then scaled by the reference-to-physical element scale, since we compute the integral on a reference element.
    LinearAlgebra::NumericalVector integrandInitialSolutionOnRefElement(const Base::Element *ptrElement, const double &startTime, const Geometry::PointReference &pRef)
    {
        std::size_t numOfBasisFunctions = ptrElement->getNrOfBasisFunctions();
        
        LinearAlgebra::NumericalVector integrand(numOfVariables_ * numOfBasisFunctions);
        
        Geometry::PointPhysical pPhys = ptrElement->referenceToPhysical(pRef);
        
        LinearAlgebra::NumericalVector initialSolution(getInitialSolution(pPhys, startTime));
        
        std::size_t iVB; // Index for both variable and basis function.
        for (std::size_t iV = 0; iV < numOfVariables_; iV++)
        {
            for (std::size_t iB = 0; iB < numOfBasisFunctions; iB++)
            {
                iVB = ptrElement->convertToSingleIndex(iB, iV);
                integrand(iVB) = ptrElement->basisFunction(iB, pRef) * initialSolution(iV);
                if (iV > 0)
                {
                    integrand(iVB) *= getCInv(pPhys);
                }
            }
        }
        
        // Scale with the reference-to-physical element ratio.
        Geometry::Jacobian jac = ptrElement->calcJacobian(pRef);
        integrand *= jac.determinant();
        
        return integrand;
    }
    
    /// \brief Compute the integrand for the stiffness matrix for the reference element.
    /// \details The integrand for the reference element is the same as the physical element, but scaled with the reference-to-physical element scale, which is the determinant of the jacobian of the reference-to-physical element mapping.
    LinearAlgebra::Matrix integrandStiffnessMatrixOnRefElement(const Base::Element *ptrElement, const Geometry::PointReference &pRef)
    {
        std::size_t numOfBasisFunctions = ptrElement->getNrOfBasisFunctions();
        LinearAlgebra::Matrix integrand(numOfVariables_ * numOfBasisFunctions, numOfVariables_ * numOfBasisFunctions);
        LinearAlgebra::NumericalVector gradientBasisFunction(DIM_);
        double valueTestFunction;
        
        std::size_t iVB, jVB; // Indices for both variable and basisfunction
        for (std::size_t jB = 0; jB < numOfBasisFunctions; jB++)
        {
            gradientBasisFunction = ptrElement->basisFunctionDeriv(jB, pRef);
            for (std::size_t iB = 0; iB < numOfBasisFunctions; iB++)
            {
                valueTestFunction = ptrElement->basisFunction(iB, pRef);
                
                iVB = ptrElement->convertToSingleIndex(iB, 0);
                for (std::size_t jD = 0; jD < DIM_; jD++) // Index for the derivatives
                {
                    jVB = ptrElement->convertToSingleIndex(jB, jD + 1);
                    integrand(iVB, jVB) = gradientBasisFunction(jD) * valueTestFunction;
                    
                }
                
                jVB = ptrElement->convertToSingleIndex(jB, 0);
                for (std::size_t iD = 0; iD < DIM_; iD++) // Index for the derivatives
                {
                    iVB = ptrElement->convertToSingleIndex(iB, iD + 1);
                    integrand(iVB, jVB) = gradientBasisFunction(iD) * valueTestFunction;
                }
            }
        }
        
        // Scale with the reference-to-physical element ratio.
        Geometry::Jacobian jac = ptrElement->calcJacobian(pRef);
        integrand *= jac.determinant();
        
        return integrand;
    }
    
    /// \brief Compute the integrand for the stiffness matrix for the reference face corresponding to an internal face.
    /// \details The integrand for the reference face is the same as the physical face, but scaled with the reference-to-physical face scale. This face scale is absorbed in the normal vector, since it is relatively cheap to compute the normal vector with a length (L2-norm) equal to the reference-to-physical face scale.
    LinearAlgebra::Matrix integrandStiffnessMatrixOnRefFace(const Base::Face *ptrFace, const Geometry::PointReference &pRef, const Base::Side &iSide, const Base::Side &jSide)
    {
        std::size_t numOfSolutionBasisFunctions = ptrFace->getPtrElement(jSide)->getNrOfBasisFunctions();
        std::size_t numOfTestBasisFunctions = ptrFace->getPtrElement(iSide)->getNrOfBasisFunctions();
        
        LinearAlgebra::Matrix integrand(numOfTestBasisFunctions * numOfVariables_, numOfSolutionBasisFunctions * numOfVariables_);
        
        double valueTestFunction;
        double valueBasisFunction;
        
        // Compute normal vector, with size of the ref-to-phys face scale, pointing outward of the left element.
        LinearAlgebra::NumericalVector normal = ptrFace->getNormalVector(pRef);
        if (jSide == Base::Side::RIGHT)
        {
            normal *= -1;
        };
        
        // Compute the integrand
        std::size_t iVB, jVB; // Indices for both variable and basisfunction.
        for (std::size_t jB = 0; jB < numOfSolutionBasisFunctions; jB++) // iB and jB are indices for the basis function.
        {
            valueBasisFunction = ptrFace->basisFunction(jSide, jB, pRef);
            for (std::size_t iB = 0; iB < numOfTestBasisFunctions; iB++)
            {
                valueTestFunction = ptrFace->basisFunction(iSide, iB, pRef);
                
                iVB = ptrFace->getPtrElement(iSide)->convertToSingleIndex(iB, 0);
                for (std::size_t jD = 0; jD < DIM_; jD++) // index for the direction
                {
                    jVB = ptrFace->getPtrElement(jSide)->convertToSingleIndex(jB, jD + 1);
                    integrand(iVB, jVB) = -0.5 * normal(jD) * valueBasisFunction * valueTestFunction;
                    
                }
                
                jVB = ptrFace->getPtrElement(jSide)->convertToSingleIndex(jB, 0);
                for (std::size_t iD = 0; iD < DIM_; iD++) // index for the direction
                {
                    iVB = ptrFace->getPtrElement(iSide)->convertToSingleIndex(iB, iD + 1);
                    integrand(iVB, jVB) = -0.5 * normal(iD) * valueBasisFunction * valueTestFunction;
                }
            }
        }
        
        return integrand;
    }
    
    /// \brief Compute the integrand for the right hand side for the reference element.
    /// \details The integrand for the reference element is the same as the physical element, but scaled with the reference-to-physical element scale, which is the determinant of the jacobian of the reference-to-physical element mapping.
    LinearAlgebra::NumericalVector integrandRightHandSideOnRefElement(const Base::Element *ptrElement, const double &time, const Geometry::PointReference &pRef, const LinearAlgebra::NumericalVector &solutionCoefficients)
    {
        std::size_t numOfBasisFunctions = ptrElement->getNrOfBasisFunctions();
        
        LinearAlgebra::NumericalVector integrand(numOfVariables_ * numOfBasisFunctions);
        
        LinearAlgebra::NumericalVector gradientBasisFunction(DIM_);
        LinearAlgebra::NumericalVector gradientScalarFunction(DIM_);
        double divergenceVectorFunction = 0;
        
        // Compute the gradient of the scalar function and the divergence of the vector function.
        std::size_t jVB; // Index for both basis function and variable
        for (std::size_t jD = 0; jD < DIM_; jD++) // Index for derivatives
        {
            gradientScalarFunction(jD) = 0;
        }
        for (std::size_t jB = 0; jB < numOfBasisFunctions; jB++) // Index for basis functions
        {
            gradientBasisFunction = ptrElement->basisFunctionDeriv(jB, pRef);
            
            for (std::size_t jD = 0; jD < DIM_; jD++) // Index for derivatives
            {
                jVB = ptrElement->convertToSingleIndex(jB, 0);
                gradientScalarFunction(jD) += gradientBasisFunction(jD) * solutionCoefficients(jVB);
                
                jVB = ptrElement->convertToSingleIndex(jB, jD + 1);
                divergenceVectorFunction += gradientBasisFunction(jD) * solutionCoefficients(jVB);
            }
        }
        
        // Compute integrand on the physical element.
        std::size_t iVB; // Index for both basis function and variable
        for (std::size_t iB = 0; iB < numOfBasisFunctions; iB++) // Index for basis function
        {
            iVB = ptrElement->convertToSingleIndex(iB, 0);
            integrand(iVB) = ptrElement->basisFunction(iB, pRef) * divergenceVectorFunction;
            
            for (std::size_t iD = 0; iD < DIM_; iD++) // Index for derivative
            {
                iVB = ptrElement->convertToSingleIndex(iB, iD + 1);
                integrand(iVB) = ptrElement->basisFunction(iB, pRef) * gradientScalarFunction(iD);
            }
        }
        
        // Scale with the reference-to-physical element ratio.
        Geometry::Jacobian jac = ptrElement->calcJacobian(pRef);
        integrand *= jac.determinant();
        
        return integrand;
    }
    
    /// \brief Compute the integrand for the right hand side for the reference face corresponding to an internal face.
    /// \details The integrand for the reference face is the same as the physical face, but scaled with the reference-to-physical face scale. This face scale is absorbed in the normal vector, since it is relatively cheap to compute the normal vector with a length (L2-norm) equal to the reference-to-physical face scale.
    LinearAlgebra::NumericalVector integrandRightHandSideOnRefFace(const Base::Face *ptrFace, const double &time, const Geometry::PointReference &pRef, const Base::Side &iSide, const LinearAlgebra::NumericalVector &solutionCoefficientsLeft, const LinearAlgebra::NumericalVector &solutionCoefficientsRight)
    {
        std::size_t numOfTestBasisFunctions = ptrFace->getPtrElement(iSide)->getNrOfBasisFunctions();
        std::size_t numOfSolutionBasisFunctionsLeft = ptrFace->getPtrElementLeft()->getNrOfBasisFunctions();
        std::size_t numOfSolutionBasisFunctionsRight = ptrFace->getPtrElementRight()->getNrOfBasisFunctions();
        
        LinearAlgebra::NumericalVector integrand(numOfVariables_ * numOfTestBasisFunctions);
        
        LinearAlgebra::NumericalVector numericalSolutionLeft(numOfVariables_);
        LinearAlgebra::NumericalVector numericalSolutionRight(numOfVariables_);
        
        // Compute the numerical solution at the given point at the left and right side.
        std::size_t jVB; // Index for both variable and basis function.
        for (std::size_t jV = 0; jV < numOfVariables_; jV++)
        {
            numericalSolutionLeft(jV) = 0;
            numericalSolutionRight(jV) = 0;
            for (std::size_t jB = 0; jB < numOfSolutionBasisFunctionsLeft; jB++)
            {
                jVB = ptrFace->getPtrElementLeft()->convertToSingleIndex(jB, jV);
                numericalSolutionLeft(jV) += ptrFace->basisFunction(Base::Side::LEFT, jB, pRef) * solutionCoefficientsLeft(jVB);
            }
            for (std::size_t jB = 0; jB < numOfSolutionBasisFunctionsRight; jB++)
            {
                jVB = ptrFace->getPtrElementRight()->convertToSingleIndex(jB, jV);
                numericalSolutionRight(jV) += ptrFace->basisFunction(Base::Side::RIGHT, jB, pRef) * solutionCoefficientsRight(jVB);
            }
        }
        
        // Compute normal vector, with size of the ref-to-phys face scale, pointing outward of the left element.
        LinearAlgebra::NumericalVector normal = ptrFace->getNormalVector(pRef);
        
        // Compute the jump of the scalar function and the vector function.
        LinearAlgebra::NumericalVector jumpScalarFunction(DIM_);
        double jumpVectorFunction = 0;
        for (std::size_t jD = 0; jD < DIM_; jD++)
        {
            jumpScalarFunction(jD) = normal(jD) * (numericalSolutionLeft(0) - numericalSolutionRight(0));
            jumpVectorFunction += normal(jD) * (numericalSolutionLeft(jD + 1) - numericalSolutionRight(jD + 1));
        }
        
        // Compute integrand on the reference element.
        std::size_t iVB; // Index for both variable and basis function.
        for (std::size_t iB = 0; iB < numOfTestBasisFunctions; iB++)
        {
            iVB = ptrFace->getPtrElement(iSide)->convertToSingleIndex(iB, 0);
            integrand(iVB) = -0.5 * ptrFace->basisFunction(iSide, iB, pRef) * jumpVectorFunction;
            
            for (std::size_t iD = 0; iD < DIM_; iD++) // Index for direction
            {
                iVB = ptrFace->getPtrElement(iSide)->convertToSingleIndex(iB, iD + 1);
                integrand(iVB) = -0.5 * ptrFace->basisFunction(iSide, iB, pRef) * jumpScalarFunction(iD);
            }
        }
        return integrand;
    }
    
    /// \brief Compute the integrand for the reference element for computing the energy-norm of the error.
    /// \details The integrand for the reference element is the same as the physical element, but scaled with the reference-to-physical element scale, which is the determinant of the jacobian of the reference-to-physical element mapping.
    LinearAlgebra::NumericalVector integrandErrorOnRefElement
    (
     const Base::Element *ptrElement,
     const double &time,
     const Geometry::PointReference &pRef,
     const LinearAlgebra::NumericalVector &solutionCoefficients
     )
    {
        std::size_t numOfBasisFunctions = ptrElement->getNrOfBasisFunctions();
        
        LinearAlgebra::NumericalVector integrand(1);
        integrand(0) = 0;
        
        Geometry::PointPhysical pPhys = ptrElement->referenceToPhysical(pRef);
        
        LinearAlgebra::NumericalVector realSolution(getRealSolution(pPhys, time));
        LinearAlgebra::NumericalVector numericalSolution(numOfVariables_);
        
        std::size_t jVB; // Index for both variable and basis function.
        for (std::size_t jV = 0; jV < numOfVariables_; jV++)
        {
            numericalSolution(jV) = 0;
            for (std::size_t jB = 0; jB < numOfBasisFunctions; jB++)
            {
                jVB = ptrElement->convertToSingleIndex(jB, jV);
                numericalSolution(jV) += ptrElement->basisFunction(jB, pRef) * solutionCoefficients(jVB);
            }
            if (jV > 0)
            {
                integrand(0) += std::pow(numericalSolution(jV) - realSolution(jV), 2);
            }
            else
            {
                integrand(0) += getCInv(pPhys) * std::pow(numericalSolution(jV) - realSolution(jV), 2);
            }
        }
        
        // Scale with the reference-to-physical element ratio.
        Geometry::Jacobian jac = ptrElement->calcJacobian(pRef);
        integrand *= jac.determinant();
        
        return integrand;
    }
    
    /// \brief Compute the mass matrix for a single element.
    LinearAlgebra::Matrix computeMassMatrixAtElement(const Base::Element *ptrElement) override
    {
        std::function<LinearAlgebra::Matrix(const Geometry::PointReference &)> integrandFunction = [=](const Geometry::PointReference & pRef) -> LinearAlgebra::Matrix{ return this -> integrandMassMatrixOnRefElement(ptrElement, pRef);};
        
        return elementIntegrator_.referenceElementIntegral(ptrElement->getGaussQuadratureRule(), integrandFunction);
    }
    
    /// \brief Solve the mass matrix equations for a single element.
    /// \details Solve the equation \f$ Mu = r \f$ for \f$ u \f$ for a single element, where \f$ r \f$ is the right-hand sid and \f$ M \f$ is the mass matrix. The input is the right hand side here called 'solutionCoefficients' and the result is returned in this same vector.
    void solveMassMatrixEquationsAtElement(const Base::Element *ptrElement, LinearAlgebra::NumericalVector &solutionCoefficients) override
    {
        LinearAlgebra::Matrix massMatrix(computeMassMatrixAtElement(ptrElement));
        massMatrix.solve(solutionCoefficients);
    }
    
    /// \brief Integrate the initial solution for a single element.
    LinearAlgebra::NumericalVector integrateInitialSolutionAtElement(const Base::Element * ptrElement, const double startTime, const std::size_t orderTimeDerivative) override
    {
        // Define the integrand function for the the initial solution integral.
        std::function<LinearAlgebra::NumericalVector(const Geometry::PointReference &)> integrandFunction = [=](const Geometry::PointReference & pRef) -> LinearAlgebra::NumericalVector { return this -> integrandInitialSolutionOnRefElement(ptrElement, startTime, pRef);};
        
        return elementIntegrator_.referenceElementIntegral(ptrElement->getGaussQuadratureRule(), integrandFunction);
    }
    
    /// \brief Compute the stiffness matrix corresponding to an element.
    LinearAlgebra::Matrix computeStiffnessMatrixAtElement(const Base::Element *ptrElement) override
    {
        // Define the integrand function for the stiffness matrix for the element.
        std::function<LinearAlgebra::Matrix(const Geometry::PointReference &)> integrandFunction = [=](const Geometry::PointReference & pRef) -> LinearAlgebra::Matrix
        {   return this->integrandStiffnessMatrixOnRefElement(ptrElement, pRef);};
        
        return elementIntegrator_.referenceElementIntegral(ptrElement->getGaussQuadratureRule(), integrandFunction);
    }
    
    /// \brief Compute the stiffness matrix corresponding to a face.
    Base::FaceMatrix computeStiffnessMatrixAtFace(const Base::Face *ptrFace) override
    {
        std::size_t numOfBasisFunctionsLeft = 0;
        std::size_t numOfBasisFunctionsRight = 0;
        
        numOfBasisFunctionsLeft = ptrFace->getPtrElementLeft()->getNrOfBasisFunctions();
        if (ptrFace->isInternal())
        {
            numOfBasisFunctionsRight = ptrFace->getPtrElementRight()->getNrOfBasisFunctions();
        }
        
        std::vector<Base::Side> allSides; // Vector with all sides of the face.
        allSides.push_back(Base::Side::LEFT);
        if (ptrFace->isInternal())
        {
            allSides.push_back(Base::Side::RIGHT);
        }
        
        Base::FaceMatrix stiffnessFaceMatrix(numOfBasisFunctionsLeft * numOfVariables_, numOfBasisFunctionsRight * numOfVariables_);
        
        for (Base::Side iSide : allSides)
        {
            for (Base::Side jSide : allSides)
            {
                // Define the integrand function for the stiffness matrix for the face.
                std::function<LinearAlgebra::Matrix(const Geometry::PointReference &)> integrandFunction = [=](const Geometry::PointReference & pRef) -> LinearAlgebra::Matrix
                {   return this->integrandStiffnessMatrixOnRefFace(ptrFace, pRef, iSide, jSide);};
                
                LinearAlgebra::Matrix stiffnessMatrix(faceIntegrator_.referenceFaceIntegral(ptrFace->getGaussQuadratureRule(), integrandFunction));
                
                stiffnessFaceMatrix.setElementMatrix(stiffnessMatrix, iSide, jSide);
            }
        }
        return stiffnessFaceMatrix;
    }
    
    /// \brief Integrate the energy of the error on a single element.
    /// \details The error is defined as error = realSolution - numericalSolution. The energy of the vector (u, s0, s1) is defined as u^2 + c^{-1} * |s|^2.
    LinearAlgebra::NumericalVector integrateErrorAtElement(const Base::Element *ptrElement, LinearAlgebra::NumericalVector &solutionCoefficients, double time) override
    {
        // Define the integrand function for the error energy.
        std::function<LinearAlgebra::NumericalVector(const Geometry::PointReference &)> integrandFunction = [=](const Geometry::PointReference & pRef) -> LinearAlgebra::NumericalVector
        {
            return this->integrandErrorOnRefElement(ptrElement, time, pRef, solutionCoefficients);
        };
        
        return elementIntegrator_.referenceElementIntegral(ptrElement->getGaussQuadratureRule(), integrandFunction);
    }
    
    /// \brief Compute the right-hand side corresponding to an element
    LinearAlgebra::NumericalVector computeRightHandSideAtElement(const Base::Element *ptrElement, LinearAlgebra::NumericalVector &solutionCoefficients, const double time, const std::size_t orderTimeDerivative = 0) override
    {
        // Define the integrand function for the right hand side for the reference element.
        std::function<LinearAlgebra::NumericalVector(const Geometry::PointReference &)> integrandFunction = [=](const Geometry::PointReference & pRef) -> LinearAlgebra::NumericalVector
        {   return this->integrandRightHandSideOnRefElement(ptrElement, time, pRef, solutionCoefficients);};
        
        return elementIntegrator_.referenceElementIntegral(ptrElement->getGaussQuadratureRule(), integrandFunction);
    }
    
    /// \brief Compute the right-hand side corresponding to a face
    LinearAlgebra::NumericalVector computeRightHandSideAtFace
    (
     const Base::Face *ptrFace,
     const Base::Side side,
     LinearAlgebra::NumericalVector &solutionCoefficientsLeft,
     LinearAlgebra::NumericalVector &solutionCoefficientsRight,
     const double time,
     const std::size_t orderTimeDerivative = 0
     ) override
    {
        // Define the integrand function for the right hand side for the reference face.
        std::function<LinearAlgebra::NumericalVector(const Geometry::PointReference &)> integrandFunction = [=](const Geometry::PointReference &pRef) -> LinearAlgebra::NumericalVector
        {   return this->integrandRightHandSideOnRefFace(ptrFace, time, pRef, side, solutionCoefficientsLeft, solutionCoefficientsRight);};
        
        return faceIntegrator_.referenceFaceIntegral(ptrFace->getGaussQuadratureRule(), integrandFunction);
    }
    
private:
    /// Dimension of the domain
    const std::size_t DIM_;
    
    /// Number of variables
    const std::size_t numOfVariables_;
    
    /// Integrator for the elements
    Integration::ElementIntegral elementIntegrator_;

    /// Integrator for the faces
    Integration::FaceIntegral faceIntegrator_;

    /// Material parameter c^{-1}
    double cInv_;
};

auto& numOfElements = Base::register_argument<std::size_t>('n', "numElems", "number of elements per dimension", true);
auto& polynomialOrder = Base::register_argument<std::size_t>('p', "order", "polynomial order of the solution", true);

// Parse options seem to do weird things. Flag 'h' also communicates with Petsc. Several other flags are also already defined in the library, because they are defined in class HpgemUISimplified.
auto& numOfOutputFrames = Base::register_argument<std::size_t>('O', "numOfOutputFrames", "Number of frames to output", false, 1);
auto& startTime = Base::register_argument<double>('S', "startTime", "start time of the simulation", false, 0.0);
auto& endTime = Base::register_argument<double>('T', "endTime", "end time of the simulation", false, 1.0);
auto& dt = Base::register_argument<double>('d', "timeStepSize", "time step of the simulation", false, 0.01);

int main(int argc, char **argv)
{
    Base::parse_options(argc, argv);
    try
    {
        // Set parameters for the PDE.
        const std::size_t dimension = 2;
        const Base::MeshType meshType = Base::MeshType::TRIANGULAR;
        const Base::ButcherTableau * const ptrButcherTableau = Base::AllTimeIntegrators::Instance().getRule(4, 4);
        const bool useMatrixStorage = true;
        const double c = 1.0;
        
        std::string variableString;
        if (dimension == 2)
        {
            variableString = "v,s0,s1";
        }
        else if (dimension == 3)
        {
            variableString = "v,s0,s1,s2";
        }
        
        // Compute parameters for PDE
        const std::size_t numOfVariables = dimension + 1;
        
        
        // Create problem solver 'test', that can solve the acoustic wave equations.
        ExampleMultipleVariableProblem test(dimension, numOfVariables, polynomialOrder.getValue(), ptrButcherTableau, useMatrixStorage);
        
        // Create the mesh
        test.createMesh(numOfElements.getValue(), meshType);
        
        // Set the material parameter
        test.setMaterialParameter(c);
        
        // Set the names for the output file
        test.setOutputNames("output","acousticWave","acousticWave",variableString);
        
        // Solve the problem over time interval [startTime,endTime].
        test.solve(startTime.getValue(), endTime.getValue(), dt.getValue(), numOfOutputFrames.getValue(), true);
        return 0;
    }
    catch (const char* e)
    {
        std::cout << e;
    }
}

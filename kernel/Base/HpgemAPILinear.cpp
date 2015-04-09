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

#include "HpgemAPILinear.h"

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

#include "Logger.h"

namespace Base
{
    /// \param[in] dimension Dimension of the domain
    /// \param[in] numOfVariables Number of variables in the PDE
    /// \param[in] polynomialOrder Polynomial order of the basis functions
    /// \param[in] butcherTableau A butcherTableau used to solve the PDE with a Runge-Kutta method.
    /// \param[in] numOfTimeLevels Number of time levels. If a butcherTableau is set and the number of time levels is too low, this will be corrected automatically.
    /// \param[in] useSourceTerm Boolean to indicate if there is a source term.
    /// \param[in] useSourceTerm Boolean to indicate if there is a source term at the domain boundary.
    HpgemAPILinear::HpgemAPILinear
    (
     const std::size_t dimension,
     const std::size_t numOfVariables,
     const std::size_t polynomialOrder,
     const Base::ButcherTableau * const ptrButcherTableau,
     const std::size_t numOfTimeLevels,
     const bool useSourceTerm,
     const bool useSourceTermAtBoundary
     ) :
    HpgemAPISimplified(dimension, numOfVariables, polynomialOrder, ptrButcherTableau, numOfTimeLevels),
    useSourceTerm_(useSourceTerm),
    useSourceTermAtBoundary_(useSourceTermAtBoundary),
    massMatrixID_(1),
    stiffnessElementMatrixID_(0),
    stiffnessFaceMatrixID_(0)
    {
    }
    
    void HpgemAPILinear::createMesh(const std::size_t numOfElementsPerDirection, const Base::MeshType meshType)
    {
        const Base::RectangularMeshDescriptor description = createMeshDescription(numOfElementsPerDirection);
        
        // Set the number of Element/Face Matrices/Vectors.
        std::size_t numOfElementMatrices = 2;   // Mass matrix and stiffness matrix
        std::size_t numOfElementVectors = 0;
        std::size_t numOfFaceMatrices = 1;      // Stiffness matrix
        std::size_t numOfFaceVectors = 0;

        // Create mesh and set basis functions.
        addMesh(description, meshType, numOfElementMatrices, numOfElementVectors, numOfFaceMatrices, numOfFaceVectors);
        meshes_[0]->useDefaultDGBasisFunctions();
        
        std::size_t nElements = meshes_[0]->getNumberOfElements();
        logger(VERBOSE, "Total number of elements: %", nElements);
    }
    
    void HpgemAPILinear::createMassMatrices()
    {
        for (Base::Element *ptrElement : meshes_[0]->getElementsList())
        {
            LinearAlgebra::Matrix massMatrix(computeMassMatrixAtElement(ptrElement));
            ptrElement->setElementMatrix(massMatrix, massMatrixID_);
        }
    }
    
    /// \details Solve the equation \f$ Mu = r \f$ for \f$ u \f$, where \f$ r \f$ is the right-hand sid and \f$ M \f$ is the mass matrix.
    void HpgemAPILinear::solveMassMatrixEquations(const std::size_t timeLevel)
    {
        for (Base::Element *ptrElement : meshes_[0]->getElementsList())
        {
            LinearAlgebra::NumericalVector &solutionCoefficients(ptrElement->getTimeLevelDataVector(timeLevel));
            
            const LinearAlgebra::Matrix &massMatrix(ptrElement->getElementMatrix(massMatrixID_));
            massMatrix.solve(solutionCoefficients);
        }
        
        synchronize(timeLevel);
    }
    
    LinearAlgebra::Matrix HpgemAPILinear::computeStiffnessMatrixAtElement(Base::Element *ptrElement)
    {
        // Define a function for the integrand of the stiffness matrix at the element.
        std::function<LinearAlgebra::Matrix(const Base::Element *, const Geometry::PointReference &)> integrandFunction = [=](const Base::Element *ptrElement, const Geometry::PointReference &pRef) -> LinearAlgebra::Matrix
        {return this->computeIntegrandStiffnessMatrixAtElement(ptrElement, pRef);};
        
        return elementIntegrator_.integrate(ptrElement, integrandFunction);
    }
    
    Base::FaceMatrix HpgemAPILinear::computeStiffnessMatrixAtFace(Base::Face *ptrFace)
    {
        // Define a function for the integrand of the stiffness matrix at a face.
        std::function<Base::FaceMatrix(const Base::Face *, const LinearAlgebra::NumericalVector &, const Geometry::PointReference &)> integrandFunction = [=](const Base::Face *ptrFace, const LinearAlgebra::NumericalVector &normal, const Geometry::PointReference &pRef) -> Base::FaceMatrix
        {return this->computeIntegrandStiffnessMatrixAtFace(ptrFace, normal, pRef);};
        
        return faceIntegrator_.integrate(ptrFace, integrandFunction);
    }
    
    LinearAlgebra::NumericalVector HpgemAPILinear::integrateSourceTermAtFace(Base::Face *ptrFace, const double time, const std::size_t orderTimeDerivative)
    {
        // Define a function for the integrand of the stiffness matrix at a element.
        std::function<LinearAlgebra::NumericalVector(const Base::Face *, const LinearAlgebra::NumericalVector &, const Geometry::PointReference &)> integrandFunction = [=](const Base::Face *ptrFace, const LinearAlgebra::NumericalVector &normal, const Geometry::PointReference &pRef) -> LinearAlgebra::NumericalVector
        {return this->computeIntegrandSourceTermAtFace(ptrFace, normal, pRef, time, orderTimeDerivative);};
        
        return faceIntegrator_.integrate(ptrFace, integrandFunction);
    }
    
    /// \details By default, the standard L2 inner product with the source term is computed.
    LinearAlgebra::NumericalVector HpgemAPILinear::integrateSourceTermAtElement(Base::Element * ptrElement, const double time, const std::size_t orderTimeDerivative)
    {
        // Get number of basis functions
        std::size_t numOfBasisFunctions = ptrElement->getNrOfBasisFunctions();
        
        // Declare integral source term
        LinearAlgebra::NumericalVector integralSourceTerm(numOfBasisFunctions * configData_->numberOfUnknowns_);
        
        // Declare integrand
        LinearAlgebra::NumericalVector integrandSourceTerm(numOfBasisFunctions * configData_->numberOfUnknowns_);
        
        // Get quadrature rule and number of points.
        const QuadratureRules::GaussQuadratureRule *ptrQdrRule = ptrElement->getGaussQuadratureRule();
        std::size_t numOfQuadPoints = ptrQdrRule->nrOfPoints();
        
        // For each quadrature point, compute the value of the product of the
        // test function and the source term, then add it with the correct weight to the integral solution.
        for (std::size_t pQuad = 0; pQuad < numOfQuadPoints; ++pQuad)
        {
            Geometry::PointReference pRef = ptrQdrRule->getPoint(pQuad);
            Geometry::PointPhysical pPhys = ptrElement->referenceToPhysical(pRef);
            
            Geometry::Jacobian jac = ptrElement->calcJacobian(pRef);
            
            LinearAlgebra::NumericalVector sourceTerm = getSourceTerm(pPhys, time, orderTimeDerivative);
            
            for(std::size_t iB = 0; iB < numOfBasisFunctions; ++iB)
            {
                double valueBasisFunction = ptrElement->basisFunction(iB, pRef);
                
                for(std::size_t iV = 0; iV < configData_->numberOfUnknowns_; ++iV)
                {
                    std::size_t iVB = ptrElement->convertToSingleIndex(iB,iV);
                    
                    integrandSourceTerm(iVB) = sourceTerm(iV) * valueBasisFunction;
                }
            }
            integralSourceTerm.axpy((ptrQdrRule->weight(pQuad)) * std::abs(jac.determinant()), integrandSourceTerm);
        }
        
        return integralSourceTerm;
    }
    
    void HpgemAPILinear::createStiffnessMatrices()
    {
        logger(INFO, "- Creating stiffness matrices for the elements.");
        for (Base::Element *ptrElement : meshes_[0]->getElementsList())
        {
            LinearAlgebra::Matrix stiffnessMatrix(computeStiffnessMatrixAtElement(ptrElement));
            // std::cout << "-- Stiffness matrix element:\n" << stiffnessMatrix << "\n";
            ptrElement->setElementMatrix(stiffnessMatrix, stiffnessElementMatrixID_);
        }
        
        logger(INFO, "- Creating stiffness matrices for the faces.");
        for (Base::Face *ptrFace : meshes_[0]->getFacesList())
        {
            Base::FaceMatrix stiffnessFaceMatrix(computeStiffnessMatrixAtFace(ptrFace));
            if(!ptrFace->isInternal())
            {
                logger.assert(stiffnessFaceMatrix.getNrOfDegreesOfFreedom(Base::Side::RIGHT) == 0,"The number of degrees of freedom corresonding to the right side of a boundary face should be 0, but is here %.", stiffnessFaceMatrix.getNrOfDegreesOfFreedom(Base::Side::RIGHT));
            }
            // std::cout << "-- Stiffness matrix face: \n";
            // std::cout << "--- Stiffness submatrix face:\n" << stiffnessFaceMatrix.getElementMatrix(Base::Side::LEFT, Base::Side::LEFT) << "\n";
            // std::cout << "--- Stiffness submatrix face:\n" << stiffnessFaceMatrix.getElementMatrix(Base::Side::LEFT, Base::Side::RIGHT) << "\n";
            // std::cout << "--- Stiffness submatrix face:\n" << stiffnessFaceMatrix.getElementMatrix(Base::Side::RIGHT, Base::Side::LEFT) << "\n";
            // std::cout << "--- Stiffness submatrix face:\n" << stiffnessFaceMatrix.getElementMatrix(Base::Side::RIGHT, Base::Side::RIGHT) << "\n";
            
            ptrFace->setFaceMatrix(stiffnessFaceMatrix, stiffnessFaceMatrixID_);
        }
    }
    
    /// \details Make sure timeLevelIn is different from timeLevelResult.
    void HpgemAPILinear::multiplyStiffnessMatrices(const std::size_t timeLevelIn, const std::size_t timeLevelResult)
    {
        // Multiply the stiffness matrices corresponding to the elements.
        for (Base::Element *ptrElement : meshes_[0]->getElementsList())
        {
            LinearAlgebra::NumericalVector &solutionCoefficients(ptrElement->getTimeLevelDataVector(timeLevelIn));
            LinearAlgebra::NumericalVector &solutionCoefficientsNew(ptrElement->getTimeLevelDataVector(timeLevelResult));
            
            solutionCoefficientsNew = ptrElement->getElementMatrix(stiffnessElementMatrixID_) * solutionCoefficients;
        }
        
        // Multiply the stiffness matrices corresponding to the faces.
        for (Base::Face *ptrFace : meshes_[0]->getFacesList())
        {
            if(ptrFace->isInternal())
            {
                LinearAlgebra::NumericalVector &solutionCoefficientsLeft(ptrFace->getPtrElementLeft()->getTimeLevelDataVector(timeLevelIn));
                LinearAlgebra::NumericalVector &solutionCoefficientsRight(ptrFace->getPtrElementRight()->getTimeLevelDataVector(timeLevelIn));
                LinearAlgebra::NumericalVector &solutionCoefficientsLeftNew(ptrFace->getPtrElementLeft()->getTimeLevelDataVector(timeLevelResult));
                LinearAlgebra::NumericalVector &solutionCoefficientsRightNew(ptrFace->getPtrElementRight()->getTimeLevelDataVector(timeLevelResult));
                
                const Base::FaceMatrix &stiffnessFaceMatrix = ptrFace->getFaceMatrix(stiffnessFaceMatrixID_);
                
                solutionCoefficientsLeftNew += stiffnessFaceMatrix.getElementMatrix(Base::Side::LEFT, Base::Side::LEFT) * solutionCoefficientsLeft + stiffnessFaceMatrix.getElementMatrix(Base::Side::LEFT, Base::Side::RIGHT) * solutionCoefficientsRight;
                solutionCoefficientsRightNew += stiffnessFaceMatrix.getElementMatrix(Base::Side::RIGHT, Base::Side::LEFT) * solutionCoefficientsLeft + stiffnessFaceMatrix.getElementMatrix(Base::Side::RIGHT, Base::Side::RIGHT) * solutionCoefficientsRight;
            }
            else
            {
                LinearAlgebra::NumericalVector &solutionCoefficients(ptrFace->getPtrElementLeft()->getTimeLevelDataVector(timeLevelIn));
                LinearAlgebra::NumericalVector &solutionCoefficientsNew(ptrFace->getPtrElementLeft()->getTimeLevelDataVector(timeLevelResult));
                
                const LinearAlgebra::Matrix &stiffnessMatrix = ptrFace->getFaceMatrix(stiffnessFaceMatrixID_).getElementMatrix(Base::Side::LEFT, Base::Side::LEFT);
                
                solutionCoefficientsNew += stiffnessMatrix * solutionCoefficients;
            }
        }
        
        synchronize(timeLevelResult);
    }
    
    /// \details Make sure timeLevelsIn are different from timeLevelResult.
    void HpgemAPILinear::multiplyStiffnessMatrices(const std::vector<std::size_t> timeLevelsIn, const std::vector<double> coefficientsTimeLevels, const std::size_t timeLevelResult)
    {
        // Multiply the stiffness matrices corresponding to the elements.
        for (Base::Element *ptrElement : meshes_[0]->getElementsList())
        {
            LinearAlgebra::NumericalVector solutionCoefficients(getSolutionCoefficients(ptrElement, timeLevelsIn, coefficientsTimeLevels));
            LinearAlgebra::NumericalVector &solutionCoefficientsNew(ptrElement->getTimeLevelDataVector(timeLevelResult));
            
            solutionCoefficientsNew = ptrElement->getElementMatrix(stiffnessElementMatrixID_) * solutionCoefficients;
        }
        
        // Multiply the stiffness matrices corresponding to the faces.
        for (Base::Face *ptrFace : meshes_[0]->getFacesList())
        {
            if(ptrFace->isInternal())
            {
                LinearAlgebra::NumericalVector solutionCoefficientsLeft(getSolutionCoefficients(ptrFace->getPtrElementLeft(), timeLevelsIn, coefficientsTimeLevels));
                LinearAlgebra::NumericalVector solutionCoefficientsRight(getSolutionCoefficients(ptrFace->getPtrElementRight(), timeLevelsIn, coefficientsTimeLevels));
                LinearAlgebra::NumericalVector &solutionCoefficientsLeftNew(ptrFace->getPtrElementLeft()->getTimeLevelDataVector(timeLevelResult));
                LinearAlgebra::NumericalVector &solutionCoefficientsRightNew(ptrFace->getPtrElementRight()->getTimeLevelDataVector(timeLevelResult));
                
                const Base::FaceMatrix &stiffnessFaceMatrix = ptrFace->getFaceMatrix(stiffnessFaceMatrixID_);
                
                solutionCoefficientsLeftNew += stiffnessFaceMatrix.getElementMatrix(Base::Side::LEFT, Base::Side::LEFT) * solutionCoefficientsLeft + stiffnessFaceMatrix.getElementMatrix(Base::Side::LEFT, Base::Side::RIGHT) * solutionCoefficientsRight;
                solutionCoefficientsRightNew += stiffnessFaceMatrix.getElementMatrix(Base::Side::RIGHT, Base::Side::LEFT) * solutionCoefficientsLeft + stiffnessFaceMatrix.getElementMatrix(Base::Side::RIGHT, Base::Side::RIGHT) * solutionCoefficientsRight;
            }
            else
            {
                LinearAlgebra::NumericalVector solutionCoefficients(getSolutionCoefficients(ptrFace->getPtrElementLeft(), timeLevelsIn, coefficientsTimeLevels));
                LinearAlgebra::NumericalVector &solutionCoefficientsNew(ptrFace->getPtrElementLeft()->getTimeLevelDataVector(timeLevelResult));
                
                const LinearAlgebra::Matrix &stiffnessMatrix = ptrFace->getFaceMatrix(stiffnessFaceMatrixID_).getElementMatrix(Base::Side::LEFT, Base::Side::LEFT);
                
                solutionCoefficientsNew += stiffnessMatrix * solutionCoefficients;
            }
        }
        
        synchronize(timeLevelResult);
    }
    
    void HpgemAPILinear::addSourceTerm(const std::size_t timeLevelResult, const double time, const std::size_t orderTimeDerivative)
    {
        
        // Add the source terms corresponding to the elements.
        for (Base::Element *ptrElement : meshes_[0]->getElementsList())
        {
            LinearAlgebra::NumericalVector &solutionCoefficientsNew(ptrElement->getTimeLevelDataVector(timeLevelResult));
            
            solutionCoefficientsNew += integrateSourceTermAtElement(ptrElement, time, orderTimeDerivative);
        }
        synchronize(timeLevelResult);
    }
    
    void HpgemAPILinear::addSourceTermAtBoundary(const std::size_t timeLevelResult, const double time, const std::size_t orderTimeDerivative)
    {
        // Add the source terms corresponding to the faces at the boundary
        for (Base::Face *ptrFace : meshes_[0]->getFacesList())
        {
            if(!ptrFace->isInternal())
            {
                LinearAlgebra::NumericalVector &solutionCoefficientsNew(ptrFace->getPtrElementLeft()->getTimeLevelDataVector(timeLevelResult));
                
                solutionCoefficientsNew += integrateSourceTermAtFace(ptrFace, time, orderTimeDerivative);
            }
        }
        synchronize(timeLevelResult);
    }
    
    /// \details Make sure timeLevelIn is different from timeLevelResult.
    void HpgemAPILinear::computeRightHandSide(const std::size_t timeLevelIn, const std::size_t timeLevelResult, const double time)
    {
        multiplyStiffnessMatrices(timeLevelIn, timeLevelResult);
        if(useSourceTerm_)
        {
            addSourceTerm(timeLevelResult, time, 0);
        }
        if(useSourceTermAtBoundary_)
        {
            addSourceTermAtBoundary(timeLevelResult, time, 0);
        }
    }
    
    /// \details Make sure timeLevelResult is different from the timeLevelsIn.
    void HpgemAPILinear::computeRightHandSide(const std::vector<std::size_t> timeLevelsIn, const std::vector<double> coefficientsTimeLevels, const std::size_t timeLevelResult, const double time)
    {
        multiplyStiffnessMatrices(timeLevelsIn, coefficientsTimeLevels, timeLevelResult);
        if(useSourceTerm_)
        {
            addSourceTerm(timeLevelResult, time, 0);
        }
        if(useSourceTermAtBoundary_)
        {
            addSourceTermAtBoundary(timeLevelResult, time, 0);
        }
    }
    
    void HpgemAPILinear::tasksBeforeSolving()
    {
        logger(INFO, "Computing the mass matrices.");
        createMassMatrices();
        logger(INFO, "Computing stiffness matrices.");
        createStiffnessMatrices();
    }
}



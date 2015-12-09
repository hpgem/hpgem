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
    /// \param[in] numberOfVariables Number of variables in the PDE
    /// \param[in] polynomialOrder Polynomial order of the basis functions
    /// \param[in] ptrButcherTableau A butcherTableau used to solve the PDE with a Runge-Kutta method.
    /// \param[in] numberOfTimeLevels Number of time levels.
    /// \param[in] useSourceTerm Boolean to indicate if there is a source term.
    /// \param[in] useSourceTermAtBoundary Boolean to indicate if there is a source term at the domain boundary.
    template<std::size_t DIM>
    HpgemAPILinear<DIM>::HpgemAPILinear
    (
     const std::size_t numberOfVariables,
     const std::size_t polynomialOrder,
     const TimeIntegration::ButcherTableau * const ptrButcherTableau,
     const std::size_t numberOfTimeLevels,
     const bool useSourceTerm,
     const bool useSourceTermAtBoundary
     ) :
    HpgemAPISimplified<DIM>(numberOfVariables, polynomialOrder, ptrButcherTableau, numberOfTimeLevels),
    useSourceTerm_(useSourceTerm),
    useSourceTermAtBoundary_(useSourceTermAtBoundary),
    massMatrixID_(1),
    stiffnessElementMatrixID_(0),
    stiffnessFaceMatrixID_(0)
    {
    }
    
    /// \param[in] dimension Dimension of the domain
    /// \param[in] numberOfVariables Number of variables in the PDE
    /// \param[in] polynomialOrder Polynomial order of the basis functions
    /// \param[in] ptrButcherTableau A butcherTableau used to solve the PDE with a Runge-Kutta method.
    /// \param[in] numberOfTimeLevels Number of time levels.
    /// \param[in] useSourceTerm Boolean to indicate if there is a source term.
    /// \param[in] useSourceTermAtBoundary Boolean to indicate if there is a source term at the domain boundary.
    template<std::size_t DIM>
    HpgemAPILinear<DIM>::HpgemAPILinear
    (
     const std::size_t numberOfVariables,
     const std::size_t polynomialOrder,
     const std::size_t globalNummberOfTimeIntegrationVectors,
     const std::size_t numberOfTimeLevels,
     const bool useSourceTerm,
     const bool useSourceTermAtBoundary
     ) :
    HpgemAPISimplified<DIM>(numberOfVariables, polynomialOrder, globalNummberOfTimeIntegrationVectors, numberOfTimeLevels),
    useSourceTerm_(useSourceTerm),
    useSourceTermAtBoundary_(useSourceTermAtBoundary),
    massMatrixID_(1),
    stiffnessElementMatrixID_(0),
    stiffnessFaceMatrixID_(0)
    {
    }

    template<std::size_t DIM>
    void HpgemAPILinear<DIM>::createMesh(const std::size_t numberOfElementsPerDirection, const Base::MeshType meshType)
    {
        const Base::RectangularMeshDescriptor<DIM> description = this->createMeshDescription(numberOfElementsPerDirection);
        
        // Set the number of Element/Face Matrices/Vectors.
        std::size_t numberOfElementMatrices = 2;   // Mass matrix and stiffness matrix
        std::size_t numberOfElementVectors = 0;
        std::size_t numberOfFaceMatrices = 1;      // Stiffness matrix
        std::size_t numberOfFaceVectors = 0;

        // Create mesh and set basis functions.
        this->addMesh(description, meshType, numberOfElementMatrices, numberOfElementVectors, numberOfFaceMatrices, numberOfFaceVectors);
        this->meshes_[0]->useDefaultDGBasisFunctions();
        
        // Set the number of time integration vectors according to the size of the Butcher tableau.
        this->setNumberOfTimeIntegrationVectorsGlobally(this->globalNumberOfTimeIntegrationVectors_);
        
        // Plot info about the mesh
        std::size_t nElements = this->meshes_[0]->getNumberOfElements();
        logger(VERBOSE, "Total number of elements: %", nElements);
    }

    template<std::size_t DIM>
    void HpgemAPILinear<DIM>::createMassMatrices()
    {
        for (Base::Element *ptrElement : this->meshes_[0]->getElementsList())
        {
            LinearAlgebra::MiddleSizeMatrix massMatrix(this->computeMassMatrixAtElement(ptrElement));
            logger(VERBOSE, "--Mass matrix element:\n %", massMatrix);
            ptrElement->setElementMatrix(massMatrix, massMatrixID_);
        }
    }
    
    /// \details Solve the equation \f$ Mu = r \f$ for \f$ u \f$, where \f$ r \f$ is the right-hand sid and \f$ M \f$ is the mass matrix.
    template<std::size_t DIM>
    void HpgemAPILinear<DIM>::solveMassMatrixEquations(const std::size_t timeIntegrationVectorId)
    {
        for (Base::Element *ptrElement : this->meshes_[0]->getElementsList())
        {
            LinearAlgebra::MiddleSizeVector &functionCoefficients(ptrElement->getTimeIntegrationVector(timeIntegrationVectorId));
            
            const LinearAlgebra::MiddleSizeMatrix &massMatrix(ptrElement->getElementMatrix(massMatrixID_));
            massMatrix.solve(functionCoefficients);
        }
        
        this->synchronize(timeIntegrationVectorId);
    }

    template<std::size_t DIM>
    LinearAlgebra::MiddleSizeMatrix HpgemAPILinear<DIM>::computeStiffnessMatrixAtElement(Base::Element *ptrElement)
    {
        // Define a function for the integrand of the stiffness matrix at the element.
        std::function<LinearAlgebra::MiddleSizeMatrix(Base::PhysicalElement<DIM>&)> integrandFunction = [=](Base::PhysicalElement<DIM> &element) -> LinearAlgebra::MiddleSizeMatrix
        {return this->computeIntegrandStiffnessMatrixAtElement(element);};
        
        return this->elementIntegrator_.integrate(ptrElement, integrandFunction);
    }

    template<std::size_t DIM>
    Base::FaceMatrix HpgemAPILinear<DIM>::computeStiffnessMatrixAtFace(Base::Face *ptrFace)
    {
        // Define a function for the integrand of the stiffness matrix at a face.
        std::function<Base::FaceMatrix(Base::PhysicalFace<DIM> &)> integrandFunction = [=](Base::PhysicalFace<DIM> &face) -> Base::FaceMatrix
        {return this->computeIntegrandStiffnessMatrixAtFace(face);};
        
        return this->faceIntegrator_.integrate(ptrFace, integrandFunction);
    }

    template<std::size_t DIM>
    LinearAlgebra::MiddleSizeVector HpgemAPILinear<DIM>::integrateSourceTermAtFace(Base::Face *ptrFace, const double time, const std::size_t orderTimeDerivative)
    {
        // Define a function for the integrand of the stiffness matrix at a element.
        std::function<LinearAlgebra::MiddleSizeVector(Base::PhysicalFace<DIM>&)> integrandFunction = [=](Base::PhysicalFace<DIM> &face) -> LinearAlgebra::MiddleSizeVector
        {return this->computeIntegrandSourceTermAtFace(face, time, orderTimeDerivative);};
        
        return this->faceIntegrator_.integrate(ptrFace, integrandFunction);
    }
    
    template<std::size_t DIM>
    LinearAlgebra::MiddleSizeVector HpgemAPILinear<DIM>::computeIntegrandSourceTermAtElement(Base::PhysicalElement<DIM> &element, const double time, const std::size_t orderTimeDerivative)
    {
        // Get a reference to the result vector.
        LinearAlgebra::MiddleSizeVector &integrand = element.getResultVector();
        
        // Get the physical point.
        PointPhysicalT pPhys = element.getPointPhysical();
        
        // Compute the source term.
        LinearAlgebra::MiddleSizeVector sourceTerm = getSourceTerm(pPhys, time, orderTimeDerivative);
        
        // Get the number of basis functions.
        const std::size_t numberOfBasisFunctions = element.getElement()->getNrOfBasisFunctions();
        
        // Compute the product of the source term and all test functions.
        std::size_t iVB, jVB; // indices for both variable and basis function.
        for (std::size_t iV = 0; iV < this->configData_->numberOfUnknowns_; iV++)
        {
            for (std::size_t iB = 0; iB < numberOfBasisFunctions; iB++)
            {
                iVB = element.getElement()->convertToSingleIndex(iB, iV);
                integrand(iVB) = element.basisFunction(iB) * sourceTerm(iV);
            }
        }
        
        return integrand;
    }
    
    /// \details By default, the standard L2 inner product with the source term is computed.
    /// \todo please use Integration::ElementIntegral::integrate() for integration over elements
    template<std::size_t DIM>
    LinearAlgebra::MiddleSizeVector HpgemAPILinear<DIM>::integrateSourceTermAtElement(Base::Element * ptrElement, const double time, const std::size_t orderTimeDerivative)
    {
        // Define the integrand function for the the source term.
        std::function<LinearAlgebra::MiddleSizeVector(Base::PhysicalElement<DIM>&)> integrandFunction = [=](Base::PhysicalElement<DIM>& element) -> LinearAlgebra::MiddleSizeVector { return this -> computeIntegrandSourceTermAtElement(element, time, orderTimeDerivative);};
        
        return this->elementIntegrator_.integrate(ptrElement, integrandFunction);
    }

    template<std::size_t DIM>
    void HpgemAPILinear<DIM>::createStiffnessMatrices()
    {
        logger(INFO, "- Creating stiffness matrices for the elements.");
        for (Base::Element *ptrElement : this->meshes_[0]->getElementsList())
        {
            LinearAlgebra::MiddleSizeMatrix stiffnessMatrix(computeStiffnessMatrixAtElement(ptrElement));
            //std::cout << "-- Stiffness matrix element:\n" << stiffnessMatrix << "\n";
            ptrElement->setElementMatrix(stiffnessMatrix, stiffnessElementMatrixID_);
        }
        
        logger(INFO, "- Creating stiffness matrices for the faces.");
        for (Base::Face *ptrFace : this->meshes_[0]->getFacesList())
        {
            Base::FaceMatrix stiffnessFaceMatrix(computeStiffnessMatrixAtFace(ptrFace));
            if(!ptrFace->isInternal())
            {
                logger.assert(stiffnessFaceMatrix.getNumberOfDegreesOfFreedom(Base::Side::RIGHT) == 0,"The number of degrees of freedom corresonding to the right side of a boundary face should be 0, but is here %.", stiffnessFaceMatrix.getNumberOfDegreesOfFreedom(Base::Side::RIGHT));
            }
            //std::cout << "-- Stiffness matrix face: \n";
            //std::cout << "--- Stiffness submatrix face:\n" << stiffnessFaceMatrix.getElementMatrix(Base::Side::LEFT, Base::Side::LEFT) << "\n";
            //std::cout << "--- Stiffness submatrix face:\n" << stiffnessFaceMatrix.getElementMatrix(Base::Side::LEFT, Base::Side::RIGHT) << "\n";
            //std::cout << "--- Stiffness submatrix face:\n" << stiffnessFaceMatrix.getElementMatrix(Base::Side::RIGHT, Base::Side::LEFT) << "\n";
            //std::cout << "--- Stiffness submatrix face:\n" << stiffnessFaceMatrix.getElementMatrix(Base::Side::RIGHT, Base::Side::RIGHT) << "\n";
            
            ptrFace->setFaceMatrix(stiffnessFaceMatrix, stiffnessFaceMatrixID_);
        }
    }
    
    /// \details Make sure inputVectorId is different from resultVectorId.
    template<std::size_t DIM>
    void HpgemAPILinear<DIM>::multiplyStiffnessMatrices(const std::size_t inputVectorId, const std::size_t resultVectorId)
    {
        // Multiply the stiffness matrices corresponding to the elements.
        for (Base::Element *ptrElement : this->meshes_[0]->getElementsList())
        {
            LinearAlgebra::MiddleSizeVector &inputFunctionCoefficients(ptrElement->getTimeIntegrationVector(inputVectorId));
            LinearAlgebra::MiddleSizeVector &resultFunctionCoefficients(ptrElement->getTimeIntegrationVector(resultVectorId));
            
            resultFunctionCoefficients = ptrElement->getElementMatrix(stiffnessElementMatrixID_) * inputFunctionCoefficients;
        }
        
        // Multiply the stiffness matrices corresponding to the faces.
        for (Base::Face *ptrFace : this->meshes_[0]->getFacesList())
        {
            if(ptrFace->isInternal())
            {
                LinearAlgebra::MiddleSizeVector &inputFunctionCoefficientsLeft(ptrFace->getPtrElementLeft()->getTimeIntegrationVector(inputVectorId));
                LinearAlgebra::MiddleSizeVector &inputFunctionCoefficientsRight(ptrFace->getPtrElementRight()->getTimeIntegrationVector(inputVectorId));
                LinearAlgebra::MiddleSizeVector &resultFunctionCoefficientsLeft(ptrFace->getPtrElementLeft()->getTimeIntegrationVector(resultVectorId));
                LinearAlgebra::MiddleSizeVector &resultFunctionCoefficientsRight(ptrFace->getPtrElementRight()->getTimeIntegrationVector(resultVectorId));
                
                const Base::FaceMatrix &stiffnessFaceMatrix = ptrFace->getFaceMatrix(stiffnessFaceMatrixID_);
                
                resultFunctionCoefficientsLeft += stiffnessFaceMatrix.getElementMatrix(Base::Side::LEFT, Base::Side::LEFT) * inputFunctionCoefficientsLeft + stiffnessFaceMatrix.getElementMatrix(Base::Side::LEFT, Base::Side::RIGHT) * inputFunctionCoefficientsRight;
                resultFunctionCoefficientsRight += stiffnessFaceMatrix.getElementMatrix(Base::Side::RIGHT, Base::Side::LEFT) * inputFunctionCoefficientsLeft + stiffnessFaceMatrix.getElementMatrix(Base::Side::RIGHT, Base::Side::RIGHT) * inputFunctionCoefficientsRight;
            }
            else
            {
                LinearAlgebra::MiddleSizeVector &inputFunctionCoefficients(ptrFace->getPtrElementLeft()->getTimeIntegrationVector(inputVectorId));
                LinearAlgebra::MiddleSizeVector &resultFunctionCoefficients(ptrFace->getPtrElementLeft()->getTimeIntegrationVector(resultVectorId));
                
                const LinearAlgebra::MiddleSizeMatrix &stiffnessMatrix = ptrFace->getFaceMatrix(stiffnessFaceMatrixID_).getElementMatrix(Base::Side::LEFT, Base::Side::LEFT);
                
                resultFunctionCoefficients += stiffnessMatrix * inputFunctionCoefficients;
            }
        }
        
        this->synchronize(resultVectorId);
    }
    
    /// \details Make sure inputVectorIds are different from resultVectorId.
    template<std::size_t DIM>
    void HpgemAPILinear<DIM>::multiplyStiffnessMatrices(const std::vector<std::size_t> inputVectorIds, const std::vector<double> coefficientsInputVectors, const std::size_t resultVectorId)
    {
        // Multiply the stiffness matrices corresponding to the elements.
        for (Base::Element *ptrElement : this->meshes_[0]->getElementsList())
        {
            LinearAlgebra::MiddleSizeVector inputFunctionCoefficients(this->getLinearCombinationOfVectors(ptrElement, inputVectorIds, coefficientsInputVectors));
            LinearAlgebra::MiddleSizeVector &resultFunctionCoefficients(ptrElement->getTimeIntegrationVector(resultVectorId));
            
            resultFunctionCoefficients = ptrElement->getElementMatrix(stiffnessElementMatrixID_) * inputFunctionCoefficients;
        }
        
        // Multiply the stiffness matrices corresponding to the faces.
        for (Base::Face *ptrFace : this->meshes_[0]->getFacesList())
        {
            if(ptrFace->isInternal())
            {
                LinearAlgebra::MiddleSizeVector inputFunctionCoefficientsLeft(this->getLinearCombinationOfVectors(ptrFace->getPtrElementLeft(), inputVectorIds, coefficientsInputVectors));
                LinearAlgebra::MiddleSizeVector inputFunctionCoefficientsRight(this->getLinearCombinationOfVectors(ptrFace->getPtrElementRight(), inputVectorIds, coefficientsInputVectors));
                LinearAlgebra::MiddleSizeVector &resultFunctionCoefficientsLeft(ptrFace->getPtrElementLeft()->getTimeIntegrationVector(resultVectorId));
                LinearAlgebra::MiddleSizeVector &resultFunctionCoefficientsRight(ptrFace->getPtrElementRight()->getTimeIntegrationVector(resultVectorId));
                
                const Base::FaceMatrix &stiffnessFaceMatrix = ptrFace->getFaceMatrix(stiffnessFaceMatrixID_);
                
                resultFunctionCoefficientsLeft += stiffnessFaceMatrix.getElementMatrix(Base::Side::LEFT, Base::Side::LEFT) * inputFunctionCoefficientsLeft;
                resultFunctionCoefficientsLeft += stiffnessFaceMatrix.getElementMatrix(Base::Side::LEFT, Base::Side::RIGHT) * inputFunctionCoefficientsRight;
                resultFunctionCoefficientsRight += stiffnessFaceMatrix.getElementMatrix(Base::Side::RIGHT, Base::Side::LEFT) * inputFunctionCoefficientsLeft;
                resultFunctionCoefficientsRight += stiffnessFaceMatrix.getElementMatrix(Base::Side::RIGHT, Base::Side::RIGHT) * inputFunctionCoefficientsRight;
            }
            else
            {
                LinearAlgebra::MiddleSizeVector inputFunctionCoefficients(this->getLinearCombinationOfVectors(ptrFace->getPtrElementLeft(), inputVectorIds, coefficientsInputVectors));
                LinearAlgebra::MiddleSizeVector &resultFunctionCoefficients(ptrFace->getPtrElementLeft()->getTimeIntegrationVector(resultVectorId));
                
                const LinearAlgebra::MiddleSizeMatrix &stiffnessMatrix = ptrFace->getFaceMatrix(stiffnessFaceMatrixID_).getElementMatrix(Base::Side::LEFT, Base::Side::LEFT);
                
                resultFunctionCoefficients += stiffnessMatrix * inputFunctionCoefficients;
            }
        }
        
        this->synchronize(resultVectorId);
    }

    template<std::size_t DIM>
    void HpgemAPILinear<DIM>::addSourceTerm(const std::size_t resultVectorId, const double time, const std::size_t orderTimeDerivative)
    {
        
        // Add the source terms corresponding to the elements.
        for (Base::Element *ptrElement : this->meshes_[0]->getElementsList())
        {
            LinearAlgebra::MiddleSizeVector &resultFunctionCoefficientsNew(ptrElement->getTimeIntegrationVector(resultVectorId));
            
            resultFunctionCoefficientsNew += integrateSourceTermAtElement(ptrElement, time, orderTimeDerivative);
        }
        this->synchronize(resultVectorId);
    }

    template<std::size_t DIM>
    void HpgemAPILinear<DIM>::addSourceTermAtBoundary(const std::size_t resultVectorId, const double time, const std::size_t orderTimeDerivative)
    {
        // Add the source terms corresponding to the faces at the boundary
        for (Base::Face *ptrFace : this->meshes_[0]->getFacesList())
        {
            if(!ptrFace->isInternal())
            {
                LinearAlgebra::MiddleSizeVector &resultFunctionCoefficientsNew(ptrFace->getPtrElementLeft()->getTimeIntegrationVector(resultVectorId));
                
                resultFunctionCoefficientsNew += integrateSourceTermAtFace(ptrFace, time, orderTimeDerivative);
            }
        }
        this->synchronize(resultVectorId);
    }
    
    /// \details Make sure inputVectorId is different from resultVectorId.
    template<std::size_t DIM>
    void HpgemAPILinear<DIM>::computeRightHandSide(const std::size_t inputVectorId, const std::size_t resultVectorId, const double time)
    {
        multiplyStiffnessMatrices(inputVectorId, resultVectorId);
        if(useSourceTerm_)
        {
            addSourceTerm(resultVectorId, time, 0);
        }
        if(useSourceTermAtBoundary_)
        {
            addSourceTermAtBoundary(resultVectorId, time, 0);
        }
    }
    
    /// \details Make sure resultVectorId is different from the inputVectorIds.
    template<std::size_t DIM>
    void HpgemAPILinear<DIM>::computeRightHandSide(const std::vector<std::size_t> inputVectorIds, const std::vector<double> coefficientsInputVectors, const std::size_t resultVectorId, const double time)
    {
        multiplyStiffnessMatrices(inputVectorIds, coefficientsInputVectors, resultVectorId);
        if(useSourceTerm_)
        {
            addSourceTerm(resultVectorId, time, 0);
        }
        if(useSourceTermAtBoundary_)
        {
            addSourceTermAtBoundary(resultVectorId, time, 0);
        }
    }

    template<std::size_t DIM>
    void HpgemAPILinear<DIM>::tasksBeforeSolving()
    {
        logger(INFO, "Computing the mass matrices.");
        createMassMatrices();
        logger(INFO, "Computing stiffness matrices.");
        createStiffnessMatrices();
    }
}



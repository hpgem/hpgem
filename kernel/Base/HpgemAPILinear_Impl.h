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
    /// \param[in] ptrButcherTableau A butcherTableau used to solve the PDE with a Runge-Kutta method.
    /// \param[in] numberOfTimeLevels Number of time levels.
    /// \param[in] useSourceTerm Boolean to indicate if there is a source term.
    /// \param[in] useSourceTermAtBoundary Boolean to indicate if there is a source term at the domain boundary.
    template<std::size_t DIM>
    HpgemAPILinear<DIM>::HpgemAPILinear
    (
     const std::size_t numOfVariables,
     const std::size_t polynomialOrder,
     const TimeIntegration::ButcherTableau * const ptrButcherTableau,
     const std::size_t numberOfTimeLevels,
     const bool useSourceTerm,
     const bool useSourceTermAtBoundary
     ) :
    HpgemAPISimplified<DIM>(numOfVariables, polynomialOrder, ptrButcherTableau, numberOfTimeLevels),
    useSourceTerm_(useSourceTerm),
    useSourceTermAtBoundary_(useSourceTermAtBoundary),
    massMatrixID_(1),
    stiffnessElementMatrixID_(0),
    stiffnessFaceMatrixID_(0)
    {
    }
    
    /// \param[in] dimension Dimension of the domain
    /// \param[in] numOfVariables Number of variables in the PDE
    /// \param[in] polynomialOrder Polynomial order of the basis functions
    /// \param[in] ptrButcherTableau A butcherTableau used to solve the PDE with a Runge-Kutta method.
    /// \param[in] numberOfTimeLevels Number of time levels.
    /// \param[in] useSourceTerm Boolean to indicate if there is a source term.
    /// \param[in] useSourceTermAtBoundary Boolean to indicate if there is a source term at the domain boundary.
    template<std::size_t DIM>
    HpgemAPILinear<DIM>::HpgemAPILinear
    (
     const std::size_t numOfVariables,
     const std::size_t polynomialOrder,
     const std::size_t globalNummberOfTimeIntegrationVectors,
     const std::size_t numberOfTimeLevels,
     const bool useSourceTerm,
     const bool useSourceTermAtBoundary
     ) :
    HpgemAPISimplified<DIM>(numOfVariables, polynomialOrder, globalNummberOfTimeIntegrationVectors, numberOfTimeLevels),
    useSourceTerm_(useSourceTerm),
    useSourceTermAtBoundary_(useSourceTermAtBoundary),
    massMatrixID_(1),
    stiffnessElementMatrixID_(0),
    stiffnessFaceMatrixID_(0)
    {
    }

    template<std::size_t DIM>
    void HpgemAPILinear<DIM>::createMesh(const std::size_t numOfElementsPerDirection, const Base::MeshType meshType)
    {
        const Base::RectangularMeshDescriptor<DIM> description = this->createMeshDescription(numOfElementsPerDirection);
        
        // Set the number of Element/Face Matrices/Vectors.
        std::size_t numOfElementMatrices = 2;   // Mass matrix and stiffness matrix
        std::size_t numOfElementVectors = 0;
        std::size_t numOfFaceMatrices = 1;      // Stiffness matrix
        std::size_t numOfFaceVectors = 0;

        // Create mesh and set basis functions.
        this->addMesh(description, meshType, numOfElementMatrices, numOfElementVectors, numOfFaceMatrices, numOfFaceVectors);
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
            //std::cout << "--Mass matrix element:\n" << massMatrix << "\n";
            ptrElement->setElementMatrix(massMatrix, massMatrixID_);
        }
    }
    
    /// \details Solve the equation \f$ Mu = r \f$ for \f$ u \f$, where \f$ r \f$ is the right-hand sid and \f$ M \f$ is the mass matrix.
    template<std::size_t DIM>
    void HpgemAPILinear<DIM>::solveMassMatrixEquations(const std::size_t timeIntegrationVectorId)
    {
        for (Base::Element *ptrElement : this->meshes_[0]->getElementsList())
        {
            LinearAlgebra::MiddleSizeVector &solutionCoefficients(ptrElement->getTimeIntegrationVector(timeIntegrationVectorId));
            
            const LinearAlgebra::MiddleSizeMatrix &massMatrix(ptrElement->getElementMatrix(massMatrixID_));
            massMatrix.solve(solutionCoefficients);
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
    
    /// \details By default, the standard L2 inner product with the source term is computed.
    /// \todo please use Integration::ElementIntegral::integrate() for integration over elements
    template<std::size_t DIM>
    LinearAlgebra::MiddleSizeVector HpgemAPILinear<DIM>::integrateSourceTermAtElement(Base::Element * ptrElement, const double time, const std::size_t orderTimeDerivative)
    {
        // Get number of basis functions
        std::size_t numOfBasisFunctions = ptrElement->getNumberOfBasisFunctions();
        
        // Declare integral source term
        LinearAlgebra::MiddleSizeVector integralSourceTerm(numOfBasisFunctions * this->configData_->numberOfUnknowns_);
        
        // Declare integrand
        LinearAlgebra::MiddleSizeVector integrandSourceTerm(numOfBasisFunctions * this->configData_->numberOfUnknowns_);
        
        // Get quadrature rule and number of points.
        const QuadratureRules::GaussQuadratureRule *ptrQdrRule = ptrElement->getGaussQuadratureRule();
        std::size_t numOfQuadPoints = ptrQdrRule->getNumberOfPoints();
        
        // For each quadrature point, compute the value of the product of the
        // test function and the source term, then add it with the correct weight to the integral solution.
        for (std::size_t pQuad = 0; pQuad < numOfQuadPoints; ++pQuad)
        {
            const Geometry::PointReference<DIM>& pRef = ptrQdrRule->getPoint(pQuad);
            Geometry::PointPhysical<DIM> pPhys = ptrElement->referenceToPhysical(pRef);
            
            Geometry::Jacobian<DIM, DIM> jac = ptrElement->calcJacobian(pRef);
            
            LinearAlgebra::MiddleSizeVector sourceTerm = getSourceTerm(pPhys, time, orderTimeDerivative);
            
            for(std::size_t iB = 0; iB < numOfBasisFunctions; ++iB)
            {
                double valueBasisFunction = ptrElement->basisFunction(iB, pRef);
                
                for(std::size_t iV = 0; iV < this->configData_->numberOfUnknowns_; ++iV)
                {
                    std::size_t iVB = ptrElement->convertToSingleIndex(iB,iV);
                    
                    integrandSourceTerm(iVB) = sourceTerm(iV) * valueBasisFunction;
                }
            }
            integralSourceTerm.axpy((ptrQdrRule->weight(pQuad)) * std::abs(jac.determinant()), integrandSourceTerm);
        }
        
        return integralSourceTerm;
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
            LinearAlgebra::MiddleSizeVector &solutionCoefficients(ptrElement->getTimeIntegrationVector(inputVectorId));
            LinearAlgebra::MiddleSizeVector &solutionCoefficientsNew(ptrElement->getTimeIntegrationVector(resultVectorId));
            
            solutionCoefficientsNew = ptrElement->getElementMatrix(stiffnessElementMatrixID_) * solutionCoefficients;
        }
        
        // Multiply the stiffness matrices corresponding to the faces.
        for (Base::Face *ptrFace : this->meshes_[0]->getFacesList())
        {
            if(ptrFace->isInternal())
            {
                LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeft(ptrFace->getPtrElementLeft()->getTimeIntegrationVector(inputVectorId));
                LinearAlgebra::MiddleSizeVector &solutionCoefficientsRight(ptrFace->getPtrElementRight()->getTimeIntegrationVector(inputVectorId));
                LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeftNew(ptrFace->getPtrElementLeft()->getTimeIntegrationVector(resultVectorId));
                LinearAlgebra::MiddleSizeVector &solutionCoefficientsRightNew(ptrFace->getPtrElementRight()->getTimeIntegrationVector(resultVectorId));
                
                const Base::FaceMatrix &stiffnessFaceMatrix = ptrFace->getFaceMatrix(stiffnessFaceMatrixID_);
                
                solutionCoefficientsLeftNew += stiffnessFaceMatrix.getElementMatrix(Base::Side::LEFT, Base::Side::LEFT) * solutionCoefficientsLeft + stiffnessFaceMatrix.getElementMatrix(Base::Side::LEFT, Base::Side::RIGHT) * solutionCoefficientsRight;
                solutionCoefficientsRightNew += stiffnessFaceMatrix.getElementMatrix(Base::Side::RIGHT, Base::Side::LEFT) * solutionCoefficientsLeft + stiffnessFaceMatrix.getElementMatrix(Base::Side::RIGHT, Base::Side::RIGHT) * solutionCoefficientsRight;
            }
            else
            {
                LinearAlgebra::MiddleSizeVector &solutionCoefficients(ptrFace->getPtrElementLeft()->getTimeIntegrationVector(inputVectorId));
                LinearAlgebra::MiddleSizeVector &solutionCoefficientsNew(ptrFace->getPtrElementLeft()->getTimeIntegrationVector(resultVectorId));
                
                const LinearAlgebra::MiddleSizeMatrix &stiffnessMatrix = ptrFace->getFaceMatrix(stiffnessFaceMatrixID_).getElementMatrix(Base::Side::LEFT, Base::Side::LEFT);
                
                solutionCoefficientsNew += stiffnessMatrix * solutionCoefficients;
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
            LinearAlgebra::MiddleSizeVector solutionCoefficients(this->getLinearCombinationOfVectors(ptrElement, inputVectorIds, coefficientsInputVectors));
            LinearAlgebra::MiddleSizeVector &solutionCoefficientsNew(ptrElement->getTimeIntegrationVector(resultVectorId));
            
            solutionCoefficientsNew = ptrElement->getElementMatrix(stiffnessElementMatrixID_) * solutionCoefficients;
        }
        
        // Multiply the stiffness matrices corresponding to the faces.
        for (Base::Face *ptrFace : this->meshes_[0]->getFacesList())
        {
            if(ptrFace->isInternal())
            {
                LinearAlgebra::MiddleSizeVector solutionCoefficientsLeft(this->getLinearCombinationOfVectors(ptrFace->getPtrElementLeft(), inputVectorIds, coefficientsInputVectors));
                LinearAlgebra::MiddleSizeVector solutionCoefficientsRight(this->getLinearCombinationOfVectors(ptrFace->getPtrElementRight(), inputVectorIds, coefficientsInputVectors));
                LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeftNew(ptrFace->getPtrElementLeft()->getTimeIntegrationVector(resultVectorId));
                LinearAlgebra::MiddleSizeVector &solutionCoefficientsRightNew(ptrFace->getPtrElementRight()->getTimeIntegrationVector(resultVectorId));
                
                const Base::FaceMatrix &stiffnessFaceMatrix = ptrFace->getFaceMatrix(stiffnessFaceMatrixID_);
                
                solutionCoefficientsLeftNew += stiffnessFaceMatrix.getElementMatrix(Base::Side::LEFT, Base::Side::LEFT) * solutionCoefficientsLeft;
                solutionCoefficientsLeftNew += stiffnessFaceMatrix.getElementMatrix(Base::Side::LEFT, Base::Side::RIGHT) * solutionCoefficientsRight;
                solutionCoefficientsRightNew += stiffnessFaceMatrix.getElementMatrix(Base::Side::RIGHT, Base::Side::LEFT) * solutionCoefficientsLeft;
                solutionCoefficientsRightNew += stiffnessFaceMatrix.getElementMatrix(Base::Side::RIGHT, Base::Side::RIGHT) * solutionCoefficientsRight;
            }
            else
            {
                LinearAlgebra::MiddleSizeVector solutionCoefficients(this->getLinearCombinationOfVectors(ptrFace->getPtrElementLeft(), inputVectorIds, coefficientsInputVectors));
                LinearAlgebra::MiddleSizeVector &solutionCoefficientsNew(ptrFace->getPtrElementLeft()->getTimeIntegrationVector(resultVectorId));
                
                const LinearAlgebra::MiddleSizeMatrix &stiffnessMatrix = ptrFace->getFaceMatrix(stiffnessFaceMatrixID_).getElementMatrix(Base::Side::LEFT, Base::Side::LEFT);
                
                solutionCoefficientsNew += stiffnessMatrix * solutionCoefficients;
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
            LinearAlgebra::MiddleSizeVector &solutionCoefficientsNew(ptrElement->getTimeIntegrationVector(resultVectorId));
            
            solutionCoefficientsNew += integrateSourceTermAtElement(ptrElement, time, orderTimeDerivative);
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
                LinearAlgebra::MiddleSizeVector &solutionCoefficientsNew(ptrFace->getPtrElementLeft()->getTimeIntegrationVector(resultVectorId));
                
                solutionCoefficientsNew += integrateSourceTermAtFace(ptrFace, time, orderTimeDerivative);
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



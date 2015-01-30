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

#include "Base/CommandLineOptions.hpp"
#include "Base/ConfigurationData.hpp"
#include "Base/Element.hpp"
#include "Base/Face.hpp"
#include "Base/HpgemUI.hpp"
#include "Base/RectangularMeshDescriptor.hpp"
#include "Base/TimeIntegration/AllTimeIntegrators.hpp"
#include "Geometry/PointReference.hpp"
#include "Integration/ElementIntegral.hpp"
#include "Output/TecplotDiscontinuousSolutionWriter.hpp"
#include "Output/TecplotSingleElementWriter.hpp"
#include "Utilities/BasisFunctions1DH1ConformingLine.hpp"
#include "Utilities/BasisFunctions2DH1ConformingSquare.hpp"
#include "Utilities/BasisFunctions2DH1ConformingTriangle.hpp"
#include "Utilities/BasisFunctions3DH1ConformingCube.hpp"
#include "Utilities/BasisFunctions3DH1ConformingTetrahedron.hpp"


/// \brief Class to demonstrate how a system of multiple PDE's can be solved.
/** \details 
This class illustrates how to solve a system of PDE's. The example can be used for solving linear PDE's as well as non-linear PDE's with a time dependent mass matrices and right hand side. Currently only periodic boundary conditions are used in this example.
 
The class can solve the scalar wave equation in 2D or 3D written as a first order system of PDE's. The original scalar wave equation is given by \f[ \partial_t^2 u = \nabla \cdot c \nabla u \f], where \f$u \f$ is the scalar variable and \f$c \f$ a material parameter corresponding to the velocity with which waves can propagate.
 
For a first order scheme we define the scalar function \f$ v := \partial_t u \f$ and the vector function \f$ s := c \nabla u \f$. We can then obtain the equations \f[ \partial_t v = \nabla \cdot s \f] and \f[ c^{-1} \partial_t s = \nabla v \f]
 
We define a new vector function \f$w = [w_0, w_1, w_2] = [v, s_0, s_1]\f$. We can then rewrite the system of PDE's as follows: \f[ \partial_t w_0 = \partial_i w_{i+1} \f] summing over i = 0 .. DIM and \f[ c^{-1} \partial_t w_{i+1} = \partial_i w_0 \f] for i = 0 .. DIM.
 
 A boolean useMatrixStorage can be set to true if you want to store all mass and stiffness matrices. This only works for linear problems, but might reduce computational costs.
 */

class ExampleMultipleVariableProblem : public Base::HpgemUI, public Output::TecplotSingleElementWriter
{
public:
    /// \param[in] DIM Dimension of the domain
    /// \param[in] n Number of elements per direction
    /// \param[in] p Polynomial order of the basis functions
    /// \param[in] ptrButcherTableau Pointer to a Butcher Tableau used to solve the time integration
    /// \param[in] useMatrixStorage Boolean to indicate if element and face matrices for the PDE should be stored
    ExampleMultipleVariableProblem
    (
        const std::size_t DIM,
        const std::size_t n,
        const std::size_t p,
        const Base::ButcherTableau *ptrButcherTableau,
        const bool useMatrixStorage = false
    ) :
    HpgemUI(new Base::GlobalData, new Base::ConfigurationData(DIM, DIM + 1, p, ptrButcherTableau->numStages() + 1)),
    DIM_(DIM),
    numOfVariables_(DIM + 1),
    n_(n),
    p_(p),
    ptrButcherTableau_(ptrButcherTableau),
    T_(0.0),
    cInv_(1.0),
    solutionTimeLevel_(0),
    useMatrixStorage_(useMatrixStorage),
    stiffnessMatrixID_(0),
    massMatrixID_(1)
    {
	}

	/// \brief Create the mesh.
    void createMesh(Base::MeshType meshType = Base::MeshType::RECTANGULAR)
    {
        // Create the domain. In this case the domain is the square [0,1]x[0,1].
        Base::RectangularMeshDescriptor description(DIM_);
    	for(int i=0;i<DIM_;++i)
        {
    		description.bottomLeft_[i]=0;
    		description.topRight_[i]=1;
    		description.numElementsInDIM_[i]=n_;
            description.boundaryConditions_[i]=Base::RectangularMeshDescriptor::PERIODIC;
        }
        
        // Set the number of Element/Face Matrices/Vectors.
        std::size_t numOfElementMatrices;
        std::size_t numOfElementVectors;
        std::size_t numOfFaceMatrices;
        std::size_t numOfFaceVectors;
        if(useMatrixStorage_)
        {
            numOfElementMatrices = 2;   // Mass matrix and element matrx
            numOfElementVectors = 1;    // Source term (not used in this example)
            numOfFaceMatrices = 1;      // Mass matrix
            numOfFaceVectors = 1;       // Boundary conditions (not used in this example)
        }
        else
        {
            numOfElementMatrices = 0;
            numOfElementVectors = 0;
            numOfFaceMatrices = 0;
            numOfFaceVectors = 0;
        }
        
        // Create mesh and set basis functions.
        if(DIM_ == 2)
        {
            if(meshType == Base::MeshType::TRIANGULAR)
            {
                addMesh(description, Base::TRIANGULAR, numOfElementMatrices, numOfElementVectors, numOfFaceMatrices, numOfFaceVectors);
                meshes_[0]->setDefaultBasisFunctionSet(Utilities::createDGBasisFunctionSet2DH1Triangle(p_));
            }
            else if(meshType == Base::MeshType::RECTANGULAR)
            {
                addMesh(description, Base::RECTANGULAR, numOfElementMatrices, numOfElementVectors, numOfFaceMatrices, numOfFaceVectors);
                meshes_[0]->setDefaultBasisFunctionSet(Utilities::createDGBasisFunctionSet2DH1Square(p_));
            }
        }
        else if(DIM_ == 3)
        {
            if(meshType == Base::MeshType::TRIANGULAR)
            {
                addMesh(description, Base::TRIANGULAR, numOfElementMatrices, numOfElementVectors, numOfFaceMatrices, numOfFaceVectors);
                meshes_[0]->setDefaultBasisFunctionSet(Utilities::createDGBasisFunctionSet3DH1Tetrahedron(p_));
            }
            else if(meshType == Base::MeshType::RECTANGULAR)
            {
                addMesh(description, Base::RECTANGULAR, numOfElementMatrices, numOfElementVectors, numOfFaceMatrices, numOfFaceVectors);
                meshes_[0]->setDefaultBasisFunctionSet(Utilities::createDGBasisFunctionSet3DH1Cube(p_));
            }
        }
        
        std::size_t nElements = meshes_[0]->getNumberOfElements();
        std::cout << "Total number of elements: " << nElements << "\n";
    }
    
    /// \brief Set the material parameter.
    /// \param[in] c Material parameter corresponding to the speed with which waves can propagate.
    void setMaterialParameter(const double c)
    {
        cInv_ = 1.0 / c;
    }
    
    /// \brief Compute the integral on a reference element.
    template <typename IntegrandType>
    IntegrandType referenceElementIntegral
    (
        const Base::Element *ptrElement,
        const std::size_t &time,
        std::function<IntegrandType(const Base::Element *, const std::size_t &, const Geometry::PointReference &)> integrandFunction
    )
    {
        const QuadratureRules::GaussQuadratureRule *ptrQdrRule = ptrElement->getGaussQuadratureRule();
        
        std::size_t numOfPoints = ptrQdrRule->nrOfPoints();
        std::size_t iPoint = 0; // Index for the quadrature points.
        Geometry::PointReference pRef(DIM_);
        
        ptrQdrRule->getPoint(iPoint,pRef);
        IntegrandType integral(ptrQdrRule->weight(iPoint) * integrandFunction(ptrElement, time, pRef));
        for(iPoint = 1; iPoint < numOfPoints; iPoint++)
        {
            ptrQdrRule->getPoint(iPoint,pRef);
            integral.axpy(ptrQdrRule->weight(iPoint), integrandFunction(ptrElement, time, pRef));
        }
        return integral;
    }
    
    /// \brief Compute the integral on a reference element.
    template <typename IntegrandType>
    IntegrandType referenceElementIntegral
    (
        const Base::Element *ptrElement,
        const std::size_t &time,
        const LinearAlgebra::NumericalVector &solutionCoefficients,
        std::function<IntegrandType (const Base::Element *, const std::size_t &, const Geometry::PointReference &, const LinearAlgebra::NumericalVector &)> integrandFunction
    )
    {
        const QuadratureRules::GaussQuadratureRule *ptrQdrRule = ptrElement->getGaussQuadratureRule();
        
        std::size_t numOfPoints = ptrQdrRule->nrOfPoints();
        std::size_t iPoint = 0; // Index for the quadrature points.
        Geometry::PointReference pRef(DIM_);
        
        ptrQdrRule->getPoint(iPoint,pRef);
        IntegrandType integral(ptrQdrRule->weight(iPoint) * integrandFunction(ptrElement, time, pRef, solutionCoefficients));
        for(iPoint = 1; iPoint < numOfPoints; iPoint++)
        {
            ptrQdrRule->getPoint(iPoint,pRef);
            integral.axpy(ptrQdrRule->weight(iPoint), integrandFunction(ptrElement, time, pRef, solutionCoefficients));
        }
        return integral;
    }
    
    /// \brief Compute the integral on a reference face that corresponds to an internal face.
    template <typename IntegrandType>
    IntegrandType referenceFaceIntegral
    (
     const Base::Face *ptrFace,
     const std::size_t &time,
     const Base::Side &iSide,
     const Base::Side &jSide,
     std::function<IntegrandType (const Base::Face *, const std::size_t &, const Geometry::PointReference &, const Base::Side &, const Base::Side &)> integrandFunction
    )
    {
        const QuadratureRules::GaussQuadratureRule *ptrQdrRule = ptrFace->getGaussQuadratureRule();
        
        std::size_t numOfPoints = ptrQdrRule->nrOfPoints();
        std::size_t iPoint = 0; // Index for the quadrature points.
        Geometry::PointReference pRef(DIM_);
        
        ptrQdrRule->getPoint(iPoint,pRef);
        IntegrandType integral(ptrQdrRule->weight(iPoint) * integrandFunction(ptrFace, time, pRef, iSide, jSide));
        for(iPoint = 1; iPoint < numOfPoints; iPoint++)
        {
            ptrQdrRule->getPoint(iPoint,pRef);
            integral.axpy(ptrQdrRule->weight(iPoint), integrandFunction(ptrFace, time, pRef, iSide, jSide));
        }
        return integral;
    }
    
    /// \brief Compute the integral on a reference face that corresponds to an internal face.
    template <typename IntegrandType>
    IntegrandType referenceFaceIntegral
    (
     const Base::Face *ptrFace,
     const std::size_t &time,
     const Base::Side &iSide,
     const LinearAlgebra::NumericalVector &solutionCoefficientsLeft,
     const LinearAlgebra::NumericalVector &solutionCoefficientsRight,
     std::function<IntegrandType (const Base::Face *, const std::size_t &, const Geometry::PointReference &, const Base::Side &, const LinearAlgebra::NumericalVector &, const LinearAlgebra::NumericalVector &)> integrandFunction
    )
    {
        const QuadratureRules::GaussQuadratureRule *ptrQdrRule = ptrFace->getGaussQuadratureRule();
        
        std::size_t numOfPoints = ptrQdrRule->nrOfPoints();
        std::size_t iPoint = 0; // Index for the quadrature points.
        Geometry::PointReference pRef(DIM_);
        
        ptrQdrRule->getPoint(iPoint,pRef);
        IntegrandType integral(ptrQdrRule->weight(iPoint) * integrandFunction(ptrFace, time, pRef, iSide, solutionCoefficientsLeft, solutionCoefficientsRight));
        for(iPoint = 1; iPoint < numOfPoints; iPoint++)
        {
            ptrQdrRule->getPoint(iPoint,pRef);
            integral.axpy(ptrQdrRule->weight(iPoint), integrandFunction(ptrFace, time, pRef, iSide, solutionCoefficientsLeft, solutionCoefficientsRight));
        }
        return integral;
    }
    
    /// \brief Get the material parameter c^{-1} at a given physical point.
    double getCInv(const Geometry::PointPhysical &pPhys)
    {
        return cInv_;
    }
    
    /// \brief Compute the source term at a given physical point.
    double getSourceTerm(const double &time, const PointPhysicalT &pPhys)
    {
        return 0.0;
    }
    
    /// \brief Compute the real solution at a given point in space and time.
    LinearAlgebra::NumericalVector getRealSolution(const double &time, const PointPhysicalT &pPhys)
    {
        LinearAlgebra::NumericalVector realSolution(numOfVariables_);
        double c = sqrt(1.0 / cInv_); // Wave velocity.
        
        double xWave = 0;   // The inner product of the physical point and some direction vector.
        for(std::size_t iD = 0; iD < DIM_; iD++) // Index for the dimension.
        {
            xWave += pPhys[iD];
        }
        
        for(std::size_t iV = 0; iV < numOfVariables_; iV++) // iV is the index for the variable.
        {
            realSolution(iV) = 2*M_PI*cos(2*M_PI*(xWave - sqrt(DIM_)*c*time));
            if(iV == 0)
            {
                realSolution(iV) *= -sqrt(DIM_) * c;
            }
            else
            {
                realSolution(iV) /= cInv_;
            }
        }
        
        /*
        LinearAlgebra::NumericalVector realSolution(numOfVariables_);
        for(std::size_t iV = 0; iV < numOfVariables_; iV++) // iV is the index for the variable.
        {
            realSolution(iV) = iV + 1;
        }
         */
        
        return realSolution;
    }
    
    /// \brief Compute the integrand for the mass matrix for the reference element.
    /// \details The integrand for the reference element is the same as the physical element, but scaled with the reference-to-physical element scale, which is the determinant of the jacobian of the reference-to-physical element mapping.
    LinearAlgebra::Matrix integrandMassMatrixOnRefElement
    (
        const Base::Element *ptrElement,
        const std::size_t &time,
        const Geometry::PointReference &pRef
    )
    {
        std::size_t numOfBasisFunctions = ptrElement->getNrOfBasisFunctions();
        LinearAlgebra::Matrix integrand(numOfVariables_ * numOfBasisFunctions,
                                     numOfVariables_ * numOfBasisFunctions );
        Geometry::PointPhysical pPhys(DIM_);
        ptrElement->referenceToPhysical(pRef, pPhys);
        
        std::size_t iVB, jVB;   // indices for both variable and basis function.
        for(std::size_t iV = 0; iV < numOfVariables_; iV++)
        {
            for(std::size_t iB = 0; iB < numOfBasisFunctions; iB++)
            {
                for(std::size_t jB = 0; jB < numOfBasisFunctions; jB++)
                {
                    iVB = ptrElement->convertToSingleIndex(iB, iV);
                    jVB = ptrElement->convertToSingleIndex(jB, iV);
                    integrand(iVB,jVB) = ptrElement->basisFunction(iB, pRef) * ptrElement->basisFunction(jB, pRef);
                    if(iV > 0)
                    {
                        integrand(iVB, jVB) *= getCInv(pPhys);
                    }
                }
            }
        }
        
        // Scale with the reference-to-physical element ratio.
        Geometry::Jacobian jac(DIM_, DIM_);
        ptrElement->calcJacobian(pRef,jac);
        integrand *= jac.determinant();
        
        return integrand;
    }
    
    /// \brief Compute the integrand for the reference element for obtaining the initial solution.
    /// \details The integrand for the initial solution is the exact solution at time 0 multiplied by a test function. The integrand is then scaled by the reference-to-physical element scale, since we compute the integral on a reference element.
    LinearAlgebra::NumericalVector integrandInitialSolutionOnRefElement
    (
        const Base::Element *ptrElement,
        const std::size_t &time0,
        const Geometry::PointReference &pRef
    )
    {
        std::size_t numOfBasisFunctions = ptrElement->getNrOfBasisFunctions();
        
        LinearAlgebra::NumericalVector integrand(numOfVariables_ * numOfBasisFunctions);
        
        Geometry::PointPhysical pPhys(DIM_);
        ptrElement->referenceToPhysical(pRef, pPhys);
        
        LinearAlgebra::NumericalVector initialSolution(getRealSolution(time0, pPhys));
        
        std::size_t iVB;   // Index for both variable and basis function.
        for(std::size_t iV = 0; iV < numOfVariables_; iV++)
        {
            for(std::size_t iB = 0; iB < numOfBasisFunctions; iB++)
            {
                iVB = ptrElement->convertToSingleIndex(iB, iV);
                integrand(iVB) = ptrElement->basisFunction(iB, pRef) * initialSolution(iV);
                if(iV > 0)
                {
                    integrand(iVB) *= getCInv(pPhys);
                }
            }
        }
        
        // Scale with the reference-to-physical element ratio.
        Geometry::Jacobian jac(DIM_, DIM_);
        ptrElement->calcJacobian(pRef,jac);
        integrand *= jac.determinant();
        
        return integrand;
    }
    
    /// \brief Compute the integrand for the stiffness matrix for the reference element.
    /// \details The integrand for the reference element is the same as the physical element, but scaled with the reference-to-physical element scale, which is the determinant of the jacobian of the reference-to-physical element mapping.
    LinearAlgebra::Matrix integrandStiffnessMatrixOnRefElement
    (
        const Base::Element *ptrElement,
        const std::size_t &time,
        const Geometry::PointReference &pRef
    )
    {
        std::size_t numOfBasisFunctions = ptrElement->getNrOfBasisFunctions();
        LinearAlgebra::Matrix integrand(numOfVariables_ * numOfBasisFunctions, numOfVariables_ * numOfBasisFunctions);
        LinearAlgebra::NumericalVector gradientBasisFunction(DIM_);
        double valueTestFunction;
        
        std::size_t iVB, jVB;   // Indices for both variable and basisfunction
        for(std::size_t jB = 0; jB < numOfBasisFunctions; jB++)
        {
            ptrElement->basisFunctionDeriv(jB, pRef, gradientBasisFunction);
            for(std::size_t iB = 0; iB < numOfBasisFunctions; iB++)
            {
                valueTestFunction = ptrElement->basisFunction(iB, pRef);
                
                iVB = ptrElement->convertToSingleIndex(iB, 0);
                for(std::size_t jD = 0; jD < DIM_; jD++)    // Index for the derivatives
                {
                    jVB = ptrElement->convertToSingleIndex(jB, jD + 1);
                    integrand(iVB, jVB) = gradientBasisFunction(jD) * valueTestFunction;
                    
                }
                
                jVB = ptrElement->convertToSingleIndex(jB, 0);
                for(std::size_t iD = 0; iD < DIM_; iD++)    // Index for the derivatives
                {
                    iVB = ptrElement->convertToSingleIndex(iB, iD + 1);
                    integrand(iVB, jVB) = gradientBasisFunction(iD) * valueTestFunction;
                }
            }
        }
        
        // Scale with the reference-to-physical element ratio.
        Geometry::Jacobian jac(DIM_, DIM_);
        ptrElement->calcJacobian(pRef,jac);
        integrand *= jac.determinant();
        
        return integrand;
    }
    
    /// \brief Compute the integrand for the stiffness matrix for the reference face corresponding to an internal face.
    /// \details The integrand for the reference face is the same as the physical face, but scaled with the reference-to-physical face scale. This face scale is absorbed in the normal vector, since it is relatively cheap to compute the normal vector with a length (L2-norm) equal to the reference-to-physical face scale.
    LinearAlgebra::Matrix integrandStiffnessMatrixOnRefFace
    (
     const Base::Face *ptrFace,
     const std::size_t &time,
     const Geometry::PointReference &pRef,
     const Base::Side &iSide,
     const Base::Side &jSide
    )
    {
        std::size_t numOfSolutionBasisFunctions = ptrFace->getPtrElement(jSide)->getNrOfBasisFunctions();
        std::size_t numOfTestBasisFunctions = ptrFace->getPtrElement(iSide)->getNrOfBasisFunctions();
        
        LinearAlgebra::Matrix integrand(numOfTestBasisFunctions * numOfVariables_, numOfSolutionBasisFunctions * numOfVariables_);
        
        double valueTestFunction;
        double valueBasisFunction;
        
        // Compute normal vector, with size of the ref-to-phys face scale, pointing outward of the left element.
        LinearAlgebra::NumericalVector normal(DIM_);
        ptrFace->getNormalVector(pRef, normal);
        if(jSide == Base::Side::RIGHT) {normal *= -1;};
        
        // Compute the integrand
        std::size_t iVB, jVB;    // Indices for both variable and basisfunction.
        for(std::size_t jB = 0; jB < numOfSolutionBasisFunctions; jB++) // iB and jB are indices for the basis function.
        {
            valueBasisFunction = ptrFace->basisFunction(jSide, jB, pRef);
            for(std::size_t iB = 0; iB < numOfTestBasisFunctions; iB++)
            {
                valueTestFunction = ptrFace->basisFunction(iSide, iB, pRef);
                
                iVB = ptrFace->getPtrElement(iSide)->convertToSingleIndex(iB, 0);
                for(std::size_t jD = 0; jD < DIM_; jD++)    // index for the direction
                {
                    jVB = ptrFace->getPtrElement(jSide)->convertToSingleIndex(jB, jD + 1);
                    integrand(iVB, jVB) = -0.5 * normal(jD) * valueBasisFunction * valueTestFunction;
                    
                }
                
                jVB = ptrFace->getPtrElement(jSide)->convertToSingleIndex(jB, 0);
                for(std::size_t iD = 0; iD < DIM_; iD++)    // index for the direction
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
    LinearAlgebra::NumericalVector integrandRightHandSideOnRefElement
    (
        const Base::Element *ptrElement,
        const std::size_t &time,
        const Geometry::PointReference &pRef,
        const LinearAlgebra::NumericalVector &solutionCoefficients
    )
    {
        std::size_t numOfBasisFunctions = ptrElement->getNrOfBasisFunctions();
        
        LinearAlgebra::NumericalVector integrand(numOfVariables_ * numOfBasisFunctions);
        
        LinearAlgebra::NumericalVector gradientBasisFunction(DIM_);
        LinearAlgebra::NumericalVector gradientScalarFunction(DIM_);
        double divergenceVectorFunction = 0;
        
        // Compute the gradient of the scalar function and the divergence of the vector function.
        std::size_t jVB;    // Index for both basis function and variable
        for(std::size_t jD = 0; jD < DIM_; jD++)    // Index for derivatives
        {
            gradientScalarFunction(jD) = 0;
        }
        for(std::size_t jB = 0; jB < numOfBasisFunctions; jB++) // Index for basis functions
        {
            ptrElement->basisFunctionDeriv(jB, pRef, gradientBasisFunction);
            
            for(std::size_t jD = 0; jD < DIM_; jD++)    // Index for derivatives
            {
                jVB = ptrElement->convertToSingleIndex(jB, 0);
                gradientScalarFunction(jD) += gradientBasisFunction(jD) * solutionCoefficients(jVB);
                
                jVB = ptrElement->convertToSingleIndex(jB, jD + 1);
                divergenceVectorFunction += gradientBasisFunction(jD) * solutionCoefficients(jVB);
            }
        }
        
        // Compute integrand on the physical element.
        std::size_t iVB;    // Index for both basis function and variable
        for(std::size_t iB = 0; iB < numOfBasisFunctions; iB++) // Index for basis function
        {
            iVB = ptrElement->convertToSingleIndex(iB, 0);
            integrand(iVB) = ptrElement->basisFunction(iB, pRef) * divergenceVectorFunction;
            
            for(std::size_t iD = 0; iD < DIM_; iD++)    // Index for derivative
            {
                iVB = ptrElement->convertToSingleIndex(iB, iD + 1);
                integrand(iVB) = ptrElement->basisFunction(iB, pRef) * gradientScalarFunction(iD);
            }
        }
        
        // Scale with the reference-to-physical element ratio.
        Geometry::Jacobian jac(DIM_, DIM_);
        ptrElement->calcJacobian(pRef,jac);
        integrand *= jac.determinant();
        
        return integrand;
    }

    /// \brief Compute the integrand for the right hand side for the reference face corresponding to an internal face.
    /// \details The integrand for the reference face is the same as the physical face, but scaled with the reference-to-physical face scale. This face scale is absorbed in the normal vector, since it is relatively cheap to compute the normal vector with a length (L2-norm) equal to the reference-to-physical face scale.
    LinearAlgebra::NumericalVector integrandRightHandSideOnRefFace
    (
        const Base::Face *ptrFace,
        const std::size_t &time,
        const Geometry::PointReference &pRef,
        const Base::Side &iSide,
        const LinearAlgebra::NumericalVector &solutionCoefficientsLeft,
        const LinearAlgebra::NumericalVector &solutionCoefficientsRight
    )
    {
        std::size_t numOfTestBasisFunctions = ptrFace->getPtrElement(iSide)->getNrOfBasisFunctions();
        std::size_t numOfSolutionBasisFunctionsLeft = ptrFace->getPtrElementLeft()->getNrOfBasisFunctions();
        std::size_t numOfSolutionBasisFunctionsRight = ptrFace->getPtrElementRight()->getNrOfBasisFunctions();
        
        LinearAlgebra::NumericalVector integrand(numOfVariables_ * numOfTestBasisFunctions);
        
        LinearAlgebra::NumericalVector numericalSolutionLeft(numOfVariables_);
        LinearAlgebra::NumericalVector numericalSolutionRight(numOfVariables_);
        
        // Compute the numerical solution at the given point at the left and right side.
        std::size_t jVB;   // Index for both variable and basis function.
        for(std::size_t jV = 0; jV < numOfVariables_; jV++)
        {
            numericalSolutionLeft(jV) = 0;
            numericalSolutionRight(jV) = 0;
            for(std::size_t jB = 0; jB < numOfSolutionBasisFunctionsLeft; jB++)
            {
                jVB = ptrFace->getPtrElementLeft()->convertToSingleIndex(jB, jV);
                numericalSolutionLeft(jV) += ptrFace->basisFunction(Base::Side::LEFT, jB, pRef) * solutionCoefficientsLeft(jVB);
            }
            for(std::size_t jB = 0; jB < numOfSolutionBasisFunctionsRight; jB++)
            {
                jVB = ptrFace->getPtrElementRight()->convertToSingleIndex(jB, jV);
                numericalSolutionRight(jV) += ptrFace->basisFunction(Base::Side::RIGHT, jB, pRef) * solutionCoefficientsRight(jVB);
            }
        }
        
        // Compute normal vector, with size of the ref-to-phys face scale, pointing outward of the left element.
        LinearAlgebra::NumericalVector normal(DIM_);
        ptrFace->getNormalVector(pRef, normal);
        
        // Compute the jump of the scalar function and the vector function.
        LinearAlgebra::NumericalVector jumpScalarFunction(DIM_);
        double jumpVectorFunction = 0;
        for(std::size_t jD = 0; jD < DIM_; jD++)
        {
            jumpScalarFunction(jD) = normal(jD) * (numericalSolutionLeft(0) - numericalSolutionRight(0));
            jumpVectorFunction += normal(jD) * (numericalSolutionLeft(jD + 1) - numericalSolutionRight(jD + 1));
        }
        
        // Compute integrand on the reference element.
        std::size_t iVB;    // Index for both variable and basis function.
        for(std::size_t iB = 0; iB < numOfTestBasisFunctions; iB++)
        {
            iVB = ptrFace->getPtrElement(iSide)->convertToSingleIndex(iB, 0);
            integrand(iVB) = - 0.5 * ptrFace->basisFunction(iSide, iB, pRef) * jumpVectorFunction;
            
            for(std::size_t iD = 0; iD < DIM_; iD++)    // Index for direction
            {
                iVB = ptrFace->getPtrElement(iSide)->convertToSingleIndex(iB, iD + 1);
                integrand(iVB) = - 0.5 * ptrFace->basisFunction(iSide, iB, pRef) * jumpScalarFunction(iD);
            }
        }
        return integrand;
    }
    
    /// \brief Compute the integrand for the reference element for computing the energy-norm of the error.
    /// \details The integrand for the reference element is the same as the physical element, but scaled with the reference-to-physical element scale, which is the determinant of the jacobian of the reference-to-physical element mapping.
    LinearAlgebra::NumericalVector integrandErrorOnRefElement
    (
        const Base::Element *ptrElement,
        const std::size_t &time,
        const Geometry::PointReference &pRef,
        const LinearAlgebra::NumericalVector &solutionCoefficients
    )
    {
        std::size_t numOfBasisFunctions = ptrElement->getNrOfBasisFunctions();
        
        LinearAlgebra::NumericalVector integrand(1);
        integrand(0) = 0;
        
        Geometry::PointPhysical pPhys(DIM_);
        ptrElement->referenceToPhysical(pRef, pPhys);
        
        LinearAlgebra::NumericalVector realSolution(getRealSolution(time, pPhys));
        LinearAlgebra::NumericalVector numericalSolution(numOfVariables_);
        
        std::size_t jVB;   // Index for both variable and basis function.
        for(std::size_t jV = 0; jV < numOfVariables_; jV++)
        {
            numericalSolution(jV) = 0;
            for(std::size_t jB = 0; jB < numOfBasisFunctions; jB++)
            {
                jVB = ptrElement->convertToSingleIndex(jB, jV);
                numericalSolution(jV) += ptrElement->basisFunction(jB, pRef) * solutionCoefficients(jVB);
            }
            if(jV > 0)
            {
                integrand(0) += pow(numericalSolution(jV) - realSolution(jV), 2);
            }
            else
            {
                integrand(0) += getCInv(pPhys) * pow(numericalSolution(jV) - realSolution(jV), 2);
            }
        }
        
        // Scale with the reference-to-physical element ratio.
        Geometry::Jacobian jac(DIM_, DIM_);
        ptrElement->calcJacobian(pRef,jac);
        integrand *= jac.determinant();
        
        return integrand;
    }
    
    
    /// \brief Compute and store the mass matrices.
    void createMassMatrices(const double time)
    {
        LinearAlgebra::Matrix massMatrix;
        
        std::size_t numOfBasisFunctions;
        
        // Define the integrand function for the mass matrix (I do not know an easier way to define this).
        std::function<LinearAlgebra::Matrix (const Base::Element *, const std::size_t &, const Geometry::PointReference &)> integrandFunction = [=](const Base::Element * ptrElement, const std::size_t & time, const Geometry::PointReference & pRef) -> LinearAlgebra::Matrix{return this->integrandMassMatrixOnRefElement(ptrElement, time, pRef);};
        
        for(Base::Element *ptrElement : meshes_[0]->getElementsList())
        {
            numOfBasisFunctions = ptrElement->getNrOfBasisFunctions();
            
            massMatrix.resize(numOfVariables_ * numOfBasisFunctions, numOfVariables_ * numOfBasisFunctions);
            massMatrix = referenceElementIntegral(ptrElement, time, integrandFunction);
            
            // std::cout << "Mass matrix: " << massMatrix << "\n";
            ptrElement->setElementMatrix(massMatrix, massMatrixID_);
        }
    }
    
    /// \brief Set the initial numerical solution (w at t=0).
    void setInitialSolution(const double time0)
    {
        // The function createMassMatrices() must have been called before this function when storing mass matrices.
        LinearAlgebra::NumericalVector initialSolution;
        LinearAlgebra::Matrix massMatrix;
        
        std::size_t numOfBasisFunctions;
        
        // Define the integrand function for the the initial solution integral.
        std::function<LinearAlgebra::NumericalVector (const Base::Element *, const std::size_t &, const Geometry::PointReference &)> integrandFunctionInitialSolution = [=](const Base::Element * ptrElement, const std::size_t & time, const Geometry::PointReference & pRef) -> LinearAlgebra::NumericalVector{return this->integrandInitialSolutionOnRefElement(ptrElement, time, pRef);};
        
        // Define the integrand function for the mass matrix.
        std::function<LinearAlgebra::Matrix (const Base::Element *, const std::size_t &, const Geometry::PointReference &)> integrandFunctionMassMatrix = [=](const Base::Element * ptrElement, const std::size_t & time, const Geometry::PointReference & pRef) -> LinearAlgebra::Matrix{return this->integrandMassMatrixOnRefElement(ptrElement, time, pRef);};
        
        for(Base::Element *ptrElement : meshes_[0]->getElementsList())
        {
            numOfBasisFunctions = ptrElement->getNrOfBasisFunctions();
            
            massMatrix.resize(numOfVariables_ * numOfBasisFunctions, numOfVariables_ * numOfBasisFunctions);
            initialSolution.resize(numOfVariables_ * numOfBasisFunctions);
            if(useMatrixStorage_)
            {
                ptrElement->getElementMatrix(massMatrix, massMatrixID_);
            }
            else
            {
                massMatrix = referenceElementIntegral(ptrElement, time0, integrandFunctionMassMatrix);
            }
            initialSolution = referenceElementIntegral(ptrElement, time0, integrandFunctionInitialSolution);
            massMatrix.solve(initialSolution);
            
            ptrElement->setTimeLevelDataVector(solutionTimeLevel_, initialSolution);
        }
    }
    
    /// \brief Compute and store matrices and vectors for computing the right hand side.
    void createMatricesAndVectorsRightHandSide(const double time)
    {
        LinearAlgebra::Matrix stiffnessMatrix;
        
        std::size_t numOfBasisFunctions;
        
        // Define the integrand function for the stiffness matrix for the element.
        std::function<LinearAlgebra::Matrix (const Base::Element *, const std::size_t &, const Geometry::PointReference &)> integrandFunctionElement = [=](const Base::Element * ptrElement, const std::size_t & time, const Geometry::PointReference & pRef) -> LinearAlgebra::Matrix{return this->integrandStiffnessMatrixOnRefElement(ptrElement, time, pRef);};
        
        for(Base::Element *ptrElement : meshes_[0]->getElementsList())
        {
            numOfBasisFunctions = ptrElement->getNrOfBasisFunctions();
            
            stiffnessMatrix.resize(numOfVariables_ * numOfBasisFunctions, numOfVariables_ * numOfBasisFunctions);
            stiffnessMatrix = referenceElementIntegral(ptrElement, time, integrandFunctionElement);
            // std::cout << "Stiffness matrix for element: " << stiffnessMatrix << "\n";
            ptrElement->setElementMatrix(stiffnessMatrix, stiffnessMatrixID_);
        }
        
        Base::FaceMatrix stiffnessFaceMatrix;
        
        std::size_t numOfBasisFunctionsLeft;
        std::size_t numOfBasisFunctionsRight;
        std::size_t numOfSolutionBasisFunctions;
        std::size_t numOfTestBasisFunctions;
        
        std::vector<Base::Side> allSides;
        allSides.push_back(Base::Side::LEFT);
        allSides.push_back(Base::Side::RIGHT);
        
        // Define the integrand function for the stiffness matrix for the face.
        std::function<LinearAlgebra::Matrix (const Base::Face *, const std::size_t &, const Geometry::PointReference &, const Base::Side &, const Base::Side &)> integrandFunctionFace = [=](const Base::Face * ptrFace, const std::size_t & time, const Geometry::PointReference & pRef, const Base::Side iSide, const Base::Side jSide) -> LinearAlgebra::Matrix{return this->integrandStiffnessMatrixOnRefFace(ptrFace, time, pRef, iSide, jSide);};
        
        for(Base::Face *ptrFace : meshes_[0]->getFacesList())
        {
            if(ptrFace->isInternal())
            {
                numOfBasisFunctionsLeft = ptrFace->getPtrElementLeft()->getNrOfBasisFunctions();
                numOfBasisFunctionsRight = ptrFace->getPtrElementRight()->getNrOfBasisFunctions();
                
                stiffnessFaceMatrix.resize(numOfBasisFunctionsLeft * numOfVariables_, numOfBasisFunctionsRight * numOfVariables_);
                
                // std::cout << "Compute stiffness matrix for a face.\n";
                for(Base::Side iSide : allSides)
                {
                    if(iSide == Base::Side::LEFT) {numOfTestBasisFunctions = numOfBasisFunctionsLeft;}
                    else {numOfTestBasisFunctions = numOfBasisFunctionsRight;}
                    
                    for(Base::Side jSide : allSides)
                    {
                        if(jSide == Base::Side::LEFT) {numOfSolutionBasisFunctions = numOfBasisFunctionsLeft;}
                        else {numOfSolutionBasisFunctions = numOfBasisFunctionsRight;}
                        
                        stiffnessMatrix.resize(numOfTestBasisFunctions, numOfSolutionBasisFunctions);
                        stiffnessMatrix = referenceFaceIntegral(ptrFace, time, iSide, jSide, integrandFunctionFace);
                        
                        // std::cout << "Partial stiffness matrix for face: " << stiffnessMatrix << "\n";
                        stiffnessFaceMatrix.setElementMatrix(stiffnessMatrix, iSide, jSide);
                    }
                }
                ptrFace->setFaceMatrix(stiffnessFaceMatrix, stiffnessMatrixID_);
            }
            else
            {
                std::cout << "Warning: there are internal faces.\n";
            }
        }
    }
    
    /// \brief Compute the energy-norm of the error at final time T_.
    /// \param[in] time Time corresponding to the current solution. If no input is given the time will be set to final time T_.
    /// \details The error is defined as error = realSolution - numericalSolution. The energy of the vector (u, s1, s2) is defined as u^2 + c^{-1} * |s|^2. The total energy is the integral over the energy over the entire domain (this energy is conserved in this numerical scheme on a discrete level). The energy-norm is the square root of the total energy.
    double computeEnergyNormError(double time = -1.0)
    {
        if(time < 0) time = T_;
        
        LinearAlgebra::NumericalVector totalEnergy(1);
        totalEnergy(0) = 0;
        
        LinearAlgebra::NumericalVector solutionCoefficients;
        
        std::size_t numOfBasisFunctions;
        
        // Define the integrand function for the error energy.
        std::function<LinearAlgebra::NumericalVector (const Base::Element *, const std::size_t &, const Geometry::PointReference &, const LinearAlgebra::NumericalVector &)> integrandFunction = [=](const Base::Element *ptrElement, const std::size_t &time, const Geometry::PointReference & pRef, const LinearAlgebra::NumericalVector &solutionCoefficients) -> LinearAlgebra::NumericalVector{return this->integrandErrorOnRefElement(ptrElement, time, pRef, solutionCoefficients);};
        
        for(Base::Element *ptrElement : meshes_[0]->getElementsList())
        {
            numOfBasisFunctions = ptrElement->getNrOfBasisFunctions();
            
            solutionCoefficients.resize(numOfVariables_ * numOfBasisFunctions);
            solutionCoefficients = ptrElement->getTimeLevelDataVector(solutionTimeLevel_);
            
            totalEnergy += referenceElementIntegral(ptrElement, time, solutionCoefficients, integrandFunction);
        }
        if(totalEnergy(0) >= 0)
        {
            return sqrt(totalEnergy(0));
        }
        else
        {
            std::cout << "Warning: the total energy of the error is negative.\n";
            return sqrt(-totalEnergy(0));
        }
    }
    
    /// \brief Write output to a tecplot file.
    void writeToTecplotFile(const ElementT *ptrElement, const  PointReferenceT &pRef, std::ostream &out)
    {
        LinearAlgebra::NumericalVector solution(numOfVariables_);
        ptrElement->getSolution(solutionTimeLevel_, pRef, solution);
        
        std::size_t iV = 0; // Index for the variable
        out << solution(iV);
        for(iV = 1; iV < numOfVariables_; iV++)
        {
            out << " " << solution(iV);
        }
    }
    
    /// \brief Compute one time step.
    void computeOneTimeStep(double &time, const double dt)
    {
        std:size_t numOfStages = ptrButcherTableau_->numStages();
        std::size_t numOfBasisFunctions;
        std::size_t numOfBasisFunctionsLeft;
        std::size_t numOfBasisFunctionsRight;
        
        LinearAlgebra::NumericalVector solutionCoefficients;
        LinearAlgebra::NumericalVector solutionCoefficientsLeft;
        LinearAlgebra::NumericalVector solutionCoefficientsRight;
        LinearAlgebra::NumericalVector solutionCoefficientsNew;
        LinearAlgebra::NumericalVector solutionCoefficientsLeftNew;
        LinearAlgebra::NumericalVector solutionCoefficientsRightNew;
        
        LinearAlgebra::Matrix massMatrix;
        LinearAlgebra::Matrix stiffnessMatrix;
        Base::FaceMatrix stiffnessFaceMatrix;
        
        // LinearAlgebra::NumericalVector testValuesIntermediate(numOfStages); // For testing
        
        // Define the integrand function for the mass matrix.
        std::function<LinearAlgebra::Matrix (const Base::Element *, const std::size_t &, const Geometry::PointReference &)> integrandFunctionMassMatrix = [=](const Base::Element * ptrElement, const std::size_t & time, const Geometry::PointReference & pRef) -> LinearAlgebra::Matrix{return this->integrandMassMatrixOnRefElement(ptrElement, time, pRef);};
        
        // Define the integrand function for the right hand side for the reference element.
        std::function<LinearAlgebra::NumericalVector (const Base::Element *, const std::size_t &, const Geometry::PointReference &, const LinearAlgebra::NumericalVector &)> integrandFunctionRefElement = [=](const Base::Element *ptrElement, const std::size_t &time, const Geometry::PointReference & pRef, const LinearAlgebra::NumericalVector &solutionCoefficients) -> LinearAlgebra::NumericalVector{return this->integrandRightHandSideOnRefElement(ptrElement, time, pRef, solutionCoefficients);};
        
        // Define the integrand function for the right hand side for the reference face corresponding to an internal face.
        std::function<LinearAlgebra::NumericalVector (const Base::Face *, const std::size_t &, const Geometry::PointReference &, const Base::Side &, const LinearAlgebra::NumericalVector &, const LinearAlgebra::NumericalVector &)> integrandFunctionRefFaceInternal = [=](const Base::Face *ptrFace, const std::size_t &time, const Geometry::PointReference &pRef, const Base::Side &iSide, const LinearAlgebra::NumericalVector &solutionCoefficientsLeft, const LinearAlgebra::NumericalVector &solutionCoefficientsRight) -> LinearAlgebra::NumericalVector{return this->integrandRightHandSideOnRefFace(ptrFace, time, pRef, iSide, solutionCoefficientsLeft, solutionCoefficientsRight);};
        
        double stageTime; // Time corresponding to the intermediate stages of the Butcher tableau.
        
        // Compute intermediate stages for the time integration.
        for(std::size_t iStage = 0; iStage < numOfStages; iStage++)
        {
            stageTime = time + ptrButcherTableau_->c(iStage) * dt;
            
            // Apply the right hand side corresponding to integration on the elements.
            // std::cout << "Apply the right hand side corresponding to integration on the elements.\n";
            for(Base::Element *ptrElement : meshes_[0]->getElementsList())
            {
                numOfBasisFunctions = ptrElement->getNrOfBasisFunctions();
                
                solutionCoefficients.resize(numOfVariables_ * numOfBasisFunctions);
                solutionCoefficientsNew.resize(numOfVariables_ * numOfBasisFunctions);
                
                solutionCoefficients = ptrElement->getTimeLevelDataVector(solutionTimeLevel_);
                for(std::size_t jStage = 0; jStage < iStage; jStage++)
                {
                    solutionCoefficients += dt * ptrElement->getTimeLevelDataVector(jStage + 1) * ptrButcherTableau_->a(iStage, jStage);
                }
                
                if(useMatrixStorage_)
                {
                    stiffnessMatrix.resize(numOfBasisFunctions * numOfVariables_, numOfBasisFunctions * numOfVariables_);
                    ptrElement->getElementMatrix(stiffnessMatrix, stiffnessMatrixID_);
                    
                    solutionCoefficientsNew = stiffnessMatrix * solutionCoefficients;
                }
                else
                {
                    solutionCoefficientsNew = referenceElementIntegral(ptrElement, stageTime, solutionCoefficients, integrandFunctionRefElement);
                }
                ptrElement->setTimeLevelDataVector(iStage + 1, solutionCoefficientsNew);
            }
            
            /* // Compute test values
            testValuesIntermediate(iStage) = testValue_;
            for(std::size_t jStage = 0; jStage < iStage; jStage++)
            {
                testValuesIntermediate(iStage) += dt * testValuesIntermediate(jStage) * ptrButcherTableau_->a(iStage, jStage);
            }
             */
            
            // Apply the right hand side corresponding to integration on the faces.
            // std::cout << "Apply the right hand side corresponding to integration on the faces.\n";
            for(Base::Face *ptrFace : meshes_[0]->getFacesList())
            {
                if(ptrFace->isInternal())
                {
                    numOfBasisFunctionsLeft = ptrFace->getPtrElementLeft()->getNrOfBasisFunctions();
                    numOfBasisFunctionsRight = ptrFace->getPtrElementRight()->getNrOfBasisFunctions();
                    
                    solutionCoefficientsLeft.resize(numOfVariables_ * numOfBasisFunctionsLeft);
                    solutionCoefficientsLeftNew.resize(numOfVariables_ * numOfBasisFunctionsLeft);
                    solutionCoefficientsRight.resize(numOfVariables_ * numOfBasisFunctionsRight);
                    solutionCoefficientsRightNew.resize(numOfVariables_ * numOfBasisFunctionsRight);
                    
                    solutionCoefficientsLeft = ptrFace->getPtrElementLeft()->getTimeLevelDataVector(solutionTimeLevel_);
                    solutionCoefficientsRight = ptrFace->getPtrElementRight()->getTimeLevelDataVector(solutionTimeLevel_);
                    for(std::size_t jStage = 0; jStage < iStage; jStage++)
                    {
                        solutionCoefficientsLeft += dt * ptrFace->getPtrElementLeft()->getTimeLevelDataVector(jStage + 1) * ptrButcherTableau_->a(iStage, jStage);
                        solutionCoefficientsRight += dt * ptrFace->getPtrElementRight()->getTimeLevelDataVector(jStage + 1) * ptrButcherTableau_->a(iStage, jStage);
                    }
                    
                    solutionCoefficientsLeftNew = ptrFace->getPtrElementLeft()->getTimeLevelDataVector(iStage + 1);
                    solutionCoefficientsRightNew = ptrFace->getPtrElementRight()->getTimeLevelDataVector(iStage + 1);
                    if(useMatrixStorage_)
                    {
                        stiffnessFaceMatrix.resize(numOfBasisFunctionsLeft * numOfVariables_, numOfBasisFunctionsRight * numOfVariables_);
                        stiffnessFaceMatrix = ptrFace->getFaceMatrix(stiffnessMatrixID_);
                        
                        solutionCoefficientsLeftNew += stiffnessFaceMatrix.getElementMatrix(Base::Side::LEFT, Base::Side::LEFT) * solutionCoefficientsLeft + stiffnessFaceMatrix.getElementMatrix(Base::Side::LEFT, Base::Side::RIGHT) * solutionCoefficientsRight;
                        solutionCoefficientsRightNew += stiffnessFaceMatrix.getElementMatrix(Base::Side::RIGHT, Base::Side::LEFT) * solutionCoefficientsLeft + stiffnessFaceMatrix.getElementMatrix(Base::Side::RIGHT, Base::Side::RIGHT) * solutionCoefficientsRight;
                    }
                    else
                    {
                        solutionCoefficientsLeftNew += referenceFaceIntegral(ptrFace, stageTime, Base::Side::LEFT, solutionCoefficientsLeft, solutionCoefficientsRight, integrandFunctionRefFaceInternal);
                        solutionCoefficientsRightNew += referenceFaceIntegral(ptrFace, stageTime, Base::Side::RIGHT, solutionCoefficientsLeft, solutionCoefficientsRight, integrandFunctionRefFaceInternal);
                    }
                    ptrFace->getPtrElementLeft()->setTimeLevelDataVector(iStage + 1, solutionCoefficientsLeftNew);
                    ptrFace->getPtrElementRight()->setTimeLevelDataVector(iStage + 1, solutionCoefficientsRightNew);
                }
                else
                {
                    std::cout << "Warning: there are internal faces.\n";
                }
            }
            
            // Solve the mass matrix systems of equations.
            // std::cout << "Solve the mass matrix systems of equations.\n";
            for(Base::Element *ptrElement : meshes_[0]->getElementsList())
            {
                numOfBasisFunctions = ptrElement->getNrOfBasisFunctions();
                
                solutionCoefficientsNew.resize(numOfBasisFunctions);
                massMatrix.resize(numOfVariables_ * numOfBasisFunctions, numOfVariables_ * numOfBasisFunctions);
                
                if(useMatrixStorage_)
                {
                    ptrElement->getElementMatrix(massMatrix, massMatrixID_);
                }
                else
                {
                    massMatrix = referenceElementIntegral(ptrElement, stageTime, integrandFunctionMassMatrix);
                }
                solutionCoefficientsNew = ptrElement->getTimeLevelDataVector(iStage + 1);
                massMatrix.solve(solutionCoefficientsNew);
                
                ptrElement->setTimeLevelDataVector(iStage + 1, solutionCoefficientsNew);
            }
        }
        
        // Update the solution.
        for(Base::Element *ptrElement : meshes_[0]->getElementsList())
        {
            numOfBasisFunctions = ptrElement->getNrOfBasisFunctions();
            solutionCoefficientsNew.resize(numOfVariables_ * numOfBasisFunctions);
            
            solutionCoefficientsNew = ptrElement->getTimeLevelDataVector(solutionTimeLevel_);
            for(std::size_t jStage = 0; jStage < numOfStages; jStage++)
            {
                solutionCoefficientsNew += dt * ptrElement->getTimeLevelDataVector(jStage + 1) * ptrButcherTableau_->b(jStage);
            }
            ptrElement->setTimeLevelDataVector(solutionTimeLevel_, solutionCoefficientsNew);
        }
        
        /* // Compute test values
        double testValueNew = testValue_;
        for(std::size_t jStage = 0; jStage < numOfStages; jStage++)
        {
            testValueNew += dt * testValuesIntermediate(jStage) * ptrButcherTableau_->b(jStage);
        }
        
        // std::cout << "Test value at intermediate stages: " << testValuesIntermediate << "\n";
        testValue_ = testValueNew;
         */
        
        // Update the time.
        time += dt;
    }

    /// \brief Solve the PDE over the time domain [0,T].
    /// \param[in] T Final time. The solution is solved over the interval [0,T].
    /// \param[in] dt Time step size.
    /// \param[in] numOutputFrames Number of times the solution is written to an output file.
    bool solve(const double T, double dt, const std::size_t numOfOutputFrames = 1)
    {
        // Compute parameters for time integration
        T_ = T;
        std::size_t numOfTimeSteps = T / dt;
        std::size_t numOfTimeStepsForOutput;
        if(numOfOutputFrames > 0)
        {
            // Round off to above such that the number of time steps is a multiple of the number of output frames.
            numOfTimeSteps += (numOfOutputFrames - (numOfTimeSteps % numOfOutputFrames)) % numOfOutputFrames;
            
            // Recompute dt.
            dt = T / numOfTimeSteps;
            
            // Compute the number of timesteps after which to create an output frame.
            numOfTimeStepsForOutput = (std::size_t) numOfTimeSteps / numOfOutputFrames;
        }
        
        // Create output file.
        std::ofstream outFile("output.dat");
        std::string dimensionsToWrite;
        std::string variableString;
        if(DIM_ == 2)
        {
            dimensionsToWrite = "01";
            variableString = "v,s1,s2";
        }
        else if(DIM_ == 3)
        {
            dimensionsToWrite = "012";
            variableString = "v,s1,s2,s3";
        }
        Output::TecplotDiscontinuousSolutionWriter writeFunc(outFile, "test", dimensionsToWrite, variableString);
        
        // Set the initial time.
        double time = 0.0;
        
        // Create and store matrices and vectors for the PDE.
        if(useMatrixStorage_)
        {
            std::cout << "Computing the mass matrices.\n";
            createMassMatrices(time);
            std::cout << "Computing stiffness matrices.\n";
            createMatricesAndVectorsRightHandSide(time);
        }
        
        // Set the initial numerical solution.
        testValue_ = 1.0;
        std::cout << "Computing and interpolating the initial solution.\n";
        setInitialSolution(time);
        writeFunc.write(meshes_[0],"discontinuous solution", false, this, time);
        
        // Solve the system of PDE's.
        std::cout << "Solving the system of PDE's.\n";
        std::cout << "dt: " << dt << ".\n";
        std::cout << "Total number of time steps: " << numOfTimeSteps << ".\n";
        for(std::size_t iT = 1; iT <= numOfTimeSteps; iT++)
        {
            computeOneTimeStep(time, dt);
            
            if(iT % numOfTimeStepsForOutput == 0)
            {
                writeFunc.write(meshes_[0],"discontinuous solution", false, this, time);
                // std::cout << "Test value: " << testValue_ << "\n";
            }
            if(iT % 10 == 0)
            {
                std::cout << iT << " time steps computed.\n";
            }
        }
        
        return true;
    }

private:
    /// Dimension of the domain
    std::size_t DIM_;
    
    /// Number of variables
    const std::size_t numOfVariables_;

    /// Number of elements per cardinal direction
    const std::size_t n_;
 
    /// Polynomial order of the basis functions
    const std::size_t p_;
    
    /// Butcher tableau for time integraion. The integration method is assumed to be explicit.
    const Base::ButcherTableau *ptrButcherTableau_;
    
    /// Final time. The PDE is solved over the time interval [0,T].
    double T_;
    
    /// Material parameter c^{-1}
    double cInv_;
    
    /// Index to indicate where the coefficients for the solution are stored
    std::size_t solutionTimeLevel_;
    
    /// Boolean to indicate if matrices (e.g. mass matrix, stiffness matrix etc) should be stored.
    const bool useMatrixStorage_;
    
    /// Index to indicate where the mass matrix is stored
    const std::size_t massMatrixID_;
    
    /// Index to indicate where the stiffness matrix is stored
    const std::size_t stiffnessMatrixID_;
    
    /// Value for testing and debugging
    double testValue_;
};


auto& n = Base::register_argument<std::size_t>('n', "numElems", "number of elements per dimension", true);
auto& p = Base::register_argument<std::size_t>('p', "order", "polynomial order of the solution", true);

// Parse options seem to do weird things. Flag 'h' also communicates with Petsc. Several other flags are also already defined in the library, because they are defined in class HpgemUISimplified.
auto& numOfOutputFrames = Base::register_argument<std::size_t>('O', "numOfOutputFrames", "Number of frames to output", false, 1);
auto& T = Base::register_argument<double>('T', "finalTime", "end time of the simulation", false, 1.0);
auto& dt = Base::register_argument<double>('d', "timeStepSize", "time step of the simulation", false, 0.01);

int main(int argc, char **argv)
{
    Base::parse_options(argc, argv);
	try
    {
		// Set parameters for the PDE.
        const std::size_t DIM = 2;
        const Base::MeshType meshType = Base::MeshType::RECTANGULAR;
        Base::ButcherTableau *ptrButcherTableau = Base::AllTimeIntegrators::Instance().getRule(4, 4);
        const bool useMatrixStorage = true;
		const double c = 1.0;
        
        // Create problem solver 'test'.
		ExampleMultipleVariableProblem test(DIM, n.getValue(), p.getValue(), ptrButcherTableau, useMatrixStorage);
        
        // Create the mesh
		test.createMesh(meshType);
        
        // Set the material parameter
        test.setMaterialParameter(c);
        
        // Solve the problem over time interval [0,T].
		test.solve(T.getValue(), dt.getValue(), numOfOutputFrames.getValue());
        
        // Compute the error.
		double error = test.computeEnergyNormError();
        std::cout << "Energy norm of the error: " << error << ".\n";
        
		return 0;
	}
	catch(const char* e)
    {
		std::cout<<e;
	}
}
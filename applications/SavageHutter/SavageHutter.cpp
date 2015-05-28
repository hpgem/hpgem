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

#include "SavageHutter.h"
#include <cmath>

using LinearAlgebra::MiddleSizeVector;

/// \param[in] dimension Dimension of the domain
/// \param[in] numberOfVariables Number of variables in the PDE
/// \param[in] polynomialOrder Polynomial order of the basis functions
/// \param[in] useMatrixStorage Boolean to indicate if element and face matrices for the PDE should be stored
/// \param[in] ptrButcherTableau Pointer to a Butcher Tableau used to do the time integration with a Runge-Kutta scheme. By default this is a RK4 scheme.
SavageHutter::SavageHutter
(
 const std::size_t dimension,
 const std::size_t numOfVariables,
 const std::size_t polynomialOrder,
 const Base::ButcherTableau * const ptrButcherTableau,
 const std::size_t numTimeSteps
 ) :
HpgemAPISimplified(dimension, numOfVariables, polynomialOrder, ptrButcherTableau),
DIM_(dimension), numOfVariables_(numOfVariables), numTimeSteps_(numTimeSteps), timeStepCounter(0)
{
    rhsComputer_.numOfVariables_ = numOfVariables;
    rhsComputer_.DIM_ = dimension;
    rhsComputer_.epsilon_ = 1.0;
    rhsComputer_.theta_ = M_PI / 6; //radians
}

Base::RectangularMeshDescriptor<DIM> SavageHutter::createMeshDescription(const std::size_t numOfElementPerDirection)
{
    // Create the domain. In this case the domain is the square [0,1]^DIM and periodic.
    Base::RectangularMeshDescriptor<DIM> description;
    for (std::size_t i = 0; i < DIM_; ++i)
    {
        description.bottomLeft_[i] = 0;
        description.topRight_[i] = 1;
        description.numElementsInDIM_[i] = numOfElementPerDirection;
        description.boundaryConditions_[i] = Base::BoundaryType::SOLID_WALL;
    }    
    return description;
}

/// \brief Compute the initial solution at a given point in space and time.
LinearAlgebra::MiddleSizeVector SavageHutter::getInitialSolution(const PointPhysicalT &pPhys, const double &startTime, const std::size_t orderTimeDerivative)
{
    LinearAlgebra::MiddleSizeVector initialSolution(numOfVariables_);
    if (pPhys[0] < M_PI/8)
    {
        initialSolution(0) = rhsComputer_.getInflowBC()(0) +  0.5*(std::cos(8*pPhys[0]) - 1);
    }
    else
    {
        initialSolution(0) = rhsComputer_.getInflowBC()(0) - 1.0;
    }
    initialSolution(1) = 2.5 * initialSolution(0);
    return initialSolution;
}

/// \details The integrand for the initial solution is the exact solution at time 0 multiplied by a test function. 
/// The integrand is then scaled by the reference-to-physical element scale, 
/// since we compute the integral on a reference element.
LinearAlgebra::MiddleSizeVector SavageHutter::integrandInitialSolutionOnElement
(const Base::Element *ptrElement, const double &startTime, const PointReferenceT &pRef)
{
    const std::size_t numOfBasisFunctions = ptrElement->getNrOfBasisFunctions();
    
    LinearAlgebra::MiddleSizeVector integrand(numOfVariables_ * numOfBasisFunctions);
    
    const PointPhysicalT pPhys = ptrElement->referenceToPhysical(pRef);
    
    const LinearAlgebra::MiddleSizeVector initialSolution(getInitialSolution(pPhys, startTime));
    
    std::size_t iVB; // Index for both variable and basis function.
    for (std::size_t iV = 0; iV < numOfVariables_; iV++)
    {
        for (std::size_t iB = 0; iB < numOfBasisFunctions; iB++)
        {
            iVB = ptrElement->convertToSingleIndex(iB, iV);
            integrand(iVB) = ptrElement->basisFunction(iB, pRef) * initialSolution(iV);
        }
    }
    
    
    return integrand;
}

/*********************Integrate over elements and faces************************/

LinearAlgebra::MiddleSizeVector SavageHutter::integrateInitialSolutionAtElement(Base::Element * ptrElement, const double startTime, const std::size_t orderTimeDerivative)
{
    // Define the integrand function for the the initial solution integral.
    const std::function<LinearAlgebra::MiddleSizeVector(const Base::Element *, const PointReferenceT &)> integrandFunction
    = [=](const Base::Element *elt, const PointReferenceT & pRef) -> LinearAlgebra::MiddleSizeVector
    { 
        return this -> integrandInitialSolutionOnElement(elt, startTime, pRef);
    };
    
    const LinearAlgebra::MiddleSizeVector solution = elementIntegrator_.integrate(ptrElement, integrandFunction, ptrElement->getGaussQuadratureRule());
    return solution;
}

LinearAlgebra::MiddleSizeVector SavageHutter::computeRightHandSideAtElement(Base::Element *ptrElement, LinearAlgebra::MiddleSizeVector &solutionCoefficients, const double time)
{
    // Define the integrand function for the right hand side for the reference element.
    const std::function<LinearAlgebra::MiddleSizeVector(const Base::Element*, const PointReferenceT &)> integrandFunction = [=](const Base::Element* elt, const PointReferenceT & pRef) -> LinearAlgebra::MiddleSizeVector
    {   
        return rhsComputer_.integrandRightHandSideOnElement(elt, time, pRef, solutionCoefficients);
    };
    
    return elementIntegrator_.integrate(ptrElement, integrandFunction, ptrElement->getGaussQuadratureRule());
}

LinearAlgebra::MiddleSizeVector SavageHutter::computeRightHandSideAtFace
(
 Base::Face *ptrFace,
 const Base::Side side,
 LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeft,
 LinearAlgebra::MiddleSizeVector &solutionCoefficientsRight,
 const double time
 )
{
    //Faster for 1D: 
    /*const std::function<LinearAlgebra::NumericalVector(const Geometry::PointReference &)> integrandFunction = [=](const Geometry::PointReference & pRef) -> LinearAlgebra::NumericalVector
    {   
        return rhsComputer_.integrandRightHandSideOnRefFace(ptrFace, side, ptrFace->getNormalVector(pRef), pRef, solutionCoefficientsLeft, solutionCoefficientsRight);
    };
    
    return faceIntegrator_.referenceFaceIntegral(ptrFace->getGaussQuadratureRule(), integrandFunction);*/
    
    //What I want is below. However, getPtrElement(side) does not seem to work in that case.
    const std::function<LinearAlgebra::MiddleSizeVector(const Base::Face *, const LinearAlgebra::SmallVector<DIM> &, const PointReferenceOnFaceT &)> integrandFunction =
    [=](const Base::Face * face, const LinearAlgebra::SmallVector<DIM> normal, const PointReferenceOnFaceT & pRef) -> LinearAlgebra::MiddleSizeVector
    {   
        return rhsComputer_.integrandRightHandSideOnRefFace(face, side, normal, pRef, solutionCoefficientsLeft, solutionCoefficientsRight);
    };
    
    return faceIntegrator_.integrate(ptrFace, integrandFunction, ptrFace->getGaussQuadratureRule());
}

LinearAlgebra::MiddleSizeVector SavageHutter::computeRightHandSideAtFace
        (
         Base::Face *ptrFace,
         LinearAlgebra::MiddleSizeVector &solutionCoefficients,
         const double time
         )
{
    const std::function<LinearAlgebra::MiddleSizeVector(const Base::Face *, const LinearAlgebra::SmallVector<DIM> &, const PointReferenceOnFaceT &)> integrandFunction =
    [=](const Base::Face *face, const LinearAlgebra::SmallVector<DIM> normal, const PointReferenceOnFaceT & pRef) -> LinearAlgebra::MiddleSizeVector
    {   
        return rhsComputer_.integrandRightHandSideOnRefFace(face, normal, pRef, solutionCoefficients);
    };
    
    return faceIntegrator_.integrate(ptrFace, integrandFunction, ptrFace->getGaussQuadratureRule());
}

/******************************Limiting****************************************/
void SavageHutter::computeOneTimeStep(double &time, const double dt)
    {
        std::size_t numOfStages = ptrButcherTableau_->getNumStages();
        
        // Compute intermediate Runge-Kutta stages
        for (std::size_t iStage = 0; iStage < numOfStages; iStage++)
        {
            double stageTime = time + ptrButcherTableau_->getC(iStage) * dt;
            
            std::vector<std::size_t> timeLevelsIn;
            std::vector<double> coefficientsTimeLevels;
            
            timeLevelsIn.push_back(solutionTimeLevel_);
            coefficientsTimeLevels.push_back(1);
            for (std::size_t jStage = 0; jStage < iStage; jStage++)
            {
                timeLevelsIn.push_back(intermediateTimeLevels_[jStage]);
                coefficientsTimeLevels.push_back(dt * ptrButcherTableau_->getA(iStage, jStage));
            }
            
            computeTimeDerivative(timeLevelsIn, coefficientsTimeLevels, intermediateTimeLevels_[iStage], stageTime);
        }
        
        // Update the solution
        for (std::size_t jStage = 0; jStage < numOfStages; jStage++)
        {
            scaleAndAddTimeLevel(solutionTimeLevel_, intermediateTimeLevels_[jStage], dt * ptrButcherTableau_->getB(jStage));
        }
        
        limitSolution();
        
        // Update the time.
        time += dt;
    }

void SavageHutter::limitSolution()
{
    for (const Base::Element *element : meshes_[0]->getElementsList())
    {
        if (useLimitierForElement(element))
        {
            logger(INFO, "Element % needs limiting!", element->getID());
        }
    }
}

///\return Returns true if the element needs limiting and false if it does not need limiting.
///\details In here, the discontinuity detector of Krivodonova et. al (2004) is
///implemented to determine whether or not the given element needs limiting. 
bool SavageHutter::useLimitierForElement(const Base::Element *element)
{
    LinearAlgebra::MiddleSizeVector totalIntegral(numOfVariables_);
    LinearAlgebra::MiddleSizeVector numericalSolution(numOfVariables_);
    
    for (const Base::Face *face : element->getFacesList())
    {
        const PointReferenceT& pRefForInflowTest = face->getReferenceGeometry()->getCenter();
        Base::Side sideOfElement = (face->getPtrElementLeft() == element) ? Base::Side::LEFT : Base::Side::RIGHT;
        LinearAlgebra::MiddleSizeVector normal = face->getNormalVector(pRefForInflowTest);
        if (sideOfElement == Base::Side::LEFT)
        {
            numericalSolution = face->getPtrElement(sideOfElement)->getSolution(0, face->mapRefFaceToRefElemL(pRefForInflowTest));
        }
        else
        {
            numericalSolution = face->getPtrElement(sideOfElement)->getSolution(0, face->mapRefFaceToRefElemR(pRefForInflowTest));
            normal *= -1;
        }
        //determine the direction of the flow across this face
        LinearAlgebra::MiddleSizeVector velocity = computeVelocity(numericalSolution);
        
        //if this is an inflow face, integrate the difference in {h,hu} over the face
        if (velocity * normal < 0)
        {
            logger(DEBUG, "Face % is an inflow face for element %.", face->getID(), element->getID());
            LinearAlgebra::MiddleSizeVector numericalSolutionOther(numOfVariables_);
            if (face->isInternal())
            {
                if (sideOfElement == Base::Side::LEFT)
                {
                    numericalSolutionOther = face->getPtrElement(Base::Side::RIGHT)->getSolution(0, face->mapRefFaceToRefElemR(pRefForInflowTest));
                }
                else
                {
                    numericalSolutionOther = face->getPtrElement(Base::Side::LEFT)->getSolution(0, face->mapRefFaceToRefElemL(pRefForInflowTest));
            }
            }
            else
            {
                logger.assert(normal(0) < 0, "This should be an outflow boundary!");
                numericalSolutionOther = rhsComputer_.getInflowBC();
            }
                
            ///\todo make this multi-dimensional by integrating over the face.
            totalIntegral += numericalSolution - numericalSolutionOther;            
        }            
    }
    
    logger(DEBUG, "Integral over all inflow boundaries of difference for element %: %", element->getID(), totalIntegral);    
    
    //divide the integral by h^{(p+1)/2}, the norm of {u,uh} and the size of the face
    std::size_t p = configData_->polynomialOrder_;
    ///\todo check if this definition of h is reasonable for other geometries than lines
    const PointReferenceT& center = element->getReferenceGeometry()->getCenter();
    double dx = std::pow(2. * std::abs(element->calcJacobian(center).determinant()) , 1./DIM_);
    logger(DEBUG, "grid size: %", dx);
    totalIntegral /= std::pow(dx, (p+1.)/2);
    LinearAlgebra::MiddleSizeVector average = computeNormOfAverageOfSolutionInElement(element);
    for (std::size_t i = 0; i < numOfVariables_; ++i)
    {
        totalIntegral(i) /= average(i);
    }
    ///\todo divide by the size of the faces for DIM > 1
    
    //compare this with 1: if >1, discontinuous, if <1 smooth
    return (Base::L2Norm(totalIntegral) > 1);
}

LinearAlgebra::MiddleSizeVector SavageHutter::computeVelocity(LinearAlgebra::MiddleSizeVector numericalSolution)
{
    return LinearAlgebra::MiddleSizeVector({numericalSolution(1)/numericalSolution(0)});
}

///\todo check if this is indeed what the code says
//this should work as an approximation, but more accurate is integral over the element of the solution divided by the volume of the element
LinearAlgebra::MiddleSizeVector SavageHutter::computeNormOfAverageOfSolutionInElement(const Base::Element* element)
{
    LinearAlgebra::MiddleSizeVector average(2);
    for(std::size_t i = 0; i < element->getGaussQuadratureRule()->nrOfPoints(); ++i)
    {
        const PointReferenceT& quadraturePoint = element->getGaussQuadratureRule()->getPoint(i);
        average += element->getSolution(0, quadraturePoint);
    }
    
    logger(DEBUG, "Average over element %: %", element->getID(), average / element->getGaussQuadratureRule()->nrOfPoints());
    return average / element->getGaussQuadratureRule()->nrOfPoints();
}


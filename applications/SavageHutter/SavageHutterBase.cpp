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

#include "SavageHutterBase.h"
#include "HelperFunctions.h"
#include "MeshMoverContraction.h"


SavageHutterBase::SavageHutterBase(const SHConstructorStruct& inputValues) :
HpgemAPISimplified(inputValues.numOfVariables, inputValues.polyOrder, inputValues.ptrButcherTableau, inputValues.ptrButcherTableau->getNumStages() + 2),
    numOfVariables_(inputValues.numOfVariables), dryLimit_(1e-5), temporaryTimeLevel_(inputValues.ptrButcherTableau->getNumStages() + 1), time_(0)
{
    createMesh(inputValues.numElements, inputValues.meshType);
    if (DIM == 2)
    {
        //Set up the move of the mesh; note that the mesh mover gets deleted in the mesh manipulator
        const MeshMoverContraction* meshMover = new MeshMoverContraction;
        initialiseMeshMover(meshMover, 0);
        meshes_[0]->move();
        
        for (Base::Element *element : meshes_[0]->getElementsList())
        {
            element->getReferenceToPhysicalMap()->reinit();
        }
    }
}

Base::RectangularMeshDescriptor<DIM> SavageHutterBase::createMeshDescription(const std::size_t numOfElementPerDirection)
{
    // Create the domain. In this case the domain is the square [0,1]^DIM and periodic.
    Base::RectangularMeshDescriptor<DIM> description;
    for (std::size_t i = 0; i < DIM; ++i)
    {
        description.bottomLeft_[i] = 0;
        description.topRight_[i] = 1;
        description.numElementsInDIM_[i] = 10;
        description.boundaryConditions_[i] = Base::BoundaryType::SOLID_WALL;
    }
    description.topRight_[0] = 5;
    description.numElementsInDIM_[0] = numOfElementPerDirection;
    return description;
}

/*********************Integrate over elements and faces************************/

LinearAlgebra::MiddleSizeVector SavageHutterBase::computeRightHandSideAtElement(Base::Element *ptrElement, LinearAlgebra::MiddleSizeVector &solutionCoefficients, const double time)
{
    // Define the integrand function for the right hand side for the reference element.
    const std::function < LinearAlgebra::MiddleSizeVector(Base::PhysicalElement<DIM>&) > integrandFunction = [ = ](Base::PhysicalElement<DIM>& element) -> LinearAlgebra::MiddleSizeVector
    {   
        return rhsComputer_->integrandRightHandSideOnElement(element, time, solutionCoefficients);
    };

    logger(DEBUG, "element integral on element %: %", ptrElement->getID(), elementIntegrator_.integrate(ptrElement, integrandFunction, ptrElement->getGaussQuadratureRule()));
    return elementIntegrator_.integrate(ptrElement, integrandFunction, ptrElement->getGaussQuadratureRule());
}

LinearAlgebra::MiddleSizeVector SavageHutterBase::computeRightHandSideAtFace
(
 Base::Face *ptrFace,
 const Base::Side side,
 LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeft,
 LinearAlgebra::MiddleSizeVector &solutionCoefficientsRight,
 const double time
 )
{
    const std::function < LinearAlgebra::MiddleSizeVector(Base::PhysicalFace<DIM>&) > integrandFunction = [ = ](Base::PhysicalFace<DIM>& face) -> LinearAlgebra::MiddleSizeVector
    {   
        return rhsComputer_->integrandRightHandSideOnRefFace(face, side, solutionCoefficientsLeft, solutionCoefficientsRight);
    };

    if (ptrFace->getPtrElementRight()->getID() == 11  || ptrFace->getPtrElementLeft()->getID() == 11)
    {
        logger(DEBUG, "face integral on internal face %: %", ptrFace->getID(), faceIntegrator_.integrate(ptrFace, integrandFunction));
        logger(DEBUG, "elements for face %: %, %", ptrFace->getID(), ptrFace->getPtrElementLeft()->getID(), ptrFace->getPtrElementRight()->getID());
    }
    return faceIntegrator_.integrate(ptrFace, integrandFunction);

}

LinearAlgebra::MiddleSizeVector SavageHutterBase::computeRightHandSideAtFace
(
 Base::Face *ptrFace,
 LinearAlgebra::MiddleSizeVector &solutionCoefficients,
 const double time
 )
{
    const std::function<LinearAlgebra::MiddleSizeVector(Base::PhysicalFace<DIM>&)> integrandFunction =
        [ = ](Base::PhysicalFace<DIM>& face) -> LinearAlgebra::MiddleSizeVector
        {
            return rhsComputer_->integrandRightHandSideOnRefFace(face, solutionCoefficients);
        };
        
        
    if (ptrFace->getPtrElementLeft()->getID() == 11)
    {
        logger(DEBUG, "face integral on boundary face %: %", ptrFace->getID(), faceIntegrator_.integrate(ptrFace, integrandFunction));
    }
    return faceIntegrator_.integrate(ptrFace, integrandFunction, ptrFace->getGaussQuadratureRule());
}

/******************************Limiting****************************************/

void SavageHutterBase::computeOneTimeStep(double &time, const double dt)
{
    std::size_t numOfStages = ptrButcherTableau_->getNumStages();

    // Compute intermediate Runge-Kutta stages
    for (std::size_t iStage = 0; iStage < numOfStages; iStage++)
    {
        double stageTime = time + ptrButcherTableau_->getC(iStage) * dt;

        std::vector<std::size_t> timeLevelsIn;
        std::vector<double> coefficientsTimeLevels;

        timeLevelsIn.push_back(solutionVectorId_);
        coefficientsTimeLevels.push_back(1);
        for (std::size_t jStage = 0; jStage < iStage; jStage++)
        {
            timeLevelsIn.push_back(auxiliaryVectorIds_[jStage]);
            coefficientsTimeLevels.push_back(dt * ptrButcherTableau_->getA(iStage, jStage));
        }

        computeTimeDerivative(timeLevelsIn, coefficientsTimeLevels, auxiliaryVectorIds_[iStage], stageTime);
    }

    // Update the solution
    for (std::size_t jStage = 0; jStage < numOfStages; jStage++)
    {
        scaleAndAddVector(solutionVectorId_, auxiliaryVectorIds_[jStage], dt * ptrButcherTableau_->getB(jStage));
    }    

    limitSolutionOuterLoop();

    // Update the time.
    time += dt;
    time_ = time;
    setInflowBC(time);
    logger(DEBUG, "time: %",time_);
}

/// \details Make sure timeLevelResult is different from the timeLevelsIn.
    void SavageHutterBase::computeRightHandSide(const std::vector<std::size_t> timeLevelsIn, const std::vector<double> coefficientsTimeLevels, const std::size_t timeLevelResult, const double time)
    {
        // Apply the right hand side corresponding to integration on the elements.
        for (Base::Element *ptrElement : meshes_[0]->getElementsList())
        {
            LinearAlgebra::MiddleSizeVector solutionCoefficients = getLinearCombinationOfVectors(ptrElement, timeLevelsIn, coefficientsTimeLevels);
            LinearAlgebra::MiddleSizeVector &solutionCoefficientsNew = ptrElement->getTimeIntegrationVector(timeLevelResult);
            logger(DEBUG, "Before: %", solutionCoefficients);
            heightLimiter_->limit(ptrElement, solutionCoefficients);
            ptrElement->setTimeIntegrationVector(temporaryTimeLevel_, solutionCoefficients);
            logger(DEBUG, "After: %", solutionCoefficients);
            
            solutionCoefficientsNew = computeRightHandSideAtElement(ptrElement,  solutionCoefficients, time);
        }
        
        // Apply the right hand side corresponding to integration on the faces.
        for (Base::Face *ptrFace : meshes_[0]->getFacesList())
        {
            if(ptrFace->isInternal())
            {
                LinearAlgebra::MiddleSizeVector solutionCoefficientsLeft = ptrFace->getPtrElementLeft()->getTimeIntegrationVector(temporaryTimeLevel_);
                LinearAlgebra::MiddleSizeVector solutionCoefficientsRight = ptrFace->getPtrElementRight()->getTimeIntegrationVector(temporaryTimeLevel_);
                LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeftNew(ptrFace->getPtrElementLeft()->getTimeIntegrationVector(timeLevelResult));
                LinearAlgebra::MiddleSizeVector &solutionCoefficientsRightNew(ptrFace->getPtrElementRight()->getTimeIntegrationVector(timeLevelResult));
                
                solutionCoefficientsLeftNew += computeRightHandSideAtFace(ptrFace, Base::Side::LEFT, solutionCoefficientsLeft, solutionCoefficientsRight, time);
                solutionCoefficientsRightNew += computeRightHandSideAtFace(ptrFace, Base::Side::RIGHT, solutionCoefficientsLeft, solutionCoefficientsRight, time);
            }
            else
            {
                LinearAlgebra::MiddleSizeVector solutionCoefficients = ptrFace->getPtrElementLeft()->getTimeIntegrationVector(temporaryTimeLevel_);
                LinearAlgebra::MiddleSizeVector &solutionCoefficientsNew(ptrFace->getPtrElementLeft()->getTimeIntegrationVector(timeLevelResult));
                
                    solutionCoefficientsNew += computeRightHandSideAtFace(ptrFace, solutionCoefficients, time);
                }
            }
        
        synchronize(timeLevelResult);
    }

void SavageHutterBase::limitSolutionOuterLoop()
{
    
    for (Base::Element *element : meshes_[0]->getElementsList())
    {
        //don't use the slope limiter if the water height is adapted with the non-negativity limiter
        const double minimum = getMinimumHeight(element);
        if (minimum < dryLimit_)
        {
            LinearAlgebra::MiddleSizeVector solutionCoefficients = element->getTimeIntegrationVector(0);
            heightLimiter_->limit(element, solutionCoefficients);
            element->setTimeIntegrationVector(0, solutionCoefficients);
        }
        else
        {
            //only limit the slope when there is a slope, so not when the solution exists of piecewise constants.
            if (element->getNumberOfBasisFunctions() > 1)
            {
                slopeLimiter_->limitSlope(element);
            }
        }
    }
}

///\details Compute the minimum height by checking the vertices of the element and the Gauss quadrature points in the element.
///While this does not guarantee to give the best result, it gives a good estimate.
double SavageHutterBase::getMinimumHeight(const Base::Element* element)
{
    const PointReferenceT &pRefL = element->getReferenceGeometry()->getReferenceNodeCoordinate(0);
    const LinearAlgebra::MiddleSizeVector &solutionCoefficients = element->getTimeIntegrationVector(0);
    const std::size_t numOfVariables = element->getNumberOfUnknowns();
    const double solutionLeft = Helpers::getSolution<DIM>(element, solutionCoefficients, pRefL, numOfVariables)(0);    
    double minimum = solutionLeft;
    for (std::size_t iPoint = 1; iPoint < element->getReferenceGeometry()->getNumberOfNodes(); ++iPoint)
    {
        const PointReferenceT &pRef = element->getReferenceGeometry()->getReferenceNodeCoordinate(iPoint);
        minimum = std::min(minimum, Helpers::getSolution<DIM>(element, solutionCoefficients, pRef, numOfVariables)(0));
    }
    for (std::size_t p = 0; p < element->getGaussQuadratureRule()->nrOfPoints(); ++p)
    {
        const PointReferenceT& pRef = element->getGaussQuadratureRule()->getPoint(p);
        minimum = std::min(minimum, Helpers::getSolution<DIM>(element, solutionCoefficients, pRef, numOfVariables)(0));
    }
    logger(DEBUG, "Minimum in element %: %", element->getID(), minimum);
    return minimum;
}

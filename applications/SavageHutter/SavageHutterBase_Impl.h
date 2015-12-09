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


#include "HelperFunctions.h"
#include "MeshMoverContraction.h"
#include "Base/CommandLineOptions.h"
#include "Base/MpiContainer.h"
#include "SavageHutterBase.h"

template<std::size_t DIM>
SavageHutterBase<DIM>::SavageHutterBase(std::size_t numberOfVariables, std::size_t polyOrder) :
    Base::HpgemAPISimplified<DIM>(numberOfVariables, polyOrder, TimeIntegration::AllTimeIntegrators::Instance().getRule(1, 1), 0, true),
    numberOfVariables_(numberOfVariables),
    dryLimit_(1e-5),
    time_(0)
{
    
}

///\details Show the number of time steps that have been computed on the console.
template <std::size_t DIM>
void SavageHutterBase<DIM>::showProgress(const double time, const std::size_t timeStepID)
{
    if (timeStepID % 1000 == 0)
    {
        logger(INFO, "% time steps computed.", timeStepID);
    }
}

/*********************Integrate over elements and faces************************/
template <std::size_t DIM>
LinearAlgebra::MiddleSizeVector SavageHutterBase<DIM>::computeRightHandSideAtElement(Base::Element *ptrElement, LinearAlgebra::MiddleSizeVector &solutionCoefficients, const double time)
{
    // Define the integrand function for the right hand side for the reference element.
    const std::function < LinearAlgebra::MiddleSizeVector(Base::PhysicalElement<DIM>&) > integrandFunction = [ = ](Base::PhysicalElement<DIM>& element) -> LinearAlgebra::MiddleSizeVector
    {   
        return integrandRightHandSideOnElement(element, time, solutionCoefficients);
    };

    logger(DEBUG, "element integral on element %: %", ptrElement->getID(), this->elementIntegrator_.integrate(ptrElement, integrandFunction, ptrElement->getGaussQuadratureRule()));
    return this->elementIntegrator_.integrate(ptrElement, integrandFunction, ptrElement->getGaussQuadratureRule());
}

template <std::size_t DIM>
LinearAlgebra::MiddleSizeVector SavageHutterBase<DIM>::computeRightHandSideAtFace
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
        return integrandRightHandSideOnRefFace(face, side, solutionCoefficientsLeft, solutionCoefficientsRight);
    };

    return this->faceIntegrator_.integrate(ptrFace, integrandFunction);

}

template <std::size_t DIM>
LinearAlgebra::MiddleSizeVector SavageHutterBase<DIM>::computeRightHandSideAtFace
(
 Base::Face *ptrFace,
 LinearAlgebra::MiddleSizeVector &solutionCoefficients,
 const double time
 )
{
    const std::function<LinearAlgebra::MiddleSizeVector(Base::PhysicalFace<DIM>&)> integrandFunction =
        [ = ](Base::PhysicalFace<DIM>& face) -> LinearAlgebra::MiddleSizeVector
        {
            return integrandRightHandSideOnRefFace(face, solutionCoefficients, time);
        };
        
    return this->faceIntegrator_.integrate(ptrFace, integrandFunction, ptrFace->getGaussQuadratureRule());
}

template <std::size_t DIM>
std::pair<LinearAlgebra::MiddleSizeVector, LinearAlgebra::MiddleSizeVector> SavageHutterBase<DIM>::computeBothRightHandSidesAtFace
(
        Base::Face *ptrFace,
        LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeft,
        LinearAlgebra::MiddleSizeVector &solutionCoefficientsRight,
        const double time
        )
{
    // Define the integrand function for the right hand side for the face.
    std::function < std::pair<LinearAlgebra::MiddleSizeVector, LinearAlgebra::MiddleSizeVector>(Base::PhysicalFace<DIM>&) > integrandFunction = [ = ](Base::PhysicalFace<DIM>& face) -> std::pair<LinearAlgebra::MiddleSizeVector, LinearAlgebra::MiddleSizeVector>
    {
        return this->integrandsAtFace(face, time, solutionCoefficientsLeft, solutionCoefficientsRight);
    };
    return this->faceIntegrator_.integratePair(ptrFace, integrandFunction);
}



/******************************Limiting****************************************/
template <std::size_t DIM>
void SavageHutterBase<DIM>::computeOneTimeStep(double &time, const double dt)
{
    std::size_t numberOfStages = this->ptrButcherTableau_->getNumStages();

    // Compute intermediate Runge-Kutta stages
    ///Currently, the limiting is only possible with the explicit Euler method.
    ///\todo think of something for height limiter with higher order RK methods.
    for (std::size_t iStage = 0; iStage < numberOfStages; iStage++)
    {
        double stageTime = time + this->ptrButcherTableau_->getC(iStage) * dt;

        std::vector<std::size_t> timeLevelsIn;
        std::vector<double> coefficientsTimeLevels;

        timeLevelsIn.push_back(this->solutionVectorId_);
        coefficientsTimeLevels.push_back(1);
        for (std::size_t jStage = 0; jStage < iStage; jStage++)
        {
            timeLevelsIn.push_back(this->auxiliaryVectorIds_[jStage]);
            coefficientsTimeLevels.push_back(dt * this->ptrButcherTableau_->getA(iStage, jStage));
        }

        this->computeTimeDerivative(timeLevelsIn, coefficientsTimeLevels, this->auxiliaryVectorIds_[iStage], stageTime);
    }

    // Update the solution
    for (std::size_t jStage = 0; jStage < numberOfStages; jStage++)
    {
        this->scaleAndAddVector(this->solutionVectorId_, this->auxiliaryVectorIds_[jStage], dt * this->ptrButcherTableau_->getB(jStage));
    }

    tasksAfterTimeStep();

    // Update the time.
    time += dt;
    this->time_ = time;
    setInflowBC(time);
}

template <std::size_t DIM>
void SavageHutterBase<DIM>::limitSolutionOuterLoop()
{

    this->synchronize(0);
    for (Base::Element *element : this->meshes_[0]->getElementsList())
    {
        //don't use the slope limiter if the water height is adapted with the non-negativity limiter
        const double minimum = getMinimumHeight(element);
        logger(DEBUG, "minimum: %", minimum);
        if (minimum < dryLimit_)
        {
            logger(DEBUG, "I should limit the height now!");
            LinearAlgebra::MiddleSizeVector &solutionCoefficients = element->getTimeIntegrationVector(0);
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
    this->synchronize(0);
}

///\details Compute the minimum height by checking the vertices of the element and the Gauss quadrature points in the element.
///While this does not guarantee to give the best result, it gives a good estimate.
template <std::size_t DIM>
const double SavageHutterBase<DIM>::getMinimumHeight(const Base::Element* element)
{
    const Geometry::PointReference<DIM> &pRefL = element->getReferenceGeometry()->getReferenceNodeCoordinate(0);
    const LinearAlgebra::MiddleSizeVector &solutionCoefficients = element->getTimeIntegrationVector(0);
    const double solutionLeft = element->getSolution(0, pRefL)(0);    
    double minimum = solutionLeft;
    for (std::size_t iPoint = 1; iPoint < element->getReferenceGeometry()->getNumberOfNodes(); ++iPoint)
    {
        const Geometry::PointReference<DIM> &pRef = element->getReferenceGeometry()->getReferenceNodeCoordinate(iPoint);
        minimum = std::min(minimum, element->getSolution(0, pRef)(0));
    }
    for (std::size_t p = 0; p < element->getGaussQuadratureRule()->getNumberOfPoints(); ++p)
    {
        const Geometry::PointReference<DIM>& pRef = element->getGaussQuadratureRule()->getPoint(p);
        minimum = std::min(minimum, element->getSolution(0, pRef)(0));
    }
    return minimum;
}


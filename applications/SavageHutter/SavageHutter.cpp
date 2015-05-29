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
#include "Geometry/PhysicalGeometry.h"
#include <cmath>

using LinearAlgebra::NumericalVector;

/// \param[in] dimension Dimension of the domain
/// \param[in] numberOfVariables Number of variables in the PDE
/// \param[in] polynomialOrder Polynomial order of the basis functions
/// \param[in] useMatrixStorage Boolean to indicate if element and face matrices for the PDE should be stored
/// \param[in] ptrButcherTableau Pointer to a Butcher Tableau used to do the time integration with a Runge-Kutta scheme. By default this is a RK4 scheme.
/// TODO: construct grid here and assign a struct LimiterData to each element?
SavageHutter::SavageHutter
(
 const std::size_t dimension,
 const std::size_t numOfVariables,
 const std::size_t polynomialOrder,
 const Base::ButcherTableau * const ptrButcherTableau) :
HpgemAPISimplified(dimension, numOfVariables, polynomialOrder, ptrButcherTableau),
DIM_(dimension), numOfVariables_(numOfVariables)
{
    rhsComputer_.numOfVariables_ = numOfVariables;
    rhsComputer_.DIM_ = dimension;
    rhsComputer_.epsilon_ = 1.0;
    rhsComputer_.theta_ = 0; //M_PI / 6; //radians
}

SavageHutter::SavageHutter(const SHConstructorStruct& inputValues) :
HpgemAPISimplified(inputValues.dimension, inputValues.numOfVariables, inputValues.polyOrder, inputValues.ptrButcherTableau),
DIM_(inputValues.dimension), numOfVariables_(inputValues.numOfVariables)
{
    rhsComputer_.numOfVariables_ = inputValues.numOfVariables;
    rhsComputer_.DIM_ = inputValues.dimension;
    rhsComputer_.epsilon_ = 1.0;
    rhsComputer_.theta_ = 0; //M_PI / 6; //radians
    createMesh(inputValues.numElements, inputValues.meshType);
    for (Base::Element* element : meshes_[0]->getElementsList())
    {
        element->setUserData(new LimiterData());
    }
}

Base::RectangularMeshDescriptor SavageHutter::createMeshDescription(const std::size_t numOfElementPerDirection)
{
    // Create the domain. In this case the domain is the square [0,1]^DIM and periodic.
    Base::RectangularMeshDescriptor description(DIM_);
    for (std::size_t i = 0; i < DIM_; ++i)
    {
        description.bottomLeft_[i] = 0;
        description.topRight_[i] = 1;
        description.numElementsInDIM_[i] = numOfElementPerDirection;
        description.boundaryConditions_[i] = Base::BoundaryType::PERIODIC;
    }
    return description;
}

/// \brief Compute the initial solution at a given point in space and time.
LinearAlgebra::NumericalVector SavageHutter::getInitialSolution(const PointPhysicalT &pPhys, const double &startTime, const std::size_t orderTimeDerivative)
{
    LinearAlgebra::NumericalVector initialSolution(numOfVariables_);
    if (pPhys[0] > 0.4 && pPhys[0] < 0.6)
    {
        initialSolution(0) = 1;
    }
    else
    {
        initialSolution(0) = 1;
    }
    initialSolution(1) = 0;
    return initialSolution;
}

/// \details The integrand for the initial solution is the exact solution at time 0 multiplied by a test function. 
/// The integrand is then scaled by the reference-to-physical element scale, 
/// since we compute the integral on a reference element.
LinearAlgebra::NumericalVector SavageHutter::integrandInitialSolutionOnElement
(const Base::Element *ptrElement, const double &startTime, const Geometry::PointReference &pRef)
{
    const std::size_t numOfBasisFunctions = ptrElement->getNrOfBasisFunctions();

    LinearAlgebra::NumericalVector integrand(numOfVariables_ * numOfBasisFunctions);

    const Geometry::PointPhysical pPhys = ptrElement->referenceToPhysical(pRef);

    const LinearAlgebra::NumericalVector initialSolution(getInitialSolution(pPhys, startTime));

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

LinearAlgebra::NumericalVector SavageHutter::integrateInitialSolutionAtElement(Base::Element * ptrElement, const double startTime, const std::size_t orderTimeDerivative)
{
    // Define the integrand function for the the initial solution integral.
    const std::function < LinearAlgebra::NumericalVector(const Base::Element *, const Geometry::PointReference &) > integrandFunction
        = [ = ](const Base::Element *elt, const Geometry::PointReference & pRef) -> LinearAlgebra::NumericalVector
        {
            return this -> integrandInitialSolutionOnElement(elt, startTime, pRef);
        };

    const LinearAlgebra::NumericalVector solution = elementIntegrator_.integrate(ptrElement, integrandFunction, ptrElement->getGaussQuadratureRule());
    return solution;
}

LinearAlgebra::NumericalVector SavageHutter::computeRightHandSideAtElement(Base::Element *ptrElement, LinearAlgebra::NumericalVector &solutionCoefficients, const double time)
{
    // Define the integrand function for the right hand side for the reference element.
    const std::function < LinearAlgebra::NumericalVector(const Base::Element*, const Geometry::PointReference &) > integrandFunction = [ = ](const Base::Element* elt, const Geometry::PointReference & pRef) -> LinearAlgebra::NumericalVector
    {
        return rhsComputer_.integrandRightHandSideOnElement(elt, time, pRef, solutionCoefficients);
    };

    return elementIntegrator_.integrate(ptrElement, integrandFunction, ptrElement->getGaussQuadratureRule());
}
LinearAlgebra::NumericalVector SavageHutter::computeRightHandSideAtFace
(
 Base::Face *ptrFace,
 const Base::Side side,
 LinearAlgebra::NumericalVector &solutionCoefficientsLeft,
 LinearAlgebra::NumericalVector &solutionCoefficientsRight,
 const double time
 )
{
    //Faster for 1D: 
    const std::function < LinearAlgebra::NumericalVector(const Geometry::PointReference &) > integrandFunction = [ = ](const Geometry::PointReference & pRef) -> LinearAlgebra::NumericalVector
    {
        return rhsComputer_.integrandRightHandSideOnRefFace(ptrFace, side, ptrFace->getNormalVector(pRef), pRef, solutionCoefficientsLeft, solutionCoefficientsRight);
    };

    return faceIntegrator_.referenceFaceIntegral(ptrFace->getGaussQuadratureRule(), integrandFunction);

    //Desirable syntax
    /*const std::function<LinearAlgebra::NumericalVector(const Base::Face *, const LinearAlgebra::NumericalVector &, const Geometry::PointReference &)> integrandFunction = 
    [=](const Base::Face * face, const LinearAlgebra::NumericalVector normal, const Geometry::PointReference & pRef) -> LinearAlgebra::NumericalVector
    {   
        return rhsComputer_.integrandRightHandSideOnRefFace(face, side, normal, pRef, solutionCoefficientsLeft, solutionCoefficientsRight);
    };
    
    return faceIntegrator_.integrate(ptrFace, integrandFunction, ptrFace->getGaussQuadratureRule());*/
}

LinearAlgebra::NumericalVector SavageHutter::computeRightHandSideAtFace
(
 Base::Face *ptrFace,
 LinearAlgebra::NumericalVector &solutionCoefficients,
 const double time
 )
{
    const std::function < LinearAlgebra::NumericalVector(const Base::Face *, const LinearAlgebra::NumericalVector &, const Geometry::PointReference &) > integrandFunction =
        [ = ](const Base::Face *face, const LinearAlgebra::NumericalVector normal, const Geometry::PointReference & pRef) -> LinearAlgebra::NumericalVector
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
    for (Base::Element *element : meshes_[0]->getElementsList())
    {
        useLimitierForElement(element);
    }
}

///\return Returns true if the element needs limiting and false if it does not need limiting.
///\details In here, the discontinuity detector of Krivodonova et. al (2004) is
///implemented to determine whether or not the given element needs limiting. 
void SavageHutter::useLimitierForElement(Base::Element *element)
{
    LinearAlgebra::NumericalVector totalIntegral(numOfVariables_);
    LinearAlgebra::NumericalVector numericalSolution(numOfVariables_);

    //For every face of this element, check the size of the jump
    for (const Base::Face *face : element->getFacesList())
    {
        //Compute the numerical solution and the (outward) normal of the face
        const Geometry::PointReference& pRefForInflowTest = face->getReferenceGeometry()->getCenter();
        Base::Side sideOfElement = (face->getPtrElementLeft() == element) ? Base::Side::LEFT : Base::Side::RIGHT;
        LinearAlgebra::NumericalVector normal = face->getNormalVector(pRefForInflowTest);
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
        LinearAlgebra::NumericalVector velocity = computeVelocity(numericalSolution);

        //if this is an inflow face, integrate the difference in {h,hu} over the face
        if (velocity * normal < -1e-14)
        {
            logger(DEBUG, "Face % is an inflow face for element %.", face->getID(), element->getID());
            LinearAlgebra::NumericalVector numericalSolutionOther(numOfVariables_);
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
                //logger.assert(normal(0) < 0, "This should be an outflow boundary!");
                numericalSolutionOther = rhsComputer_.getInflowBC();
            }

            ///\todo make this multi-dimensional by integrating over the face.
            totalIntegral += numericalSolution - numericalSolutionOther;
        }
    }

    logger(DEBUG, "Integral over all inflow boundaries of difference for element %: %", element->getID(), totalIntegral);

    //divide the integral by dx^{(p+1)/2}, the norm of {u,uh} and the size of the face
    std::size_t p = configData_->polynomialOrder_;
    ///\todo check if this definition of dx is reasonable for other geometries than lines
    const double dx = std::pow(2. * std::abs(element->calcJacobian(element->getReferenceGeometry()->getCenter()).determinant()), 1. / DIM_);
    logger(DEBUG, "grid size: %", dx);
    totalIntegral /= std::pow(dx, (p + 1.) / 2);

    LinearAlgebra::NumericalVector average = computeNormOfAverageOfSolutionInElement(element);
    for (std::size_t i = 0; i < numOfVariables_; ++i)
    {
        if (average(i) > 1e-10)
            totalIntegral(i) /= average(i);
    }
    ///\todo divide by the size of the faces for DIM > 1

    //compare this with 1: if >1, discontinuous, if <1 smooth
    LimiterData* lData = static_cast<LimiterData*> (element->getUserData());
    for (std::size_t i = 0; i < numOfVariables_; ++i)
    {
        if (totalIntegral(i) > 1)
        {
            logger(DEBUG, "Element % with variable % will be limited.", element->getID(), i);
            lData->isLimited[i] = false;
            //if (element->getFace(0)->isInternal() && element->getFace(1)->isInternal())
            limitWithMinMod(element, i);
            /*else // for now, just take the mean of the first element if it needs limiting.
            {
                lData->valLeft[i] = element->getSolution(0,element->getReferenceGeometry()->getCenter())(i);
                lData->valRight[i] = element->getSolution(0,element->getReferenceGeometry()->getCenter())(i);
            }*/
        }
        else
        {
            lData->isLimited[i] = false;
        }
    }
}

LinearAlgebra::NumericalVector SavageHutter::computeVelocity(LinearAlgebra::NumericalVector numericalSolution)
{
    return LinearAlgebra::NumericalVector({numericalSolution(1) / numericalSolution(0)});
}

///\todo check if this is indeed what the paper says
LinearAlgebra::NumericalVector SavageHutter::computeNormOfAverageOfSolutionInElement(const Base::Element* element)
{
    LinearAlgebra::NumericalVector average(2);
    for (std::size_t i = 0; i < element->getGaussQuadratureRule()->nrOfPoints(); ++i)
    {
        average += element->getSolution(0, element->getGaussQuadratureRule()->getPoint(i));
    }

    logger(DEBUG, "Average over element %: %", element->getID(), average / element->getGaussQuadratureRule()->nrOfPoints());
    return average / element->getGaussQuadratureRule()->nrOfPoints();
}

/*void SavageHutter::limitWithMinMod(Base::Element* element, const std::size_t iVar)
{
    LimiterData* ld = static_cast<LimiterData*> (element->getUserData());
    logger.assert(ld->isLimited[iVar] == true, "Called limiter on variable that should not be limited.");

    const PointReferenceT &pRef = element->getReferenceGeometry()->getCenter();
    const PointReferenceT &pRefFace = element->getFace(0)->getReferenceGeometry()->getCenter();
    const PointReferenceT &pRefL = element->getFace(1)->mapRefFaceToRefElemR(pRefFace); //TODO: do this with getNodeCoordinates()
    const PointReferenceT &pRefR = element->getFace(1)->mapRefFaceToRefElemL(pRefFace);
    logger.assert(0 == pRef[0] && -1 == pRefL[0] && 1 == pRefR[0], "Coordinates to evaluate basis functions are wrong");

    const double u0 = element->getSolution(0, pRef)(iVar);
    const double uPlus = element->getSolution(0, pRefR)(iVar) - u0;
    const double uMinus = u0 - element->getSolution(0, pRefL)(iVar);
    if (element->getNrOfBasisFunctions() < 2) //constants or linear polynomials
    {
        logger.assert(std::abs(uPlus - uMinus) < 1e-10, "Linear polynomial is not straight?");
    }
    //this does not work for first boundary, probably also not for triangular mesh.
    //maybe write getNeighbour(face)?
    const Base::Element * const elemL = element->getFace(0)->getPtrElementLeft();
    const Base::Element * const elemR = element->getFace(1)->getPtrElementRight();

    //check if left is indeed of the left side of right
    logger.assert((elemL->getPhysicalGeometry()->getLocalNodeCoordinates(0)) < (elemR->getPhysicalGeometry()->getLocalNodeCoordinates(0)), "elements left/right in wrong order");

    const double uElemR = elemR->getSolution(0, pRef)(iVar);
    const double uElemL = elemL->getSolution(0, pRef)(iVar);

    logger(DEBUG, "u0: %, uPlus: %, : %, : %", u0, uPlus, uElemR - u0, u0 - uElemL);
    if ((sign(uPlus) == sign(uElemR - u0)) && sign(uPlus) == sign(u0 - uElemL))
    {
        double minAbs = std::min(std::abs(uPlus), std::abs(uElemR - u0));
        minAbs = std::min(minAbs, std::abs(u0 - uElemL));
        logger.assert(minAbs >= 0, "absolute value is smaller than 0");
        ld->valRight[iVar] = u0 + sign(uPlus) * minAbs;
    }
    else
    {
        ld->valRight[iVar] = u0;
    }

    if (sign(uMinus) == sign(uElemR - u0) && sign(uMinus) == sign(u0 - uElemL))
    {
        double minAbs = std::min(std::abs(uMinus), std::abs(uElemR - u0));
        minAbs = std::min(minAbs, std::abs(u0 - uElemL));
        logger.assert(minAbs >= 0, "absolute value is smaller than 0");
        ld->valLeft[iVar] = u0 - sign(uMinus) * minAbs;
    }
    else
    {
        ld->valLeft[iVar] = u0;
    }

    //std::cout << element->getID() << " " << iVar << " " << u0 << " " << ld->valLeft[iVar] << "  " << ld->valRight[iVar] << std::endl;
}*/

void SavageHutter::limitWithMinMod(Base::Element* element, const std::size_t iVar)
{
    LimiterData* ld = static_cast<LimiterData*> (element->getUserData());
    logger.assert(ld->isLimited[iVar] == true, "Called limiter on variable that should not be limited.");
    
    const PointReferenceT &pRef = element->getReferenceGeometry()->getCenter();
    logger.assert(0 == pRef[0], "xi != 0");
    
    //this does not work for first boundary, probably also not for triangular mesh.
    //maybe write getNeighbour(face)?
    const Base::Element * const elemL = element->getFace(0)->getPtrElementLeft();
    const Base::Element * const elemR = element->getFace(1)->getPtrElementRight();

    //check if left is indeed of the left side of right
    logger.assert((elemL->getPhysicalGeometry()->getLocalNodeCoordinates(0)) < (elemR->getPhysicalGeometry()->getLocalNodeCoordinates(0)), "elements left/right in wrong order");

    const double u0 = element->getSolution(0, pRef)(iVar);
    const double uElemR = elemR->getSolution(0, pRef)(iVar);
    const double uElemL = elemL->getSolution(0, pRef)(iVar);
    
    double slope =  0;
    if (sign(uElemR - u0) == sign(u0 - uElemL))
    {
        slope = sign(uElemR - u0) * std::min(std::abs(uElemR - u0), std::abs(u0 - uElemL));
    }
    //replace coordinates with "u0 + slope * x" coordinates
}
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

using LinearAlgebra::MiddleSizeVector;

/// \param[in] dimension Dimension of the domain
/// \param[in] numberOfVariables Number of variables in the PDE
/// \param[in] polynomialOrder Polynomial order of the basis functions
/// \param[in] useMatrixStorage Boolean to indicate if element and face matrices for the PDE should be stored
/// \param[in] ptrButcherTableau Pointer to a Butcher Tableau used to do the time integration with a Runge-Kutta scheme. By default this is a RK4 scheme.
SavageHutter::SavageHutter
(
 const std::size_t numOfVariables,
 const std::size_t polynomialOrder,
 const Base::ButcherTableau * const ptrButcherTableau) :
HpgemAPISimplified(numOfVariables, polynomialOrder, ptrButcherTableau),
numOfVariables_(numOfVariables), minH_(1e-5)
{
    rhsComputer_.numOfVariables_ = numOfVariables;
    rhsComputer_.epsilon_ = 1.0;
    rhsComputer_.chuteAngle_ = 0; //radians
}

SavageHutter::SavageHutter(const SHConstructorStruct& inputValues) :
HpgemAPISimplified(inputValues.numOfVariables, inputValues.polyOrder, inputValues.ptrButcherTableau),
numOfVariables_(inputValues.numOfVariables), minH_(1e-5)
{
    rhsComputer_.numOfVariables_ = inputValues.numOfVariables;
    rhsComputer_.epsilon_ = 1.0;
    rhsComputer_.chuteAngle_ = 0; //radians
    createMesh(inputValues.numElements, inputValues.meshType);
    
    for (Base::Element* elt : meshes_[0]->getElementsList())
    {
        elt->setUserData(new DryFlag);
    }
}

Base::RectangularMeshDescriptor<DIM> SavageHutter::createMeshDescription(const std::size_t numOfElementPerDirection)
{
    // Create the domain. In this case the domain is the square [0,1]^DIM and periodic.
    Base::RectangularMeshDescriptor<DIM> description;
    for (std::size_t i = 0; i < DIM; ++i)
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
    
    if (pPhys[0] < 0.5 && pPhys[0] > 0)
    {
        initialSolution(0) =  .1;
    }
    else
    {
        initialSolution(0) = 0.;
    }
    initialSolution(1) = 0;
    return initialSolution;
}
/*********************Integrate over elements and faces************************/

LinearAlgebra::MiddleSizeVector SavageHutter::computeRightHandSideAtElement(Base::Element *ptrElement, LinearAlgebra::MiddleSizeVector &solutionCoefficients, const double time)
{
    // Define the integrand function for the right hand side for the reference element.
    const std::function < LinearAlgebra::MiddleSizeVector(Base::PhysicalElement<DIM>&) > integrandFunction = [ = ](Base::PhysicalElement<DIM>& element) -> LinearAlgebra::MiddleSizeVector
    {   
        return rhsComputer_.integrandRightHandSideOnElement(element, time, solutionCoefficients);
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
    const std::function < LinearAlgebra::MiddleSizeVector(Base::PhysicalFace<DIM>&) > integrandFunction = [ = ](Base::PhysicalFace<DIM>& face) -> LinearAlgebra::MiddleSizeVector
    {   
        return rhsComputer_.integrandRightHandSideOnRefFace(face, side, solutionCoefficientsLeft, solutionCoefficientsRight);
    };

    return faceIntegrator_.integrate(ptrFace, integrandFunction);

}

LinearAlgebra::MiddleSizeVector SavageHutter::computeRightHandSideAtFace
(
 Base::Face *ptrFace,
 LinearAlgebra::MiddleSizeVector &solutionCoefficients,
 const double time
 )
{
    const std::function<LinearAlgebra::MiddleSizeVector(Base::PhysicalFace<DIM>&)> integrandFunction =
        [ = ](Base::PhysicalFace<DIM>& face) -> LinearAlgebra::MiddleSizeVector
        {
            return rhsComputer_.integrandRightHandSideOnRefFace(face, solutionCoefficients);
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
        //don't use the slope limiter if the water height is adapted with the non-negativity limiter
        const double minimum = getMinimumHeight(element);
        if (minimum < minH_)
        {
            changeHeight(element, minimum);
            changeDischarge(element);
        }
        else
        {
            if (element->getNrOfBasisFunctions() > 1)
            {
                useLimiterForElement(element);
            }
        }
        isDryElement(element);
    }
}

///\return Returns true if the element needs limiting and false if it does not need limiting.
///\details In here, the discontinuity detector of Krivodonova et. al (2004) is
///implemented to determine whether or not the given element needs limiting. 
void SavageHutter::useLimiterForElement(Base::Element *element)
{
    LinearAlgebra::MiddleSizeVector totalIntegral(numOfVariables_);
    LinearAlgebra::MiddleSizeVector numericalSolution(numOfVariables_);

    //For every face of this element, check the size of the jump
    for (const Base::Face *face : element->getFacesList())
    {
        //Compute the numerical solution and the (outward) normal of the face
        const Geometry::PointReference<0>& pRefForInflowTest = face->getReferenceGeometry()->getCenter();
        Base::Side sideOfElement = (face->getPtrElementLeft() == element) ? Base::Side::LEFT : Base::Side::RIGHT;
        LinearAlgebra::SmallVector<1> normal = face->getNormalVector(pRefForInflowTest);
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
        LinearAlgebra::SmallVector<1> velocity = computeVelocity(numericalSolution);
        //if this is an inflow face, integrate the difference in {h,hu} over the face
        if (velocity * normal < -1e-14)
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
                numericalSolutionOther = rhsComputer_.getInflowBC();
            }

            ///\todo make this multi-dimensional by integrating over the face.
            totalIntegral += numericalSolution - numericalSolutionOther;
        }
        
    }

    logger(DEBUG, "Integral over all inflow boundaries of difference for element %: %", element->getID(), totalIntegral);

    //divide the integral by dx^{(p+1)/2}, the norm of {h,hu} and the size of the face
    std::size_t p = configData_->polynomialOrder_;
        ///\todo make dx the maximum of the edge-lengths
    const PointReferenceT& center = element->getReferenceGeometry()->getCenter();
    const double dx = 2. * std::abs(element->calcJacobian(center).determinant());
    logger(DEBUG, "grid size: %", dx);
    totalIntegral /= std::pow(dx, (p + 1.) / 2);

    LinearAlgebra::MiddleSizeVector average = computeAverageOfSolution(element);
    for (std::size_t i = 0; i < numOfVariables_; ++i)
    {
        if (average(i) > 1e-10)
            totalIntegral(i) /= std::abs(average(i));
    }
    ///\todo divide by the size of the faces for DIM > 1

    //compare this with 1: if >1, discontinuous, if <1 smooth
    for (std::size_t i = 0; i < numOfVariables_; ++i)
    {
        if (std::abs(totalIntegral(i)) > 1)
        {
            logger(DEBUG, "Element % with variable % will be limited.", element->getID(), i);
            limitWithMinMod(element, i); //needs to be adapted for boundaries!
        }
    }
}

LinearAlgebra::SmallVector<1> SavageHutter::computeVelocity(LinearAlgebra::MiddleSizeVector numericalSolution)
{
    if (numericalSolution(0) < 1e-5)
    {
        return LinearAlgebra::SmallVector<1>({0.});
    }
    
    return LinearAlgebra::SmallVector<1>({numericalSolution(1) / numericalSolution(0)});
}

LinearAlgebra::MiddleSizeVector SavageHutter::computeAverageOfSolution(Base::Element* element)
{
    const std::function <LinearAlgebra::MiddleSizeVector(Base::PhysicalElement<DIM>&)> integrandFunction = 
    [ = ](Base::PhysicalElement<DIM>& elt) -> LinearAlgebra::MiddleSizeVector
    {
        LinearAlgebra::MiddleSizeVector solution = elt.getSolution();
        logger(DEBUG, "Solution at quadrature point: %", elt.getSolution());
        return solution;
    };
    LinearAlgebra::MiddleSizeVector average = (elementIntegrator_.integrate(element, integrandFunction, element->getGaussQuadratureRule()));
    //\todo: generalize to more than 1D
    PointPhysicalT p0 = element->getPhysicalGeometry()->getLocalNodeCoordinates(0);
    PointPhysicalT p1 = element->getPhysicalGeometry()->getLocalNodeCoordinates(1);
    average /= Base::L2Norm(p1-p0);
    
    logger(DEBUG, "Average over element %: %", element->getID(), average);
    logger.assert_always(average(0) > -1e-16, "Average water height negative on "
    "element %! (%), u: %",element->getID(),  average, average(1)/average(0));
    return average;
}

void SavageHutter::limitWithMinMod(Base::Element* element, const std::size_t iVar)
{    
    const Geometry::PointReference<0> &pRefFace = element->getFace(0)->getReferenceGeometry()->getCenter();
    const PointReferenceT &pRefL = element->getReferenceGeometry()->getReferenceNodeCoordinate(0); 
    const PointReferenceT &pRefR = element->getReferenceGeometry()->getReferenceNodeCoordinate(1);
    logger.assert(std::abs(-1  - pRefL[0]) < 1e-10, "xi_L != -1");    
    logger.assert(std::abs(1 - pRefR[0]) < 1e-10, "xi_R != 1");
    
    //this does not work for first boundary, probably also not for triangular mesh.
    //maybe write getNeighbour(face)?
    const Base::Element *elemL;
    const Base::Element *elemR;
    if (element->getFace(0)->isInternal())
        elemL  = element->getFace(0)->getPtrElementLeft();
    else
        return; //for now, just don't use the limiter here...
    if (element->getFace(1)->isInternal())
        elemR = element->getFace(1)->getPtrElementRight();
    else
        return; //for now, just don't use the limiter here...

    //check if left is indeed of the left side of right //TODO: this may be false for general meshes (or may even examine the same node twice)
    logger.assert((PointPhysicalT(elemL->getPhysicalGeometry()->getLocalNodeCoordinates(0)))[0] < PointPhysicalT(elemR->getPhysicalGeometry()->getLocalNodeCoordinates(0))[0], "elements left/right in wrong order");

    
    const double uPlus = element->getSolution(0, pRefR)(iVar);
    const double uMinus = element->getSolution(0, pRefL)(iVar);
    //TVB term
    const double M = 1;
    if (std::abs(uPlus - uMinus) < M * (2*element->calcJacobian(pRefL).determinant())* (2*element->calcJacobian(pRefL).determinant()))
    {
        return;
    }
    const double u0 = computeAverageOfSolution(element)(iVar);
    const double uElemR = computeAverageOfSolution(const_cast<Base::Element*>(elemR))(iVar);
    const double uElemL = computeAverageOfSolution(const_cast<Base::Element*>(elemL))(iVar);
    logger(DEBUG, "coefficients: %", element->getTimeLevelData(0, iVar));
    logger(DEBUG, "uPlus: %, basis function vals: %, %", uPlus , element->basisFunction(0,pRefR), element->basisFunction(1,pRefR));
    logger(DEBUG, "uMinus: %, basis function vals: %, %", uMinus, element->basisFunction(0,pRefL), element->basisFunction(1,pRefL));
    logger(DEBUG, "%, %, %",(uPlus - uMinus), uElemR - u0, u0 - uElemL);
    logger(DEBUG, "u0: %, average of u: %, (uPlus+uMinus)/2: %", u0, computeAverageOfSolution(element)[iVar], (uPlus + uMinus)/2);
    double slope =  0;
    if (sign(uElemR - u0) == sign(u0 - uElemL) && sign(uElemR - u0) == sign(uPlus - uMinus))
    {
        slope = sign(uElemR - u0) * std::min(std::abs(uPlus - uMinus), std::min(std::abs(uElemR - u0), std::abs(u0 - uElemL)));
    }
    
    if (std::abs(slope - std::abs(uPlus - uMinus)) > 1e-16 )
    {
        //replace coefficients with "u0 + slope/2 * xi" coefficients
        std::function<double(const PointReferenceT&)> newFun = [=] (const PointReferenceT& pRef) {return u0 + slope/2*pRef[0];};
        LinearAlgebra::MiddleSizeVector newCoeffs = projectOnBasisFuns(element, newFun);
        element->setTimeLevelData(0, iVar, newCoeffs);
    }
}

double SavageHutter::getMinimumHeight(const Base::Element* element)
{
    const Geometry::PointReference<DIM-1> &pRefFace = element->getFace(0)->getReferenceGeometry()->getCenter();
    const PointReferenceT &pRefL = element->getReferenceGeometry()->getReferenceNodeCoordinate(0); 
    const PointReferenceT &pRefR = element->getReferenceGeometry()->getReferenceNodeCoordinate(1);
    
    double minimum = std::min(element->getSolution(0,pRefL)(0), element->getSolution(0,pRefR)(0));
    for (std::size_t p = 0; p < element->getGaussQuadratureRule()->nrOfPoints(); ++p)
    {
        const PointReferenceT& pRef = element->getGaussQuadratureRule()->getPoint(p);
        minimum = std::min(minimum, element->getSolution(0,pRef)(0));
    }
    logger(DEBUG, "Minimum in element %: %", element->getID(), minimum);
    return minimum;
}

double SavageHutter::getMaximumHeight(const Base::Element* element)
{
    const Geometry::PointReference<DIM-1> &pRefFace = element->getFace(0)->getReferenceGeometry()->getCenter();
    const PointReferenceT &pRefL = element->getReferenceGeometry()->getReferenceNodeCoordinate(0); 
    const PointReferenceT &pRefR = element->getReferenceGeometry()->getReferenceNodeCoordinate(1);
    
    double maximum = std::max(element->getSolution(0,pRefL)(0), element->getSolution(0,pRefR)(0));
    for (std::size_t p = 0; p < element->getGaussQuadratureRule()->nrOfPoints(); ++p)
    {
        const PointReferenceT& pRef = element->getGaussQuadratureRule()->getPoint(p);
        maximum = std::max(maximum, element->getSolution(0,pRef)(0));
    }
    logger(DEBUG, "Maximum in element %: %", element->getID(), maximum);
    return maximum;
}

///Adapt the heigth as given in Bunya et. al. (2009) 
///\param[in] minimum the (negative) minimum height in the given element
void SavageHutter::changeHeight( Base::Element* element, double minimum)
{
    logger(DEBUG, "minimum before adaption in element %: %",element->getID(), minimum);
    const double averageH = computeAverageOfSolution(element)(0);
    const double averageHU = computeAverageOfSolution(element)(1);
    if (averageH < minH_)
    {
        //solution is constant with value average
        LinearAlgebra::MiddleSizeVector newCoeffsH = projectOnBasisFuns(element, [=](const PointReferenceT& pRef){return averageH;});
        element->setTimeLevelData(0,0,newCoeffsH);
        return;
    }
    //squeeze around average, such that the minimum height is at least minH_
    const double theta = 0;//std::min(1.0, (averageH - minH_)/(averageH - minimum));
    logger(DEBUG, "Theta: %", theta);
    std::function<double(const PointReferenceT&)> newFun = [=] (const PointReferenceT& pRef){return theta*(element->getSolution(0,pRef)(0) - averageH) + averageH;};
    LinearAlgebra::MiddleSizeVector newCoeffsH = projectOnBasisFuns(element, newFun);
    element->setTimeLevelData(0, 0, newCoeffsH);
    logger(DEBUG, "minimum after adaption: % (average %)", getMinimumHeight(element), averageH);
}

void SavageHutter::changeDischarge(Base::Element* element)
{
    
    const PointReferenceT &pRefL = element->getReferenceGeometry()->getReferenceNodeCoordinate(0); 
    const PointReferenceT &pRefR = element->getReferenceGeometry()->getReferenceNodeCoordinate(1);
    if (element->getSolution(0,pRefL)(0) <= minH_ && element->getSolution(0,pRefR)(0) <= minH_ )
    {
        LinearAlgebra::MiddleSizeVector newCoeffsHU = projectOnBasisFuns(element, [=](const PointReferenceT& pRef){return 0;});
        element->setTimeLevelData(0,1,newCoeffsHU);
        return;
    }
    if (element->getSolution(0,pRefL)(0) > minH_ && element->getSolution(0,pRefR)(0) > minH_ )
    {
        return;
    }
    double huL = element->getSolution(0,pRefL)(1);
    double huR = element->getSolution(0,pRefR)(1);
    double average = (huL + huR) / 2;
    double slope = 0;//(huL + huR)/2;

    if (element->getSolution(0, pRefR)(0) <= minH_)
    {
        slope *= -1;
    }
    std::function<double(const PointReferenceT&) > newFun = [ = ] (const PointReferenceT & pRef){return average + slope * pRef[0];};
    LinearAlgebra::MiddleSizeVector newCoeffs = projectOnBasisFuns(element, newFun);
    element->setTimeLevelData(0, 1, newCoeffs);
    logger(DEBUG, "Values of hu: [%,%]", element->getSolution(0,pRefL)(1), element->getSolution(0,pRefR)(1));
}

LinearAlgebra::MiddleSizeVector SavageHutter::projectOnBasisFuns(Base::Element *elt, std::function<double(const PointReferenceT&)> myFun)
{
    const std::size_t numBasisFuns = elt->getNrOfBasisFunctions();
    LinearAlgebra::MiddleSizeVector projection(numBasisFuns);
    for (std::size_t i = 0; i < numBasisFuns; ++i)
    {
        const std::function < double(Base::PhysicalElement<DIM>&) > integrandFunction = [ = ](Base::PhysicalElement<DIM>& element) -> double
        {   
            return myFun(element.getPointReference())*element.basisFunction(i);
        };
        double val = elementIntegrator_.integrate(elt, integrandFunction, elt->getGaussQuadratureRule());
        projection[i] = val;
    }
    LinearAlgebra::MiddleSizeMatrix massMatrix(numBasisFuns, numBasisFuns);
    for (std::size_t i = 0; i < numBasisFuns; ++i)
    {
        for (std::size_t j = 0; j < numBasisFuns; ++j)
        {
            const std::function < double(Base::PhysicalElement<DIM>&) > massFun = [ = ](Base::PhysicalElement<DIM>& element) -> double
            {   
                return element.basisFunction(j)*element.basisFunction(i);
            };
            massMatrix(i,j) = elementIntegrator_.integrate(elt, massFun);
        }
    }
    massMatrix.solve(projection);
    return projection;
}

//Note that this only holds for a flat bottom!
void SavageHutter::isDryElement(Base::Element* elt)
{
    bool flag = static_cast<DryFlag*>(elt->getUserData())->isDry;
    logger(DEBUG, "isDry: %", flag);
    double averageH = computeAverageOfSolution(elt)(0);
    flag = (averageH > minH_);
}
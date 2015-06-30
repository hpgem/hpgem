
#include "SavageHutterBase.h"


SavageHutterBase::SavageHutterBase(const SHConstructorStruct& inputValues) :
HpgemAPISimplified(inputValues.numOfVariables, inputValues.polyOrder, inputValues.ptrButcherTableau),
    numOfVariables_(inputValues.numOfVariables), dryLimit_(1e-5)
{
    createMesh(inputValues.numElements, inputValues.meshType);
}

Base::RectangularMeshDescriptor<DIM> SavageHutterBase::createMeshDescription(const std::size_t numOfElementPerDirection)
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

/*********************Integrate over elements and faces************************/

LinearAlgebra::MiddleSizeVector SavageHutterBase::computeRightHandSideAtElement(Base::Element *ptrElement, LinearAlgebra::MiddleSizeVector &solutionCoefficients, const double time)
{
    // Define the integrand function for the right hand side for the reference element.
    const std::function < LinearAlgebra::MiddleSizeVector(Base::PhysicalElement<DIM>&) > integrandFunction = [ = ](Base::PhysicalElement<DIM>& element) -> LinearAlgebra::MiddleSizeVector
    {   
        return rhsComputer_->integrandRightHandSideOnElement(element, time, solutionCoefficients);
    };

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
    //Faster for 1D: 
    const std::function < LinearAlgebra::MiddleSizeVector(Base::PhysicalFace<DIM>&) > integrandFunction = [ = ](Base::PhysicalFace<DIM>& face) -> LinearAlgebra::MiddleSizeVector
    {   
        return rhsComputer_->integrandRightHandSideOnRefFace(face, side, solutionCoefficientsLeft, solutionCoefficientsRight);
    };

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

void SavageHutterBase::limitSolution()
{
    for (Base::Element *element : meshes_[0]->getElementsList())
    {
        //don't use the slope limiter if the water height is adapted with the non-negativity limiter
        const double minimum = getMinimumHeight(element);
        if (minimum < dryLimit_)
        {
            heightLimiter_->limit(element);
        }
        else
        {
            //only limit the slope when there is a slope, so not when the solution exists of piecewise constants.
            if (element->getNrOfBasisFunctions() > 1)
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
    ///\todo generalize for more than 1D
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

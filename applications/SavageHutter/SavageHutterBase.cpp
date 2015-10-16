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
#include "Base/CommandLineOptions.h"
#include "Base/MpiContainer.h"


SavageHutterBase::SavageHutterBase(const SHConstructorStruct& inputValues) :
HpgemAPISimplified(inputValues.numberOfVariables, inputValues.polyOrder, inputValues.ptrButcherTableau),
    numberOfVariables_(inputValues.numberOfVariables), dryLimit_(1e-5), time_(0)
{
    createMesh(inputValues.numberOfElements, inputValues.meshType);
    if (DIM == 2)
    {
        //Set up the move of the mesh; note that the mesh mover gets deleted in the mesh manipulator
        const MeshMoverContraction* meshMover = new MeshMoverContraction;
        initialiseMeshMover(meshMover, 0);
        meshes_[0]->move();
        
        for (Base::Element *element : meshes_[0]->getElementsList(Base::IteratorType::GLOBAL))
        {
            element->getReferenceToPhysicalMap()->reinit();
        }
    }
}

Base::RectangularMeshDescriptor<DIM> SavageHutterBase::createMeshDescription(const std::size_t numberOfElementPerDirection)
{
    // Create the domain. In this case the domain is the square [0,1]^DIM and periodic.
    Base::RectangularMeshDescriptor<DIM> description;
    for (std::size_t i = 0; i < DIM; ++i)
    {
        description.bottomLeft_[i] = 0;
        description.topRight_[i] = 1;
        description.numElementsInDIM_[i] = 20;
        description.boundaryConditions_[i] = Base::BoundaryType::SOLID_WALL;
    }
    description.topRight_[0] = 2;
    description.numElementsInDIM_[0] = numberOfElementPerDirection;
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
            return rhsComputer_->integrandRightHandSideOnRefFace(face, solutionCoefficients, time);
        };
        
    return faceIntegrator_.integrate(ptrFace, integrandFunction, ptrFace->getGaussQuadratureRule());
}

/******************************Limiting****************************************/

void SavageHutterBase::computeOneTimeStep(double &time, const double dt)
{
    std::size_t numberOfStages = ptrButcherTableau_->getNumStages();

    // Compute intermediate Runge-Kutta stages
    ///\todo think of something for height limiter with higher order RK methods
    for (std::size_t iStage = 0; iStage < numberOfStages; iStage++)
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
    for (std::size_t jStage = 0; jStage < numberOfStages; jStage++)
    {
        scaleAndAddVector(solutionVectorId_, auxiliaryVectorIds_[jStage], dt * ptrButcherTableau_->getB(jStage));
    }    

    limitSolutionOuterLoop();

    // Update the time.
    time += dt;
    time_ = time;
    setInflowBC(time);
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
    const std::size_t numberOfVariables = element->getNumberOfUnknowns();
    const double solutionLeft = Helpers::getSolution<DIM>(element, solutionCoefficients, pRefL, numberOfVariables)(0);    
    double minimum = solutionLeft;
    for (std::size_t iPoint = 1; iPoint < element->getReferenceGeometry()->getNumberOfNodes(); ++iPoint)
    {
        const PointReferenceT &pRef = element->getReferenceGeometry()->getReferenceNodeCoordinate(iPoint);
        minimum = std::min(minimum, Helpers::getSolution<DIM>(element, solutionCoefficients, pRef, numberOfVariables)(0));
    }
    for (std::size_t p = 0; p < element->getGaussQuadratureRule()->getNumberOfPoints(); ++p)
    {
        const PointReferenceT& pRef = element->getGaussQuadratureRule()->getPoint(p);
        minimum = std::min(minimum, Helpers::getSolution<DIM>(element, solutionCoefficients, pRef, numberOfVariables)(0));
    }
    return minimum;
}

///\details function that computes the width-average of the solution by simply adding
/// the values of all elements and then divide by the number of nodes that have been
/// added. Finally, each point is multiplied by the width of the chute at that point.
// Sorry for the ugliness...
std::vector<std::pair<double, LinearAlgebra::MiddleSizeVector>> SavageHutterBase::widthAverage()
{
    
    //Since we're using a rectangular grid, we can get the number of nodes in x direction and y direction
    //The number of elements in x direction is given by the user in the commandline.
    extern Base::CommandLineOption<std::size_t>& numOfElements;
    const std::size_t nodesInXDirection = numOfElements.getValue() + 1;
    const std::size_t elementsInYDirection = meshes_[0]->getNumberOfElements(Base::IteratorType::GLOBAL) / (nodesInXDirection - 1);
    logger(DEBUG, "elements in y direction: % ", elementsInYDirection);
    
    //make xs
    ///\todo insert length of the domain here automatically instead of hardcoded
    const double dx = 11./(nodesInXDirection - 1);
    std::vector<std::pair<double, LinearAlgebra::MiddleSizeVector>> totals;
    for(std::size_t i = 0; i < nodesInXDirection; ++i)
    {
        totals.push_back(std::make_pair(i*dx, LinearAlgebra::MiddleSizeVector(numberOfVariables_)));
    }
    
    //add all values at a certain x-coordinate. To do that, first check if this value
    //of x is already in the vector. If not, make a pair of this x-value and the value of the variables
    //if there was already an entry for this x, add the value of the current point
    for (Base::Element* element : meshes_[0]->getElementsList())
    {
        const Geometry::ReferenceGeometry *referenceElement = element->getReferenceGeometry();
        for (std::size_t i = 0; i < referenceElement->getNumberOfNodes(); ++i)
        {
            const Geometry::PointReference<DIM> &nodeReference = referenceElement->getReferenceNodeCoordinate(i);
            const LinearAlgebra::MiddleSizeVector& value = element->getSolution(0, nodeReference);
            logger(DEBUG, "value: %", value);
            const Geometry::PointPhysical<DIM> node = element->referenceToPhysical(nodeReference);
            
            const auto xPosInVector = std::find_if(totals.begin(), totals.end(), [=](const std::pair<double, LinearAlgebra::MiddleSizeVector> current){return std::abs(current.first - (node)[0]) < 1e-10;});
            if (xPosInVector == totals.end())
            {
                totals.push_back(std::make_pair(node[0], value));
                logger(WARN, "x = % not found", node[0]);
            }
            else
            {
                (*xPosInVector).second += value;
            }
        }
        
    }
#ifdef HPGEM_USE_MPI
    int world_rank = Base::MPIContainer::Instance().getProcessorID();
    auto& comm = Base::MPIContainer::Instance().getComm();

    //split the pairs, since MPI can't send over a vector of pairs
    std::vector<LinearAlgebra::MiddleSizeVector> solutions;
    solutions.reserve(totals.size());
    for (std::pair<double, LinearAlgebra::MiddleSizeVector> p : totals)
    {
        solutions.push_back(p.second);
    }

    
    std::vector<LinearAlgebra::MiddleSizeVector> globalSolutions(solutions.size());
    for(std::size_t i = 0; i < solutions.size(); ++i)
    {
        LinearAlgebra::MiddleSizeVector v = solutions[i];
        comm.Reduce(v.data(), globalSolutions[i].data(), v.size(), Base::Detail::toMPIType(*v.data()), MPI::SUM, 0);
    }

    if(world_rank == 0)
    {
        logger.assert(totals.size() == globalSolutions.size(), "wrong size");
        for(std::size_t i = 0; i < globalSolutions.size(); ++i)
        {
            totals[i].second = globalSolutions[i];
        }
 #endif      
    //divide by the number of times a value for the given x-point is added, which
    //is 2 times the number of elements in y-direction for boundary nodes and 
    //4 times the number of elements in y-direction otherwise
    (*totals.begin()).second *= 2;
    (totals.back()).second *= 2;
    for (std::pair<double, LinearAlgebra::MiddleSizeVector> &val : totals)
    {
        val.second /= 4*elementsInYDirection;
        logger(INFO, "x: %, average: %", val.first, val.second);        
    }
    //just make sure we did not forget any points or that points that are the same have been put in different rows
    logger.assert(totals.size() == nodesInXDirection, "wrong number of points in vector of width-averaging");
#ifdef HPGEM_USE_MPI
    }
#endif     
    return totals;
}
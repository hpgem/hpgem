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
#include "SavageHutterRHS2D.h"
#include "HelperFunctions.h"

MiddleSizeVector SavageHutterRHS2D::integrandRightHandSideOnElement
(Base::PhysicalElement<DIM>& element, const double &time, const MiddleSizeVector &solutionCoefficients)
{
    const std::size_t numBasisFuncs = element.getElement()->getNumberOfBasisFunctions();

    MiddleSizeVector& integrand = element.getResultVector(); //just to have the correct length    
    const PointPhysicalT& pPhys = element.getPointPhysical();
    const PointReferenceT& pRef = element.getPointReference();
    const MiddleSizeVector numericalSolution = Helpers::getSolution<DIM>(element.getElement(), solutionCoefficients, pRef, numOfVariables_);
    const MiddleSizeVector physicalFlux = computePhysicalFlux(numericalSolution);
    const MiddleSizeVector source = computeSourceTerm(numericalSolution, pPhys, time);

    // Compute integrand on the physical element.
    std::size_t iVB; // Index for both basis function and variable
    for (std::size_t iB = 0; iB < numBasisFuncs; iB++) // Index for basis function
    {
        for (std::size_t iV = 0; iV < numOfVariables_; ++iV)
        {
            iVB = element.getElement()->convertToSingleIndex(iB, iV);
            integrand(iVB) = physicalFlux(2 * iV) * element.basisFunctionDeriv(iB)(0);
            integrand(iVB) += physicalFlux(2 * iV + 1) * element.basisFunctionDeriv(iB)(1);
            integrand(iVB) += source(iV) * element.basisFunction(iB);
        }
    }
    return integrand;
}

/// \details The integrand for the reference face is the same as the physical face, but scaled with the reference-to-physical face scale. This face scale is absorbed in the normal vector, since it is relatively cheap to compute the normal vector with a length (L2-norm) equal to the reference-to-physical face scale.
MiddleSizeVector SavageHutterRHS2D::integrandRightHandSideOnRefFace
(Base::PhysicalFace<DIM>& face, const Base::Side &iSide, const MiddleSizeVector &solutionCoefficientsLeft, const MiddleSizeVector &solutionCoefficientsRight)
{
    const std::size_t numTestBasisFuncs = face.getFace()->getPtrElement(iSide)->getNumberOfBasisFunctions();

    //compute numerical solution at the left side and right side of this face
    const PointReferenceOnFaceT& pRef = face.getPointReference();
    const PointReferenceT& pRefL = face.getFace()->mapRefFaceToRefElemL(pRef);
    const PointReferenceT& pRefR = face.getFace()->mapRefFaceToRefElemR(pRef);
    MiddleSizeVector solutionLeft = Helpers::getSolution<DIM>(face.getFace()->getPtrElementLeft(), solutionCoefficientsLeft, pRefL, numOfVariables_);
    MiddleSizeVector solutionRight = Helpers::getSolution<DIM>(face.getFace()->getPtrElementRight(), solutionCoefficientsRight, pRefR, numOfVariables_);

    MiddleSizeVector flux = hllcFlux(solutionLeft, solutionRight, face.getUnitNormalVector());

    if (iSide == Base::Side::RIGHT) //the normal is defined for the left element
    {
        flux *= -1;
    }

    MiddleSizeVector& integrand = face.getResultVector(iSide); // Integrand value based on n number of testbasisfunctions from element corresponding to side iSide

    for (std::size_t iFun = 0; iFun < numTestBasisFuncs; ++iFun)
    {
        for (std::size_t iVar = 0; iVar < numOfVariables_; ++iVar)
        {
            std::size_t iVarFun = face.getFace()->getPtrElement(iSide)->convertToSingleIndex(iFun, iVar);
            integrand(iVarFun) = -flux(iVar) * face.basisFunction(iSide, iFun);
        }
    }
    return integrand;
}

MiddleSizeVector SavageHutterRHS2D::integrandRightHandSideOnRefFace
(
 Base::PhysicalFace<DIM>& face,
 const MiddleSizeVector &solutionCoefficients,
 const double time
 )
{
    double normalX = face.getUnitNormalVector()[0];
    double normalY = face.getUnitNormalVector()[1];
    const std::size_t numBasisFuncs = face.getFace()->getNumberOfBasisFunctions();

    const PointReferenceOnFaceT& pRef = face.getPointReference();
    //note that at the boundary, the element is the left element by definition
    const PointReferenceT& pRefL = face.getFace()->mapRefFaceToRefElemL(pRef);
    MiddleSizeVector solution = Helpers::getSolution<DIM>(face.getFace()->getPtrElementLeft(), solutionCoefficients, pRefL, numOfVariables_);

    MiddleSizeVector flux(3);

    //outflow
    if (normalX > 0 && std::abs(normalY) < 1e-10)
    {
        logger.assert_always( solution[1]/solution[0] - std::sqrt(epsilon_ * std::cos(chuteAngle_) * solution[0]) > 0, "subcritical outflow");
        flux = hllcFlux(solution, solution, face.getUnitNormalVector());
    }
    else if (std::abs(normalY) < 1e-16) //inflow
    {
        flux = hllcFlux(inflowBC_, inflowBC_, face.getUnitNormalVector());
    }
    else //solid wall
    {
        LinearAlgebra::SmallVector<DIM> velocity({solution[1] / solution[0], solution[2] / solution[0]});
        LinearAlgebra::SmallVector<DIM> velocityReflected = velocity - 2 * (velocity * face.getUnitNormalVector()) * face.getUnitNormalVector();
        const MiddleSizeVector &reflection = MiddleSizeVector({solution[0], velocityReflected[0] * solution[0], velocityReflected[1] * solution[0]});
        flux = hllcFlux(solution, reflection, face.getUnitNormalVector());
    }

    MiddleSizeVector integrand(numOfVariables_ * numBasisFuncs);

    for (std::size_t iFun = 0; iFun < numBasisFuncs; ++iFun)
    {
        for (std::size_t iVar = 0; iVar < numOfVariables_; ++iVar)
        {
            std::size_t iVarFun = face.getFace()->getPtrElementLeft()->convertToSingleIndex(iFun, iVar);
            integrand(iVarFun) = -flux(iVar) * face.basisFunction(iFun);
        }
    }

    return integrand;
}

MiddleSizeVector SavageHutterRHS2D::computePhysicalFlux(const MiddleSizeVector &numericalSolution)
{
    const double h = numericalSolution(0);
    logger.assert(h > -1e-16, "Negative height (%)", h);
    double hu = numericalSolution(1);
    double hv = numericalSolution(2);
    double u = 0;
    double v = 0;
    if (h > minH_)
    {
        u = hu / h;
        v = hv / h;
    }
    MiddleSizeVector flux(6);
    flux(0) = hu;
    flux(1) = hv;
    flux(2) = hu * u + epsilon_ / 2 * std::cos(chuteAngle_) * h * h;
    flux(3) = hu*v;
    flux(4) = hu*v;
    flux(5) = hv * v + epsilon_ / 2 * std::cos(chuteAngle_) * h * h;
    logger(DEBUG, "flux values: %, ", flux);
    return flux;
}

MiddleSizeVector SavageHutterRHS2D::computeSourceTerm(const MiddleSizeVector& numericalSolution, const PointPhysicalT& pPhys, const double time)
{
    logger.assert(chuteAngle_ < M_PI / 2, "Angle must be in radians, not degrees!");
    const double h = numericalSolution(0);
    const double hu = numericalSolution(1);
    const double hv = numericalSolution(2);
    
    double u = 0;
    double v = 0;
    if (h > minH_)
    {
        u = hu / h;
        v = hv / h;
    }

    double uNormalized = 0;
    double vNormalized = 0;
    if (Base::L2Norm({u, v}) > 1e-16)
    {
        uNormalized = u / Base::L2Norm({u, v});
        vNormalized = v / Base::L2Norm({u, v});
    }
    const double mu = computeFriction(numericalSolution);
    const double sourceX = h * std::sin(chuteAngle_) - h * mu * uNormalized * std::cos(chuteAngle_);
    const double sourceY = -h * mu * vNormalized * std::cos(chuteAngle_);
    return MiddleSizeVector({0, sourceX, sourceY});
}

MiddleSizeVector SavageHutterRHS2D::localLaxFriedrichsFlux(const MiddleSizeVector& numericalSolutionLeft, const MiddleSizeVector& numericalSolutionRight,const LinearAlgebra::SmallVector<DIM>& normal)
{
    logger.assert(std::abs(Base::L2Norm(normal) - 1) < 1e-16, "LLF flux needs a unit normal vector");
    const double hLeft = numericalSolutionLeft[0];
    const double hRight = numericalSolutionRight[0];
    double uLeft = 0;
    double vLeft = 0;
    if (hLeft > minH_)
    {
        uLeft = numericalSolutionLeft(1) / hLeft;
        vLeft = numericalSolutionLeft(2) / hLeft;
    }

    double uRight = 0;
    double vRight = 0;
    if (hRight > minH_)
    {
        uRight = numericalSolutionRight(1) / hRight;
        vRight = numericalSolutionRight(2) / hRight;
    }

    //take the maximum of |u+-sqrt(epsilon cos(theta) h)| and |v+-sqrt(epsilon cos(theta) h)| for left and right side
    const double eigenSpeed1 = std::max(std::abs(uLeft + std::sqrt(epsilon_ * std::cos(chuteAngle_) * hLeft)), std::abs(uLeft - std::sqrt(epsilon_ * std::cos(chuteAngle_) * hLeft)));
    const double eigenSpeed2 = std::max(std::abs(vLeft + std::sqrt(epsilon_ * std::cos(chuteAngle_) * hLeft)), std::abs(vLeft - std::sqrt(epsilon_ * std::cos(chuteAngle_) * hLeft)));
    const double eigenSpeed3 = std::max(std::abs(uRight + std::sqrt(epsilon_ * std::cos(chuteAngle_) * hRight)), std::abs(uRight - std::sqrt(epsilon_ * std::cos(chuteAngle_) * hRight)));
    const double eigenSpeed4 = std::max(std::abs(vRight + std::sqrt(epsilon_ * std::cos(chuteAngle_) * hRight)), std::abs(vRight - std::sqrt(epsilon_ * std::cos(chuteAngle_) * hRight)));

    //alpha is the maximum eigenvalue of the system
    double alpha = std::max(eigenSpeed1, eigenSpeed2);
    alpha = std::max(alpha, eigenSpeed3);
    alpha = std::max(alpha, eigenSpeed4);

    MiddleSizeVector fluxLeft = computePhysicalFlux(numericalSolutionLeft);
    MiddleSizeVector fluxRight = computePhysicalFlux(numericalSolutionRight);
    MiddleSizeVector fluxNormalLeft(numOfVariables_);
    MiddleSizeVector fluxNormalRight(numOfVariables_);
    for (std::size_t i = 0; i < numOfVariables_; ++i)
    {
        fluxNormalLeft(i) = fluxLeft(2 * i) * normal(0) + fluxLeft(2 * i + 1) * normal(1);
        fluxNormalRight(i) = fluxRight(2 * i) * normal(0) + fluxRight(2 * i + 1) * normal(1);
    }

    const MiddleSizeVector numericalFlux = 0.5 *
        (fluxNormalLeft + fluxNormalRight - alpha * (numericalSolutionRight - numericalSolutionLeft));

    return numericalFlux;
}

///The HLLC flux is an approximate Riemann solver
MiddleSizeVector SavageHutterRHS2D::hllcFlux(const MiddleSizeVector& numericalSolutionLeft, const MiddleSizeVector& numericalSolutionRight, const LinearAlgebra::SmallVector<DIM>& normal)
{
    logger.assert(std::abs(Base::L2Norm(normal)-1) < 1e-14, "hllc flux needs unit normal vector");
    const double nx = normal[0];
    const double ny = normal[1];
    const double hLeft = numericalSolutionLeft[0];
    const double hRight = numericalSolutionRight[0];
    double uLeft = 0;
    double vLeft = 0;
    if (hLeft > minH_)
    {
        uLeft = numericalSolutionLeft(1) / hLeft;
        vLeft = numericalSolutionLeft(2) / hLeft;
    }
    double uRight = 0;
    double vRight = 0;
    if (hRight > minH_)
    {
        uRight = numericalSolutionRight(1) / hRight;
        vRight = numericalSolutionRight(2) / hRight;
    }
    MiddleSizeVector fluxLeft = computePhysicalFlux(numericalSolutionLeft);
    MiddleSizeVector fluxRight = computePhysicalFlux(numericalSolutionRight);
    MiddleSizeVector fluxNormalLeft(numOfVariables_);
    MiddleSizeVector fluxNormalRight(numOfVariables_);
    for (std::size_t i = 0; i < numOfVariables_; ++i)
    {
        fluxNormalLeft(i) = fluxLeft(2 * i) * nx + fluxLeft(2 * i + 1) * ny;
        fluxNormalRight(i) = fluxRight(2 * i) * nx + fluxRight(2 * i + 1) * ny;
    }
    double normalSpeedLeft = uLeft * nx + vLeft * ny;
    double normalSpeedRight = uRight * nx + vRight * ny;
    double phaseSpeedLeft = std::sqrt(epsilon_ * std::cos(chuteAngle_) * hLeft);
    double phaseSpeedRight = std::sqrt(epsilon_ * std::cos(chuteAngle_) * hRight);
    double sl=std::min(normalSpeedLeft-phaseSpeedLeft, normalSpeedRight-phaseSpeedRight);
    double sr=std::max(normalSpeedLeft+phaseSpeedLeft, normalSpeedRight+phaseSpeedRight);
    if (sl >= 0)
    {
        return fluxNormalLeft;
    }
    if (sr <= 0)
    {
        return fluxNormalRight;
    }
	return ((sr*fluxNormalLeft-sl*fluxNormalRight+sl*sr*(numericalSolutionRight-numericalSolutionLeft))/(sr-sl));
}

double SavageHutterRHS2D::computeFriction(const MiddleSizeVector& numericalSolution)
{
    const double delta1 = 17.518/180 * M_PI;
    const double delta2 = 29.712/180 * M_PI;
    const double h = numericalSolution[0];
    if (h < 1e-10)
        return std::tan(delta1);
    const double u = numericalSolution[1] / h;
    const double v = numericalSolution[2] / h;
    const double froude = std::sqrt(u * u + v * v) / std::sqrt(epsilon_ * std::cos(chuteAngle_) * h);
    const double A = 5.29;
    const double beta = 0.189;
    const double gamma = -0.080;
    const double d = 2;

    //return std::tan(45./180*M_PI);
    return std::tan(delta1) + (std::tan(delta2) - std::tan(delta1)) / (beta * h / (A*d * (froude - gamma)) + 1);
}

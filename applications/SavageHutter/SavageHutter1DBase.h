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

#ifndef SAVAGEHUTTER1DBASE_H
#define	SAVAGEHUTTER1DBASE_H

#include "SavageHutterBase.h"


class SavageHutter1DBase : public SavageHutterBase<1>
{
public:
    using PointPhysicalT = Geometry::PointPhysical<1>;
    using PointReferenceT = Geometry::PointReference<1>;
    using PointReferenceOnFaceT = Geometry::PointReference<0>;
    
    SavageHutter1DBase(std::size_t numberOfVariables, std::size_t polyOrder) 
    : SavageHutterBase(numberOfVariables, polyOrder)
    {
        
    }
    
    virtual ~SavageHutter1DBase()
    {
        
    }
    
protected:
        
    /// \brief Create a description of the domain
    Base::RectangularMeshDescriptor<1> createMeshDescription(const std::size_t numOfElementPerDirection, const double endOfDomain, const Base::BoundaryType boundary);
    
    /// \brief Function to compute the integrand for the right hand side for the reference element.
    const LinearAlgebra::MiddleSizeVector integrandRightHandSideOnElement
    (
        Base::PhysicalElement<1> &element,
        const double &time,
        const LinearAlgebra::MiddleSizeVector &solutionCoefficients
        ) override final;

    /// \brief Function to compute the integrand for the right hand side for the reference face corresponding to a boundary face.
    const LinearAlgebra::MiddleSizeVector integrandRightHandSideOnRefFace
    (
        Base::PhysicalFace<1> &face,
        const LinearAlgebra::MiddleSizeVector &solutionCoefficients,
        const double time
        ) override final;

    /// \brief Function to compute the integrand for the right hand side for the reference face corresponding to an internal face.
    const LinearAlgebra::MiddleSizeVector integrandRightHandSideOnRefFace
    (
        Base::PhysicalFace<1> &face,
        const Base::Side &iSide,
        const LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeft,
        const LinearAlgebra::MiddleSizeVector &solutionCoefficientsRight
        ) override final;

    std::pair<LinearAlgebra::MiddleSizeVector,LinearAlgebra::MiddleSizeVector> integrandsAtFace(
	Base::PhysicalFace<1> &face,
	const double &time,
	const LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeft,
	const LinearAlgebra::MiddleSizeVector &solutionCoefficientsRight);
    
    ///\brief Compute the friction as in Weinhart et. al. (2012)
    ///\todo make the friction depend on h and u (or F) instead of numericalSolution
    virtual double computeFriction(const LinearAlgebra::MiddleSizeVector &numericalSolution);
    
    ///\brief Compute friction with exponential friction law
    double computeFrictionExponential(const LinearAlgebra::MiddleSizeVector &numericalSolution);
    
    ///\brief Compute the friction with \mu = \tan frictionAngle
    double computeFrictionCoulomb(const LinearAlgebra::MiddleSizeVector &numericalSolution, const double frictionAngle);
    
private:
    virtual LinearAlgebra::MiddleSizeVector computePhysicalFlux(const LinearAlgebra::MiddleSizeVector &numericalSolution, const PointPhysicalT& pPhys) = 0;
    virtual LinearAlgebra::MiddleSizeVector computeSourceTerm(const LinearAlgebra::MiddleSizeVector &numericalSolution, const PointPhysicalT &pPhys, const double time) = 0;
    
    ///\brief Compute the local Lax-Friedrichs flux for the two given numerical solutions across a face.
    LinearAlgebra::MiddleSizeVector localLaxFriedrichsFlux(const LinearAlgebra::MiddleSizeVector &numericalSolutionLeft, const LinearAlgebra::MiddleSizeVector &NumericalSolutionRight, Base::PhysicalFace<1> &face);
    
    ///\brief Compute the HLLC flux for the two given numerical solutions across a face.
    virtual LinearAlgebra::MiddleSizeVector hllcFlux(const LinearAlgebra::MiddleSizeVector &numericalSolutionLeft, const LinearAlgebra::MiddleSizeVector &NumericalSolutionRight, const double normal, Base::PhysicalFace<1> &face);
    
};

#endif	/* SAVAGEHUTTER1DBASE_H */


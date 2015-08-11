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

#ifndef SAVAGEHUTTERRIGHTHANDSIDECOMPUTER_H
#define	SAVAGEHUTTERRIGHTHANDSIDECOMPUTER_H

#include "Base/Element.h"
#include "Base/Face.h"
#include "RightHandSideComputer.h"


using LinearAlgebra::MiddleSizeVector;

class SavageHutterRightHandSideComputer : public RightHandSideComputer
{
    using PointPhysicalT = Geometry::PointPhysical<DIM>;
    using PointReferenceT = Geometry::PointReference<DIM>;
    using PointReferenceOnFaceT = Geometry::PointReference<DIM - 1 >;

public:

    SavageHutterRightHandSideComputer(const std::size_t numOfVariables, const double epsilon, const double chuteAngle, const MiddleSizeVector inflowBC) :
    RightHandSideComputer(numOfVariables), epsilon_(epsilon), chuteAngle_(chuteAngle), inflowBC_(inflowBC), minH_(1e-10), alpha_(1) { } //simple shear: 4./3, plug flow: 1, Bagnold: 5./4

    /// \brief Compute the integrand for the right hand side for the reference element.
    MiddleSizeVector integrandRightHandSideOnElement
    (
        Base::PhysicalElement<DIM>& element,
        const double &time,
        const MiddleSizeVector &solutionCoefficients
        ) override final;

    /// \brief Compute the integrand for the right hand side for the reference face corresponding to a boundary face.
    MiddleSizeVector integrandRightHandSideOnRefFace
    (
        Base::PhysicalFace<DIM>& face,
        const MiddleSizeVector &solutionCoefficients
        ) override final;

    /// \brief Compute the integrand for the right hand side for the reference face corresponding to an internal face.
    /// Note that a face in 1D is a point.
    MiddleSizeVector integrandRightHandSideOnRefFace
    (
        Base::PhysicalFace<DIM>& face,
        const Base::Side &iSide,
        const MiddleSizeVector &solutionCoefficientsLeft,
        const MiddleSizeVector &solutionCoefficientsRight
        ) override final;
    
    void setInflowBC(MiddleSizeVector inflowBC) override final
    {
        inflowBC_ = inflowBC;
    }

private:
    MiddleSizeVector computePhysicalFlux(const MiddleSizeVector &numericalSolution);
    MiddleSizeVector computeSourceTerm(const MiddleSizeVector &numericalSolution, const PointPhysicalT &pPhys, const double time);
    MiddleSizeVector localLaxFriedrichsFlux(const MiddleSizeVector &numericalSolutionLeft, const MiddleSizeVector &NumericalSolutionRight);
    double computeFriction(const MiddleSizeVector &numericalSolution);
    
    double computeFrictionExponential(const MiddleSizeVector &numericalSolution);
    double computeFrictionCoulomb(const MiddleSizeVector &numericalSolution);

    double epsilon_;
    double chuteAngle_; //in radians
    MiddleSizeVector inflowBC_;
    double minH_; //below this height, don't divide by it, but set u to 0
    double alpha_;
};

#endif	/* SAVAGEHUTTERRIGHTHANDSIDECOMPUTER_H */


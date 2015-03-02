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

//------------------------------------------------------------------------------
// File: GaussQuadratureRulesForHypercube.cpp 
// Implementation of Gauss quadrature rules for reference hypercube.
// Lars Pesch, Fri Mar  3 12:59:11 CET 2006
//----
// Modified from original file allGaussQuadratureRules.hpp
// by M.T. Julianto, Wed Feb 25 10:45:06 UTC 2013
//---------------------------------------------------------------------------
// System includes and names imported from them:
#include <cmath>
//---------------------------------------------------------------------------
// Package includes:
#include "Integration/GlobalNamespaceIntegration.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRulesForHypercube.hpp"
#include "Geometry/ReferenceHypercube.hpp"
#include "Geometry/PointReference.hpp"
using Geometry::ReferenceHypercube;

//---------------------------------------------------------------------------
namespace QuadratureRules
{
//---------------------------------------------------------------------------
    std::string
    Cn4_1_1::getName() const
    {
        return name_;
    }

    unsigned int
    Cn4_1_1::order() const
    {
        return 1;
    }

    unsigned int
    Cn4_1_1::dimension() const
    {
        return 4;
    }

    unsigned int
    Cn4_1_1::nrOfPoints() const
    {
        return 1;
    }

    double
    Cn4_1_1::weight(unsigned int i) const
    {
        if (i < 1)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }

    void
    Cn4_1_1::getPoint(unsigned int i, PointReferenceT& p) const
    {
        if (i < 1)
            p=gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }

    Cn4_1_1::ReferenceGeometryT*
    Cn4_1_1::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    Cn4_1_1::Cn4_1_1():
        name_("Cn4_1_1"),
        refGeoPtr_(&ReferenceHypercube::Instance()),gp_(1,4)
    {
        weight_[0] = 16.0;
        gp_[0][0] = 0.0;
        gp_[0][1] = 0.0;
        gp_[0][2] = 0.0;
        gp_[0][3] = 0.0;

    }

    Cn4_1_1::~Cn4_1_1()
    {
    
    }

//---------------------------------------------------------------------------
    std::string
    Cn4_3_4::getName() const
    {
        return name_;
    }

    unsigned int
    Cn4_3_4::order() const
    {
        return 3;
    }

    unsigned int
    Cn4_3_4::dimension() const
    {
        return 4;
    }

    unsigned int
    Cn4_3_4::nrOfPoints() const
    {
        return 16;
    }

    double
    Cn4_3_4::weight(unsigned int i) const
    {
        if (i < 16)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }

    void
    Cn4_3_4::getPoint(unsigned int i, PointReferenceT& p) const
    {
        if (i < 16)
            p=gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }

    Cn4_3_4::ReferenceGeometryT*
    Cn4_3_4::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    Cn4_3_4::Cn4_3_4():
        name_("Cn4_3_4"),
        refGeoPtr_(&ReferenceHypercube::Instance()),gp_(16,4)
    {
        weight_[0] = 1.0;
        gp_[0][0] = -std::sqrt(3.0) / 3.0;
        gp_[0][1] = -std::sqrt(3.0) / 3.0;
        gp_[0][2] = -std::sqrt(3.0) / 3.0;
        gp_[0][3] = -std::sqrt(3.0) / 3.0;

        weight_[1] = 1.0;
        gp_[1][0] = +std::sqrt(3.0) / 3.0;
        gp_[1][1] = -std::sqrt(3.0) / 3.0;
        gp_[1][2] = -std::sqrt(3.0) / 3.0;
        gp_[1][3] = -std::sqrt(3.0) / 3.0;

        weight_[2] = 1.0;
        gp_[2][0] = -std::sqrt(3.0) / 3.0;
        gp_[2][1] = +std::sqrt(3.0) / 3.0;
        gp_[2][2] = -std::sqrt(3.0) / 3.0;
        gp_[2][3] = -std::sqrt(3.0) / 3.0;

        weight_[3] = 1.0;
        gp_[3][0] = +std::sqrt(3.0) / 3.0;
        gp_[3][1] = +std::sqrt(3.0) / 3.0;
        gp_[3][2] = -std::sqrt(3.0) / 3.0;
        gp_[3][3] = -std::sqrt(3.0) / 3.0;

        weight_[4] = 1.0;
        gp_[4][0] = -std::sqrt(3.0) / 3.0;
        gp_[4][1] = -std::sqrt(3.0) / 3.0;
        gp_[4][2] = +std::sqrt(3.0) / 3.0;
        gp_[4][3] = -std::sqrt(3.0) / 3.0;

        weight_[5] = 1.0;
        gp_[5][0] = +std::sqrt(3.0) / 3.0;
        gp_[5][1] = -std::sqrt(3.0) / 3.0;
        gp_[5][2] = +std::sqrt(3.0) / 3.0;
        gp_[5][3] = -std::sqrt(3.0) / 3.0;

        weight_[6] = 1.0;
        gp_[6][0] = -std::sqrt(3.0) / 3.0;
        gp_[6][1] = +std::sqrt(3.0) / 3.0;
        gp_[6][2] = +std::sqrt(3.0) / 3.0;
        gp_[6][3] = -std::sqrt(3.0) / 3.0;

        weight_[7] = 1.0;
        gp_[7][0] = +std::sqrt(3.0) / 3.0;
        gp_[7][1] = +std::sqrt(3.0) / 3.0;
        gp_[7][2] = +std::sqrt(3.0) / 3.0;
        gp_[7][3] = -std::sqrt(3.0) / 3.0;

        weight_[8] = 1.0;
        gp_[8][0] = -std::sqrt(3.0) / 3.0;
        gp_[8][1] = -std::sqrt(3.0) / 3.0;
        gp_[8][2] = -std::sqrt(3.0) / 3.0;
        gp_[8][3] = +std::sqrt(3.0) / 3.0;

        weight_[9] = 1.0;
        gp_[9][0] = +std::sqrt(3.0) / 3.0;
        gp_[9][1] = -std::sqrt(3.0) / 3.0;
        gp_[9][2] = -std::sqrt(3.0) / 3.0;
        gp_[9][3] = +std::sqrt(3.0) / 3.0;

        weight_[10] = 1.0;
        gp_[10][0] = -std::sqrt(3.0) / 3.0;
        gp_[10][1] = +std::sqrt(3.0) / 3.0;
        gp_[10][2] = -std::sqrt(3.0) / 3.0;
        gp_[10][3] = +std::sqrt(3.0) / 3.0;

        weight_[11] = 1.0;
        gp_[11][0] = +std::sqrt(3.0) / 3.0;
        gp_[11][1] = +std::sqrt(3.0) / 3.0;
        gp_[11][2] = -std::sqrt(3.0) / 3.0;
        gp_[11][3] = +std::sqrt(3.0) / 3.0;

        weight_[12] = 1.0;
        gp_[12][0] = -std::sqrt(3.0) / 3.0;
        gp_[12][1] = -std::sqrt(3.0) / 3.0;
        gp_[12][2] = +std::sqrt(3.0) / 3.0;
        gp_[12][3] = +std::sqrt(3.0) / 3.0;

        weight_[13] = 1.0;
        gp_[13][0] = +std::sqrt(3.0) / 3.0;
        gp_[13][1] = -std::sqrt(3.0) / 3.0;
        gp_[13][2] = +std::sqrt(3.0) / 3.0;
        gp_[13][3] = +std::sqrt(3.0) / 3.0;

        weight_[14] = 1.0;
        gp_[14][0] = -std::sqrt(3.0) / 3.0;
        gp_[14][1] = +std::sqrt(3.0) / 3.0;
        gp_[14][2] = +std::sqrt(3.0) / 3.0;
        gp_[14][3] = +std::sqrt(3.0) / 3.0;

        weight_[15] = 1.0;
        gp_[15][0] = +std::sqrt(3.0) / 3.0;
        gp_[15][1] = +std::sqrt(3.0) / 3.0;
        gp_[15][2] = +std::sqrt(3.0) / 3.0;
        gp_[15][3] = +std::sqrt(3.0) / 3.0;

    }

    Cn4_3_4::~Cn4_3_4()
    {
    }
//---------------------------------------------------------------------------
} // close namespace IntegrationRules

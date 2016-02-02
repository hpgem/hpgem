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

#include <cmath>
#include "Integration/QuadratureRules/GaussQuadratureRulesForHypercube.h"
#include "Geometry/ReferenceHypercube.h"
#include "Geometry/PointReference.h"

//---------------------------------------------------------------------------
namespace QuadratureRules
    {
    //---------------------------------------------------------------------------

    std::string
    Cn4_1_1::getName() const
    {
        return name_;
    }

    std::size_t
    Cn4_1_1::order() const
    {
        return 1;
    }

    std::size_t
    Cn4_1_1::dimension() const
    {
        return 4;
    }

    std::size_t
    Cn4_1_1::getNumberOfPoints() const
    {
        return 1;
    }

    double
    Cn4_1_1::weight(std::size_t i) const
    {
        logger.assert(i < 1, "%::weight - wrong index!", name_);
        return weight_[i];
    }

    const Geometry::PointReferenceBase&
    Cn4_1_1::getPoint(std::size_t i) const
    {
        logger.assert(i < 1, "%::getPoint - wrong index!", name_);
        return gp_[i];
    }

    Geometry::ReferenceGeometry*
    Cn4_1_1::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    Cn4_1_1::Cn4_1_1()
    : name_("Cn4_1_1"), refGeoPtr_(&Geometry::ReferenceHypercube::Instance()), gp_(0)
    {
        weight_[0] = 16.0;
        gp_.push_back({0.0, 0.0, 0.0, 0.0});

    }

    //---------------------------------------------------------------------------

    std::string
    Cn4_3_16::getName() const
    {
        return name_;
    }

    std::size_t
    Cn4_3_16::order() const
    {
        return 3;
    }

    std::size_t
    Cn4_3_16::dimension() const
    {
        return 4;
    }

    std::size_t
    Cn4_3_16::getNumberOfPoints() const
    {
        return 16;
    }

    double
    Cn4_3_16::weight(std::size_t i) const
    {
        logger.assert(i < 16, "%::weight - wrong index!", name_);
        return weight_[i];
    }

    const Geometry::PointReferenceBase&
    Cn4_3_16::getPoint(std::size_t i) const
    {
        logger.assert(i < 16, "%::getPoint - wrong index!", name_);
        return gp_[i];
    }

    Geometry::ReferenceGeometry*
    Cn4_3_16::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    Cn4_3_16::Cn4_3_16()
    : name_("Cn4_3_4"), refGeoPtr_(&Geometry::ReferenceHypercube::Instance()), gp_(0)
    {
        weight_[0] = 1.0;
        gp_.push_back({-std::sqrt(3.0) / 3.0, -std::sqrt(3.0) / 3.0, -std::sqrt(3.0) / 3.0, -std::sqrt(3.0) / 3.0});

        weight_[1] = 1.0;
        gp_.push_back({+std::sqrt(3.0) / 3.0, -std::sqrt(3.0) / 3.0, -std::sqrt(3.0) / 3.0, -std::sqrt(3.0) / 3.0});

        weight_[2] = 1.0;
        gp_.push_back({-std::sqrt(3.0) / 3.0, +std::sqrt(3.0) / 3.0, -std::sqrt(3.0) / 3.0, -std::sqrt(3.0) / 3.0});

        weight_[3] = 1.0;
        gp_.push_back({+std::sqrt(3.0) / 3.0, +std::sqrt(3.0) / 3.0, -std::sqrt(3.0) / 3.0, -std::sqrt(3.0) / 3.0});

        weight_[4] = 1.0;
        gp_.push_back({-std::sqrt(3.0) / 3.0, -std::sqrt(3.0) / 3.0, +std::sqrt(3.0) / 3.0, -std::sqrt(3.0) / 3.0});

        weight_[5] = 1.0;
        gp_.push_back({+std::sqrt(3.0) / 3.0, -std::sqrt(3.0) / 3.0, +std::sqrt(3.0) / 3.0, -std::sqrt(3.0) / 3.0});

        weight_[6] = 1.0;
        gp_.push_back({-std::sqrt(3.0) / 3.0, +std::sqrt(3.0) / 3.0, +std::sqrt(3.0) / 3.0, -std::sqrt(3.0) / 3.0});

        weight_[7] = 1.0;
        gp_.push_back({+std::sqrt(3.0) / 3.0, +std::sqrt(3.0) / 3.0, +std::sqrt(3.0) / 3.0, -std::sqrt(3.0) / 3.0});

        weight_[8] = 1.0;
        gp_.push_back({-std::sqrt(3.0) / 3.0, -std::sqrt(3.0) / 3.0, -std::sqrt(3.0) / 3.0, +std::sqrt(3.0) / 3.0});

        weight_[9] = 1.0;
        gp_.push_back({+std::sqrt(3.0) / 3.0, -std::sqrt(3.0) / 3.0, -std::sqrt(3.0) / 3.0, +std::sqrt(3.0) / 3.0});

        weight_[10] = 1.0;
        gp_.push_back({-std::sqrt(3.0) / 3.0, +std::sqrt(3.0) / 3.0, -std::sqrt(3.0) / 3.0, +std::sqrt(3.0) / 3.0});

        weight_[11] = 1.0;
        gp_.push_back({+std::sqrt(3.0) / 3.0, +std::sqrt(3.0) / 3.0, -std::sqrt(3.0) / 3.0, +std::sqrt(3.0) / 3.0});

        weight_[12] = 1.0;
        gp_.push_back({-std::sqrt(3.0) / 3.0, -std::sqrt(3.0) / 3.0, +std::sqrt(3.0) / 3.0, +std::sqrt(3.0) / 3.0});

        weight_[13] = 1.0;
        gp_.push_back({+std::sqrt(3.0) / 3.0, -std::sqrt(3.0) / 3.0, +std::sqrt(3.0) / 3.0, +std::sqrt(3.0) / 3.0});

        weight_[14] = 1.0;
        gp_.push_back({-std::sqrt(3.0) / 3.0, +std::sqrt(3.0) / 3.0, +std::sqrt(3.0) / 3.0, +std::sqrt(3.0) / 3.0});

        weight_[15] = 1.0;
        gp_.push_back({+std::sqrt(3.0) / 3.0, +std::sqrt(3.0) / 3.0, +std::sqrt(3.0) / 3.0, +std::sqrt(3.0) / 3.0});

    }
    //---------------------------------------------------------------------------
    }// close namespace IntegrationRules

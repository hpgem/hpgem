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
//---------------------------------------------------------------------------
// System includes and names imported from them:
#include <cmath>
//---------------------------------------------------------------------------
// Package includes:
#include "Integration/QuadratureRules/GaussQuadratureRulesForLine.h"

#include "Geometry/ReferenceLine.h"
#include "Geometry/PointReference.h"
using Geometry::ReferenceLine;

//---------------------------------------------------------------------------
namespace QuadratureRules
    {
    //---------------------------------------------------------------------------

    std::string
    Cn1_1_1::getName() const
    {
        return name_;
    }

    std::size_t
    Cn1_1_1::order() const
    {
        return 1;
    }

    std::size_t
    Cn1_1_1::dimension() const
    {
        return 1;
    }

    std::size_t
    Cn1_1_1::nrOfPoints() const
    {
        return 1;
    }

    double
    Cn1_1_1::weight(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }

    const Geometry::PointReference&
    Cn1_1_1::getPoint(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }

    Cn1_1_1::ReferenceGeometryT*
    Cn1_1_1::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    Cn1_1_1::Cn1_1_1()
    : name_("Cn1_1_1"), refGeoPtr_(&ReferenceLine::Instance()), gp_(0)
    {
        weight_[0] = 2.0;
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({0.0}));

    }

    //---------------------------------------------------------------------------

    std::string
    Cn1_3_2::getName() const
    {
        return name_;
    }

    std::size_t
    Cn1_3_2::order() const
    {
        return 3;
    }

    std::size_t
    Cn1_3_2::dimension() const
    {
        return 1;
    }

    std::size_t
    Cn1_3_2::nrOfPoints() const
    {
        return 2;
    }

    double
    Cn1_3_2::weight(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }

    const Geometry::PointReference&
    Cn1_3_2::getPoint(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }

    Cn1_3_2::ReferenceGeometryT*
    Cn1_3_2::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    Cn1_3_2::Cn1_3_2()
    : name_("Cn1_3_4"), refGeoPtr_(&ReferenceLine::Instance()), gp_(0)
    {
        weight_[0] = 1.0;
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({-std::sqrt(3.0) / 3.0}));

        weight_[1] = 1.0;
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({+std::sqrt(3.0) / 3.0}));

    }

    //---------------------------------------------------------------------------

    std::string
    Cn1_5_3::getName() const
    {
        return name_;
    }

    std::size_t
    Cn1_5_3::order() const
    {
        return 5;
    }

    std::size_t
    Cn1_5_3::dimension() const
    {
        return 1;
    }

    std::size_t
    Cn1_5_3::nrOfPoints() const
    {
        return 3;
    }

    double
    Cn1_5_3::weight(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }

    const Geometry::PointReference&
    Cn1_5_3::getPoint(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }

    Cn1_5_3::ReferenceGeometryT*
    Cn1_5_3::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    Cn1_5_3::Cn1_5_3()
    : name_("Cn1_5_9"), refGeoPtr_(&ReferenceLine::Instance()), gp_(0)
    {
        weight_[0] = 5. / 9.;
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({-std::sqrt(3.0 / 5.0)}));

        weight_[1] = 8. / 9.;
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({0.0}));

        weight_[2] = 5. / 9.;
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({+std::sqrt(3.0 / 5.0)}));

    }

    //---------------------------------------------------------------------------

    std::string
    C1_7_4::getName() const
    {
        return name_;
    }

    std::size_t
    C1_7_4::order() const
    {
        return 7;
    }

    std::size_t
    C1_7_4::dimension() const
    {
        return 1;
    }

    std::size_t
    C1_7_4::nrOfPoints() const
    {
        return 4;
    }

    double
    C1_7_4::weight(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }

    const Geometry::PointReference&
    C1_7_4::getPoint(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }

    C1_7_4::ReferenceGeometryT*
    C1_7_4::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    C1_7_4::C1_7_4()
    : name_("C1_7_x"), refGeoPtr_(&ReferenceLine::Instance()), gp_(0)
    {
        weight_[0] = (0.347854845137453);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({-0.861136311594053}));

        weight_[1] = (0.652145154862546);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({-0.339981043584856}));

        weight_[2] = (0.652145154862546);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({+0.339981043584856}));

        weight_[3] = (0.347854845137453);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({+0.861136311594053}));

    }

    //---------------------------------------------------------------------------

    std::string
    C1_9_5::getName() const
    {
        return name_;
    }

    std::size_t
    C1_9_5::order() const
    {
        return 9;
    }

    std::size_t
    C1_9_5::dimension() const
    {
        return 1;
    }

    std::size_t
    C1_9_5::nrOfPoints() const
    {
        return 5;
    }

    double
    C1_9_5::weight(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReference&
    C1_9_5::getPoint(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }

    C1_9_5::ReferenceGeometryT*
    C1_9_5::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    C1_9_5::C1_9_5()
    : name_("C1_9_25"), refGeoPtr_(&ReferenceLine::Instance()), gp_(0)
    {
        weight_[0] = (0.236926885056189);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({-0.906179845938663}));

        weight_[1] = (0.478628670499366);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({-0.538469310105683}));

        weight_[2] = (0.56888888888888888888888888);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({0.0}));

        weight_[3] = (0.478628670499366);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({0.538469310105683}));

        weight_[4] = (0.236926885056189);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({0.906179845938663}));

    }

    //---------------------------------------------------------------------------

    std::string
    C1_11_6::getName() const
    {
        return name_;
    }

    std::size_t
    C1_11_6::order() const
    {
        return 11;
    }

    std::size_t
    C1_11_6::dimension() const
    {
        return 1;
    }

    std::size_t
    C1_11_6::nrOfPoints() const
    {
        return 6;
    }

    double
    C1_11_6::weight(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }

    const Geometry::PointReference&
    C1_11_6::getPoint(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }

    C1_11_6::ReferenceGeometryT*
    C1_11_6::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    C1_11_6::C1_11_6()
    : name_("C1_11_36"), refGeoPtr_(&ReferenceLine::Instance()), gp_(0)
    {
        weight_[0] = (0.171324492379170);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({-0.932469514203152}));

        weight_[1] = (0.360761573048138);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({-0.661209386466264}));

        weight_[2] = (0.467913934572691);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({-0.238619186083196}));

        weight_[3] = (0.467913934572691);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({0.238619186083196}));

        weight_[4] = (0.360761573048138);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({0.661209386466264}));

        weight_[5] = (0.171324492379170);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({0.932469514203152}));

    }

    //---------------------------------------------------------------------------
    }// close namespace IntegrationRules

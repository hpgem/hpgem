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
#include "Integration/QuadratureRules/GaussQuadratureRulesForSquare.h"
#include "Geometry/ReferenceSquare.h"
#include "GaussQuadratureRulesForLine.h"
#include "Geometry/PointReference.h"

//---------------------------------------------------------------------------
namespace QuadratureRules
{
//---------------------------------------------------------------------------
    std::string Cn2_1_1::getName() const
    {
        return name_;
    }
    
    std::size_t Cn2_1_1::order() const
    {
        return 1;
    }
    
    std::size_t Cn2_1_1::dimension() const
    {
        return 2;
    }
    
    std::size_t Cn2_1_1::getNumberOfPoints() const
    {
        return 1;
    }
    
    double Cn2_1_1::weight(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReferenceBase&
    Cn2_1_1::getPoint(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }
    
    Geometry::ReferenceGeometry*
    Cn2_1_1::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    Cn2_1_1::Cn2_1_1()
            : name_("Cn2_1_1"), refGeoPtr_(&Geometry::ReferenceSquare::Instance()), gp_(0)
    {
        weight_[0] = (2.0) * (2.0);
        gp_.push_back(Geometry::PointReferenceFactory<2>::instance()->makePoint({0.0, 0.0}));
        
    }
    
//---------------------------------------------------------------------------
    std::string Cn2_3_4::getName() const
    {
        return name_;
    }
    
    std::size_t Cn2_3_4::order() const
    {
        return 3;
    }
    
    std::size_t Cn2_3_4::dimension() const
    {
        return 2;
    }
    
    std::size_t Cn2_3_4::getNumberOfPoints() const
    {
        return 4;
    }
    
    double Cn2_3_4::weight(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReferenceBase&
    Cn2_3_4::getPoint(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }
    
    Geometry::ReferenceGeometry*
    Cn2_3_4::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    Cn2_3_4::Cn2_3_4()
            : name_("Cn2_3_4"), refGeoPtr_(&Geometry::ReferenceSquare::Instance()), gp_(0)
    {
        weight_[0] = (1.0) * (1.0);
        gp_.push_back(Geometry::PointReferenceFactory<2>::instance()->makePoint({-std::sqrt(3.0) / 3.0, -std::sqrt(3.0) / 3.0}));
        
        weight_[1] = (1.0) * (1.0);
        gp_.push_back(Geometry::PointReferenceFactory<2>::instance()->makePoint({+std::sqrt(3.0) / 3.0, -std::sqrt(3.0) / 3.0}));
        
        weight_[2] = (1.0) * (1.0);
        gp_.push_back(Geometry::PointReferenceFactory<2>::instance()->makePoint({-std::sqrt(3.0) / 3.0, +std::sqrt(3.0) / 3.0}));
        
        weight_[3] = (1.0) * (1.0);
        gp_.push_back(Geometry::PointReferenceFactory<2>::instance()->makePoint({+std::sqrt(3.0) / 3.0, +std::sqrt(3.0) / 3.0}));
        
    }
    
//---------------------------------------------------------------------------
    std::string Cn2_5_9::getName() const
    {
        return name_;
    }
    
    std::size_t Cn2_5_9::order() const
    {
        return 5;
    }
    
    std::size_t Cn2_5_9::dimension() const
    {
        return 2;
    }
    
    std::size_t Cn2_5_9::getNumberOfPoints() const
    {
        return 9;
    }
    
    double Cn2_5_9::weight(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReferenceBase&
    Cn2_5_9::getPoint(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }
    
    Geometry::ReferenceGeometry*
    Cn2_5_9::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    Cn2_5_9::Cn2_5_9()
            : name_("Cn2_5_9"), refGeoPtr_(&Geometry::ReferenceSquare::Instance()), gp_(0)
    {
        weight_[0] = (5. / 9.) * (5. / 9.);
        gp_.push_back(Geometry::PointReferenceFactory<2>::instance()->makePoint({-std::sqrt(3.0 / 5.0), -std::sqrt(3.0 / 5.0)}));
        
        weight_[1] = (8. / 9.) * (5. / 9.);
        gp_.push_back(Geometry::PointReferenceFactory<2>::instance()->makePoint({0.0, -std::sqrt(3.0 / 5.0)}));
        
        weight_[2] = (5. / 9.) * (5. / 9.);
        gp_.push_back(Geometry::PointReferenceFactory<2>::instance()->makePoint({+std::sqrt(3.0 / 5.0), -std::sqrt(3.0 / 5.0)}));
        
        weight_[3] = (5. / 9.) * (8. / 9.);
        gp_.push_back(Geometry::PointReferenceFactory<2>::instance()->makePoint({-std::sqrt(3.0 / 5.0), 0.0}));
        
        weight_[4] = (8. / 9.) * (8. / 9.);
        gp_.push_back(Geometry::PointReferenceFactory<2>::instance()->makePoint({0.0, 0.0}));
        
        weight_[5] = (5. / 9.) * (8. / 9.);
        gp_.push_back(Geometry::PointReferenceFactory<2>::instance()->makePoint({+std::sqrt(3.0 / 5.0), 0.0}));
        
        weight_[6] = (5. / 9.) * (5. / 9.);
        gp_.push_back(Geometry::PointReferenceFactory<2>::instance()->makePoint({-std::sqrt(3.0 / 5.0), +std::sqrt(3.0 / 5.0)}));
        
        weight_[7] = (8. / 9.) * (5. / 9.);
        gp_.push_back(Geometry::PointReferenceFactory<2>::instance()->makePoint({0.0, +std::sqrt(3.0 / 5.0)}));
        
        weight_[8] = (5. / 9.) * (5. / 9.);
        gp_.push_back(Geometry::PointReferenceFactory<2>::instance()->makePoint({+std::sqrt(3.0 / 5.0), +std::sqrt(3.0 / 5.0)}));
        
    }
    
//---------------------------------------------------------------------------
    std::string C2_7_4::getName() const
    {
        return name_;
    }
    
    std::size_t C2_7_4::order() const
    {
        return 7;
    }
    
    std::size_t C2_7_4::dimension() const
    {
        return 2;
    }
    
    std::size_t C2_7_4::getNumberOfPoints() const
    {
        return 16;
    }
    
    double C2_7_4::weight(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReferenceBase&
    C2_7_4::getPoint(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }
    
    Geometry::ReferenceGeometry*
    C2_7_4::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    C2_7_4::C2_7_4()
            : name_("C2_7_4"), refGeoPtr_(&Geometry::ReferenceSquare::Instance()), gp_(0)
    {
        weight_[0] = (59. + 6. * std::sqrt(30.)) / 216.;
        gp_.push_back(Geometry::PointReferenceFactory<2>::instance()->makePoint({+std::sqrt((15. - 2. * std::sqrt(30.)) / 35.), +std::sqrt((15. - 2. * std::sqrt(30.)) / 35.)}));
        
        weight_[1] = (59. + 6. * std::sqrt(30.)) / 216.;
        gp_.push_back(Geometry::PointReferenceFactory<2>::instance()->makePoint({-std::sqrt((15. - 2. * std::sqrt(30.)) / 35.), +std::sqrt((15. - 2. * std::sqrt(30.)) / 35.)}));
        
        weight_[2] = (59. + 6. * std::sqrt(30.)) / 216.;
        gp_.push_back(Geometry::PointReferenceFactory<2>::instance()->makePoint({+std::sqrt((15. - 2. * std::sqrt(30.)) / 35.), -std::sqrt((15. - 2. * std::sqrt(30.)) / 35.)}));
        
        weight_[3] = (59. + 6. * std::sqrt(30.)) / 216.;
        gp_.push_back(Geometry::PointReferenceFactory<2>::instance()->makePoint({-std::sqrt((15. - 2. * std::sqrt(30.)) / 35.), -std::sqrt((15. - 2. * std::sqrt(30.)) / 35.)}));
        
        weight_[4] = (59. - 6. * std::sqrt(30.)) / 216.;
        gp_.push_back(Geometry::PointReferenceFactory<2>::instance()->makePoint({+std::sqrt((15. + 2. * std::sqrt(30.)) / 35.), +std::sqrt((15. + 2. * std::sqrt(30.)) / 35.)}));
        
        weight_[5] = (59. - 6. * std::sqrt(30.)) / 216.;
        gp_.push_back(Geometry::PointReferenceFactory<2>::instance()->makePoint({-std::sqrt((15. + 2. * std::sqrt(30.)) / 35.), +std::sqrt((15. + 2. * std::sqrt(30.)) / 35.)}));
        
        weight_[6] = (59. - 6. * std::sqrt(30.)) / 216.;
        gp_.push_back(Geometry::PointReferenceFactory<2>::instance()->makePoint({+std::sqrt((15. + 2. * std::sqrt(30.)) / 35.), -std::sqrt((15. + 2. * std::sqrt(30.)) / 35.)}));
        
        weight_[7] = (59. - 6. * std::sqrt(30.)) / 216.;
        gp_.push_back(Geometry::PointReferenceFactory<2>::instance()->makePoint({-std::sqrt((15. + 2. * std::sqrt(30.)) / 35.), -std::sqrt((15. + 2. * std::sqrt(30.)) / 35.)}));
        
        weight_[8] = 49. / 216.;
        gp_.push_back(Geometry::PointReferenceFactory<2>::instance()->makePoint({+std::sqrt((15. - 2. * std::sqrt(30.)) / 35.), +std::sqrt((15. + 2. * std::sqrt(30.)) / 35.)}));
        
        weight_[9] = 49. / 216.;
        gp_.push_back(Geometry::PointReferenceFactory<2>::instance()->makePoint({-std::sqrt((15. - 2. * std::sqrt(30.)) / 35.), +std::sqrt((15. + 2. * std::sqrt(30.)) / 35.)}));
        
        weight_[10] = 49. / 216.;
        gp_.push_back(Geometry::PointReferenceFactory<2>::instance()->makePoint({+std::sqrt((15. - 2. * std::sqrt(30.)) / 35.), -std::sqrt((15. + 2. * std::sqrt(30.)) / 35.)}));
        
        weight_[11] = 49. / 216.;
        gp_.push_back(Geometry::PointReferenceFactory<2>::instance()->makePoint({-std::sqrt((15. - 2. * std::sqrt(30.)) / 35.), -std::sqrt((15. + 2. * std::sqrt(30.)) / 35.)}));
        
        weight_[12] = 49. / 216.;
        gp_.push_back(Geometry::PointReferenceFactory<2>::instance()->makePoint({+std::sqrt((15. + 2. * std::sqrt(30.)) / 35.), +std::sqrt((15. - 2. * std::sqrt(30.)) / 35.)}));
        
        weight_[13] = 49. / 216.;
        gp_.push_back(Geometry::PointReferenceFactory<2>::instance()->makePoint({-std::sqrt((15. + 2. * std::sqrt(30.)) / 35.), +std::sqrt((15. - 2. * std::sqrt(30.)) / 35.)}));
        
        weight_[14] = 49. / 216.;
        gp_.push_back(Geometry::PointReferenceFactory<2>::instance()->makePoint({+std::sqrt((15. + 2. * std::sqrt(30.)) / 35.), -std::sqrt((15. - 2. * std::sqrt(30.)) / 35.)}));
        
        weight_[15] = 49. / 216.;
        gp_.push_back(Geometry::PointReferenceFactory<2>::instance()->makePoint({-std::sqrt((15. + 2. * std::sqrt(30.)) / 35.), -std::sqrt((15. - 2. * std::sqrt(30.)) / 35.)}));
        
    }
    
//---------------------------------------------------------------------------
    std::string C2_9_5::getName() const
    {
        return name_;
    }
    
    std::size_t C2_9_5::order() const
    {
        return 9;
    }
    
    std::size_t C2_9_5::dimension() const
    {
        return 2;
    }
    
    std::size_t C2_9_5::getNumberOfPoints() const
    {
        return 25;
    }
    
    double C2_9_5::weight(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReferenceBase&
    C2_9_5::getPoint(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }
    
    Geometry::ReferenceGeometry*
    C2_9_5::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    C2_9_5::C2_9_5()
            : name_("C2_9_5"), refGeoPtr_(&Geometry::ReferenceSquare::Instance()), gp_(0)
    {
        std::size_t position(0);
        C1_9_5& ruleForLine = C1_9_5::Instance();
        for (std::size_t i = 0; i < ruleForLine.getNumberOfPoints(); ++i)
        {
            for (std::size_t j = 0; j < ruleForLine.getNumberOfPoints(); ++j)
            {
                weight_[position] = ruleForLine.weight(i) * ruleForLine.weight(j);
                gp_.push_back(Geometry::PointReferenceFactory<2>::instance()->makePoint({static_cast<const Geometry::PointReference<1>&>(ruleForLine.getPoint(i))[0], static_cast<const Geometry::PointReference<1>&>(ruleForLine.getPoint(j))[0]}));
                ++position;
            }
        }
        
    }
    
//---------------------------------------------------------------------------
    std::string C2_11_6::getName() const
    {
        return name_;
    }
    
    std::size_t C2_11_6::order() const
    {
        return 11;
    }
    
    std::size_t C2_11_6::dimension() const
    {
        return 2;
    }
    
    std::size_t C2_11_6::getNumberOfPoints() const
    {
        return 36;
    }
    
    double C2_11_6::weight(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReferenceBase&
    C2_11_6::getPoint(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }
    
    Geometry::ReferenceGeometry*
    C2_11_6::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    C2_11_6::C2_11_6()
            : name_("C2_11_6"), refGeoPtr_(&Geometry::ReferenceSquare::Instance()), gp_(0)
    {
        std::size_t position(0);
        C1_11_6& ruleForLine = C1_11_6::Instance();
        for (std::size_t i = 0; i < ruleForLine.getNumberOfPoints(); ++i)
        {
            for (std::size_t j = 0; j < ruleForLine.getNumberOfPoints(); ++j)
            {
                weight_[position] = ruleForLine.weight(i) * ruleForLine.weight(j);
                gp_.push_back(Geometry::PointReferenceFactory<2>::instance()->makePoint({static_cast<const Geometry::PointReference<1>&>(ruleForLine.getPoint(i))[0], static_cast<const Geometry::PointReference<1>&>(ruleForLine.getPoint(j))[0]}));
                ++position;
            }
        }
        
    }

//---------------------------------------------------------------------------
}// close namespace IntegrationRules

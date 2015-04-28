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
#include "Integration/QuadratureRules/GaussQuadratureRulesForPyramid.h"
#include "Geometry/ReferencePyramid.h"
#include "Geometry/PointReference.h"
using Geometry::ReferencePyramid;

//---------------------------------------------------------------------------
namespace QuadratureRules
{
//---------------------------------------------------------------------------
    std::string Pyramid_1_4::getName() const
    {
        return name_;
    }
    
    std::size_t Pyramid_1_4::order() const
    {
        return 1;
    }
    
    std::size_t Pyramid_1_4::dimension() const
    {
        return 3;
    }
    
    std::size_t Pyramid_1_4::nrOfPoints() const
    {
        return 4;
    }
    
    double Pyramid_1_4::weight(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReference&
    Pyramid_1_4::getPoint(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }
    
    Pyramid_1_4::ReferenceGeometryT*
    Pyramid_1_4::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    Pyramid_1_4::Pyramid_1_4()
            : name_("Pyramid_1_1"), refGeoPtr_(&ReferencePyramid::Instance()), gp_(0)
    {
        weight_[0] = ((2.0) * (2.0)) * (0.1108884156);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(0.0) * (1. - (0.4850054945e-1)), (0.0) * (1. - (0.4850054945e-1)), 0.4850054945e-1}));
        
        weight_[1] = ((2.0) * (2.0)) * (0.1434587898);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(0.0) * (1. - (0.2386007376)), (0.0) * (1. - (0.2386007376)), 0.2386007376}));
        
        weight_[2] = ((2.0) * (2.0)) * (0.6863388717e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(0.0) * (1. - (0.5170472951)), (0.0) * (1. - (0.5170472951)), 0.5170472951}));
        
        weight_[3] = ((2.0) * (2.0)) * (0.1035224075e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(0.0) * (1. - (0.7958514179)), (0.0) * (1. - (0.7958514179)), 0.7958514179}));
        
    }
    
//---------------------------------------------------------------------------
    std::string Pyramid_3_16::getName() const
    {
        return name_;
    }
    
    std::size_t Pyramid_3_16::order() const
    {
        return 3;
    }
    
    std::size_t Pyramid_3_16::dimension() const
    {
        return 3;
    }
    
    std::size_t Pyramid_3_16::nrOfPoints() const
    {
        return 16;
    }
    
    double Pyramid_3_16::weight(std::size_t i) const    
    {
        logger.assert(i < nrOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReference&
    Pyramid_3_16::getPoint(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }
    
    Pyramid_3_16::ReferenceGeometryT*
    Pyramid_3_16::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    Pyramid_3_16::Pyramid_3_16()
            : name_("Pyramid_3_1"), refGeoPtr_(&ReferencePyramid::Instance()), gp_(0)
    {
        weight_[0] = ((1.0) * (1.0)) * (0.1108884156);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-std::sqrt(3.0) / 3.0) * (1. - (0.4850054945e-1)), (-std::sqrt(3.0) / 3.0) * (1. - (0.4850054945e-1)), 0.4850054945e-1}));
        
        weight_[1] = ((1.0) * (1.0)) * (0.1108884156);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+std::sqrt(3.0) / 3.0) * (1. - (0.4850054945e-1)), (-std::sqrt(3.0) / 3.0) * (1. - (0.4850054945e-1)), 0.4850054945e-1}));
        
        weight_[2] = ((1.0) * (1.0)) * (0.1108884156);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-std::sqrt(3.0) / 3.0) * (1. - (0.4850054945e-1)), (+std::sqrt(3.0) / 3.0) * (1. - (0.4850054945e-1)), 0.4850054945e-1}));
        
        weight_[3] = ((1.0) * (1.0)) * (0.1108884156);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+std::sqrt(3.0) / 3.0) * (1. - (0.4850054945e-1)), (+std::sqrt(3.0) / 3.0) * (1. - (0.4850054945e-1)), 0.4850054945e-1}));
        
        weight_[4] = ((1.0) * (1.0)) * (0.1434587898);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-std::sqrt(3.0) / 3.0) * (1. - (0.2386007376)), (-std::sqrt(3.0) / 3.0) * (1. - (0.2386007376)), 0.2386007376}));
        
        weight_[5] = ((1.0) * (1.0)) * (0.1434587898);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+std::sqrt(3.0) / 3.0) * (1. - (0.2386007376)), (-std::sqrt(3.0) / 3.0) * (1. - (0.2386007376)), 0.2386007376}));
        
        weight_[6] = ((1.0) * (1.0)) * (0.1434587898);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-std::sqrt(3.0) / 3.0) * (1. - (0.2386007376)), (+std::sqrt(3.0) / 3.0) * (1. - (0.2386007376)), 0.2386007376}));
        
        weight_[7] = ((1.0) * (1.0)) * (0.1434587898);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+std::sqrt(3.0) / 3.0) * (1. - (0.2386007376)), (+std::sqrt(3.0) / 3.0) * (1. - (0.2386007376)), 0.2386007376}));
        
        weight_[8] = ((1.0) * (1.0)) * (0.6863388717e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-std::sqrt(3.0) / 3.0) * (1. - (0.5170472951)), (-std::sqrt(3.0) / 3.0) * (1. - (0.5170472951)), 0.5170472951}));
        
        weight_[9] = ((1.0) * (1.0)) * (0.6863388717e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+std::sqrt(3.0) / 3.0) * (1. - (0.5170472951)), (-std::sqrt(3.0) / 3.0) * (1. - (0.5170472951)), 0.5170472951}));
        
        weight_[10] = ((1.0) * (1.0)) * (0.6863388717e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-std::sqrt(3.0) / 3.0) * (1. - (0.5170472951)), (+std::sqrt(3.0) / 3.0) * (1. - (0.5170472951)), 0.5170472951}));
        
        weight_[11] = ((1.0) * (1.0)) * (0.6863388717e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+std::sqrt(3.0) / 3.0) * (1. - (0.5170472951)), (+std::sqrt(3.0) / 3.0) * (1. - (0.5170472951)), 0.5170472951}));
        
        weight_[12] = ((1.0) * (1.0)) * (0.1035224075e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-std::sqrt(3.0) / 3.0) * (1. - (0.7958514179)), (-std::sqrt(3.0) / 3.0) * (1. - (0.7958514179)), 0.7958514179}));
        
        weight_[13] = ((1.0) * (1.0)) * (0.1035224075e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+std::sqrt(3.0) / 3.0) * (1. - (0.7958514179)), (-std::sqrt(3.0) / 3.0) * (1. - (0.7958514179)), 0.7958514179}));
        
        weight_[14] = ((1.0) * (1.0)) * (0.1035224075e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-std::sqrt(3.0) / 3.0) * (1. - (0.7958514179)), (+std::sqrt(3.0) / 3.0) * (1. - (0.7958514179)), 0.7958514179}));
        
        weight_[15] = ((1.0) * (1.0)) * (0.1035224075e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+std::sqrt(3.0) / 3.0) * (1. - (0.7958514179)), (+std::sqrt(3.0) / 3.0) * (1. - (0.7958514179)), 0.7958514179}));
        
    }
    
//---------------------------------------------------------------------------
    std::string Pyramid_5_36::getName() const
    {
        return name_;
    }
    
    std::size_t Pyramid_5_36::order() const
    {
        return 5;
    }
    
    std::size_t Pyramid_5_36::dimension() const
    {
        return 3;
    }
    
    std::size_t Pyramid_5_36::nrOfPoints() const
    {
        return 36;
    }
    
    double Pyramid_5_36::weight(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReference&
    Pyramid_5_36::getPoint(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }
    
    Pyramid_5_36::ReferenceGeometryT*
    Pyramid_5_36::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    Pyramid_5_36::Pyramid_5_36()
            : name_("Pyramid_5_1"), refGeoPtr_(&ReferencePyramid::Instance()), gp_(0)
    {
        weight_[0] = ((5. / 9.) * (5. / 9.)) * (0.1108884156);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-std::sqrt(3.0 / 5.0)) * (1. - (0.4850054945e-1)), (-std::sqrt(3.0 / 5.0)) * (1. - (0.4850054945e-1)), 0.4850054945e-1}));
        
        weight_[1] = ((8. / 9.) * (5. / 9.)) * (0.1108884156);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(0.0) * (1. - (0.4850054945e-1)), (-std::sqrt(3.0 / 5.0)) * (1. - (0.4850054945e-1)), 0.4850054945e-1}));
        
        weight_[2] = ((5. / 9.) * (5. / 9.)) * (0.1108884156);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+std::sqrt(3.0 / 5.0)) * (1. - (0.4850054945e-1)), (-std::sqrt(3.0 / 5.0)) * (1. - (0.4850054945e-1)), 0.4850054945e-1}));
        
        weight_[3] = ((5. / 9.) * (8. / 9.)) * (0.1108884156);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-std::sqrt(3.0 / 5.0)) * (1. - (0.4850054945e-1)), (0.0) * (1. - (0.4850054945e-1)), 0.4850054945e-1}));
        
        weight_[4] = ((8. / 9.) * (8. / 9.)) * (0.1108884156);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(0.0) * (1. - (0.4850054945e-1)), (0.0) * (1. - (0.4850054945e-1)), 0.4850054945e-1}));
        
        weight_[5] = ((5. / 9.) * (8. / 9.)) * (0.1108884156);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+std::sqrt(3.0 / 5.0)) * (1. - (0.4850054945e-1)), (0.0) * (1. - (0.4850054945e-1)), 0.4850054945e-1}));
        
        weight_[6] = ((5. / 9.) * (5. / 9.)) * (0.1108884156);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-std::sqrt(3.0 / 5.0)) * (1. - (0.4850054945e-1)), (+std::sqrt(3.0 / 5.0)) * (1. - (0.4850054945e-1)), 0.4850054945e-1}));
        
        weight_[7] = ((8. / 9.) * (5. / 9.)) * (0.1108884156);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(0.0) * (1. - (0.4850054945e-1)), (+std::sqrt(3.0 / 5.0)) * (1. - (0.4850054945e-1)), 0.4850054945e-1}));
        
        weight_[8] = ((5. / 9.) * (5. / 9.)) * (0.1108884156);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+std::sqrt(3.0 / 5.0)) * (1. - (0.4850054945e-1)), (+std::sqrt(3.0 / 5.0)) * (1. - (0.4850054945e-1)), 0.4850054945e-1}));
        
        weight_[9] = ((5. / 9.) * (5. / 9.)) * (0.1434587898);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-std::sqrt(3.0 / 5.0)) * (1. - (0.2386007376)), (-std::sqrt(3.0 / 5.0)) * (1. - (0.2386007376)), 0.2386007376}));
        
        weight_[10] = ((8. / 9.) * (5. / 9.)) * (0.1434587898);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(0.0) * (1. - (0.2386007376)), (-std::sqrt(3.0 / 5.0)) * (1. - (0.2386007376)), 0.2386007376}));
        
        weight_[11] = ((5. / 9.) * (5. / 9.)) * (0.1434587898);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+std::sqrt(3.0 / 5.0)) * (1. - (0.2386007376)), (-std::sqrt(3.0 / 5.0)) * (1. - (0.2386007376)), 0.2386007376}));
        
        weight_[12] = ((5. / 9.) * (8. / 9.)) * (0.1434587898);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-std::sqrt(3.0 / 5.0)) * (1. - (0.2386007376)), (0.0) * (1. - (0.2386007376)), 0.2386007376}));
        
        weight_[13] = ((8. / 9.) * (8. / 9.)) * (0.1434587898);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(0.0) * (1. - (0.2386007376)), (0.0) * (1. - (0.2386007376)), 0.2386007376}));
        
        weight_[14] = ((5. / 9.) * (8. / 9.)) * (0.1434587898);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+std::sqrt(3.0 / 5.0)) * (1. - (0.2386007376)), (0.0) * (1. - (0.2386007376)), 0.2386007376}));
        
        weight_[15] = ((5. / 9.) * (5. / 9.)) * (0.1434587898);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-std::sqrt(3.0 / 5.0)) * (1. - (0.2386007376)), (+std::sqrt(3.0 / 5.0)) * (1. - (0.2386007376)), 0.2386007376}));
        
        weight_[16] = ((8. / 9.) * (5. / 9.)) * (0.1434587898);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(0.0) * (1. - (0.2386007376)), (+std::sqrt(3.0 / 5.0)) * (1. - (0.2386007376)), 0.2386007376}));
        
        weight_[17] = ((5. / 9.) * (5. / 9.)) * (0.1434587898);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+std::sqrt(3.0 / 5.0)) * (1. - (0.2386007376)), (+std::sqrt(3.0 / 5.0)) * (1. - (0.2386007376)), 0.2386007376}));
        
        weight_[18] = ((5. / 9.) * (5. / 9.)) * (0.6863388717e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-std::sqrt(3.0 / 5.0)) * (1. - (0.5170472951)), (-std::sqrt(3.0 / 5.0)) * (1. - (0.5170472951)), 0.5170472951}));
        
        weight_[19] = ((8. / 9.) * (5. / 9.)) * (0.6863388717e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(0.0) * (1. - (0.5170472951)), (-std::sqrt(3.0 / 5.0)) * (1. - (0.5170472951)), 0.5170472951}));
        
        weight_[20] = ((5. / 9.) * (5. / 9.)) * (0.6863388717e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+std::sqrt(3.0 / 5.0)) * (1. - (0.5170472951)), (-std::sqrt(3.0 / 5.0)) * (1. - (0.5170472951)), 0.5170472951}));
        
        weight_[21] = ((5. / 9.) * (8. / 9.)) * (0.6863388717e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-std::sqrt(3.0 / 5.0)) * (1. - (0.5170472951)), (0.0) * (1. - (0.5170472951)), 0.5170472951}));
        
        weight_[22] = ((8. / 9.) * (8. / 9.)) * (0.6863388717e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(0.0) * (1. - (0.5170472951)), (0.0) * (1. - (0.5170472951)), 0.5170472951}));
        
        weight_[23] = ((5. / 9.) * (8. / 9.)) * (0.6863388717e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+std::sqrt(3.0 / 5.0)) * (1. - (0.5170472951)), (0.0) * (1. - (0.5170472951)), 0.5170472951}));
        
        weight_[24] = ((5. / 9.) * (5. / 9.)) * (0.6863388717e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-std::sqrt(3.0 / 5.0)) * (1. - (0.5170472951)), (+std::sqrt(3.0 / 5.0)) * (1. - (0.5170472951)), 0.5170472951}));
        
        weight_[25] = ((8. / 9.) * (5. / 9.)) * (0.6863388717e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(0.0) * (1. - (0.5170472951)), (+std::sqrt(3.0 / 5.0)) * (1. - (0.5170472951)), 0.5170472951}));
        
        weight_[26] = ((5. / 9.) * (5. / 9.)) * (0.6863388717e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+std::sqrt(3.0 / 5.0)) * (1. - (0.5170472951)), (+std::sqrt(3.0 / 5.0)) * (1. - (0.5170472951)), 0.5170472951}));
        
        weight_[27] = ((5. / 9.) * (5. / 9.)) * (0.1035224075e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-std::sqrt(3.0 / 5.0)) * (1. - (0.7958514179)), (-std::sqrt(3.0 / 5.0)) * (1. - (0.7958514179)), 0.7958514179}));
        
        weight_[28] = ((8. / 9.) * (5. / 9.)) * (0.1035224075e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(0.0) * (1. - (0.7958514179)), (-std::sqrt(3.0 / 5.0)) * (1. - (0.7958514179)), 0.7958514179}));
        
        weight_[29] = ((5. / 9.) * (5. / 9.)) * (0.1035224075e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+std::sqrt(3.0 / 5.0)) * (1. - (0.7958514179)), (-std::sqrt(3.0 / 5.0)) * (1. - (0.7958514179)), 0.7958514179}));
        
        weight_[30] = ((5. / 9.) * (8. / 9.)) * (0.1035224075e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-std::sqrt(3.0 / 5.0)) * (1. - (0.7958514179)), (0.0) * (1. - (0.7958514179)), 0.7958514179}));
        
        weight_[31] = ((8. / 9.) * (8. / 9.)) * (0.1035224075e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(0.0) * (1. - (0.7958514179)), (0.0) * (1. - (0.7958514179)), 0.7958514179}));
        
        weight_[32] = ((5. / 9.) * (8. / 9.)) * (0.1035224075e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+std::sqrt(3.0 / 5.0)) * (1. - (0.7958514179)), (0.0) * (1. - (0.7958514179)), 0.7958514179}));
        
        weight_[33] = ((5. / 9.) * (5. / 9.)) * (0.1035224075e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-std::sqrt(3.0 / 5.0)) * (1. - (0.7958514179)), (+std::sqrt(3.0 / 5.0)) * (1. - (0.7958514179)), 0.7958514179}));
        
        weight_[34] = ((8. / 9.) * (5. / 9.)) * (0.1035224075e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(0.0) * (1. - (0.7958514179)), (+std::sqrt(3.0 / 5.0)) * (1. - (0.7958514179)), 0.7958514179}));
        
        weight_[35] = ((5. / 9.) * (5. / 9.)) * (0.1035224075e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+std::sqrt(3.0 / 5.0)) * (1. - (0.7958514179)), (+std::sqrt(3.0 / 5.0)) * (1. - (0.7958514179)), 0.7958514179}));
        
    }
    
//---------------------------------------------------------------------------
    std::string Pyramid_7_48::getName() const
    {
        return name_;
    }
    
    std::size_t Pyramid_7_48::order() const
    {
        return 7;
    }
    
    std::size_t Pyramid_7_48::dimension() const
    {
        return 3;
    }
    
    std::size_t Pyramid_7_48::nrOfPoints() const
    {
        return 48;
    }
    
    double Pyramid_7_48::weight(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReference&
    Pyramid_7_48::getPoint(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }
    
    Pyramid_7_48::ReferenceGeometryT*
    Pyramid_7_48::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    Pyramid_7_48::Pyramid_7_48()
            : name_("Pyramid_7_1"), refGeoPtr_(&ReferencePyramid::Instance()), gp_(0)
    {
        weight_[0] = (98. / 405.) * (0.1108884156);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+(std::sqrt(6. / 7.))) * (1. - (0.4850054945e-1)), (0.) * (1. - (0.4850054945e-1)), 0.4850054945e-1}));
        
        weight_[1] = (98. / 405.) * (0.1108884156);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-(std::sqrt(6. / 7.))) * (1. - (0.4850054945e-1)), (0.) * (1. - (0.4850054945e-1)), 0.4850054945e-1}));
        
        weight_[2] = (98. / 405.) * (0.1108884156);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(0.) * (1. - (0.4850054945e-1)), (+(std::sqrt(6. / 7.))) * (1. - (0.4850054945e-1)), 0.4850054945e-1}));
        
        weight_[3] = (98. / 405.) * (0.1108884156);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(0.) * (1. - (0.4850054945e-1)), (-(std::sqrt(6. / 7.))) * (1. - (0.4850054945e-1)), 0.4850054945e-1}));
        
        weight_[4] = ((178981. + 2769. * std::sqrt(583.)) / 472230.) * (0.1108884156);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+(std::sqrt((114. - 3. * std::sqrt(583.)) / 287.))) * (1. - (0.4850054945e-1)), (+(std::sqrt((114. - 3. * std::sqrt(583.)) / 287.))) * (1. - (0.4850054945e-1)), 0.4850054945e-1}));
        
        weight_[5] = ((178981. + 2769. * std::sqrt(583.)) / 472230.) * (0.1108884156);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+(std::sqrt((114. - 3. * std::sqrt(583.)) / 287.))) * (1. - (0.4850054945e-1)), (-(std::sqrt((114. - 3. * std::sqrt(583.)) / 287.))) * (1. - (0.4850054945e-1)), 0.4850054945e-1}));
        
        weight_[6] = ((178981. + 2769. * std::sqrt(583.)) / 472230.) * (0.1108884156);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-(std::sqrt((114. - 3. * std::sqrt(583.)) / 287.))) * (1. - (0.4850054945e-1)), (+(std::sqrt((114. - 3. * std::sqrt(583.)) / 287.))) * (1. - (0.4850054945e-1)), 0.4850054945e-1}));
        
        weight_[7] = ((178981. + 2769. * std::sqrt(583.)) / 472230.) * (0.1108884156);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-(std::sqrt((114. - 3. * std::sqrt(583.)) / 287.))) * (1. - (0.4850054945e-1)), (-(std::sqrt((114. - 3. * std::sqrt(583.)) / 287.))) * (1. - (0.4850054945e-1)), 0.4850054945e-1}));
        
        weight_[8] = ((178981. - 2769. * std::sqrt(583.)) / 472230.) * (0.1108884156);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+(std::sqrt((114. + 3. * std::sqrt(583.)) / 287.))) * (1. - (0.4850054945e-1)), (+(std::sqrt((114. + 3. * std::sqrt(583.)) / 287.))) * (1. - (0.4850054945e-1)), 0.4850054945e-1}));
        
        weight_[9] = ((178981. - 2769. * std::sqrt(583.)) / 472230.) * (0.1108884156);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+(std::sqrt((114. + 3. * std::sqrt(583.)) / 287.))) * (1. - (0.4850054945e-1)), (-(std::sqrt((114. + 3. * std::sqrt(583.)) / 287.))) * (1. - (0.4850054945e-1)), 0.4850054945e-1}));
        
        weight_[10] = ((178981. - 2769. * std::sqrt(583.)) / 472230.) * (0.1108884156);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-(std::sqrt((114. + 3. * std::sqrt(583.)) / 287.))) * (1. - (0.4850054945e-1)), (+(std::sqrt((114. + 3. * std::sqrt(583.)) / 287.))) * (1. - (0.4850054945e-1)), 0.4850054945e-1}));
        
        weight_[11] = ((178981. - 2769. * std::sqrt(583.)) / 472230.) * (0.1108884156);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-(std::sqrt((114. + 3. * std::sqrt(583.)) / 287.))) * (1. - (0.4850054945e-1)), (-(std::sqrt((114. + 3. * std::sqrt(583.)) / 287.))) * (1. - (0.4850054945e-1)), 0.4850054945e-1}));
        
        weight_[12] = (98. / 405.) * (0.1434587898);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+(std::sqrt(6. / 7.))) * (1. - (0.2386007376)), (0.) * (1. - (0.2386007376)), 0.2386007376}));
        
        weight_[13] = (98. / 405.) * (0.1434587898);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-(std::sqrt(6. / 7.))) * (1. - (0.2386007376)), (0.) * (1. - (0.2386007376)), 0.2386007376}));
        
        weight_[14] = (98. / 405.) * (0.1434587898);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(0.) * (1. - (0.2386007376)), (+(std::sqrt(6. / 7.))) * (1. - (0.2386007376)), 0.2386007376}));
        
        weight_[15] = (98. / 405.) * (0.1434587898);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(0.) * (1. - (0.2386007376)), (-(std::sqrt(6. / 7.))) * (1. - (0.2386007376)), 0.2386007376}));
        
        weight_[16] = ((178981. + 2769. * std::sqrt(583.)) / 472230.) * (0.1434587898);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+(std::sqrt((114. - 3. * std::sqrt(583.)) / 287.))) * (1. - (0.2386007376)), (+(std::sqrt((114. - 3. * std::sqrt(583.)) / 287.))) * (1. - (0.2386007376)), 0.2386007376}));
        
        weight_[17] = ((178981. + 2769. * std::sqrt(583.)) / 472230.) * (0.1434587898);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+(std::sqrt((114. - 3. * std::sqrt(583.)) / 287.))) * (1. - (0.2386007376)), (-(std::sqrt((114. - 3. * std::sqrt(583.)) / 287.))) * (1. - (0.2386007376)), 0.2386007376}));
        
        weight_[18] = ((178981. + 2769. * std::sqrt(583.)) / 472230.) * (0.1434587898);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-(std::sqrt((114. - 3. * std::sqrt(583.)) / 287.))) * (1. - (0.2386007376)), (+(std::sqrt((114. - 3. * std::sqrt(583.)) / 287.))) * (1. - (0.2386007376)), 0.2386007376}));
        
        weight_[19] = ((178981. + 2769. * std::sqrt(583.)) / 472230.) * (0.1434587898);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-(std::sqrt((114. - 3. * std::sqrt(583.)) / 287.))) * (1. - (0.2386007376)), (-(std::sqrt((114. - 3. * std::sqrt(583.)) / 287.))) * (1. - (0.2386007376)), 0.2386007376}));
        
        weight_[20] = ((178981. - 2769. * std::sqrt(583.)) / 472230.) * (0.1434587898);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+(std::sqrt((114. + 3. * std::sqrt(583.)) / 287.))) * (1. - (0.2386007376)), (+(std::sqrt((114. + 3. * std::sqrt(583.)) / 287.))) * (1. - (0.2386007376)), 0.2386007376}));
        
        weight_[21] = ((178981. - 2769. * std::sqrt(583.)) / 472230.) * (0.1434587898);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+(std::sqrt((114. + 3. * std::sqrt(583.)) / 287.))) * (1. - (0.2386007376)), (-(std::sqrt((114. + 3. * std::sqrt(583.)) / 287.))) * (1. - (0.2386007376)), 0.2386007376}));
        
        weight_[22] = ((178981. - 2769. * std::sqrt(583.)) / 472230.) * (0.1434587898);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-(std::sqrt((114. + 3. * std::sqrt(583.)) / 287.))) * (1. - (0.2386007376)), (+(std::sqrt((114. + 3. * std::sqrt(583.)) / 287.))) * (1. - (0.2386007376)), 0.2386007376}));
        
        weight_[23] = ((178981. - 2769. * std::sqrt(583.)) / 472230.) * (0.1434587898);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-(std::sqrt((114. + 3. * std::sqrt(583.)) / 287.))) * (1. - (0.2386007376)), (-(std::sqrt((114. + 3. * std::sqrt(583.)) / 287.))) * (1. - (0.2386007376)), 0.2386007376}));
        
        weight_[24] = (98. / 405.) * (0.6863388717e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+(std::sqrt(6. / 7.))) * (1. - (0.5170472951)), (0.) * (1. - (0.5170472951)), 0.5170472951}));
        
        weight_[25] = (98. / 405.) * (0.6863388717e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-(std::sqrt(6. / 7.))) * (1. - (0.5170472951)), (0.) * (1. - (0.5170472951)), 0.5170472951}));
        
        weight_[26] = (98. / 405.) * (0.6863388717e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(0.) * (1. - (0.5170472951)), (+(std::sqrt(6. / 7.))) * (1. - (0.5170472951)), 0.5170472951}));
        
        weight_[27] = (98. / 405.) * (0.6863388717e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(0.) * (1. - (0.5170472951)), (-(std::sqrt(6. / 7.))) * (1. - (0.5170472951)), 0.5170472951}));
        
        weight_[28] = ((178981. + 2769. * std::sqrt(583.)) / 472230.) * (0.6863388717e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+(std::sqrt((114. - 3. * std::sqrt(583.)) / 287.))) * (1. - (0.5170472951)), (+(std::sqrt((114. - 3. * std::sqrt(583.)) / 287.))) * (1. - (0.5170472951)), 0.5170472951}));
        
        weight_[29] = ((178981. + 2769. * std::sqrt(583.)) / 472230.) * (0.6863388717e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+(std::sqrt((114. - 3. * std::sqrt(583.)) / 287.))) * (1. - (0.5170472951)), (-(std::sqrt((114. - 3. * std::sqrt(583.)) / 287.))) * (1. - (0.5170472951)), 0.5170472951}));
        
        weight_[30] = ((178981. + 2769. * std::sqrt(583.)) / 472230.) * (0.6863388717e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-(std::sqrt((114. - 3. * std::sqrt(583.)) / 287.))) * (1. - (0.5170472951)), (+(std::sqrt((114. - 3. * std::sqrt(583.)) / 287.))) * (1. - (0.5170472951)), 0.5170472951}));
        
        weight_[31] = ((178981. + 2769. * std::sqrt(583.)) / 472230.) * (0.6863388717e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-(std::sqrt((114. - 3. * std::sqrt(583.)) / 287.))) * (1. - (0.5170472951)), (-(std::sqrt((114. - 3. * std::sqrt(583.)) / 287.))) * (1. - (0.5170472951)), 0.5170472951}));
        
        weight_[32] = ((178981. - 2769. * std::sqrt(583.)) / 472230.) * (0.6863388717e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+(std::sqrt((114. + 3. * std::sqrt(583.)) / 287.))) * (1. - (0.5170472951)), (+(std::sqrt((114. + 3. * std::sqrt(583.)) / 287.))) * (1. - (0.5170472951)), 0.5170472951}));
        
        weight_[33] = ((178981. - 2769. * std::sqrt(583.)) / 472230.) * (0.6863388717e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+(std::sqrt((114. + 3. * std::sqrt(583.)) / 287.))) * (1. - (0.5170472951)), (-(std::sqrt((114. + 3. * std::sqrt(583.)) / 287.))) * (1. - (0.5170472951)), 0.5170472951}));
        
        weight_[34] = ((178981. - 2769. * std::sqrt(583.)) / 472230.) * (0.6863388717e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-(std::sqrt((114. + 3. * std::sqrt(583.)) / 287.))) * (1. - (0.5170472951)), (+(std::sqrt((114. + 3. * std::sqrt(583.)) / 287.))) * (1. - (0.5170472951)), 0.5170472951}));
        
        weight_[35] = ((178981. - 2769. * std::sqrt(583.)) / 472230.) * (0.6863388717e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-(std::sqrt((114. + 3. * std::sqrt(583.)) / 287.))) * (1. - (0.5170472951)), (-(std::sqrt((114. + 3. * std::sqrt(583.)) / 287.))) * (1. - (0.5170472951)), 0.5170472951}));
        
        weight_[36] = (98. / 405.) * (0.1035224075e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+(std::sqrt(6. / 7.))) * (1. - (0.7958514179)), (0.) * (1. - (0.7958514179)), 0.7958514179}));
        
        weight_[37] = (98. / 405.) * (0.1035224075e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-(std::sqrt(6. / 7.))) * (1. - (0.7958514179)), (0.) * (1. - (0.7958514179)), 0.7958514179}));
        
        weight_[38] = (98. / 405.) * (0.1035224075e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(0.) * (1. - (0.7958514179)), (+(std::sqrt(6. / 7.))) * (1. - (0.7958514179)), 0.7958514179}));
        
        weight_[39] = (98. / 405.) * (0.1035224075e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(0.) * (1. - (0.7958514179)), (-(std::sqrt(6. / 7.))) * (1. - (0.7958514179)), 0.7958514179}));
        
        weight_[40] = ((178981. + 2769. * std::sqrt(583.)) / 472230.) * (0.1035224075e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+(std::sqrt((114. - 3. * std::sqrt(583.)) / 287.))) * (1. - (0.7958514179)), (+(std::sqrt((114. - 3. * std::sqrt(583.)) / 287.))) * (1. - (0.7958514179)), 0.7958514179}));
        
        weight_[41] = ((178981. + 2769. * std::sqrt(583.)) / 472230.) * (0.1035224075e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+(std::sqrt((114. - 3. * std::sqrt(583.)) / 287.))) * (1. - (0.7958514179)), (-(std::sqrt((114. - 3. * std::sqrt(583.)) / 287.))) * (1. - (0.7958514179)), 0.7958514179}));
        
        weight_[42] = ((178981. + 2769. * std::sqrt(583.)) / 472230.) * (0.1035224075e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-(std::sqrt((114. - 3. * std::sqrt(583.)) / 287.))) * (1. - (0.7958514179)), (+(std::sqrt((114. - 3. * std::sqrt(583.)) / 287.))) * (1. - (0.7958514179)), 0.7958514179}));
        
        weight_[43] = ((178981. + 2769. * std::sqrt(583.)) / 472230.) * (0.1035224075e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-(std::sqrt((114. - 3. * std::sqrt(583.)) / 287.))) * (1. - (0.7958514179)), (-(std::sqrt((114. - 3. * std::sqrt(583.)) / 287.))) * (1. - (0.7958514179)), 0.7958514179}));
        
        weight_[44] = ((178981. - 2769. * std::sqrt(583.)) / 472230.) * (0.1035224075e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+(std::sqrt((114. + 3. * std::sqrt(583.)) / 287.))) * (1. - (0.7958514179)), (+(std::sqrt((114. + 3. * std::sqrt(583.)) / 287.))) * (1. - (0.7958514179)), 0.7958514179}));
        
        weight_[45] = ((178981. - 2769. * std::sqrt(583.)) / 472230.) * (0.1035224075e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(+(std::sqrt((114. + 3. * std::sqrt(583.)) / 287.))) * (1. - (0.7958514179)), (-(std::sqrt((114. + 3. * std::sqrt(583.)) / 287.))) * (1. - (0.7958514179)), 0.7958514179}));
        
        weight_[46] = ((178981. - 2769. * std::sqrt(583.)) / 472230.) * (0.1035224075e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-(std::sqrt((114. + 3. * std::sqrt(583.)) / 287.))) * (1. - (0.7958514179)), (+(std::sqrt((114. + 3. * std::sqrt(583.)) / 287.))) * (1. - (0.7958514179)), 0.7958514179}));
        
        weight_[47] = ((178981. - 2769. * std::sqrt(583.)) / 472230.) * (0.1035224075e-1);
        gp_.push_back(Geometry::PointReferenceFactory::instance()->makePoint({(-(std::sqrt((114. + 3. * std::sqrt(583.)) / 287.))) * (1. - (0.7958514179)), (-(std::sqrt((114. + 3. * std::sqrt(583.)) / 287.))) * (1. - (0.7958514179)), 0.7958514179}));
        
    }

//---------------------------------------------------------------------------
}// close namespace IntegrationRules

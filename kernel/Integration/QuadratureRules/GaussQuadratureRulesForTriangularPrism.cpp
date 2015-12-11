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
#include "Integration/QuadratureRules/GaussQuadratureRulesForTriangularPrism.h"
#include "Geometry/ReferenceTriangularPrism.h"
#include "Geometry/PointReference.h"

//---------------------------------------------------------------------------
namespace QuadratureRules
{
//---------------------------------------------------------------------------
    std::string TriPrism_1_1::getName() const
    {
        return name_;
    }
    
    std::size_t TriPrism_1_1::order() const
    {
        return 1;
    }
    
    std::size_t TriPrism_1_1::dimension() const
    {
        return 3;
    }
    
    std::size_t TriPrism_1_1::getNumberOfPoints() const
    {
        return 1;
    }
    
    double TriPrism_1_1::weight(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReferenceBase&
    TriPrism_1_1::getPoint(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }
    
    Geometry::ReferenceGeometry*
    TriPrism_1_1::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    TriPrism_1_1::TriPrism_1_1()
            : name_("TriPrism_1_1"), refGeoPtr_(&Geometry::ReferenceTriangularPrism::Instance()), gp_(0)
    {
        weight_[0] = (0.5) * (2.0);
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.0 / 3.0, 1.0 / 3.0, 0.0}));
        
    }
    
//---------------------------------------------------------------------------
    std::string TriPrism_3_8::getName() const
    {
        return name_;
    }
    
    std::size_t TriPrism_3_8::order() const
    {
        return 3;
    }
    
    std::size_t TriPrism_3_8::dimension() const
    {
        return 3;
    }
    
    std::size_t TriPrism_3_8::getNumberOfPoints() const
    {
        return 8;
    }
    
    double TriPrism_3_8::weight(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReferenceBase&
    TriPrism_3_8::getPoint(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }
    
    Geometry::ReferenceGeometry*
    TriPrism_3_8::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    TriPrism_3_8::TriPrism_3_8()
            : name_("TriPrism_3_1"), refGeoPtr_(&Geometry::ReferenceTriangularPrism::Instance()), gp_(0)
    {
        weight_[0] = (-9. / 32.) * (1.0);
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1. / 3., 1. / 3., -std::sqrt(3.0) / 3.0}));
        
        weight_[1] = (25. / 96.) * (1.0);
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1. / 5., 1. / 5., -std::sqrt(3.0) / 3.0}));
        
        weight_[2] = (25. / 96.) * (1.0);
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1. / 5., 3. / 5., -std::sqrt(3.0) / 3.0}));
        
        weight_[3] = (25. / 96.) * (1.0);
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3. / 5., 1. / 5., -std::sqrt(3.0) / 3.0}));
        
        weight_[4] = (-9. / 32.) * (1.0);
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1. / 3., 1. / 3., +std::sqrt(3.0) / 3.0}));
        
        weight_[5] = (25. / 96.) * (1.0);
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1. / 5., 1. / 5., +std::sqrt(3.0) / 3.0}));
        
        weight_[6] = (25. / 96.) * (1.0);
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1. / 5., 3. / 5., +std::sqrt(3.0) / 3.0}));
        
        weight_[7] = (25. / 96.) * (1.0);
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3. / 5., 1. / 5., +std::sqrt(3.0) / 3.0}));
        
    }
    
//---------------------------------------------------------------------------
    std::string TriPrism_5_21::getName() const
    {
        return name_;
    }
    
    std::size_t TriPrism_5_21::order() const
    {
        return 5;
    }
    
    std::size_t TriPrism_5_21::dimension() const
    {
        return 3;
    }
    
    std::size_t TriPrism_5_21::getNumberOfPoints() const
    {
        return 21;
    }
    
    double TriPrism_5_21::weight(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReferenceBase&
    TriPrism_5_21::getPoint(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }
    
    Geometry::ReferenceGeometry*
    TriPrism_5_21::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    TriPrism_5_21::TriPrism_5_21()
            : name_("TriPrism_5_1"), refGeoPtr_(&Geometry::ReferenceTriangularPrism::Instance()), gp_(0)
    {
        weight_[0] = (9. / 80.) * (5. / 9.);
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1. / 3., 1. / 3., -std::sqrt(3.0 / 5.0)}));
        
        weight_[1] = ((155. - std::sqrt(15.)) / 2400.) * (5. / 9.);
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({(6. - std::sqrt(15.)) / 21., (6. - std::sqrt(15.)) / 21., -std::sqrt(3.0 / 5.0)}));
        
        weight_[2] = ((155. - std::sqrt(15.)) / 2400.) * (5. / 9.);
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({(9. + 2. * std::sqrt(15.)) / 21., (6. - std::sqrt(15.)) / 21., -std::sqrt(3.0 / 5.0)}));
        
        weight_[3] = ((155. - std::sqrt(15.)) / 2400.) * (5. / 9.);
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({(6. - std::sqrt(15.)) / 21., (9. + 2. * std::sqrt(15.)) / 21., -std::sqrt(3.0 / 5.0)}));
        
        weight_[4] = ((155. + std::sqrt(15.)) / 2400.) * (5. / 9.);
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({(6. + std::sqrt(15.)) / 21., (6. + std::sqrt(15.)) / 21., -std::sqrt(3.0 / 5.0)}));
        
        weight_[5] = ((155. + std::sqrt(15.)) / 2400.) * (5. / 9.);
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({(9. - 2. * std::sqrt(15.)) / 21., (6. + std::sqrt(15.)) / 21., -std::sqrt(3.0 / 5.0)}));
        
        weight_[6] = ((155. + std::sqrt(15.)) / 2400.) * (5. / 9.);
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({(6. + std::sqrt(15.)) / 21., (9. - 2. * std::sqrt(15.)) / 21., -std::sqrt(3.0 / 5.0)}));
        
        weight_[7] = (9. / 80.) * (8. / 9.);
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1. / 3., 1. / 3., 0.0}));
        
        weight_[8] = ((155. - std::sqrt(15.)) / 2400.) * (8. / 9.);
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({(6. - std::sqrt(15.)) / 21., (6. - std::sqrt(15.)) / 21., 0.0}));
        
        weight_[9] = ((155. - std::sqrt(15.)) / 2400.) * (8. / 9.);
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({(9. + 2. * std::sqrt(15.)) / 21., (6. - std::sqrt(15.)) / 21., 0.0}));
        
        weight_[10] = ((155. - std::sqrt(15.)) / 2400.) * (8. / 9.);
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({(6. - std::sqrt(15.)) / 21., (9. + 2. * std::sqrt(15.)) / 21., 0.0}));
        
        weight_[11] = ((155. + std::sqrt(15.)) / 2400.) * (8. / 9.);
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({(6. + std::sqrt(15.)) / 21., (6. + std::sqrt(15.)) / 21., 0.0}));
        
        weight_[12] = ((155. + std::sqrt(15.)) / 2400.) * (8. / 9.);
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({(9. - 2. * std::sqrt(15.)) / 21., (6. + std::sqrt(15.)) / 21., 0.0}));
        
        weight_[13] = ((155. + std::sqrt(15.)) / 2400.) * (8. / 9.);
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({(6. + std::sqrt(15.)) / 21., (9. - 2. * std::sqrt(15.)) / 21., 0.0}));
        
        weight_[14] = (9. / 80.) * (5. / 9.);
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1. / 3., 1. / 3., +std::sqrt(3.0 / 5.0)}));
        
        weight_[15] = ((155. - std::sqrt(15.)) / 2400.) * (5. / 9.);
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({(6. - std::sqrt(15.)) / 21., (6. - std::sqrt(15.)) / 21., +std::sqrt(3.0 / 5.0)}));
        
        weight_[16] = ((155. - std::sqrt(15.)) / 2400.) * (5. / 9.);
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({(9. + 2. * std::sqrt(15.)) / 21., (6. - std::sqrt(15.)) / 21., +std::sqrt(3.0 / 5.0)}));
        
        weight_[17] = ((155. - std::sqrt(15.)) / 2400.) * (5. / 9.);
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({(6. - std::sqrt(15.)) / 21., (9. + 2. * std::sqrt(15.)) / 21., +std::sqrt(3.0 / 5.0)}));
        
        weight_[18] = ((155. + std::sqrt(15.)) / 2400.) * (5. / 9.);
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({(6. + std::sqrt(15.)) / 21., (6. + std::sqrt(15.)) / 21., +std::sqrt(3.0 / 5.0)}));
        
        weight_[19] = ((155. + std::sqrt(15.)) / 2400.) * (5. / 9.);
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({(9. - 2. * std::sqrt(15.)) / 21., (6. + std::sqrt(15.)) / 21., +std::sqrt(3.0 / 5.0)}));
        
        weight_[20] = ((155. + std::sqrt(15.)) / 2400.) * (5. / 9.);
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({(6. + std::sqrt(15.)) / 21., (9. - 2. * std::sqrt(15.)) / 21., +std::sqrt(3.0 / 5.0)}));
        
    }
    
//---------------------------------------------------------------------------
    std::string TriPrism_7_64::getName() const
    {
        return name_;
    }
    
    std::size_t TriPrism_7_64::order() const
    {
        return 7;
    }
    
    std::size_t TriPrism_7_64::dimension() const
    {
        return 3;
    }
    
    std::size_t TriPrism_7_64::getNumberOfPoints() const
    {
        return 64;
    }
    
    double TriPrism_7_64::weight(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReferenceBase&
    TriPrism_7_64::getPoint(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }
    
    Geometry::ReferenceGeometry*
    TriPrism_7_64::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    TriPrism_7_64::TriPrism_7_64()
            : name_("TriPrism_7_1"), refGeoPtr_(&Geometry::ReferenceTriangularPrism::Instance()), gp_(0)
    {
        weight_[0] = ((0.1739274226) * (0.1355069134)) * ((2. * 0.1739274226));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((-0.861136312 + 1.) / 2.) * (1. - (0.0571041961)), 0.0571041961, (-0.861136312)}));
        
        weight_[1] = ((0.3260725774) * (0.1355069134)) * ((2. * 0.1739274226));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((-0.3399810436 + 1.) / 2.) * (1. - (0.0571041961)), 0.0571041961, (-0.861136312)}));
        
        weight_[2] = ((0.3260725774) * (0.1355069134)) * ((2. * 0.1739274226));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((+0.3399810436 + 1.) / 2.) * (1. - (0.0571041961)), 0.0571041961, (-0.861136312)}));
        
        weight_[3] = ((0.1739274226) * (0.1355069134)) * ((2. * 0.1739274226));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((+0.861136312 + 1.) / 2.) * (1. - (0.0571041961)), 0.0571041961, (-0.861136312)}));
        
        weight_[4] = ((0.1739274226) * (0.2034645680)) * ((2. * 0.1739274226));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((-0.861136312 + 1.) / 2.) * (1. - (0.2768430136)), 0.2768430136, (-0.861136312)}));
        
        weight_[5] = ((0.3260725774) * (0.2034645680)) * ((2. * 0.1739274226));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((-0.3399810436 + 1.) / 2.) * (1. - (0.2768430136)), 0.2768430136, (-0.861136312)}));
        
        weight_[6] = ((0.3260725774) * (0.2034645680)) * ((2. * 0.1739274226));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((+0.3399810436 + 1.) / 2.) * (1. - (0.2768430136)), 0.2768430136, (-0.861136312)}));
        
        weight_[7] = ((0.1739274226) * (0.2034645680)) * ((2. * 0.1739274226));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((+0.861136312 + 1.) / 2.) * (1. - (0.2768430136)), 0.2768430136, (-0.861136312)}));
        
        weight_[8] = ((0.1739274226) * (0.1298475476)) * ((2. * 0.1739274226));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((-0.861136312 + 1.) / 2.) * (1. - (0.5835904324)), 0.5835904324, (-0.861136312)}));
        
        weight_[9] = ((0.3260725774) * (0.1298475476)) * ((2. * 0.1739274226));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((-0.3399810436 + 1.) / 2.) * (1. - (0.5835904324)), 0.5835904324, (-0.861136312)}));
        
        weight_[10] = ((0.3260725774) * (0.1298475476)) * ((2. * 0.1739274226));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((+0.3399810436 + 1.) / 2.) * (1. - (0.5835904324)), 0.5835904324, (-0.861136312)}));
        
        weight_[11] = ((0.1739274226) * (0.1298475476)) * ((2. * 0.1739274226));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((+0.861136312 + 1.) / 2.) * (1. - (0.5835904324)), 0.5835904324, (-0.861136312)}));
        
        weight_[12] = ((0.1739274226) * (0.0311809709)) * ((2. * 0.1739274226));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((-0.861136312 + 1.) / 2.) * (1. - (0.8602401357)), 0.8602401357, (-0.861136312)}));
        
        weight_[13] = ((0.3260725774) * (0.0311809709)) * ((2. * 0.1739274226));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((-0.3399810436 + 1.) / 2.) * (1. - (0.8602401357)), 0.8602401357, (-0.861136312)}));
        
        weight_[14] = ((0.3260725774) * (0.0311809709)) * ((2. * 0.1739274226));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((+0.3399810436 + 1.) / 2.) * (1. - (0.8602401357)), 0.8602401357, (-0.861136312)}));
        
        weight_[15] = ((0.1739274226) * (0.0311809709)) * ((2. * 0.1739274226));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((+0.861136312 + 1.) / 2.) * (1. - (0.8602401357)), 0.8602401357, (-0.861136312)}));
        
        weight_[16] = ((0.1739274226) * (0.1355069134)) * ((2. * 0.3260725774));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((-0.861136312 + 1.) / 2.) * (1. - (0.0571041961)), 0.0571041961, (-0.3399810436)}));
        
        weight_[17] = ((0.3260725774) * (0.1355069134)) * ((2. * 0.3260725774));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((-0.3399810436 + 1.) / 2.) * (1. - (0.0571041961)), 0.0571041961, (-0.3399810436)}));
        
        weight_[18] = ((0.3260725774) * (0.1355069134)) * ((2. * 0.3260725774));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((+0.3399810436 + 1.) / 2.) * (1. - (0.0571041961)), 0.0571041961, (-0.3399810436)}));
        
        weight_[19] = ((0.1739274226) * (0.1355069134)) * ((2. * 0.3260725774));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((+0.861136312 + 1.) / 2.) * (1. - (0.0571041961)), 0.0571041961, (-0.3399810436)}));
        
        weight_[20] = ((0.1739274226) * (0.2034645680)) * ((2. * 0.3260725774));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((-0.861136312 + 1.) / 2.) * (1. - (0.2768430136)), 0.2768430136, (-0.3399810436)}));
        
        weight_[21] = ((0.3260725774) * (0.2034645680)) * ((2. * 0.3260725774));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((-0.3399810436 + 1.) / 2.) * (1. - (0.2768430136)), 0.2768430136, (-0.3399810436)}));
        
        weight_[22] = ((0.3260725774) * (0.2034645680)) * ((2. * 0.3260725774));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((+0.3399810436 + 1.) / 2.) * (1. - (0.2768430136)), 0.2768430136, (-0.3399810436)}));
        
        weight_[23] = ((0.1739274226) * (0.2034645680)) * ((2. * 0.3260725774));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((+0.861136312 + 1.) / 2.) * (1. - (0.2768430136)), 0.2768430136, (-0.3399810436)}));
        
        weight_[24] = ((0.1739274226) * (0.1298475476)) * ((2. * 0.3260725774));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((-0.861136312 + 1.) / 2.) * (1. - (0.5835904324)), 0.5835904324, (-0.3399810436)}));
        
        weight_[25] = ((0.3260725774) * (0.1298475476)) * ((2. * 0.3260725774));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((-0.3399810436 + 1.) / 2.) * (1. - (0.5835904324)), 0.5835904324, (-0.3399810436)}));
        
        weight_[26] = ((0.3260725774) * (0.1298475476)) * ((2. * 0.3260725774));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((+0.3399810436 + 1.) / 2.) * (1. - (0.5835904324)), 0.5835904324, (-0.3399810436)}));
        
        weight_[27] = ((0.1739274226) * (0.1298475476)) * ((2. * 0.3260725774));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((+0.861136312 + 1.) / 2.) * (1. - (0.5835904324)), 0.5835904324, (-0.3399810436)}));
        
        weight_[28] = ((0.1739274226) * (0.0311809709)) * ((2. * 0.3260725774));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((-0.861136312 + 1.) / 2.) * (1. - (0.8602401357)), 0.8602401357, (-0.3399810436)}));
        
        weight_[29] = ((0.3260725774) * (0.0311809709)) * ((2. * 0.3260725774));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((-0.3399810436 + 1.) / 2.) * (1. - (0.8602401357)), 0.8602401357, (-0.3399810436)}));
        
        weight_[30] = ((0.3260725774) * (0.0311809709)) * ((2. * 0.3260725774));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((+0.3399810436 + 1.) / 2.) * (1. - (0.8602401357)), 0.8602401357, (-0.3399810436)}));
        
        weight_[31] = ((0.1739274226) * (0.0311809709)) * ((2. * 0.3260725774));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((+0.861136312 + 1.) / 2.) * (1. - (0.8602401357)), 0.8602401357, (-0.3399810436)}));
        
        weight_[32] = ((0.1739274226) * (0.1355069134)) * ((2. * 0.3260725774));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((-0.861136312 + 1.) / 2.) * (1. - (0.0571041961)), 0.0571041961, (+0.3399810436)}));
        
        weight_[33] = ((0.3260725774) * (0.1355069134)) * ((2. * 0.3260725774));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((-0.3399810436 + 1.) / 2.) * (1. - (0.0571041961)), 0.0571041961, (+0.3399810436)}));
        
        weight_[34] = ((0.3260725774) * (0.1355069134)) * ((2. * 0.3260725774));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((+0.3399810436 + 1.) / 2.) * (1. - (0.0571041961)), 0.0571041961, (+0.3399810436)}));
        
        weight_[35] = ((0.1739274226) * (0.1355069134)) * ((2. * 0.3260725774));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((+0.861136312 + 1.) / 2.) * (1. - (0.0571041961)), 0.0571041961, (+0.3399810436)}));
        
        weight_[36] = ((0.1739274226) * (0.2034645680)) * ((2. * 0.3260725774));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((-0.861136312 + 1.) / 2.) * (1. - (0.2768430136)), 0.2768430136, (+0.3399810436)}));
        
        weight_[37] = ((0.3260725774) * (0.2034645680)) * ((2. * 0.3260725774));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((-0.3399810436 + 1.) / 2.) * (1. - (0.2768430136)), 0.2768430136, (+0.3399810436)}));
        
        weight_[38] = ((0.3260725774) * (0.2034645680)) * ((2. * 0.3260725774));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((+0.3399810436 + 1.) / 2.) * (1. - (0.2768430136)), 0.2768430136, (+0.3399810436)}));
        
        weight_[39] = ((0.1739274226) * (0.2034645680)) * ((2. * 0.3260725774));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((+0.861136312 + 1.) / 2.) * (1. - (0.2768430136)), 0.2768430136, (+0.3399810436)}));
        
        weight_[40] = ((0.1739274226) * (0.1298475476)) * ((2. * 0.3260725774));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((-0.861136312 + 1.) / 2.) * (1. - (0.5835904324)), 0.5835904324, (+0.3399810436)}));
        
        weight_[41] = ((0.3260725774) * (0.1298475476)) * ((2. * 0.3260725774));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((-0.3399810436 + 1.) / 2.) * (1. - (0.5835904324)), 0.5835904324, (+0.3399810436)}));
        
        weight_[42] = ((0.3260725774) * (0.1298475476)) * ((2. * 0.3260725774));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((+0.3399810436 + 1.) / 2.) * (1. - (0.5835904324)), 0.5835904324, (+0.3399810436)}));
        
        weight_[43] = ((0.1739274226) * (0.1298475476)) * ((2. * 0.3260725774));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((+0.861136312 + 1.) / 2.) * (1. - (0.5835904324)), 0.5835904324, (+0.3399810436)}));
        
        weight_[44] = ((0.1739274226) * (0.0311809709)) * ((2. * 0.3260725774));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((-0.861136312 + 1.) / 2.) * (1. - (0.8602401357)), 0.8602401357, (+0.3399810436)}));
        
        weight_[45] = ((0.3260725774) * (0.0311809709)) * ((2. * 0.3260725774));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((-0.3399810436 + 1.) / 2.) * (1. - (0.8602401357)), 0.8602401357, (+0.3399810436)}));
        
        weight_[46] = ((0.3260725774) * (0.0311809709)) * ((2. * 0.3260725774));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((+0.3399810436 + 1.) / 2.) * (1. - (0.8602401357)), 0.8602401357, (+0.3399810436)}));
        
        weight_[47] = ((0.1739274226) * (0.0311809709)) * ((2. * 0.3260725774));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((+0.861136312 + 1.) / 2.) * (1. - (0.8602401357)), 0.8602401357, (+0.3399810436)}));
        
        weight_[48] = ((0.1739274226) * (0.1355069134)) * ((2. * 0.1739274226));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((-0.861136312 + 1.) / 2.) * (1. - (0.0571041961)), 0.0571041961, (+0.861136312)}));
        
        weight_[49] = ((0.3260725774) * (0.1355069134)) * ((2. * 0.1739274226));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((-0.3399810436 + 1.) / 2.) * (1. - (0.0571041961)), 0.0571041961, (+0.861136312)}));
        
        weight_[50] = ((0.3260725774) * (0.1355069134)) * ((2. * 0.1739274226));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((+0.3399810436 + 1.) / 2.) * (1. - (0.0571041961)), 0.0571041961, (+0.861136312)}));
        
        weight_[51] = ((0.1739274226) * (0.1355069134)) * ((2. * 0.1739274226));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((+0.861136312 + 1.) / 2.) * (1. - (0.0571041961)), 0.0571041961, (+0.861136312)}));
        
        weight_[52] = ((0.1739274226) * (0.2034645680)) * ((2. * 0.1739274226));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((-0.861136312 + 1.) / 2.) * (1. - (0.2768430136)), 0.2768430136, (+0.861136312)}));
        
        weight_[53] = ((0.3260725774) * (0.2034645680)) * ((2. * 0.1739274226));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((-0.3399810436 + 1.) / 2.) * (1. - (0.2768430136)), 0.2768430136, (+0.861136312)}));
        
        weight_[54] = ((0.3260725774) * (0.2034645680)) * ((2. * 0.1739274226));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((+0.3399810436 + 1.) / 2.) * (1. - (0.2768430136)), 0.2768430136, (+0.861136312)}));
        
        weight_[55] = ((0.1739274226) * (0.2034645680)) * ((2. * 0.1739274226));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((+0.861136312 + 1.) / 2.) * (1. - (0.2768430136)), 0.2768430136, (+0.861136312)}));
        
        weight_[56] = ((0.1739274226) * (0.1298475476)) * ((2. * 0.1739274226));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((-0.861136312 + 1.) / 2.) * (1. - (0.5835904324)), 0.5835904324, (+0.861136312)}));
        
        weight_[57] = ((0.3260725774) * (0.1298475476)) * ((2. * 0.1739274226));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((-0.3399810436 + 1.) / 2.) * (1. - (0.5835904324)), 0.5835904324, (+0.861136312)}));
        
        weight_[58] = ((0.3260725774) * (0.1298475476)) * ((2. * 0.1739274226));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((+0.3399810436 + 1.) / 2.) * (1. - (0.5835904324)), 0.5835904324, (+0.861136312)}));
        
        weight_[59] = ((0.1739274226) * (0.1298475476)) * ((2. * 0.1739274226));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((+0.861136312 + 1.) / 2.) * (1. - (0.5835904324)), 0.5835904324, (+0.861136312)}));
        
        weight_[60] = ((0.1739274226) * (0.0311809709)) * ((2. * 0.1739274226));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((-0.861136312 + 1.) / 2.) * (1. - (0.8602401357)), 0.8602401357, (+0.861136312)}));
        
        weight_[61] = ((0.3260725774) * (0.0311809709)) * ((2. * 0.1739274226));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((-0.3399810436 + 1.) / 2.) * (1. - (0.8602401357)), 0.8602401357, (+0.861136312)}));
        
        weight_[62] = ((0.3260725774) * (0.0311809709)) * ((2. * 0.1739274226));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((+0.3399810436 + 1.) / 2.) * (1. - (0.8602401357)), 0.8602401357, (+0.861136312)}));
        
        weight_[63] = ((0.1739274226) * (0.0311809709)) * ((2. * 0.1739274226));
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({((+0.861136312 + 1.) / 2.) * (1. - (0.8602401357)), 0.8602401357, (+0.861136312)}));
        
    }

//---------------------------------------------------------------------------
}// close namespace IntegrationRules

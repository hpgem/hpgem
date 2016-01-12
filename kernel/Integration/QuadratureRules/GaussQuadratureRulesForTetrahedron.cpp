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
#include "Integration/QuadratureRules/GaussQuadratureRulesForTetrahedron.h"
#include "Geometry/ReferenceTetrahedron.h"
#include "Geometry/PointReference.h"

//---------------------------------------------------------------------------
namespace QuadratureRules
{
//---------------------------------------------------------------------------
    std::string Tn3_1_1::getName() const
    {
        return name_;
    }
    
    std::size_t Tn3_1_1::order() const
    {
        return 1;
    }
    
    std::size_t Tn3_1_1::dimension() const
    {
        return 3;
    }
    
    std::size_t Tn3_1_1::getNumberOfPoints() const
    {
        return 1;
    }
    
    double Tn3_1_1::weight(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReferenceBase&
    Tn3_1_1::getPoint(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }
    
    Geometry::ReferenceGeometry*
    Tn3_1_1::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    Tn3_1_1::Tn3_1_1()
            : name_("Tn3_1_1"), refGeoPtr_(&Geometry::ReferenceTetrahedron::Instance()), gp_(0)
    {
        weight_[0] = 1.6666666666666663e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.5000000000000000e-01, 2.5000000000000000e-01, 2.5000000000000000e-01}));
        
    }
    
//---------------------------------------------------------------------------
    std::string Tn3_2_4::getName() const
    {
        return name_;
    }
    
    std::size_t Tn3_2_4::order() const
    {
        return 2;
    }
    
    std::size_t Tn3_2_4::dimension() const
    {
        return 3;
    }
    
    std::size_t Tn3_2_4::getNumberOfPoints() const
    {
        return 4;
    }
    
    double Tn3_2_4::weight(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReferenceBase&
    Tn3_2_4::getPoint(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }
    
    Geometry::ReferenceGeometry*
    Tn3_2_4::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    Tn3_2_4::Tn3_2_4()
            : name_("Tn3_2_1"), refGeoPtr_(&Geometry::ReferenceTetrahedron::Instance()), gp_(0)
    {
        weight_[0] = 4.1666666666666623e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.3819660112501048e-01, 1.3819660112501048e-01, 1.3819660112501048e-01}));
        
        weight_[1] = 4.1666666666666623e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({5.8541019662496852e-01, 1.3819660112501048e-01, 1.3819660112501048e-01}));
        
        weight_[2] = 4.1666666666666623e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.3819660112501048e-01, 5.8541019662496852e-01, 1.3819660112501048e-01}));
        
        weight_[3] = 4.1666666666666623e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.3819660112501048e-01, 1.3819660112501048e-01, 5.8541019662496852e-01}));
        
    }
    
//---------------------------------------------------------------------------
    std::string Tn3_3_5::getName() const
    {
        return name_;
    }
    
    std::size_t Tn3_3_5::order() const
    {
        return 3;
    }
    
    std::size_t Tn3_3_5::dimension() const
    {
        return 3;
    }
    
    std::size_t Tn3_3_5::getNumberOfPoints() const
    {
        return 5;
    }
    
    double Tn3_3_5::weight(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReferenceBase&
    Tn3_3_5::getPoint(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }
    
    Geometry::ReferenceGeometry*
    Tn3_3_5::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    Tn3_3_5::Tn3_3_5()
            : name_("Tn3_3_1"), refGeoPtr_(&Geometry::ReferenceTetrahedron::Instance()), gp_(0)
    {
        weight_[0] = -1.3333333333333339e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.5000000000000000e-01, 2.5000000000000000e-01, 2.5000000000000000e-01}));
        
        weight_[1] = 7.4999999999999997e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.6666666666666652e-01, 1.6666666666666652e-01, 1.6666666666666652e-01}));
        
        weight_[2] = 7.4999999999999997e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.6666666666666652e-01, 1.6666666666666652e-01, 5.0000000000000000e-01}));
        
        weight_[3] = 7.4999999999999997e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.6666666666666652e-01, 5.0000000000000000e-01, 1.6666666666666652e-01}));
        
        weight_[4] = 7.4999999999999997e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({5.0000000000000000e-01, 1.6666666666666652e-01, 1.6666666666666652e-01}));
        
    }
    
//---------------------------------------------------------------------------
    std::string Tn3_4_11::getName() const
    {
        return name_;
    }
    
    std::size_t Tn3_4_11::order() const
    {
        return 4;
    }
    
    std::size_t Tn3_4_11::dimension() const
    {
        return 3;
    }
    
    std::size_t Tn3_4_11::getNumberOfPoints() const
    {
        return 11;
    }
    
    double Tn3_4_11::weight(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReferenceBase&
    Tn3_4_11::getPoint(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }
    
    Geometry::ReferenceGeometry*
    Tn3_4_11::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    Tn3_4_11::Tn3_4_11()
            : name_("Tn3_4_1"), refGeoPtr_(&Geometry::ReferenceTetrahedron::Instance()), gp_(0)
    {
        weight_[0] = -1.3155555555555500e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.5000000000000000e-01, 2.5000000000000000e-01, 2.5000000000000000e-01}));
        
        weight_[1] = 7.6222222222222498e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.1428571428571508e-02, 7.1428571428571508e-02, 7.1428571428571508e-02}));
        
        weight_[2] = 7.6222222222222498e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.1428571428571508e-02, 7.1428571428571508e-02, 7.8571428571428548e-01}));
        
        weight_[3] = 7.6222222222222498e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.1428571428571508e-02, 7.8571428571428548e-01, 7.1428571428571508e-02}));
        
        weight_[4] = 7.6222222222222498e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.8571428571428548e-01, 7.1428571428571508e-02, 7.1428571428571508e-02}));
        
        weight_[5] = 2.4888888888888874e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.9940357616679900e-01, 3.9940357616679900e-01, 1.0059642383320100e-01}));
        
        weight_[6] = 2.4888888888888874e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.9940357616679900e-01, 1.0059642383320100e-01, 3.9940357616679900e-01}));
        
        weight_[7] = 2.4888888888888874e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.0059642383320100e-01, 3.9940357616679900e-01, 3.9940357616679900e-01}));
        
        weight_[8] = 2.4888888888888874e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.9940357616679900e-01, 1.0059642383320100e-01, 1.0059642383320100e-01}));
        
        weight_[9] = 2.4888888888888874e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.0059642383320100e-01, 3.9940357616679900e-01, 1.0059642383320100e-01}));
        
        weight_[10] = 2.4888888888888874e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.0059642383320100e-01, 1.0059642383320100e-01, 3.9940357616679900e-01}));
        
    }
    
//---------------------------------------------------------------------------
    std::string T3_5_14::getName() const
    {
        return name_;
    }
    
    std::size_t T3_5_14::order() const
    {
        return 5;
    }
    
    std::size_t T3_5_14::dimension() const
    {
        return 3;
    }
    
    std::size_t T3_5_14::getNumberOfPoints() const
    {
        return 14;
    }
    
    double T3_5_14::weight(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReferenceBase&
    T3_5_14::getPoint(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }
    
    Geometry::ReferenceGeometry*
    T3_5_14::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    T3_5_14::T3_5_14()
            : name_("T3_5_1"), refGeoPtr_(&Geometry::ReferenceTetrahedron::Instance()), gp_(0)
    {
        weight_[0] = 1.2248840519393626e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({9.2735250310891026e-02, 9.2735250310891026e-02, 9.2735250310891026e-02}));
        
        weight_[1] = 1.2248840519393626e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.2179424906732648e-01, 9.2735250310891026e-02, 9.2735250310891026e-02}));
        
        weight_[2] = 1.2248840519393626e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({9.2735250310891026e-02, 7.2179424906732648e-01, 9.2735250310891026e-02}));
        
        weight_[3] = 1.2248840519393626e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({9.2735250310891026e-02, 9.2735250310891026e-02, 7.2179424906732648e-01}));
        
        weight_[4] = 1.8781320953002625e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.1088591926330050e-01, 3.1088591926330050e-01, 3.1088591926330050e-01}));
        
        weight_[5] = 1.8781320953002625e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({6.7342242210097991e-02, 3.1088591926330050e-01, 3.1088591926330050e-01}));
        
        weight_[6] = 1.8781320953002625e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.1088591926330050e-01, 6.7342242210097991e-02, 3.1088591926330050e-01}));
        
        weight_[7] = 1.8781320953002625e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.1088591926330050e-01, 3.1088591926330050e-01, 6.7342242210097991e-02}));
        
        weight_[8] = 7.0910034628468748e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({4.5449629587435048e-01, 4.5449629587435048e-01, 4.5503704125649524e-02}));
        
        weight_[9] = 7.0910034628468748e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({4.5449629587435048e-01, 4.5503704125649524e-02, 4.5449629587435048e-01}));
        
        weight_[10] = 7.0910034628468748e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({4.5503704125649524e-02, 4.5449629587435048e-01, 4.5449629587435048e-01}));
        
        weight_[11] = 7.0910034628468748e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({4.5449629587435048e-01, 4.5503704125649524e-02, 4.5503704125649524e-02}));
        
        weight_[12] = 7.0910034628468748e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({4.5503704125649524e-02, 4.5449629587435048e-01, 4.5503704125649524e-02}));
        
        weight_[13] = 7.0910034628468748e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({4.5503704125649524e-02, 4.5503704125649524e-02, 4.5449629587435048e-01}));
        
    }
    
//---------------------------------------------------------------------------
    std::string T3_6_24::getName() const
    {
        return name_;
    }
    
    std::size_t T3_6_24::order() const
    {
        return 6;
    }
    
    std::size_t T3_6_24::dimension() const
    {
        return 3;
    }
    
    std::size_t T3_6_24::getNumberOfPoints() const
    {
        return 24;
    }
    
    double T3_6_24::weight(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReferenceBase&
    T3_6_24::getPoint(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }
    
    Geometry::ReferenceGeometry*
    T3_6_24::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    T3_6_24::T3_6_24()
            : name_("T3_6_1"), refGeoPtr_(&Geometry::ReferenceTetrahedron::Instance()), gp_(0)
    {
        weight_[0] = 6.6537917096946252e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.1460287125915201e-01, 2.1460287125915201e-01, 2.1460287125915201e-01}));
        
        weight_[1] = 6.6537917096946252e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.5619138622254398e-01, 2.1460287125915201e-01, 2.1460287125915201e-01}));
        
        weight_[2] = 6.6537917096946252e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.1460287125915201e-01, 3.5619138622254398e-01, 2.1460287125915201e-01}));
        
        weight_[3] = 6.6537917096946252e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.1460287125915201e-01, 2.1460287125915201e-01, 3.5619138622254398e-01}));
        
        weight_[4] = 1.6795351758867501e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({4.0673958534611476e-02, 4.0673958534611476e-02, 4.0673958534611476e-02}));
        
        weight_[5] = 1.6795351758867501e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({8.7797812439616596e-01, 4.0673958534611476e-02, 4.0673958534611476e-02}));
        
        weight_[6] = 1.6795351758867501e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({4.0673958534611476e-02, 8.7797812439616596e-01, 4.0673958534611476e-02}));
        
        weight_[7] = 1.6795351758867501e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({4.0673958534611476e-02, 4.0673958534611476e-02, 8.7797812439616596e-01}));
        
        weight_[8] = 9.2261969239424996e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.2233789014227554e-01, 3.2233789014227554e-01, 3.2233789014227554e-01}));
        
        weight_[9] = 9.2261969239424996e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.2986329573173490e-02, 3.2233789014227554e-01, 3.2233789014227554e-01}));
        
        weight_[10] = 9.2261969239424996e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.2233789014227554e-01, 3.2986329573173490e-02, 3.2233789014227554e-01}));
        
        weight_[11] = 9.2261969239424996e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.2233789014227554e-01, 3.2233789014227554e-01, 3.2986329573173490e-02}));
        
        weight_[12] = 8.0357142857142502e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({6.3661001875017498e-02, 6.3661001875017498e-02, 2.6967233145831604e-01}));
        
        weight_[13] = 8.0357142857142502e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({6.3661001875017498e-02, 2.6967233145831604e-01, 6.3661001875017498e-02}));
        
        weight_[14] = 8.0357142857142502e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({6.3661001875017498e-02, 6.3661001875017498e-02, 6.0300566479164897e-01}));
        
        weight_[15] = 8.0357142857142502e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({6.3661001875017498e-02, 6.0300566479164897e-01, 6.3661001875017498e-02}));
        
        weight_[16] = 8.0357142857142502e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({6.3661001875017498e-02, 2.6967233145831604e-01, 6.0300566479164897e-01}));
        
        weight_[17] = 8.0357142857142502e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({6.3661001875017498e-02, 6.0300566479164897e-01, 2.6967233145831604e-01}));
        
        weight_[18] = 8.0357142857142502e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.6967233145831604e-01, 6.3661001875017498e-02, 6.3661001875017498e-02}));
        
        weight_[19] = 8.0357142857142502e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.6967233145831604e-01, 6.3661001875017498e-02, 6.0300566479164897e-01}));
        
        weight_[20] = 8.0357142857142502e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.6967233145831604e-01, 6.0300566479164897e-01, 6.3661001875017498e-02}));
        
        weight_[21] = 8.0357142857142502e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({6.0300566479164897e-01, 6.3661001875017498e-02, 2.6967233145831604e-01}));
        
        weight_[22] = 8.0357142857142502e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({6.0300566479164897e-01, 6.3661001875017498e-02, 6.3661001875017498e-02}));
        
        weight_[23] = 8.0357142857142502e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({6.0300566479164897e-01, 2.6967233145831604e-01, 6.3661001875017498e-02}));
        
    }
    
//---------------------------------------------------------------------------
    std::string T3_7_31::getName() const
    {
        return name_;
    }
    
    std::size_t T3_7_31::order() const
    {
        return 7;
    }
    
    std::size_t T3_7_31::dimension() const
    {
        return 3;
    }
    
    std::size_t T3_7_31::getNumberOfPoints() const
    {
        return 31;
    }
    
    double T3_7_31::weight(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReferenceBase&
    T3_7_31::getPoint(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }
    
    Geometry::ReferenceGeometry*
    T3_7_31::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    T3_7_31::T3_7_31()
            : name_("T3_7_1"), refGeoPtr_(&Geometry::ReferenceTetrahedron::Instance()), gp_(0)
    {
        weight_[0] = 9.7001763668425002e-04;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({5.0000000000000000e-01, 5.0000000000000000e-01, 0.0000000000000000e+00}));
        
        weight_[1] = 9.7001763668425002e-04;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({5.0000000000000000e-01, 0.0000000000000000e+00, 5.0000000000000000e-01}));
        
        weight_[2] = 9.7001763668425002e-04;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({0.0000000000000000e+00, 5.0000000000000000e-01, 5.0000000000000000e-01}));
        
        weight_[3] = 9.7001763668425002e-04;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({0.0000000000000000e+00, 0.0000000000000000e+00, 5.0000000000000000e-01}));
        
        weight_[4] = 9.7001763668425002e-04;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({0.0000000000000000e+00, 5.0000000000000000e-01, 0.0000000000000000e+00}));
        
        weight_[5] = 9.7001763668425002e-04;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({5.0000000000000000e-01, 0.0000000000000000e+00, 0.0000000000000000e+00}));
        
        weight_[6] = 1.8264223466108873e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.5000000000000000e-01, 2.5000000000000000e-01, 2.5000000000000000e-01}));
        
        weight_[7] = 1.0599941524413625e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.8213192330317982e-02, 7.8213192330317982e-02, 7.8213192330317982e-02}));
        
        weight_[8] = 1.0599941524413625e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.8213192330317982e-02, 7.8213192330317982e-02, 7.6536042300904605e-01}));
        
        weight_[9] = 1.0599941524413625e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.8213192330317982e-02, 7.6536042300904605e-01, 7.8213192330317982e-02}));
        
        weight_[10] = 1.0599941524413625e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.6536042300904605e-01, 7.8213192330317982e-02, 7.8213192330317982e-02}));
        
        weight_[11] = -6.2517740114331879e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.2184321666390502e-01, 1.2184321666390502e-01, 1.2184321666390502e-01}));
        
        weight_[12] = -6.2517740114331879e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.2184321666390502e-01, 1.2184321666390502e-01, 6.3447035000828444e-01}));
        
        weight_[13] = -6.2517740114331879e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.2184321666390502e-01, 6.3447035000828444e-01, 1.2184321666390502e-01}));
        
        weight_[14] = -6.2517740114331879e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({6.3447035000828444e-01, 1.2184321666390502e-01, 1.2184321666390502e-01}));
        
        weight_[15] = 4.8914252630735001e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.3253916444642051e-01, 3.3253916444642051e-01, 3.3253916444642051e-01}));
        
        weight_[16] = 4.8914252630735001e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.3253916444642051e-01, 3.3253916444642051e-01, 2.3825066607380263e-03}));
        
        weight_[17] = 4.8914252630735001e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.3253916444642051e-01, 2.3825066607380263e-03, 3.3253916444642051e-01}));
        
        weight_[18] = 4.8914252630735001e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.3825066607380263e-03, 3.3253916444642051e-01, 3.3253916444642051e-01}));
        
        weight_[19] = 2.7557319223985875e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({9.9999999999999978e-02, 9.9999999999999978e-02, 2.0000000000000001e-01}));
        
        weight_[20] = 2.7557319223985875e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({9.9999999999999978e-02, 2.0000000000000001e-01, 9.9999999999999978e-02}));
        
        weight_[21] = 2.7557319223985875e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({9.9999999999999978e-02, 9.9999999999999978e-02, 5.9999999999999998e-01}));
        
        weight_[22] = 2.7557319223985875e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({9.9999999999999978e-02, 5.9999999999999998e-01, 9.9999999999999978e-02}));
        
        weight_[23] = 2.7557319223985875e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({9.9999999999999978e-02, 2.0000000000000001e-01, 5.9999999999999998e-01}));
        
        weight_[24] = 2.7557319223985875e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({9.9999999999999978e-02, 5.9999999999999998e-01, 2.0000000000000001e-01}));
        
        weight_[25] = 2.7557319223985875e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.0000000000000001e-01, 9.9999999999999978e-02, 9.9999999999999978e-02}));
        
        weight_[26] = 2.7557319223985875e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.0000000000000001e-01, 9.9999999999999978e-02, 5.9999999999999998e-01}));
        
        weight_[27] = 2.7557319223985875e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.0000000000000001e-01, 5.9999999999999998e-01, 9.9999999999999978e-02}));
        
        weight_[28] = 2.7557319223985875e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({5.9999999999999998e-01, 9.9999999999999978e-02, 2.0000000000000001e-01}));
        
        weight_[29] = 2.7557319223985875e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({5.9999999999999998e-01, 9.9999999999999978e-02, 9.9999999999999978e-02}));
        
        weight_[30] = 2.7557319223985875e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({5.9999999999999998e-01, 2.0000000000000001e-01, 9.9999999999999978e-02}));
        
    }
    
//---------------------------------------------------------------------------
    std::string T3_8_43::getName() const
    {
        return name_;
    }
    
    std::size_t T3_8_43::order() const
    {
        return 8;
    }
    
    std::size_t T3_8_43::dimension() const
    {
        return 3;
    }
    
    std::size_t T3_8_43::getNumberOfPoints() const
    {
        return 43;
    }
    
    double T3_8_43::weight(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReferenceBase&
    T3_8_43::getPoint(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }
    
    Geometry::ReferenceGeometry*
    T3_8_43::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    T3_8_43::T3_8_43()
            : name_("T3_8_1"), refGeoPtr_(&Geometry::ReferenceTetrahedron::Instance()), gp_(0)
    {
        weight_[0] = -2.0500188658639874e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.5000000000000000e-01, 2.5000000000000000e-01, 2.5000000000000000e-01}));
        
        weight_[1] = 1.4250305822866875e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.0682993161067298e-01, 2.0682993161067298e-01, 2.0682993161067298e-01}));
        
        weight_[2] = 1.4250305822866875e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.0682993161067298e-01, 2.0682993161067298e-01, 3.7951020516798051e-01}));
        
        weight_[3] = 1.4250305822866875e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.0682993161067298e-01, 3.7951020516798051e-01, 2.0682993161067298e-01}));
        
        weight_[4] = 1.4250305822866875e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.7951020516798051e-01, 2.0682993161067298e-01, 2.0682993161067298e-01}));
        
        weight_[5] = 1.9670333131338751e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({8.2103588310546483e-02, 8.2103588310546483e-02, 8.2103588310546483e-02}));
        
        weight_[6] = 1.9670333131338751e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({8.2103588310546483e-02, 8.2103588310546483e-02, 7.5368923506835994e-01}));
        
        weight_[7] = 1.9670333131338751e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({8.2103588310546483e-02, 7.5368923506835994e-01, 8.2103588310546483e-02}));
        
        weight_[8] = 1.9670333131338751e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.5368923506835994e-01, 8.2103588310546483e-02, 8.2103588310546483e-02}));
        
        weight_[9] = 1.6983410909287501e-04;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({5.7819505051979747e-03, 5.7819505051979747e-03, 5.7819505051979747e-03}));
        
        weight_[10] = 1.6983410909287501e-04;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({5.7819505051979747e-03, 5.7819505051979747e-03, 9.8265414848440602e-01}));
        
        weight_[11] = 1.6983410909287501e-04;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({5.7819505051979747e-03, 9.8265414848440602e-01, 5.7819505051979747e-03}));
        
        weight_[12] = 1.6983410909287501e-04;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({9.8265414848440602e-01, 5.7819505051979747e-03, 5.7819505051979747e-03}));
        
        weight_[13] = 4.5796838244672499e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({5.0532740018894018e-02, 5.0532740018894018e-02, 4.4946725998110598e-01}));
        
        weight_[14] = 4.5796838244672499e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({5.0532740018894018e-02, 4.4946725998110598e-01, 5.0532740018894018e-02}));
        
        weight_[15] = 4.5796838244672499e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({4.4946725998110598e-01, 5.0532740018894018e-02, 5.0532740018894018e-02}));
        
        weight_[16] = 4.5796838244672499e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({5.0532740018894018e-02, 4.4946725998110598e-01, 4.4946725998110598e-01}));
        
        weight_[17] = 4.5796838244672499e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({4.4946725998110598e-01, 5.0532740018894018e-02, 4.4946725998110598e-01}));
        
        weight_[18] = 4.5796838244672499e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({4.4946725998110598e-01, 4.4946725998110598e-01, 5.0532740018894018e-02}));
        
        weight_[19] = 5.7044858086818754e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.2906653611681099e-01, 2.2906653611681099e-01, 3.5639582788534019e-02}));
        
        weight_[20] = 5.7044858086818754e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.2906653611681099e-01, 3.5639582788534019e-02, 2.2906653611681099e-01}));
        
        weight_[21] = 5.7044858086818754e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.2906653611681099e-01, 2.2906653611681099e-01, 5.0622734497784350e-01}));
        
        weight_[22] = 5.7044858086818754e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.2906653611681099e-01, 5.0622734497784350e-01, 2.2906653611681099e-01}));
        
        weight_[23] = 5.7044858086818754e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.2906653611681099e-01, 3.5639582788534019e-02, 5.0622734497784350e-01}));
        
        weight_[24] = 5.7044858086818754e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.2906653611681099e-01, 5.0622734497784350e-01, 3.5639582788534019e-02}));
        
        weight_[25] = 5.7044858086818754e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.5639582788534019e-02, 2.2906653611681099e-01, 2.2906653611681099e-01}));
        
        weight_[26] = 5.7044858086818754e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.5639582788534019e-02, 2.2906653611681099e-01, 5.0622734497784350e-01}));
        
        weight_[27] = 5.7044858086818754e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.5639582788534019e-02, 5.0622734497784350e-01, 2.2906653611681099e-01}));
        
        weight_[28] = 5.7044858086818754e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({5.0622734497784350e-01, 2.2906653611681099e-01, 3.5639582788534019e-02}));
        
        weight_[29] = 5.7044858086818754e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({5.0622734497784350e-01, 2.2906653611681099e-01, 2.2906653611681099e-01}));
        
        weight_[30] = 5.7044858086818754e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({5.0622734497784350e-01, 3.5639582788534019e-02, 2.2906653611681099e-01}));
        
        weight_[31] = 2.1405191411621250e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.6607749553197511e-02, 3.6607749553197511e-02, 1.9048604193463348e-01}));
        
        weight_[32] = 2.1405191411621250e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.6607749553197511e-02, 1.9048604193463348e-01, 3.6607749553197511e-02}));
        
        weight_[33] = 2.1405191411621250e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.6607749553197511e-02, 3.6607749553197511e-02, 7.3629845895897150e-01}));
        
        weight_[34] = 2.1405191411621250e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.6607749553197511e-02, 7.3629845895897150e-01, 3.6607749553197511e-02}));
        
        weight_[35] = 2.1405191411621250e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.6607749553197511e-02, 1.9048604193463348e-01, 7.3629845895897150e-01}));
        
        weight_[36] = 2.1405191411621250e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.6607749553197511e-02, 7.3629845895897150e-01, 1.9048604193463348e-01}));
        
        weight_[37] = 2.1405191411621250e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.9048604193463348e-01, 3.6607749553197511e-02, 3.6607749553197511e-02}));
        
        weight_[38] = 2.1405191411621250e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.9048604193463348e-01, 3.6607749553197511e-02, 7.3629845895897150e-01}));
        
        weight_[39] = 2.1405191411621250e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.9048604193463348e-01, 7.3629845895897150e-01, 3.6607749553197511e-02}));
        
        weight_[40] = 2.1405191411621250e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.3629845895897150e-01, 3.6607749553197511e-02, 1.9048604193463348e-01}));
        
        weight_[41] = 2.1405191411621250e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.3629845895897150e-01, 3.6607749553197511e-02, 3.6607749553197511e-02}));
        
        weight_[42] = 2.1405191411621250e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.3629845895897150e-01, 1.9048604193463348e-01, 3.6607749553197511e-02}));
        
    }
    
//---------------------------------------------------------------------------
    std::string T3_9_53::getName() const
    {
        return name_;
    }
    
    std::size_t T3_9_53::order() const
    {
        return 9;
    }
    
    std::size_t T3_9_53::dimension() const
    {
        return 3;
    }
    
    std::size_t T3_9_53::getNumberOfPoints() const
    {
        return 53;
    }
    
    double T3_9_53::weight(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReferenceBase&
    T3_9_53::getPoint(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }
    
    Geometry::ReferenceGeometry*
    T3_9_53::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    T3_9_53::T3_9_53()
            : name_("T3_9_1"), refGeoPtr_(&Geometry::ReferenceTetrahedron::Instance()), gp_(0)
    {
        weight_[0] = -1.3779903832610862e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.5000000000000000e-01, 2.5000000000000000e-01, 2.5000000000000000e-01}));
        
        weight_[1] = 1.8653365690852501e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({4.8351038549736991e-02, 4.8351038549736991e-02, 4.8351038549736991e-02}));
        
        weight_[2] = 1.8653365690852501e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({4.8351038549736991e-02, 4.8351038549736991e-02, 8.5494688435079003e-01}));
        
        weight_[3] = 1.8653365690852501e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({4.8351038549736991e-02, 8.5494688435079003e-01, 4.8351038549736991e-02}));
        
        weight_[4] = 1.8653365690852501e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({8.5494688435079003e-01, 4.8351038549736991e-02, 4.8351038549736991e-02}));
        
        weight_[5] = 4.3094239694933751e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.2457928011788251e-01, 3.2457928011788251e-01, 3.2457928011788251e-01}));
        
        weight_[6] = 4.3094239694933751e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.2457928011788251e-01, 3.2457928011788251e-01, 2.6262159646352978e-02}));
        
        weight_[7] = 4.3094239694933751e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.2457928011788251e-01, 2.6262159646352978e-02, 3.2457928011788251e-01}));
        
        weight_[8] = 4.3094239694933751e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.6262159646352978e-02, 3.2457928011788251e-01, 3.2457928011788251e-01}));
        
        weight_[9] = -9.0184766481201495e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.1461654022399498e-01, 1.1461654022399498e-01, 1.1461654022399498e-01}));
        
        weight_[10] = -9.0184766481201495e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.1461654022399498e-01, 1.1461654022399498e-01, 6.5615037932801457e-01}));
        
        weight_[11] = -9.0184766481201495e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.1461654022399498e-01, 6.5615037932801457e-01, 1.1461654022399498e-01}));
        
        weight_[12] = -9.0184766481201495e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({6.5615037932801457e-01, 1.1461654022399498e-01, 1.1461654022399498e-01}));
        
        weight_[13] = 4.4672576202511499e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.2548995191151400e-01, 2.2548995191151400e-01, 2.2548995191151400e-01}));
        
        weight_[14] = 4.4672576202511499e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.2548995191151400e-01, 2.2548995191151400e-01, 3.2353014426545801e-01}));
        
        weight_[15] = 4.4672576202511499e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.2548995191151400e-01, 3.2353014426545801e-01, 2.2548995191151400e-01}));
        
        weight_[16] = 4.4672576202511499e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.2353014426545801e-01, 2.2548995191151400e-01, 2.2548995191151400e-01}));
        
        weight_[17] = 3.4700405884550749e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.3162780924687001e-01, 1.3162780924687001e-01, 8.3664701617184978e-02}));
        
        weight_[18] = 3.4700405884550749e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.3162780924687001e-01, 8.3664701617184978e-02, 1.3162780924687001e-01}));
        
        weight_[19] = 3.4700405884550749e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.3162780924687001e-01, 1.3162780924687001e-01, 6.5307967988907545e-01}));
        
        weight_[20] = 3.4700405884550749e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.3162780924687001e-01, 6.5307967988907545e-01, 1.3162780924687001e-01}));
        
        weight_[21] = 3.4700405884550749e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.3162780924687001e-01, 8.3664701617184978e-02, 6.5307967988907545e-01}));
        
        weight_[22] = 3.4700405884550749e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.3162780924687001e-01, 6.5307967988907545e-01, 8.3664701617184978e-02}));
        
        weight_[23] = 3.4700405884550749e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({8.3664701617184978e-02, 1.3162780924687001e-01, 1.3162780924687001e-01}));
        
        weight_[24] = 3.4700405884550749e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({8.3664701617184978e-02, 1.3162780924687001e-01, 6.5307967988907545e-01}));
        
        weight_[25] = 3.4700405884550749e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({8.3664701617184978e-02, 6.5307967988907545e-01, 1.3162780924687001e-01}));
        
        weight_[26] = 3.4700405884550749e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({6.5307967988907545e-01, 1.3162780924687001e-01, 8.3664701617184978e-02}));
        
        weight_[27] = 3.4700405884550749e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({6.5307967988907545e-01, 1.3162780924687001e-01, 1.3162780924687001e-01}));
        
        weight_[28] = 3.4700405884550749e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({6.5307967988907545e-01, 8.3664701617184978e-02, 1.3162780924687001e-01}));
        
        weight_[29] = 3.3525839026606252e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({4.3395146141140700e-01, 4.3395146141140700e-01, 1.0776985954942853e-01}));
        
        weight_[30] = 3.3525839026606252e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({4.3395146141140700e-01, 1.0776985954942853e-01, 4.3395146141140700e-01}));
        
        weight_[31] = 3.3525839026606252e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({4.3395146141140700e-01, 4.3395146141140700e-01, 2.4327217627758024e-02}));
        
        weight_[32] = 3.3525839026606252e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({4.3395146141140700e-01, 2.4327217627758024e-02, 4.3395146141140700e-01}));
        
        weight_[33] = 3.3525839026606252e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({4.3395146141140700e-01, 1.0776985954942853e-01, 2.4327217627758024e-02}));
        
        weight_[34] = 3.3525839026606252e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({4.3395146141140700e-01, 2.4327217627758024e-02, 1.0776985954942853e-01}));
        
        weight_[35] = 3.3525839026606252e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.0776985954942853e-01, 4.3395146141140700e-01, 4.3395146141140700e-01}));
        
        weight_[36] = 3.3525839026606252e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.0776985954942853e-01, 4.3395146141140700e-01, 2.4327217627758024e-02}));
        
        weight_[37] = 3.3525839026606252e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.0776985954942853e-01, 2.4327217627758024e-02, 4.3395146141140700e-01}));
        
        weight_[38] = 3.3525839026606252e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.4327217627758024e-02, 4.3395146141140700e-01, 1.0776985954942853e-01}));
        
        weight_[39] = 3.3525839026606252e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.4327217627758024e-02, 4.3395146141140700e-01, 4.3395146141140700e-01}));
        
        weight_[40] = 3.3525839026606252e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.4327217627758024e-02, 1.0776985954942853e-01, 4.3395146141140700e-01}));
        
        weight_[41] = 4.3162887555699998e-04;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({-1.3762773181380528e-03, -1.3762773181380528e-03, 2.7655347263680752e-01}));
        
        weight_[42] = 4.3162887555699998e-04;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({-1.3762773181380528e-03, 2.7655347263680752e-01, -1.3762773181380528e-03}));
        
        weight_[43] = 4.3162887555699998e-04;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({-1.3762773181380528e-03, -1.3762773181380528e-03, 7.2619908199946903e-01}));
        
        weight_[44] = 4.3162887555699998e-04;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({-1.3762773181380528e-03, 7.2619908199946903e-01, -1.3762773181380528e-03}));
        
        weight_[45] = 4.3162887555699998e-04;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({-1.3762773181380528e-03, 2.7655347263680752e-01, 7.2619908199946903e-01}));
        
        weight_[46] = 4.3162887555699998e-04;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({-1.3762773181380528e-03, 7.2619908199946903e-01, 2.7655347263680752e-01}));
        
        weight_[47] = 4.3162887555699998e-04;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.7655347263680752e-01, -1.3762773181380528e-03, -1.3762773181380528e-03}));
        
        weight_[48] = 4.3162887555699998e-04;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.7655347263680752e-01, -1.3762773181380528e-03, 7.2619908199946903e-01}));
        
        weight_[49] = 4.3162887555699998e-04;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.7655347263680752e-01, 7.2619908199946903e-01, -1.3762773181380528e-03}));
        
        weight_[50] = 4.3162887555699998e-04;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.2619908199946903e-01, -1.3762773181380528e-03, 2.7655347263680752e-01}));
        
        weight_[51] = 4.3162887555699998e-04;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.2619908199946903e-01, -1.3762773181380528e-03, -1.3762773181380528e-03}));
        
        weight_[52] = 4.3162887555699998e-04;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.2619908199946903e-01, 2.7655347263680752e-01, -1.3762773181380528e-03}));
        
    }
    
//---------------------------------------------------------------------------
    std::string T3_10_126::getName() const
    {
        return name_;
    }
    
    std::size_t T3_10_126::order() const
    {
        return 10;
    }
    
    std::size_t T3_10_126::dimension() const
    {
        return 3;
    }
    
    std::size_t T3_10_126::getNumberOfPoints() const
    {
        return 126;
    }
    
    double T3_10_126::weight(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::getPoint - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReferenceBase&
    T3_10_126::getPoint(std::size_t i) const
    {
        logger.assert(i < getNumberOfPoints(), "%::weight - wrong index!", name_);
        return *gp_[i];
    }
    
    Geometry::ReferenceGeometry*
    T3_10_126::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    T3_10_126::T3_10_126()
            : name_("T3_10_1"), refGeoPtr_(&Geometry::ReferenceTetrahedron::Instance()), gp_(0)
    {
        weight_[0] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.1428571428571508e-02, 7.1428571428571508e-02, 7.8571428571428548e-01}));
        
        weight_[1] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.1428571428571508e-02, 2.1428571428571452e-01, 6.4285714285714302e-01}));
        
        weight_[2] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.1428571428571508e-02, 3.5714285714285698e-01, 5.0000000000000000e-01}));
        
        weight_[3] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.1428571428571508e-02, 5.0000000000000000e-01, 3.5714285714285698e-01}));
        
        weight_[4] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.1428571428571508e-02, 6.4285714285714302e-01, 2.1428571428571452e-01}));
        
        weight_[5] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.1428571428571508e-02, 7.8571428571428548e-01, 7.1428571428571508e-02}));
        
        weight_[6] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.1428571428571452e-01, 7.1428571428571508e-02, 6.4285714285714302e-01}));
        
        weight_[7] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.1428571428571452e-01, 2.1428571428571452e-01, 5.0000000000000000e-01}));
        
        weight_[8] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.1428571428571452e-01, 3.5714285714285698e-01, 3.5714285714285698e-01}));
        
        weight_[9] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.1428571428571452e-01, 5.0000000000000000e-01, 2.1428571428571452e-01}));
        
        weight_[10] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.1428571428571452e-01, 6.4285714285714302e-01, 7.1428571428571508e-02}));
        
        weight_[11] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.5714285714285698e-01, 7.1428571428571508e-02, 5.0000000000000000e-01}));
        
        weight_[12] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.5714285714285698e-01, 2.1428571428571452e-01, 3.5714285714285698e-01}));
        
        weight_[13] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.5714285714285698e-01, 3.5714285714285698e-01, 2.1428571428571452e-01}));
        
        weight_[14] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.5714285714285698e-01, 5.0000000000000000e-01, 7.1428571428571508e-02}));
        
        weight_[15] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({5.0000000000000000e-01, 7.1428571428571508e-02, 3.5714285714285698e-01}));
        
        weight_[16] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({5.0000000000000000e-01, 2.1428571428571452e-01, 2.1428571428571452e-01}));
        
        weight_[17] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({5.0000000000000000e-01, 3.5714285714285698e-01, 7.1428571428571508e-02}));
        
        weight_[18] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({6.4285714285714302e-01, 7.1428571428571508e-02, 2.1428571428571452e-01}));
        
        weight_[19] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({6.4285714285714302e-01, 2.1428571428571452e-01, 7.1428571428571508e-02}));
        
        weight_[20] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.8571428571428548e-01, 7.1428571428571508e-02, 7.1428571428571508e-02}));
        
        weight_[21] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.1428571428571508e-02, 7.1428571428571508e-02, 6.4285714285714302e-01}));
        
        weight_[22] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.1428571428571508e-02, 2.1428571428571452e-01, 5.0000000000000000e-01}));
        
        weight_[23] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.1428571428571508e-02, 3.5714285714285698e-01, 3.5714285714285698e-01}));
        
        weight_[24] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.1428571428571508e-02, 5.0000000000000000e-01, 2.1428571428571452e-01}));
        
        weight_[25] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.1428571428571508e-02, 6.4285714285714302e-01, 7.1428571428571508e-02}));
        
        weight_[26] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.1428571428571452e-01, 7.1428571428571508e-02, 5.0000000000000000e-01}));
        
        weight_[27] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.1428571428571452e-01, 2.1428571428571452e-01, 3.5714285714285698e-01}));
        
        weight_[28] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.1428571428571452e-01, 3.5714285714285698e-01, 2.1428571428571452e-01}));
        
        weight_[29] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.1428571428571452e-01, 5.0000000000000000e-01, 7.1428571428571508e-02}));
        
        weight_[30] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.5714285714285698e-01, 7.1428571428571508e-02, 3.5714285714285698e-01}));
        
        weight_[31] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.5714285714285698e-01, 2.1428571428571452e-01, 2.1428571428571452e-01}));
        
        weight_[32] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.5714285714285698e-01, 3.5714285714285698e-01, 7.1428571428571508e-02}));
        
        weight_[33] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({5.0000000000000000e-01, 7.1428571428571508e-02, 2.1428571428571452e-01}));
        
        weight_[34] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({5.0000000000000000e-01, 2.1428571428571452e-01, 7.1428571428571508e-02}));
        
        weight_[35] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({6.4285714285714302e-01, 7.1428571428571508e-02, 7.1428571428571508e-02}));
        
        weight_[36] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.1428571428571508e-02, 7.1428571428571508e-02, 5.0000000000000000e-01}));
        
        weight_[37] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.1428571428571508e-02, 2.1428571428571452e-01, 3.5714285714285698e-01}));
        
        weight_[38] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.1428571428571508e-02, 3.5714285714285698e-01, 2.1428571428571452e-01}));
        
        weight_[39] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.1428571428571508e-02, 5.0000000000000000e-01, 7.1428571428571508e-02}));
        
        weight_[40] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.1428571428571452e-01, 7.1428571428571508e-02, 3.5714285714285698e-01}));
        
        weight_[41] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.1428571428571452e-01, 2.1428571428571452e-01, 2.1428571428571452e-01}));
        
        weight_[42] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.1428571428571452e-01, 3.5714285714285698e-01, 7.1428571428571508e-02}));
        
        weight_[43] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.5714285714285698e-01, 7.1428571428571508e-02, 2.1428571428571452e-01}));
        
        weight_[44] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.5714285714285698e-01, 2.1428571428571452e-01, 7.1428571428571508e-02}));
        
        weight_[45] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({5.0000000000000000e-01, 7.1428571428571508e-02, 7.1428571428571508e-02}));
        
        weight_[46] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.1428571428571508e-02, 7.1428571428571508e-02, 3.5714285714285698e-01}));
        
        weight_[47] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.1428571428571508e-02, 2.1428571428571452e-01, 2.1428571428571452e-01}));
        
        weight_[48] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.1428571428571508e-02, 3.5714285714285698e-01, 7.1428571428571508e-02}));
        
        weight_[49] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.1428571428571452e-01, 7.1428571428571508e-02, 2.1428571428571452e-01}));
        
        weight_[50] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.1428571428571452e-01, 2.1428571428571452e-01, 7.1428571428571508e-02}));
        
        weight_[51] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.5714285714285698e-01, 7.1428571428571508e-02, 7.1428571428571508e-02}));
        
        weight_[52] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.1428571428571508e-02, 7.1428571428571508e-02, 2.1428571428571452e-01}));
        
        weight_[53] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.1428571428571508e-02, 2.1428571428571452e-01, 7.1428571428571508e-02}));
        
        weight_[54] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.1428571428571452e-01, 7.1428571428571508e-02, 7.1428571428571508e-02}));
        
        weight_[55] = 4.5362824065080999e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.1428571428571508e-02, 7.1428571428571508e-02, 7.1428571428571508e-02}));
        
        weight_[56] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({8.3333333333333481e-02, 8.3333333333333481e-02, 7.5000000000000000e-01}));
        
        weight_[57] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({8.3333333333333481e-02, 2.5000000000000000e-01, 5.8333333333333348e-01}));
        
        weight_[58] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({8.3333333333333481e-02, 4.1666666666666652e-01, 4.1666666666666652e-01}));
        
        weight_[59] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({8.3333333333333481e-02, 5.8333333333333348e-01, 2.5000000000000000e-01}));
        
        weight_[60] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({8.3333333333333481e-02, 7.5000000000000000e-01, 8.3333333333333481e-02}));
        
        weight_[61] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.5000000000000000e-01, 8.3333333333333481e-02, 5.8333333333333348e-01}));
        
        weight_[62] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.5000000000000000e-01, 2.5000000000000000e-01, 4.1666666666666652e-01}));
        
        weight_[63] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.5000000000000000e-01, 4.1666666666666652e-01, 2.5000000000000000e-01}));
        
        weight_[64] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.5000000000000000e-01, 5.8333333333333348e-01, 8.3333333333333481e-02}));
        
        weight_[65] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({4.1666666666666652e-01, 8.3333333333333481e-02, 4.1666666666666652e-01}));
        
        weight_[66] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({4.1666666666666652e-01, 2.5000000000000000e-01, 2.5000000000000000e-01}));
        
        weight_[67] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({4.1666666666666652e-01, 4.1666666666666652e-01, 8.3333333333333481e-02}));
        
        weight_[68] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({5.8333333333333348e-01, 8.3333333333333481e-02, 2.5000000000000000e-01}));
        
        weight_[69] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({5.8333333333333348e-01, 2.5000000000000000e-01, 8.3333333333333481e-02}));
        
        weight_[70] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({7.5000000000000000e-01, 8.3333333333333481e-02, 8.3333333333333481e-02}));
        
        weight_[71] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({8.3333333333333481e-02, 8.3333333333333481e-02, 5.8333333333333348e-01}));
        
        weight_[72] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({8.3333333333333481e-02, 2.5000000000000000e-01, 4.1666666666666652e-01}));
        
        weight_[73] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({8.3333333333333481e-02, 4.1666666666666652e-01, 2.5000000000000000e-01}));
        
        weight_[74] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({8.3333333333333481e-02, 5.8333333333333348e-01, 8.3333333333333481e-02}));
        
        weight_[75] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.5000000000000000e-01, 8.3333333333333481e-02, 4.1666666666666652e-01}));
        
        weight_[76] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.5000000000000000e-01, 2.5000000000000000e-01, 2.5000000000000000e-01}));
        
        weight_[77] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.5000000000000000e-01, 4.1666666666666652e-01, 8.3333333333333481e-02}));
        
        weight_[78] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({4.1666666666666652e-01, 8.3333333333333481e-02, 2.5000000000000000e-01}));
        
        weight_[79] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({4.1666666666666652e-01, 2.5000000000000000e-01, 8.3333333333333481e-02}));
        
        weight_[80] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({5.8333333333333348e-01, 8.3333333333333481e-02, 8.3333333333333481e-02}));
        
        weight_[81] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({8.3333333333333481e-02, 8.3333333333333481e-02, 4.1666666666666652e-01}));
        
        weight_[82] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({8.3333333333333481e-02, 2.5000000000000000e-01, 2.5000000000000000e-01}));
        
        weight_[83] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({8.3333333333333481e-02, 4.1666666666666652e-01, 8.3333333333333481e-02}));
        
        weight_[84] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.5000000000000000e-01, 8.3333333333333481e-02, 2.5000000000000000e-01}));
        
        weight_[85] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.5000000000000000e-01, 2.5000000000000000e-01, 8.3333333333333481e-02}));
        
        weight_[86] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({4.1666666666666652e-01, 8.3333333333333481e-02, 8.3333333333333481e-02}));
        
        weight_[87] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({8.3333333333333481e-02, 8.3333333333333481e-02, 2.5000000000000000e-01}));
        
        weight_[88] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({8.3333333333333481e-02, 2.5000000000000000e-01, 8.3333333333333481e-02}));
        
        weight_[89] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.5000000000000000e-01, 8.3333333333333481e-02, 8.3333333333333481e-02}));
        
        weight_[90] = -1.1652347652347650e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({8.3333333333333481e-02, 8.3333333333333481e-02, 8.3333333333333481e-02}));
        
        weight_[91] = 1.0193728997982475e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({9.9999999999999978e-02, 9.9999999999999978e-02, 6.9999999999999996e-01}));
        
        weight_[92] = 1.0193728997982475e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({9.9999999999999978e-02, 2.9999999999999999e-01, 5.0000000000000000e-01}));
        
        weight_[93] = 1.0193728997982475e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({9.9999999999999978e-02, 5.0000000000000000e-01, 2.9999999999999999e-01}));
        
        weight_[94] = 1.0193728997982475e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({9.9999999999999978e-02, 6.9999999999999996e-01, 9.9999999999999978e-02}));
        
        weight_[95] = 1.0193728997982475e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.9999999999999999e-01, 9.9999999999999978e-02, 5.0000000000000000e-01}));
        
        weight_[96] = 1.0193728997982475e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.9999999999999999e-01, 2.9999999999999999e-01, 2.9999999999999999e-01}));
        
        weight_[97] = 1.0193728997982475e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.9999999999999999e-01, 5.0000000000000000e-01, 9.9999999999999978e-02}));
        
        weight_[98] = 1.0193728997982475e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({5.0000000000000000e-01, 9.9999999999999978e-02, 2.9999999999999999e-01}));
        
        weight_[99] = 1.0193728997982475e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({5.0000000000000000e-01, 2.9999999999999999e-01, 9.9999999999999978e-02}));
        
        weight_[100] = 1.0193728997982475e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({6.9999999999999996e-01, 9.9999999999999978e-02, 9.9999999999999978e-02}));
        
        weight_[101] = 1.0193728997982475e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({9.9999999999999978e-02, 9.9999999999999978e-02, 5.0000000000000000e-01}));
        
        weight_[102] = 1.0193728997982475e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({9.9999999999999978e-02, 2.9999999999999999e-01, 2.9999999999999999e-01}));
        
        weight_[103] = 1.0193728997982475e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({9.9999999999999978e-02, 5.0000000000000000e-01, 9.9999999999999978e-02}));
        
        weight_[104] = 1.0193728997982475e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.9999999999999999e-01, 9.9999999999999978e-02, 2.9999999999999999e-01}));
        
        weight_[105] = 1.0193728997982475e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.9999999999999999e-01, 2.9999999999999999e-01, 9.9999999999999978e-02}));
        
        weight_[106] = 1.0193728997982475e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({5.0000000000000000e-01, 9.9999999999999978e-02, 9.9999999999999978e-02}));
        
        weight_[107] = 1.0193728997982475e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({9.9999999999999978e-02, 9.9999999999999978e-02, 2.9999999999999999e-01}));
        
        weight_[108] = 1.0193728997982475e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({9.9999999999999978e-02, 2.9999999999999999e-01, 9.9999999999999978e-02}));
        
        weight_[109] = 1.0193728997982475e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.9999999999999999e-01, 9.9999999999999978e-02, 9.9999999999999978e-02}));
        
        weight_[110] = 1.0193728997982475e-01;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({9.9999999999999978e-02, 9.9999999999999978e-02, 9.9999999999999978e-02}));
        
        weight_[111] = -3.5025386136497250e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.2500000000000000e-01, 1.2500000000000000e-01, 6.2500000000000000e-01}));
        
        weight_[112] = -3.5025386136497250e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.2500000000000000e-01, 3.7500000000000000e-01, 3.7500000000000000e-01}));
        
        weight_[113] = -3.5025386136497250e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.2500000000000000e-01, 6.2500000000000000e-01, 1.2500000000000000e-01}));
        
        weight_[114] = -3.5025386136497250e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.7500000000000000e-01, 1.2500000000000000e-01, 3.7500000000000000e-01}));
        
        weight_[115] = -3.5025386136497250e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.7500000000000000e-01, 3.7500000000000000e-01, 1.2500000000000000e-01}));
        
        weight_[116] = -3.5025386136497250e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({6.2500000000000000e-01, 1.2500000000000000e-01, 1.2500000000000000e-01}));
        
        weight_[117] = -3.5025386136497250e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.2500000000000000e-01, 1.2500000000000000e-01, 3.7500000000000000e-01}));
        
        weight_[118] = -3.5025386136497250e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.2500000000000000e-01, 3.7500000000000000e-01, 1.2500000000000000e-01}));
        
        weight_[119] = -3.5025386136497250e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({3.7500000000000000e-01, 1.2500000000000000e-01, 1.2500000000000000e-01}));
        
        weight_[120] = -3.5025386136497250e-02;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.2500000000000000e-01, 1.2500000000000000e-01, 1.2500000000000000e-01}));
        
        weight_[121] = 4.0680803571428751e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.6666666666666652e-01, 1.6666666666666652e-01, 5.0000000000000000e-01}));
        
        weight_[122] = 4.0680803571428751e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.6666666666666652e-01, 5.0000000000000000e-01, 1.6666666666666652e-01}));
        
        weight_[123] = 4.0680803571428751e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({5.0000000000000000e-01, 1.6666666666666652e-01, 1.6666666666666652e-01}));
        
        weight_[124] = 4.0680803571428751e-03;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({1.6666666666666652e-01, 1.6666666666666652e-01, 1.6666666666666652e-01}));
        
        weight_[125] = -9.4062316284499996e-05;
        gp_.push_back(Geometry::PointReferenceFactory<3>::instance()->makePoint({2.5000000000000000e-01, 2.5000000000000000e-01, 2.5000000000000000e-01}));
        
    }
//---------------------------------------------------------------------------
}// close namespace IntegrationRules

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
using Geometry::ReferenceTetrahedron;

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
    
    std::size_t Tn3_1_1::nrOfPoints() const
    {
        return 1;
    }
    
    double Tn3_1_1::weight(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReference&
    Tn3_1_1::getPoint(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }
    
    Tn3_1_1::ReferenceGeometryT*
    Tn3_1_1::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    Tn3_1_1::Tn3_1_1()
            : name_("Tn3_1_1"), refGeoPtr_(&ReferenceTetrahedron::Instance()), gp_(0)
    {
        weight_[0] = 1.6666666666666663e-01;
        gp_.push_back(new Geometry::PointReference{2.5000000000000000e-01, 2.5000000000000000e-01, 2.5000000000000000e-01});
        
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
    
    std::size_t Tn3_2_4::nrOfPoints() const
    {
        return 4;
    }
    
    double Tn3_2_4::weight(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReference&
    Tn3_2_4::getPoint(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }
    
    Tn3_2_4::ReferenceGeometryT*
    Tn3_2_4::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    Tn3_2_4::Tn3_2_4()
            : name_("Tn3_2_1"), refGeoPtr_(&ReferenceTetrahedron::Instance()), gp_(0)
    {
        weight_[0] = 4.1666666666666623e-02;
        gp_.push_back(new Geometry::PointReference{1.3819660112501048e-01, 1.3819660112501048e-01, 1.3819660112501048e-01});
        
        weight_[1] = 4.1666666666666623e-02;
        gp_.push_back(new Geometry::PointReference{5.8541019662496852e-01, 1.3819660112501048e-01, 1.3819660112501048e-01});
        
        weight_[2] = 4.1666666666666623e-02;
        gp_.push_back(new Geometry::PointReference{1.3819660112501048e-01, 5.8541019662496852e-01, 1.3819660112501048e-01});
        
        weight_[3] = 4.1666666666666623e-02;
        gp_.push_back(new Geometry::PointReference{1.3819660112501048e-01, 1.3819660112501048e-01, 5.8541019662496852e-01});
        
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
    
    std::size_t Tn3_3_5::nrOfPoints() const
    {
        return 5;
    }
    
    double Tn3_3_5::weight(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReference&
    Tn3_3_5::getPoint(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }
    
    Tn3_3_5::ReferenceGeometryT*
    Tn3_3_5::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    Tn3_3_5::Tn3_3_5()
            : name_("Tn3_3_1"), refGeoPtr_(&ReferenceTetrahedron::Instance()), gp_(0)
    {
        weight_[0] = -1.3333333333333339e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[0])[0] = 2.5000000000000000e-01;
        (*gp_[0])[1] = 2.5000000000000000e-01;
        (*gp_[0])[2] = 2.5000000000000000e-01;
        
        weight_[1] = 7.4999999999999997e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[1])[0] = 1.6666666666666652e-01;
        (*gp_[1])[1] = 1.6666666666666652e-01;
        (*gp_[1])[2] = 1.6666666666666652e-01;
        
        weight_[2] = 7.4999999999999997e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[2])[0] = 1.6666666666666652e-01;
        (*gp_[2])[1] = 1.6666666666666652e-01;
        (*gp_[2])[2] = 5.0000000000000000e-01;
        
        weight_[3] = 7.4999999999999997e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[3])[0] = 1.6666666666666652e-01;
        (*gp_[3])[1] = 5.0000000000000000e-01;
        (*gp_[3])[2] = 1.6666666666666652e-01;
        
        weight_[4] = 7.4999999999999997e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[4])[0] = 5.0000000000000000e-01;
        (*gp_[4])[1] = 1.6666666666666652e-01;
        (*gp_[4])[2] = 1.6666666666666652e-01;
        
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
    
    std::size_t Tn3_4_11::nrOfPoints() const
    {
        return 11;
    }
    
    double Tn3_4_11::weight(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReference&
    Tn3_4_11::getPoint(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }
    
    Tn3_4_11::ReferenceGeometryT*
    Tn3_4_11::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    Tn3_4_11::Tn3_4_11()
            : name_("Tn3_4_1"), refGeoPtr_(&ReferenceTetrahedron::Instance()), gp_(0)
    {
        weight_[0] = -1.3155555555555500e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[0])[0] = 2.5000000000000000e-01;
        (*gp_[0])[1] = 2.5000000000000000e-01;
        (*gp_[0])[2] = 2.5000000000000000e-01;
        
        weight_[1] = 7.6222222222222498e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[1])[0] = 7.1428571428571508e-02;
        (*gp_[1])[1] = 7.1428571428571508e-02;
        (*gp_[1])[2] = 7.1428571428571508e-02;
        
        weight_[2] = 7.6222222222222498e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[2])[0] = 7.1428571428571508e-02;
        (*gp_[2])[1] = 7.1428571428571508e-02;
        (*gp_[2])[2] = 7.8571428571428548e-01;
        
        weight_[3] = 7.6222222222222498e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[3])[0] = 7.1428571428571508e-02;
        (*gp_[3])[1] = 7.8571428571428548e-01;
        (*gp_[3])[2] = 7.1428571428571508e-02;
        
        weight_[4] = 7.6222222222222498e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[4])[0] = 7.8571428571428548e-01;
        (*gp_[4])[1] = 7.1428571428571508e-02;
        (*gp_[4])[2] = 7.1428571428571508e-02;
        
        weight_[5] = 2.4888888888888874e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[5])[0] = 3.9940357616679900e-01;
        (*gp_[5])[1] = 3.9940357616679900e-01;
        (*gp_[5])[2] = 1.0059642383320100e-01;
        
        weight_[6] = 2.4888888888888874e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[6])[0] = 3.9940357616679900e-01;
        (*gp_[6])[1] = 1.0059642383320100e-01;
        (*gp_[6])[2] = 3.9940357616679900e-01;
        
        weight_[7] = 2.4888888888888874e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[7])[0] = 1.0059642383320100e-01;
        (*gp_[7])[1] = 3.9940357616679900e-01;
        (*gp_[7])[2] = 3.9940357616679900e-01;
        
        weight_[8] = 2.4888888888888874e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[8])[0] = 3.9940357616679900e-01;
        (*gp_[8])[1] = 1.0059642383320100e-01;
        (*gp_[8])[2] = 1.0059642383320100e-01;
        
        weight_[9] = 2.4888888888888874e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[9])[0] = 1.0059642383320100e-01;
        (*gp_[9])[1] = 3.9940357616679900e-01;
        (*gp_[9])[2] = 1.0059642383320100e-01;
        
        weight_[10] = 2.4888888888888874e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[10])[0] = 1.0059642383320100e-01;
        (*gp_[10])[1] = 1.0059642383320100e-01;
        (*gp_[10])[2] = 3.9940357616679900e-01;
        
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
    
    std::size_t T3_5_14::nrOfPoints() const
    {
        return 14;
    }
    
    double T3_5_14::weight(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReference&
    T3_5_14::getPoint(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }
    
    T3_5_14::ReferenceGeometryT*
    T3_5_14::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    T3_5_14::T3_5_14()
            : name_("T3_5_1"), refGeoPtr_(&ReferenceTetrahedron::Instance()), gp_(0)
    {
        weight_[0] = 1.2248840519393626e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[0])[0] = 9.2735250310891026e-02;
        (*gp_[0])[1] = 9.2735250310891026e-02;
        (*gp_[0])[2] = 9.2735250310891026e-02;
        
        weight_[1] = 1.2248840519393626e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[1])[0] = 7.2179424906732648e-01;
        (*gp_[1])[1] = 9.2735250310891026e-02;
        (*gp_[1])[2] = 9.2735250310891026e-02;
        
        weight_[2] = 1.2248840519393626e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[2])[0] = 9.2735250310891026e-02;
        (*gp_[2])[1] = 7.2179424906732648e-01;
        (*gp_[2])[2] = 9.2735250310891026e-02;
        
        weight_[3] = 1.2248840519393626e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[3])[0] = 9.2735250310891026e-02;
        (*gp_[3])[1] = 9.2735250310891026e-02;
        (*gp_[3])[2] = 7.2179424906732648e-01;
        
        weight_[4] = 1.8781320953002625e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[4])[0] = 3.1088591926330050e-01;
        (*gp_[4])[1] = 3.1088591926330050e-01;
        (*gp_[4])[2] = 3.1088591926330050e-01;
        
        weight_[5] = 1.8781320953002625e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[5])[0] = 6.7342242210097991e-02;
        (*gp_[5])[1] = 3.1088591926330050e-01;
        (*gp_[5])[2] = 3.1088591926330050e-01;
        
        weight_[6] = 1.8781320953002625e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[6])[0] = 3.1088591926330050e-01;
        (*gp_[6])[1] = 6.7342242210097991e-02;
        (*gp_[6])[2] = 3.1088591926330050e-01;
        
        weight_[7] = 1.8781320953002625e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[7])[0] = 3.1088591926330050e-01;
        (*gp_[7])[1] = 3.1088591926330050e-01;
        (*gp_[7])[2] = 6.7342242210097991e-02;
        
        weight_[8] = 7.0910034628468748e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[8])[0] = 4.5449629587435048e-01;
        (*gp_[8])[1] = 4.5449629587435048e-01;
        (*gp_[8])[2] = 4.5503704125649524e-02;
        
        weight_[9] = 7.0910034628468748e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[9])[0] = 4.5449629587435048e-01;
        (*gp_[9])[1] = 4.5503704125649524e-02;
        (*gp_[9])[2] = 4.5449629587435048e-01;
        
        weight_[10] = 7.0910034628468748e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[10])[0] = 4.5503704125649524e-02;
        (*gp_[10])[1] = 4.5449629587435048e-01;
        (*gp_[10])[2] = 4.5449629587435048e-01;
        
        weight_[11] = 7.0910034628468748e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[11])[0] = 4.5449629587435048e-01;
        (*gp_[11])[1] = 4.5503704125649524e-02;
        (*gp_[11])[2] = 4.5503704125649524e-02;
        
        weight_[12] = 7.0910034628468748e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[12])[0] = 4.5503704125649524e-02;
        (*gp_[12])[1] = 4.5449629587435048e-01;
        (*gp_[12])[2] = 4.5503704125649524e-02;
        
        weight_[13] = 7.0910034628468748e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[13])[0] = 4.5503704125649524e-02;
        (*gp_[13])[1] = 4.5503704125649524e-02;
        (*gp_[13])[2] = 4.5449629587435048e-01;
        
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
    
    std::size_t T3_6_24::nrOfPoints() const
    {
        return 24;
    }
    
    double T3_6_24::weight(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReference&
    T3_6_24::getPoint(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }
    
    T3_6_24::ReferenceGeometryT*
    T3_6_24::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    T3_6_24::T3_6_24()
            : name_("T3_6_1"), refGeoPtr_(&ReferenceTetrahedron::Instance()), gp_(0)
    {
        weight_[0] = 6.6537917096946252e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[0])[0] = 2.1460287125915201e-01;
        (*gp_[0])[1] = 2.1460287125915201e-01;
        (*gp_[0])[2] = 2.1460287125915201e-01;
        
        weight_[1] = 6.6537917096946252e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[1])[0] = 3.5619138622254398e-01;
        (*gp_[1])[1] = 2.1460287125915201e-01;
        (*gp_[1])[2] = 2.1460287125915201e-01;
        
        weight_[2] = 6.6537917096946252e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[2])[0] = 2.1460287125915201e-01;
        (*gp_[2])[1] = 3.5619138622254398e-01;
        (*gp_[2])[2] = 2.1460287125915201e-01;
        
        weight_[3] = 6.6537917096946252e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[3])[0] = 2.1460287125915201e-01;
        (*gp_[3])[1] = 2.1460287125915201e-01;
        (*gp_[3])[2] = 3.5619138622254398e-01;
        
        weight_[4] = 1.6795351758867501e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[4])[0] = 4.0673958534611476e-02;
        (*gp_[4])[1] = 4.0673958534611476e-02;
        (*gp_[4])[2] = 4.0673958534611476e-02;
        
        weight_[5] = 1.6795351758867501e-03;
        (*gp_[5])[0] = 8.7797812439616596e-01;
        (*gp_[5])[1] = 4.0673958534611476e-02;
        (*gp_[5])[2] = 4.0673958534611476e-02;
        
        weight_[6] = 1.6795351758867501e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[6])[0] = 4.0673958534611476e-02;
        (*gp_[6])[1] = 8.7797812439616596e-01;
        (*gp_[6])[2] = 4.0673958534611476e-02;
        
        weight_[7] = 1.6795351758867501e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[7])[0] = 4.0673958534611476e-02;
        (*gp_[7])[1] = 4.0673958534611476e-02;
        (*gp_[7])[2] = 8.7797812439616596e-01;
        
        weight_[8] = 9.2261969239424996e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[8])[0] = 3.2233789014227554e-01;
        (*gp_[8])[1] = 3.2233789014227554e-01;
        (*gp_[8])[2] = 3.2233789014227554e-01;
        
        weight_[9] = 9.2261969239424996e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[9])[0] = 3.2986329573173490e-02;
        (*gp_[9])[1] = 3.2233789014227554e-01;
        (*gp_[9])[2] = 3.2233789014227554e-01;
        
        weight_[10] = 9.2261969239424996e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[10])[0] = 3.2233789014227554e-01;
        (*gp_[10])[1] = 3.2986329573173490e-02;
        (*gp_[10])[2] = 3.2233789014227554e-01;
        
        weight_[11] = 9.2261969239424996e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[11])[0] = 3.2233789014227554e-01;
        (*gp_[11])[1] = 3.2233789014227554e-01;
        (*gp_[11])[2] = 3.2986329573173490e-02;
        
        weight_[12] = 8.0357142857142502e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[12])[0] = 6.3661001875017498e-02;
        (*gp_[12])[1] = 6.3661001875017498e-02;
        (*gp_[12])[2] = 2.6967233145831604e-01;
        
        weight_[13] = 8.0357142857142502e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[13])[0] = 6.3661001875017498e-02;
        (*gp_[13])[1] = 2.6967233145831604e-01;
        (*gp_[13])[2] = 6.3661001875017498e-02;
        
        weight_[14] = 8.0357142857142502e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[14])[0] = 6.3661001875017498e-02;
        (*gp_[14])[1] = 6.3661001875017498e-02;
        (*gp_[14])[2] = 6.0300566479164897e-01;
        
        weight_[15] = 8.0357142857142502e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[15])[0] = 6.3661001875017498e-02;
        (*gp_[15])[1] = 6.0300566479164897e-01;
        (*gp_[15])[2] = 6.3661001875017498e-02;
        
        weight_[16] = 8.0357142857142502e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[16])[0] = 6.3661001875017498e-02;
        (*gp_[16])[1] = 2.6967233145831604e-01;
        (*gp_[16])[2] = 6.0300566479164897e-01;
        
        weight_[17] = 8.0357142857142502e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[17])[0] = 6.3661001875017498e-02;
        (*gp_[17])[1] = 6.0300566479164897e-01;
        (*gp_[17])[2] = 2.6967233145831604e-01;
        
        weight_[18] = 8.0357142857142502e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[18])[0] = 2.6967233145831604e-01;
        (*gp_[18])[1] = 6.3661001875017498e-02;
        (*gp_[18])[2] = 6.3661001875017498e-02;
        
        weight_[19] = 8.0357142857142502e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[19])[0] = 2.6967233145831604e-01;
        (*gp_[19])[1] = 6.3661001875017498e-02;
        (*gp_[19])[2] = 6.0300566479164897e-01;
        
        weight_[20] = 8.0357142857142502e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[20])[0] = 2.6967233145831604e-01;
        (*gp_[20])[1] = 6.0300566479164897e-01;
        (*gp_[20])[2] = 6.3661001875017498e-02;
        
        weight_[21] = 8.0357142857142502e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[21])[0] = 6.0300566479164897e-01;
        (*gp_[21])[1] = 6.3661001875017498e-02;
        (*gp_[21])[2] = 2.6967233145831604e-01;
        
        weight_[22] = 8.0357142857142502e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[22])[0] = 6.0300566479164897e-01;
        (*gp_[22])[1] = 6.3661001875017498e-02;
        (*gp_[22])[2] = 6.3661001875017498e-02;
        
        weight_[23] = 8.0357142857142502e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[23])[0] = 6.0300566479164897e-01;
        (*gp_[23])[1] = 2.6967233145831604e-01;
        (*gp_[23])[2] = 6.3661001875017498e-02;
        
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
    
    std::size_t T3_7_31::nrOfPoints() const
    {
        return 31;
    }
    
    double T3_7_31::weight(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReference&
    T3_7_31::getPoint(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }
    
    T3_7_31::ReferenceGeometryT*
    T3_7_31::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    T3_7_31::T3_7_31()
            : name_("T3_7_1"), refGeoPtr_(&ReferenceTetrahedron::Instance()), gp_(0)
    {
        weight_[0] = 9.7001763668425002e-04;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[0])[0] = 5.0000000000000000e-01;
        (*gp_[0])[1] = 5.0000000000000000e-01;
        (*gp_[0])[2] = 0.0000000000000000e+00;
        
        weight_[1] = 9.7001763668425002e-04;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[1])[0] = 5.0000000000000000e-01;
        (*gp_[1])[1] = 0.0000000000000000e+00;
        (*gp_[1])[2] = 5.0000000000000000e-01;
        
        weight_[2] = 9.7001763668425002e-04;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[2])[0] = 0.0000000000000000e+00;
        (*gp_[2])[1] = 5.0000000000000000e-01;
        (*gp_[2])[2] = 5.0000000000000000e-01;
        
        weight_[3] = 9.7001763668425002e-04;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[3])[0] = 0.0000000000000000e+00;
        (*gp_[3])[1] = 0.0000000000000000e+00;
        (*gp_[3])[2] = 5.0000000000000000e-01;
        
        weight_[4] = 9.7001763668425002e-04;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[4])[0] = 0.0000000000000000e+00;
        (*gp_[4])[1] = 5.0000000000000000e-01;
        (*gp_[4])[2] = 0.0000000000000000e+00;
        
        weight_[5] = 9.7001763668425002e-04;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[5])[0] = 5.0000000000000000e-01;
        (*gp_[5])[1] = 0.0000000000000000e+00;
        (*gp_[5])[2] = 0.0000000000000000e+00;
        
        weight_[6] = 1.8264223466108873e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[6])[0] = 2.5000000000000000e-01;
        (*gp_[6])[1] = 2.5000000000000000e-01;
        (*gp_[6])[2] = 2.5000000000000000e-01;
        
        weight_[7] = 1.0599941524413625e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[7])[0] = 7.8213192330317982e-02;
        (*gp_[7])[1] = 7.8213192330317982e-02;
        (*gp_[7])[2] = 7.8213192330317982e-02;
        
        weight_[8] = 1.0599941524413625e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[8])[0] = 7.8213192330317982e-02;
        (*gp_[8])[1] = 7.8213192330317982e-02;
        (*gp_[8])[2] = 7.6536042300904605e-01;
        
        weight_[9] = 1.0599941524413625e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[9])[0] = 7.8213192330317982e-02;
        (*gp_[9])[1] = 7.6536042300904605e-01;
        (*gp_[9])[2] = 7.8213192330317982e-02;
        
        weight_[10] = 1.0599941524413625e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[10])[0] = 7.6536042300904605e-01;
        (*gp_[10])[1] = 7.8213192330317982e-02;
        (*gp_[10])[2] = 7.8213192330317982e-02;
        
        weight_[11] = -6.2517740114331879e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[11])[0] = 1.2184321666390502e-01;
        (*gp_[11])[1] = 1.2184321666390502e-01;
        (*gp_[11])[2] = 1.2184321666390502e-01;
        
        weight_[12] = -6.2517740114331879e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[12])[0] = 1.2184321666390502e-01;
        (*gp_[12])[1] = 1.2184321666390502e-01;
        (*gp_[12])[2] = 6.3447035000828444e-01;
        
        weight_[13] = -6.2517740114331879e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[13])[0] = 1.2184321666390502e-01;
        (*gp_[13])[1] = 6.3447035000828444e-01;
        (*gp_[13])[2] = 1.2184321666390502e-01;
        
        weight_[14] = -6.2517740114331879e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[14])[0] = 6.3447035000828444e-01;
        (*gp_[14])[1] = 1.2184321666390502e-01;
        (*gp_[14])[2] = 1.2184321666390502e-01;
        
        weight_[15] = 4.8914252630735001e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[15])[0] = 3.3253916444642051e-01;
        (*gp_[15])[1] = 3.3253916444642051e-01;
        (*gp_[15])[2] = 3.3253916444642051e-01;
        
        weight_[16] = 4.8914252630735001e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[16])[0] = 3.3253916444642051e-01;
        (*gp_[16])[1] = 3.3253916444642051e-01;
        (*gp_[16])[2] = 2.3825066607380263e-03;
        
        weight_[17] = 4.8914252630735001e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[17])[0] = 3.3253916444642051e-01;
        (*gp_[17])[1] = 2.3825066607380263e-03;
        (*gp_[17])[2] = 3.3253916444642051e-01;
        
        weight_[18] = 4.8914252630735001e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[18])[0] = 2.3825066607380263e-03;
        (*gp_[18])[1] = 3.3253916444642051e-01;
        (*gp_[18])[2] = 3.3253916444642051e-01;
        
        weight_[19] = 2.7557319223985875e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[19])[0] = 9.9999999999999978e-02;
        (*gp_[19])[1] = 9.9999999999999978e-02;
        (*gp_[19])[2] = 2.0000000000000001e-01;
        
        weight_[20] = 2.7557319223985875e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[20])[0] = 9.9999999999999978e-02;
        (*gp_[20])[1] = 2.0000000000000001e-01;
        (*gp_[20])[2] = 9.9999999999999978e-02;
        
        weight_[21] = 2.7557319223985875e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[21])[0] = 9.9999999999999978e-02;
        (*gp_[21])[1] = 9.9999999999999978e-02;
        (*gp_[21])[2] = 5.9999999999999998e-01;
        
        weight_[22] = 2.7557319223985875e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[22])[0] = 9.9999999999999978e-02;
        (*gp_[22])[1] = 5.9999999999999998e-01;
        (*gp_[22])[2] = 9.9999999999999978e-02;
        
        weight_[23] = 2.7557319223985875e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[23])[0] = 9.9999999999999978e-02;
        (*gp_[23])[1] = 2.0000000000000001e-01;
        (*gp_[23])[2] = 5.9999999999999998e-01;
        
        weight_[24] = 2.7557319223985875e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[24])[0] = 9.9999999999999978e-02;
        (*gp_[24])[1] = 5.9999999999999998e-01;
        (*gp_[24])[2] = 2.0000000000000001e-01;
        
        weight_[25] = 2.7557319223985875e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[25])[0] = 2.0000000000000001e-01;
        (*gp_[25])[1] = 9.9999999999999978e-02;
        (*gp_[25])[2] = 9.9999999999999978e-02;
        
        weight_[26] = 2.7557319223985875e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[26])[0] = 2.0000000000000001e-01;
        (*gp_[26])[1] = 9.9999999999999978e-02;
        (*gp_[26])[2] = 5.9999999999999998e-01;
        
        weight_[27] = 2.7557319223985875e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[27])[0] = 2.0000000000000001e-01;
        (*gp_[27])[1] = 5.9999999999999998e-01;
        (*gp_[27])[2] = 9.9999999999999978e-02;
        
        weight_[28] = 2.7557319223985875e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[28])[0] = 5.9999999999999998e-01;
        (*gp_[28])[1] = 9.9999999999999978e-02;
        (*gp_[28])[2] = 2.0000000000000001e-01;
        
        weight_[29] = 2.7557319223985875e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[29])[0] = 5.9999999999999998e-01;
        (*gp_[29])[1] = 9.9999999999999978e-02;
        (*gp_[29])[2] = 9.9999999999999978e-02;
        
        weight_[30] = 2.7557319223985875e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[30])[0] = 5.9999999999999998e-01;
        (*gp_[30])[1] = 2.0000000000000001e-01;
        (*gp_[30])[2] = 9.9999999999999978e-02;
        
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
    
    std::size_t T3_8_43::nrOfPoints() const
    {
        return 43;
    }
    
    double T3_8_43::weight(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReference&
    T3_8_43::getPoint(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }
    
    T3_8_43::ReferenceGeometryT*
    T3_8_43::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    T3_8_43::T3_8_43()
            : name_("T3_8_1"), refGeoPtr_(&ReferenceTetrahedron::Instance()), gp_(0)
    {
        weight_[0] = -2.0500188658639874e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[0])[0] = 2.5000000000000000e-01;
        (*gp_[0])[1] = 2.5000000000000000e-01;
        (*gp_[0])[2] = 2.5000000000000000e-01;
        
        weight_[1] = 1.4250305822866875e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[1])[0] = 2.0682993161067298e-01;
        (*gp_[1])[1] = 2.0682993161067298e-01;
        (*gp_[1])[2] = 2.0682993161067298e-01;
        
        weight_[2] = 1.4250305822866875e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[2])[0] = 2.0682993161067298e-01;
        (*gp_[2])[1] = 2.0682993161067298e-01;
        (*gp_[2])[2] = 3.7951020516798051e-01;
        
        weight_[3] = 1.4250305822866875e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[3])[0] = 2.0682993161067298e-01;
        (*gp_[3])[1] = 3.7951020516798051e-01;
        (*gp_[3])[2] = 2.0682993161067298e-01;
        
        weight_[4] = 1.4250305822866875e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[4])[0] = 3.7951020516798051e-01;
        (*gp_[4])[1] = 2.0682993161067298e-01;
        (*gp_[4])[2] = 2.0682993161067298e-01;
        
        weight_[5] = 1.9670333131338751e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[5])[0] = 8.2103588310546483e-02;
        (*gp_[5])[1] = 8.2103588310546483e-02;
        (*gp_[5])[2] = 8.2103588310546483e-02;
        
        weight_[6] = 1.9670333131338751e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[6])[0] = 8.2103588310546483e-02;
        (*gp_[6])[1] = 8.2103588310546483e-02;
        (*gp_[6])[2] = 7.5368923506835994e-01;
        
        weight_[7] = 1.9670333131338751e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[7])[0] = 8.2103588310546483e-02;
        (*gp_[7])[1] = 7.5368923506835994e-01;
        (*gp_[7])[2] = 8.2103588310546483e-02;
        
        weight_[8] = 1.9670333131338751e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[8])[0] = 7.5368923506835994e-01;
        (*gp_[8])[1] = 8.2103588310546483e-02;
        (*gp_[8])[2] = 8.2103588310546483e-02;
        
        weight_[9] = 1.6983410909287501e-04;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[9])[0] = 5.7819505051979747e-03;
        (*gp_[9])[1] = 5.7819505051979747e-03;
        (*gp_[9])[2] = 5.7819505051979747e-03;
        
        weight_[10] = 1.6983410909287501e-04;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[10])[0] = 5.7819505051979747e-03;
        (*gp_[10])[1] = 5.7819505051979747e-03;
        (*gp_[10])[2] = 9.8265414848440602e-01;
        
        weight_[11] = 1.6983410909287501e-04;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[11])[0] = 5.7819505051979747e-03;
        (*gp_[11])[1] = 9.8265414848440602e-01;
        (*gp_[11])[2] = 5.7819505051979747e-03;
        
        weight_[12] = 1.6983410909287501e-04;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[12])[0] = 9.8265414848440602e-01;
        (*gp_[12])[1] = 5.7819505051979747e-03;
        (*gp_[12])[2] = 5.7819505051979747e-03;
        
        weight_[13] = 4.5796838244672499e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[13])[0] = 5.0532740018894018e-02;
        (*gp_[13])[1] = 5.0532740018894018e-02;
        (*gp_[13])[2] = 4.4946725998110598e-01;
        
        weight_[14] = 4.5796838244672499e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[14])[0] = 5.0532740018894018e-02;
        (*gp_[14])[1] = 4.4946725998110598e-01;
        (*gp_[14])[2] = 5.0532740018894018e-02;
        
        weight_[15] = 4.5796838244672499e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[15])[0] = 4.4946725998110598e-01;
        (*gp_[15])[1] = 5.0532740018894018e-02;
        (*gp_[15])[2] = 5.0532740018894018e-02;
        
        weight_[16] = 4.5796838244672499e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[16])[0] = 5.0532740018894018e-02;
        (*gp_[16])[1] = 4.4946725998110598e-01;
        (*gp_[16])[2] = 4.4946725998110598e-01;
        
        weight_[17] = 4.5796838244672499e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[17])[0] = 4.4946725998110598e-01;
        (*gp_[17])[1] = 5.0532740018894018e-02;
        (*gp_[17])[2] = 4.4946725998110598e-01;
        
        weight_[18] = 4.5796838244672499e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[18])[0] = 4.4946725998110598e-01;
        (*gp_[18])[1] = 4.4946725998110598e-01;
        (*gp_[18])[2] = 5.0532740018894018e-02;
        
        weight_[19] = 5.7044858086818754e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[19])[0] = 2.2906653611681099e-01;
        (*gp_[19])[1] = 2.2906653611681099e-01;
        (*gp_[19])[2] = 3.5639582788534019e-02;
        
        weight_[20] = 5.7044858086818754e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[20])[0] = 2.2906653611681099e-01;
        (*gp_[20])[1] = 3.5639582788534019e-02;
        (*gp_[20])[2] = 2.2906653611681099e-01;
        
        weight_[21] = 5.7044858086818754e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[21])[0] = 2.2906653611681099e-01;
        (*gp_[21])[1] = 2.2906653611681099e-01;
        (*gp_[21])[2] = 5.0622734497784350e-01;
        
        weight_[22] = 5.7044858086818754e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[22])[0] = 2.2906653611681099e-01;
        (*gp_[22])[1] = 5.0622734497784350e-01;
        (*gp_[22])[2] = 2.2906653611681099e-01;
        
        weight_[23] = 5.7044858086818754e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[23])[0] = 2.2906653611681099e-01;
        (*gp_[23])[1] = 3.5639582788534019e-02;
        (*gp_[23])[2] = 5.0622734497784350e-01;
        
        weight_[24] = 5.7044858086818754e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[24])[0] = 2.2906653611681099e-01;
        (*gp_[24])[1] = 5.0622734497784350e-01;
        (*gp_[24])[2] = 3.5639582788534019e-02;
        
        weight_[25] = 5.7044858086818754e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[25])[0] = 3.5639582788534019e-02;
        (*gp_[25])[1] = 2.2906653611681099e-01;
        (*gp_[25])[2] = 2.2906653611681099e-01;
        
        weight_[26] = 5.7044858086818754e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[26])[0] = 3.5639582788534019e-02;
        (*gp_[26])[1] = 2.2906653611681099e-01;
        (*gp_[26])[2] = 5.0622734497784350e-01;
        
        weight_[27] = 5.7044858086818754e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[27])[0] = 3.5639582788534019e-02;
        (*gp_[27])[1] = 5.0622734497784350e-01;
        (*gp_[27])[2] = 2.2906653611681099e-01;
        
        weight_[28] = 5.7044858086818754e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[28])[0] = 5.0622734497784350e-01;
        (*gp_[28])[1] = 2.2906653611681099e-01;
        (*gp_[28])[2] = 3.5639582788534019e-02;
        
        weight_[29] = 5.7044858086818754e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[29])[0] = 5.0622734497784350e-01;
        (*gp_[29])[1] = 2.2906653611681099e-01;
        (*gp_[29])[2] = 2.2906653611681099e-01;
        
        weight_[30] = 5.7044858086818754e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[30])[0] = 5.0622734497784350e-01;
        (*gp_[30])[1] = 3.5639582788534019e-02;
        (*gp_[30])[2] = 2.2906653611681099e-01;
        
        weight_[31] = 2.1405191411621250e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[31])[0] = 3.6607749553197511e-02;
        (*gp_[31])[1] = 3.6607749553197511e-02;
        (*gp_[31])[2] = 1.9048604193463348e-01;
        
        weight_[32] = 2.1405191411621250e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[32])[0] = 3.6607749553197511e-02;
        (*gp_[32])[1] = 1.9048604193463348e-01;
        (*gp_[32])[2] = 3.6607749553197511e-02;
        
        weight_[33] = 2.1405191411621250e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[33])[0] = 3.6607749553197511e-02;
        (*gp_[33])[1] = 3.6607749553197511e-02;
        (*gp_[33])[2] = 7.3629845895897150e-01;
        
        weight_[34] = 2.1405191411621250e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[34])[0] = 3.6607749553197511e-02;
        (*gp_[34])[1] = 7.3629845895897150e-01;
        (*gp_[34])[2] = 3.6607749553197511e-02;
        
        weight_[35] = 2.1405191411621250e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[35])[0] = 3.6607749553197511e-02;
        (*gp_[35])[1] = 1.9048604193463348e-01;
        (*gp_[35])[2] = 7.3629845895897150e-01;
        
        weight_[36] = 2.1405191411621250e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[36])[0] = 3.6607749553197511e-02;
        (*gp_[36])[1] = 7.3629845895897150e-01;
        (*gp_[36])[2] = 1.9048604193463348e-01;
        
        weight_[37] = 2.1405191411621250e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[37])[0] = 1.9048604193463348e-01;
        (*gp_[37])[1] = 3.6607749553197511e-02;
        (*gp_[37])[2] = 3.6607749553197511e-02;
        
        weight_[38] = 2.1405191411621250e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[38])[0] = 1.9048604193463348e-01;
        (*gp_[38])[1] = 3.6607749553197511e-02;
        (*gp_[38])[2] = 7.3629845895897150e-01;
        
        weight_[39] = 2.1405191411621250e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[39])[0] = 1.9048604193463348e-01;
        (*gp_[39])[1] = 7.3629845895897150e-01;
        (*gp_[39])[2] = 3.6607749553197511e-02;
        
        weight_[40] = 2.1405191411621250e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[40])[0] = 7.3629845895897150e-01;
        (*gp_[40])[1] = 3.6607749553197511e-02;
        (*gp_[40])[2] = 1.9048604193463348e-01;
        
        weight_[41] = 2.1405191411621250e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[41])[0] = 7.3629845895897150e-01;
        (*gp_[41])[1] = 3.6607749553197511e-02;
        (*gp_[41])[2] = 3.6607749553197511e-02;
        
        weight_[42] = 2.1405191411621250e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[42])[0] = 7.3629845895897150e-01;
        (*gp_[42])[1] = 1.9048604193463348e-01;
        (*gp_[42])[2] = 3.6607749553197511e-02;
        
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
    
    std::size_t T3_9_53::nrOfPoints() const
    {
        return 53;
    }
    
    double T3_9_53::weight(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReference&
    T3_9_53::getPoint(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::getPoint - wrong index!", name_);
        return *gp_[i];
    }
    
    T3_9_53::ReferenceGeometryT*
    T3_9_53::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    T3_9_53::T3_9_53()
            : name_("T3_9_1"), refGeoPtr_(&ReferenceTetrahedron::Instance()), gp_(0)
    {
        weight_[0] = -1.3779903832610862e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[0])[0] = 2.5000000000000000e-01;
        (*gp_[0])[1] = 2.5000000000000000e-01;
        (*gp_[0])[2] = 2.5000000000000000e-01;
        
        weight_[1] = 1.8653365690852501e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[1])[0] = 4.8351038549736991e-02;
        (*gp_[1])[1] = 4.8351038549736991e-02;
        (*gp_[1])[2] = 4.8351038549736991e-02;
        
        weight_[2] = 1.8653365690852501e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[2])[0] = 4.8351038549736991e-02;
        (*gp_[2])[1] = 4.8351038549736991e-02;
        (*gp_[2])[2] = 8.5494688435079003e-01;
        
        weight_[3] = 1.8653365690852501e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[3])[0] = 4.8351038549736991e-02;
        (*gp_[3])[1] = 8.5494688435079003e-01;
        (*gp_[3])[2] = 4.8351038549736991e-02;
        
        weight_[4] = 1.8653365690852501e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[4])[0] = 8.5494688435079003e-01;
        (*gp_[4])[1] = 4.8351038549736991e-02;
        (*gp_[4])[2] = 4.8351038549736991e-02;
        
        weight_[5] = 4.3094239694933751e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[5])[0] = 3.2457928011788251e-01;
        (*gp_[5])[1] = 3.2457928011788251e-01;
        (*gp_[5])[2] = 3.2457928011788251e-01;
        
        weight_[6] = 4.3094239694933751e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[6])[0] = 3.2457928011788251e-01;
        (*gp_[6])[1] = 3.2457928011788251e-01;
        (*gp_[6])[2] = 2.6262159646352978e-02;
        
        weight_[7] = 4.3094239694933751e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[7])[0] = 3.2457928011788251e-01;
        (*gp_[7])[1] = 2.6262159646352978e-02;
        (*gp_[7])[2] = 3.2457928011788251e-01;
        
        weight_[8] = 4.3094239694933751e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[8])[0] = 2.6262159646352978e-02;
        (*gp_[8])[1] = 3.2457928011788251e-01;
        (*gp_[8])[2] = 3.2457928011788251e-01;
        
        weight_[9] = -9.0184766481201495e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[9])[0] = 1.1461654022399498e-01;
        (*gp_[9])[1] = 1.1461654022399498e-01;
        (*gp_[9])[2] = 1.1461654022399498e-01;
        
        weight_[10] = -9.0184766481201495e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[10])[0] = 1.1461654022399498e-01;
        (*gp_[10])[1] = 1.1461654022399498e-01;
        (*gp_[10])[2] = 6.5615037932801457e-01;
        
        weight_[11] = -9.0184766481201495e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[11])[0] = 1.1461654022399498e-01;
        (*gp_[11])[1] = 6.5615037932801457e-01;
        (*gp_[11])[2] = 1.1461654022399498e-01;
        
        weight_[12] = -9.0184766481201495e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[12])[0] = 6.5615037932801457e-01;
        (*gp_[12])[1] = 1.1461654022399498e-01;
        (*gp_[12])[2] = 1.1461654022399498e-01;
        
        weight_[13] = 4.4672576202511499e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[13])[0] = 2.2548995191151400e-01;
        (*gp_[13])[1] = 2.2548995191151400e-01;
        (*gp_[13])[2] = 2.2548995191151400e-01;
        
        weight_[14] = 4.4672576202511499e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[14])[0] = 2.2548995191151400e-01;
        (*gp_[14])[1] = 2.2548995191151400e-01;
        (*gp_[14])[2] = 3.2353014426545801e-01;
        
        weight_[15] = 4.4672576202511499e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[15])[0] = 2.2548995191151400e-01;
        (*gp_[15])[1] = 3.2353014426545801e-01;
        (*gp_[15])[2] = 2.2548995191151400e-01;
        
        weight_[16] = 4.4672576202511499e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[16])[0] = 3.2353014426545801e-01;
        (*gp_[16])[1] = 2.2548995191151400e-01;
        (*gp_[16])[2] = 2.2548995191151400e-01;
        
        weight_[17] = 3.4700405884550749e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[17])[0] = 1.3162780924687001e-01;
        (*gp_[17])[1] = 1.3162780924687001e-01;
        (*gp_[17])[2] = 8.3664701617184978e-02;
        
        weight_[18] = 3.4700405884550749e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[18])[0] = 1.3162780924687001e-01;
        (*gp_[18])[1] = 8.3664701617184978e-02;
        (*gp_[18])[2] = 1.3162780924687001e-01;
        
        weight_[19] = 3.4700405884550749e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[19])[0] = 1.3162780924687001e-01;
        (*gp_[19])[1] = 1.3162780924687001e-01;
        (*gp_[19])[2] = 6.5307967988907545e-01;
        
        weight_[20] = 3.4700405884550749e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[20])[0] = 1.3162780924687001e-01;
        (*gp_[20])[1] = 6.5307967988907545e-01;
        (*gp_[20])[2] = 1.3162780924687001e-01;
        
        weight_[21] = 3.4700405884550749e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[21])[0] = 1.3162780924687001e-01;
        (*gp_[21])[1] = 8.3664701617184978e-02;
        (*gp_[21])[2] = 6.5307967988907545e-01;
        
        weight_[22] = 3.4700405884550749e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[22])[0] = 1.3162780924687001e-01;
        (*gp_[22])[1] = 6.5307967988907545e-01;
        (*gp_[22])[2] = 8.3664701617184978e-02;
        
        weight_[23] = 3.4700405884550749e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[23])[0] = 8.3664701617184978e-02;
        (*gp_[23])[1] = 1.3162780924687001e-01;
        (*gp_[23])[2] = 1.3162780924687001e-01;
        
        weight_[24] = 3.4700405884550749e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[24])[0] = 8.3664701617184978e-02;
        (*gp_[24])[1] = 1.3162780924687001e-01;
        (*gp_[24])[2] = 6.5307967988907545e-01;
        
        weight_[25] = 3.4700405884550749e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[25])[0] = 8.3664701617184978e-02;
        (*gp_[25])[1] = 6.5307967988907545e-01;
        (*gp_[25])[2] = 1.3162780924687001e-01;
        
        weight_[26] = 3.4700405884550749e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[26])[0] = 6.5307967988907545e-01;
        (*gp_[26])[1] = 1.3162780924687001e-01;
        (*gp_[26])[2] = 8.3664701617184978e-02;
        
        weight_[27] = 3.4700405884550749e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[27])[0] = 6.5307967988907545e-01;
        (*gp_[27])[1] = 1.3162780924687001e-01;
        (*gp_[27])[2] = 1.3162780924687001e-01;
        
        weight_[28] = 3.4700405884550749e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[28])[0] = 6.5307967988907545e-01;
        (*gp_[28])[1] = 8.3664701617184978e-02;
        (*gp_[28])[2] = 1.3162780924687001e-01;
        
        weight_[29] = 3.3525839026606252e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[29])[0] = 4.3395146141140700e-01;
        (*gp_[29])[1] = 4.3395146141140700e-01;
        (*gp_[29])[2] = 1.0776985954942853e-01;
        
        weight_[30] = 3.3525839026606252e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[30])[0] = 4.3395146141140700e-01;
        (*gp_[30])[1] = 1.0776985954942853e-01;
        (*gp_[30])[2] = 4.3395146141140700e-01;
        
        weight_[31] = 3.3525839026606252e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[31])[0] = 4.3395146141140700e-01;
        (*gp_[31])[1] = 4.3395146141140700e-01;
        (*gp_[31])[2] = 2.4327217627758024e-02;
        
        weight_[32] = 3.3525839026606252e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[32])[0] = 4.3395146141140700e-01;
        (*gp_[32])[1] = 2.4327217627758024e-02;
        (*gp_[32])[2] = 4.3395146141140700e-01;
        
        weight_[33] = 3.3525839026606252e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[33])[0] = 4.3395146141140700e-01;
        (*gp_[33])[1] = 1.0776985954942853e-01;
        (*gp_[33])[2] = 2.4327217627758024e-02;
        
        weight_[34] = 3.3525839026606252e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[34])[0] = 4.3395146141140700e-01;
        (*gp_[34])[1] = 2.4327217627758024e-02;
        (*gp_[34])[2] = 1.0776985954942853e-01;
        
        weight_[35] = 3.3525839026606252e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[35])[0] = 1.0776985954942853e-01;
        (*gp_[35])[1] = 4.3395146141140700e-01;
        (*gp_[35])[2] = 4.3395146141140700e-01;
        
        weight_[36] = 3.3525839026606252e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[36])[0] = 1.0776985954942853e-01;
        (*gp_[36])[1] = 4.3395146141140700e-01;
        (*gp_[36])[2] = 2.4327217627758024e-02;
        
        weight_[37] = 3.3525839026606252e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[37])[0] = 1.0776985954942853e-01;
        (*gp_[37])[1] = 2.4327217627758024e-02;
        (*gp_[37])[2] = 4.3395146141140700e-01;
        
        weight_[38] = 3.3525839026606252e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[38])[0] = 2.4327217627758024e-02;
        (*gp_[38])[1] = 4.3395146141140700e-01;
        (*gp_[38])[2] = 1.0776985954942853e-01;
        
        weight_[39] = 3.3525839026606252e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[39])[0] = 2.4327217627758024e-02;
        (*gp_[39])[1] = 4.3395146141140700e-01;
        (*gp_[39])[2] = 4.3395146141140700e-01;
        
        weight_[40] = 3.3525839026606252e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[40])[0] = 2.4327217627758024e-02;
        (*gp_[40])[1] = 1.0776985954942853e-01;
        (*gp_[40])[2] = 4.3395146141140700e-01;
        
        weight_[41] = 4.3162887555699998e-04;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[41])[0] = -1.3762773181380528e-03;
        (*gp_[41])[1] = -1.3762773181380528e-03;
        (*gp_[41])[2] = 2.7655347263680752e-01;
        
        weight_[42] = 4.3162887555699998e-04;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[42])[0] = -1.3762773181380528e-03;
        (*gp_[42])[1] = 2.7655347263680752e-01;
        (*gp_[42])[2] = -1.3762773181380528e-03;
        
        weight_[43] = 4.3162887555699998e-04;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[43])[0] = -1.3762773181380528e-03;
        (*gp_[43])[1] = -1.3762773181380528e-03;
        (*gp_[43])[2] = 7.2619908199946903e-01;
        
        weight_[44] = 4.3162887555699998e-04;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[44])[0] = -1.3762773181380528e-03;
        (*gp_[44])[1] = 7.2619908199946903e-01;
        (*gp_[44])[2] = -1.3762773181380528e-03;
        
        weight_[45] = 4.3162887555699998e-04;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[45])[0] = -1.3762773181380528e-03;
        (*gp_[45])[1] = 2.7655347263680752e-01;
        (*gp_[45])[2] = 7.2619908199946903e-01;
        
        weight_[46] = 4.3162887555699998e-04;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[46])[0] = -1.3762773181380528e-03;
        (*gp_[46])[1] = 7.2619908199946903e-01;
        (*gp_[46])[2] = 2.7655347263680752e-01;
        
        weight_[47] = 4.3162887555699998e-04;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[47])[0] = 2.7655347263680752e-01;
        (*gp_[47])[1] = -1.3762773181380528e-03;
        (*gp_[47])[2] = -1.3762773181380528e-03;
        
        weight_[48] = 4.3162887555699998e-04;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[48])[0] = 2.7655347263680752e-01;
        (*gp_[48])[1] = -1.3762773181380528e-03;
        (*gp_[48])[2] = 7.2619908199946903e-01;
        
        weight_[49] = 4.3162887555699998e-04;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[49])[0] = 2.7655347263680752e-01;
        (*gp_[49])[1] = 7.2619908199946903e-01;
        (*gp_[49])[2] = -1.3762773181380528e-03;
        
        weight_[50] = 4.3162887555699998e-04;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[50])[0] = 7.2619908199946903e-01;
        (*gp_[50])[1] = -1.3762773181380528e-03;
        (*gp_[50])[2] = 2.7655347263680752e-01;
        
        weight_[51] = 4.3162887555699998e-04;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[51])[0] = 7.2619908199946903e-01;
        (*gp_[51])[1] = -1.3762773181380528e-03;
        (*gp_[51])[2] = -1.3762773181380528e-03;
        
        weight_[52] = 4.3162887555699998e-04;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[52])[0] = 7.2619908199946903e-01;
        (*gp_[52])[1] = 2.7655347263680752e-01;
        (*gp_[52])[2] = -1.3762773181380528e-03;
        
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
    
    std::size_t T3_10_126::nrOfPoints() const
    {
        return 126;
    }
    
    double T3_10_126::weight(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::getPoint - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReference&
    T3_10_126::getPoint(std::size_t i) const
    {
        logger.assert(i < nrOfPoints(), "%::weight - wrong index!", name_);
        return *gp_[i];
    }
    
    T3_10_126::ReferenceGeometryT*
    T3_10_126::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    T3_10_126::T3_10_126()
            : name_("T3_10_1"), refGeoPtr_(&ReferenceTetrahedron::Instance()), gp_(0)
    {
        weight_[0] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[0])[0] = 7.1428571428571508e-02;
        (*gp_[0])[1] = 7.1428571428571508e-02;
        (*gp_[0])[2] = 7.8571428571428548e-01;
        
        weight_[1] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[1])[0] = 7.1428571428571508e-02;
        (*gp_[1])[1] = 2.1428571428571452e-01;
        (*gp_[1])[2] = 6.4285714285714302e-01;
        
        weight_[2] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[2])[0] = 7.1428571428571508e-02;
        (*gp_[2])[1] = 3.5714285714285698e-01;
        (*gp_[2])[2] = 5.0000000000000000e-01;
        
        weight_[3] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[3])[0] = 7.1428571428571508e-02;
        (*gp_[3])[1] = 5.0000000000000000e-01;
        (*gp_[3])[2] = 3.5714285714285698e-01;
        
        weight_[4] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[4])[0] = 7.1428571428571508e-02;
        (*gp_[4])[1] = 6.4285714285714302e-01;
        (*gp_[4])[2] = 2.1428571428571452e-01;
        
        weight_[5] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[5])[0] = 7.1428571428571508e-02;
        (*gp_[5])[1] = 7.8571428571428548e-01;
        (*gp_[5])[2] = 7.1428571428571508e-02;
        
        weight_[6] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[6])[0] = 2.1428571428571452e-01;
        (*gp_[6])[1] = 7.1428571428571508e-02;
        (*gp_[6])[2] = 6.4285714285714302e-01;
        
        weight_[7] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[7])[0] = 2.1428571428571452e-01;
        (*gp_[7])[1] = 2.1428571428571452e-01;
        (*gp_[7])[2] = 5.0000000000000000e-01;
        
        weight_[8] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[8])[0] = 2.1428571428571452e-01;
        (*gp_[8])[1] = 3.5714285714285698e-01;
        (*gp_[8])[2] = 3.5714285714285698e-01;
        
        weight_[9] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[9])[0] = 2.1428571428571452e-01;
        (*gp_[9])[1] = 5.0000000000000000e-01;
        (*gp_[9])[2] = 2.1428571428571452e-01;
        
        weight_[10] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[10])[0] = 2.1428571428571452e-01;
        (*gp_[10])[1] = 6.4285714285714302e-01;
        (*gp_[10])[2] = 7.1428571428571508e-02;
        
        weight_[11] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[11])[0] = 3.5714285714285698e-01;
        (*gp_[11])[1] = 7.1428571428571508e-02;
        (*gp_[11])[2] = 5.0000000000000000e-01;
        
        weight_[12] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[12])[0] = 3.5714285714285698e-01;
        (*gp_[12])[1] = 2.1428571428571452e-01;
        (*gp_[12])[2] = 3.5714285714285698e-01;
        
        weight_[13] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[13])[0] = 3.5714285714285698e-01;
        (*gp_[13])[1] = 3.5714285714285698e-01;
        (*gp_[13])[2] = 2.1428571428571452e-01;
        
        weight_[14] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[14])[0] = 3.5714285714285698e-01;
        (*gp_[14])[1] = 5.0000000000000000e-01;
        (*gp_[14])[2] = 7.1428571428571508e-02;
        
        weight_[15] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[15])[0] = 5.0000000000000000e-01;
        (*gp_[15])[1] = 7.1428571428571508e-02;
        (*gp_[15])[2] = 3.5714285714285698e-01;
        
        weight_[16] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[16])[0] = 5.0000000000000000e-01;
        (*gp_[16])[1] = 2.1428571428571452e-01;
        (*gp_[16])[2] = 2.1428571428571452e-01;
        
        weight_[17] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[17])[0] = 5.0000000000000000e-01;
        (*gp_[17])[1] = 3.5714285714285698e-01;
        (*gp_[17])[2] = 7.1428571428571508e-02;
        
        weight_[18] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[18])[0] = 6.4285714285714302e-01;
        (*gp_[18])[1] = 7.1428571428571508e-02;
        (*gp_[18])[2] = 2.1428571428571452e-01;
        
        weight_[19] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[19])[0] = 6.4285714285714302e-01;
        (*gp_[19])[1] = 2.1428571428571452e-01;
        (*gp_[19])[2] = 7.1428571428571508e-02;
        
        weight_[20] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[20])[0] = 7.8571428571428548e-01;
        (*gp_[20])[1] = 7.1428571428571508e-02;
        (*gp_[20])[2] = 7.1428571428571508e-02;
        
        weight_[21] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[21])[0] = 7.1428571428571508e-02;
        (*gp_[21])[1] = 7.1428571428571508e-02;
        (*gp_[21])[2] = 6.4285714285714302e-01;
        
        weight_[22] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[22])[0] = 7.1428571428571508e-02;
        (*gp_[22])[1] = 2.1428571428571452e-01;
        (*gp_[22])[2] = 5.0000000000000000e-01;
        
        weight_[23] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[23])[0] = 7.1428571428571508e-02;
        (*gp_[23])[1] = 3.5714285714285698e-01;
        (*gp_[23])[2] = 3.5714285714285698e-01;
        
        weight_[24] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[24])[0] = 7.1428571428571508e-02;
        (*gp_[24])[1] = 5.0000000000000000e-01;
        (*gp_[24])[2] = 2.1428571428571452e-01;
        
        weight_[25] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[25])[0] = 7.1428571428571508e-02;
        (*gp_[25])[1] = 6.4285714285714302e-01;
        (*gp_[25])[2] = 7.1428571428571508e-02;
        
        weight_[26] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[26])[0] = 2.1428571428571452e-01;
        (*gp_[26])[1] = 7.1428571428571508e-02;
        (*gp_[26])[2] = 5.0000000000000000e-01;
        
        weight_[27] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[27])[0] = 2.1428571428571452e-01;
        (*gp_[27])[1] = 2.1428571428571452e-01;
        (*gp_[27])[2] = 3.5714285714285698e-01;
        
        weight_[28] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[28])[0] = 2.1428571428571452e-01;
        (*gp_[28])[1] = 3.5714285714285698e-01;
        (*gp_[28])[2] = 2.1428571428571452e-01;
        
        weight_[29] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[29])[0] = 2.1428571428571452e-01;
        (*gp_[29])[1] = 5.0000000000000000e-01;
        (*gp_[29])[2] = 7.1428571428571508e-02;
        
        weight_[30] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[30])[0] = 3.5714285714285698e-01;
        (*gp_[30])[1] = 7.1428571428571508e-02;
        (*gp_[30])[2] = 3.5714285714285698e-01;
        
        weight_[31] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[31])[0] = 3.5714285714285698e-01;
        (*gp_[31])[1] = 2.1428571428571452e-01;
        (*gp_[31])[2] = 2.1428571428571452e-01;
        
        weight_[32] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[32])[0] = 3.5714285714285698e-01;
        (*gp_[32])[1] = 3.5714285714285698e-01;
        (*gp_[32])[2] = 7.1428571428571508e-02;
        
        weight_[33] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[33])[0] = 5.0000000000000000e-01;
        (*gp_[33])[1] = 7.1428571428571508e-02;
        (*gp_[33])[2] = 2.1428571428571452e-01;
        
        weight_[34] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[34])[0] = 5.0000000000000000e-01;
        (*gp_[34])[1] = 2.1428571428571452e-01;
        (*gp_[34])[2] = 7.1428571428571508e-02;
        
        weight_[35] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[35])[0] = 6.4285714285714302e-01;
        (*gp_[35])[1] = 7.1428571428571508e-02;
        (*gp_[35])[2] = 7.1428571428571508e-02;
        
        weight_[36] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[36])[0] = 7.1428571428571508e-02;
        (*gp_[36])[1] = 7.1428571428571508e-02;
        (*gp_[36])[2] = 5.0000000000000000e-01;
        
        weight_[37] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[37])[0] = 7.1428571428571508e-02;
        (*gp_[37])[1] = 2.1428571428571452e-01;
        (*gp_[37])[2] = 3.5714285714285698e-01;
        
        weight_[38] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[38])[0] = 7.1428571428571508e-02;
        (*gp_[38])[1] = 3.5714285714285698e-01;
        (*gp_[38])[2] = 2.1428571428571452e-01;
        
        weight_[39] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[39])[0] = 7.1428571428571508e-02;
        (*gp_[39])[1] = 5.0000000000000000e-01;
        (*gp_[39])[2] = 7.1428571428571508e-02;
        
        weight_[40] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[40])[0] = 2.1428571428571452e-01;
        (*gp_[40])[1] = 7.1428571428571508e-02;
        (*gp_[40])[2] = 3.5714285714285698e-01;
        
        weight_[41] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[41])[0] = 2.1428571428571452e-01;
        (*gp_[41])[1] = 2.1428571428571452e-01;
        (*gp_[41])[2] = 2.1428571428571452e-01;
        
        weight_[42] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[42])[0] = 2.1428571428571452e-01;
        (*gp_[42])[1] = 3.5714285714285698e-01;
        (*gp_[42])[2] = 7.1428571428571508e-02;
        
        weight_[43] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[43])[0] = 3.5714285714285698e-01;
        (*gp_[43])[1] = 7.1428571428571508e-02;
        (*gp_[43])[2] = 2.1428571428571452e-01;
        
        weight_[44] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[44])[0] = 3.5714285714285698e-01;
        (*gp_[44])[1] = 2.1428571428571452e-01;
        (*gp_[44])[2] = 7.1428571428571508e-02;
        
        weight_[45] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[45])[0] = 5.0000000000000000e-01;
        (*gp_[45])[1] = 7.1428571428571508e-02;
        (*gp_[45])[2] = 7.1428571428571508e-02;
        
        weight_[46] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[46])[0] = 7.1428571428571508e-02;
        (*gp_[46])[1] = 7.1428571428571508e-02;
        (*gp_[46])[2] = 3.5714285714285698e-01;
        
        weight_[47] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[47])[0] = 7.1428571428571508e-02;
        (*gp_[47])[1] = 2.1428571428571452e-01;
        (*gp_[47])[2] = 2.1428571428571452e-01;
        
        weight_[48] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[48])[0] = 7.1428571428571508e-02;
        (*gp_[48])[1] = 3.5714285714285698e-01;
        (*gp_[48])[2] = 7.1428571428571508e-02;
        
        weight_[49] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[49])[0] = 2.1428571428571452e-01;
        (*gp_[49])[1] = 7.1428571428571508e-02;
        (*gp_[49])[2] = 2.1428571428571452e-01;
        
        weight_[50] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[50])[0] = 2.1428571428571452e-01;
        (*gp_[50])[1] = 2.1428571428571452e-01;
        (*gp_[50])[2] = 7.1428571428571508e-02;
        
        weight_[51] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[51])[0] = 3.5714285714285698e-01;
        (*gp_[51])[1] = 7.1428571428571508e-02;
        (*gp_[51])[2] = 7.1428571428571508e-02;
        
        weight_[52] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[52])[0] = 7.1428571428571508e-02;
        (*gp_[52])[1] = 7.1428571428571508e-02;
        (*gp_[52])[2] = 2.1428571428571452e-01;
        
        weight_[53] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[53])[0] = 7.1428571428571508e-02;
        (*gp_[53])[1] = 2.1428571428571452e-01;
        (*gp_[53])[2] = 7.1428571428571508e-02;
        
        weight_[54] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[54])[0] = 2.1428571428571452e-01;
        (*gp_[54])[1] = 7.1428571428571508e-02;
        (*gp_[54])[2] = 7.1428571428571508e-02;
        
        weight_[55] = 4.5362824065080999e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[55])[0] = 7.1428571428571508e-02;
        (*gp_[55])[1] = 7.1428571428571508e-02;
        (*gp_[55])[2] = 7.1428571428571508e-02;
        
        weight_[56] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[56])[0] = 8.3333333333333481e-02;
        (*gp_[56])[1] = 8.3333333333333481e-02;
        (*gp_[56])[2] = 7.5000000000000000e-01;
        
        weight_[57] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[57])[0] = 8.3333333333333481e-02;
        (*gp_[57])[1] = 2.5000000000000000e-01;
        (*gp_[57])[2] = 5.8333333333333348e-01;
        
        weight_[58] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[58])[0] = 8.3333333333333481e-02;
        (*gp_[58])[1] = 4.1666666666666652e-01;
        (*gp_[58])[2] = 4.1666666666666652e-01;
        
        weight_[59] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[59])[0] = 8.3333333333333481e-02;
        (*gp_[59])[1] = 5.8333333333333348e-01;
        (*gp_[59])[2] = 2.5000000000000000e-01;
        
        weight_[60] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[60])[0] = 8.3333333333333481e-02;
        (*gp_[60])[1] = 7.5000000000000000e-01;
        (*gp_[60])[2] = 8.3333333333333481e-02;
        
        weight_[61] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[61])[0] = 2.5000000000000000e-01;
        (*gp_[61])[1] = 8.3333333333333481e-02;
        (*gp_[61])[2] = 5.8333333333333348e-01;
        
        weight_[62] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[62])[0] = 2.5000000000000000e-01;
        (*gp_[62])[1] = 2.5000000000000000e-01;
        (*gp_[62])[2] = 4.1666666666666652e-01;
        
        weight_[63] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[63])[0] = 2.5000000000000000e-01;
        (*gp_[63])[1] = 4.1666666666666652e-01;
        (*gp_[63])[2] = 2.5000000000000000e-01;
        
        weight_[64] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[64])[0] = 2.5000000000000000e-01;
        (*gp_[64])[1] = 5.8333333333333348e-01;
        (*gp_[64])[2] = 8.3333333333333481e-02;
        
        weight_[65] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[65])[0] = 4.1666666666666652e-01;
        (*gp_[65])[1] = 8.3333333333333481e-02;
        (*gp_[65])[2] = 4.1666666666666652e-01;
        
        weight_[66] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[66])[0] = 4.1666666666666652e-01;
        (*gp_[66])[1] = 2.5000000000000000e-01;
        (*gp_[66])[2] = 2.5000000000000000e-01;
        
        weight_[67] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[67])[0] = 4.1666666666666652e-01;
        (*gp_[67])[1] = 4.1666666666666652e-01;
        (*gp_[67])[2] = 8.3333333333333481e-02;
        
        weight_[68] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[68])[0] = 5.8333333333333348e-01;
        (*gp_[68])[1] = 8.3333333333333481e-02;
        (*gp_[68])[2] = 2.5000000000000000e-01;
        
        weight_[69] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[69])[0] = 5.8333333333333348e-01;
        (*gp_[69])[1] = 2.5000000000000000e-01;
        (*gp_[69])[2] = 8.3333333333333481e-02;
        
        weight_[70] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[70])[0] = 7.5000000000000000e-01;
        (*gp_[70])[1] = 8.3333333333333481e-02;
        (*gp_[70])[2] = 8.3333333333333481e-02;
        
        weight_[71] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[71])[0] = 8.3333333333333481e-02;
        (*gp_[71])[1] = 8.3333333333333481e-02;
        (*gp_[71])[2] = 5.8333333333333348e-01;
        
        weight_[72] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[72])[0] = 8.3333333333333481e-02;
        (*gp_[72])[1] = 2.5000000000000000e-01;
        (*gp_[72])[2] = 4.1666666666666652e-01;
        
        weight_[73] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[73])[0] = 8.3333333333333481e-02;
        (*gp_[73])[1] = 4.1666666666666652e-01;
        (*gp_[73])[2] = 2.5000000000000000e-01;
        
        weight_[74] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[74])[0] = 8.3333333333333481e-02;
        (*gp_[74])[1] = 5.8333333333333348e-01;
        (*gp_[74])[2] = 8.3333333333333481e-02;
        
        weight_[75] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[75])[0] = 2.5000000000000000e-01;
        (*gp_[75])[1] = 8.3333333333333481e-02;
        (*gp_[75])[2] = 4.1666666666666652e-01;
        
        weight_[76] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[76])[0] = 2.5000000000000000e-01;
        (*gp_[76])[1] = 2.5000000000000000e-01;
        (*gp_[76])[2] = 2.5000000000000000e-01;
        
        weight_[77] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[77])[0] = 2.5000000000000000e-01;
        (*gp_[77])[1] = 4.1666666666666652e-01;
        (*gp_[77])[2] = 8.3333333333333481e-02;
        
        weight_[78] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[78])[0] = 4.1666666666666652e-01;
        (*gp_[78])[1] = 8.3333333333333481e-02;
        (*gp_[78])[2] = 2.5000000000000000e-01;
        
        weight_[79] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[79])[0] = 4.1666666666666652e-01;
        (*gp_[79])[1] = 2.5000000000000000e-01;
        (*gp_[79])[2] = 8.3333333333333481e-02;
        
        weight_[80] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[80])[0] = 5.8333333333333348e-01;
        (*gp_[80])[1] = 8.3333333333333481e-02;
        (*gp_[80])[2] = 8.3333333333333481e-02;
        
        weight_[81] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[81])[0] = 8.3333333333333481e-02;
        (*gp_[81])[1] = 8.3333333333333481e-02;
        (*gp_[81])[2] = 4.1666666666666652e-01;
        
        weight_[82] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[82])[0] = 8.3333333333333481e-02;
        (*gp_[82])[1] = 2.5000000000000000e-01;
        (*gp_[82])[2] = 2.5000000000000000e-01;
        
        weight_[83] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[83])[0] = 8.3333333333333481e-02;
        (*gp_[83])[1] = 4.1666666666666652e-01;
        (*gp_[83])[2] = 8.3333333333333481e-02;
        
        weight_[84] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[84])[0] = 2.5000000000000000e-01;
        (*gp_[84])[1] = 8.3333333333333481e-02;
        (*gp_[84])[2] = 2.5000000000000000e-01;
        
        weight_[85] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[85])[0] = 2.5000000000000000e-01;
        (*gp_[85])[1] = 2.5000000000000000e-01;
        (*gp_[85])[2] = 8.3333333333333481e-02;
        
        weight_[86] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[86])[0] = 4.1666666666666652e-01;
        (*gp_[86])[1] = 8.3333333333333481e-02;
        (*gp_[86])[2] = 8.3333333333333481e-02;
        
        weight_[87] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[87])[0] = 8.3333333333333481e-02;
        (*gp_[87])[1] = 8.3333333333333481e-02;
        (*gp_[87])[2] = 2.5000000000000000e-01;
        
        weight_[88] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[88])[0] = 8.3333333333333481e-02;
        (*gp_[88])[1] = 2.5000000000000000e-01;
        (*gp_[88])[2] = 8.3333333333333481e-02;
        
        weight_[89] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[89])[0] = 2.5000000000000000e-01;
        (*gp_[89])[1] = 8.3333333333333481e-02;
        (*gp_[89])[2] = 8.3333333333333481e-02;
        
        weight_[90] = -1.1652347652347650e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[90])[0] = 8.3333333333333481e-02;
        (*gp_[90])[1] = 8.3333333333333481e-02;
        (*gp_[90])[2] = 8.3333333333333481e-02;
        
        weight_[91] = 1.0193728997982475e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[91])[0] = 9.9999999999999978e-02;
        (*gp_[91])[1] = 9.9999999999999978e-02;
        (*gp_[91])[2] = 6.9999999999999996e-01;
        
        weight_[92] = 1.0193728997982475e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[92])[0] = 9.9999999999999978e-02;
        (*gp_[92])[1] = 2.9999999999999999e-01;
        (*gp_[92])[2] = 5.0000000000000000e-01;
        
        weight_[93] = 1.0193728997982475e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[93])[0] = 9.9999999999999978e-02;
        (*gp_[93])[1] = 5.0000000000000000e-01;
        (*gp_[93])[2] = 2.9999999999999999e-01;
        
        weight_[94] = 1.0193728997982475e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[94])[0] = 9.9999999999999978e-02;
        (*gp_[94])[1] = 6.9999999999999996e-01;
        (*gp_[94])[2] = 9.9999999999999978e-02;
        
        weight_[95] = 1.0193728997982475e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[95])[0] = 2.9999999999999999e-01;
        (*gp_[95])[1] = 9.9999999999999978e-02;
        (*gp_[95])[2] = 5.0000000000000000e-01;
        
        weight_[96] = 1.0193728997982475e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[96])[0] = 2.9999999999999999e-01;
        (*gp_[96])[1] = 2.9999999999999999e-01;
        (*gp_[96])[2] = 2.9999999999999999e-01;
        
        weight_[97] = 1.0193728997982475e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[97])[0] = 2.9999999999999999e-01;
        (*gp_[97])[1] = 5.0000000000000000e-01;
        (*gp_[97])[2] = 9.9999999999999978e-02;
        
        weight_[98] = 1.0193728997982475e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[98])[0] = 5.0000000000000000e-01;
        (*gp_[98])[1] = 9.9999999999999978e-02;
        (*gp_[98])[2] = 2.9999999999999999e-01;
        
        weight_[99] = 1.0193728997982475e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[99])[0] = 5.0000000000000000e-01;
        (*gp_[99])[1] = 2.9999999999999999e-01;
        (*gp_[99])[2] = 9.9999999999999978e-02;
        
        weight_[100] = 1.0193728997982475e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[100])[0] = 6.9999999999999996e-01;
        (*gp_[100])[1] = 9.9999999999999978e-02;
        (*gp_[100])[2] = 9.9999999999999978e-02;
        
        weight_[101] = 1.0193728997982475e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[101])[0] = 9.9999999999999978e-02;
        (*gp_[101])[1] = 9.9999999999999978e-02;
        (*gp_[101])[2] = 5.0000000000000000e-01;
        
        weight_[102] = 1.0193728997982475e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[102])[0] = 9.9999999999999978e-02;
        (*gp_[102])[1] = 2.9999999999999999e-01;
        (*gp_[102])[2] = 2.9999999999999999e-01;
        
        weight_[103] = 1.0193728997982475e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[103])[0] = 9.9999999999999978e-02;
        (*gp_[103])[1] = 5.0000000000000000e-01;
        (*gp_[103])[2] = 9.9999999999999978e-02;
        
        weight_[104] = 1.0193728997982475e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[104])[0] = 2.9999999999999999e-01;
        (*gp_[104])[1] = 9.9999999999999978e-02;
        (*gp_[104])[2] = 2.9999999999999999e-01;
        
        weight_[105] = 1.0193728997982475e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[105])[0] = 2.9999999999999999e-01;
        (*gp_[105])[1] = 2.9999999999999999e-01;
        (*gp_[105])[2] = 9.9999999999999978e-02;
        
        weight_[106] = 1.0193728997982475e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[106])[0] = 5.0000000000000000e-01;
        (*gp_[106])[1] = 9.9999999999999978e-02;
        (*gp_[106])[2] = 9.9999999999999978e-02;
        
        weight_[107] = 1.0193728997982475e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[107])[0] = 9.9999999999999978e-02;
        (*gp_[107])[1] = 9.9999999999999978e-02;
        (*gp_[107])[2] = 2.9999999999999999e-01;
        
        weight_[108] = 1.0193728997982475e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[108])[0] = 9.9999999999999978e-02;
        (*gp_[108])[1] = 2.9999999999999999e-01;
        (*gp_[108])[2] = 9.9999999999999978e-02;
        
        weight_[109] = 1.0193728997982475e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[109])[0] = 2.9999999999999999e-01;
        (*gp_[109])[1] = 9.9999999999999978e-02;
        (*gp_[109])[2] = 9.9999999999999978e-02;
        
        weight_[110] = 1.0193728997982475e-01;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[110])[0] = 9.9999999999999978e-02;
        (*gp_[110])[1] = 9.9999999999999978e-02;
        (*gp_[110])[2] = 9.9999999999999978e-02;
        
        weight_[111] = -3.5025386136497250e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[111])[0] = 1.2500000000000000e-01;
        (*gp_[111])[1] = 1.2500000000000000e-01;
        (*gp_[111])[2] = 6.2500000000000000e-01;
        
        weight_[112] = -3.5025386136497250e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[112])[0] = 1.2500000000000000e-01;
        (*gp_[112])[1] = 3.7500000000000000e-01;
        (*gp_[112])[2] = 3.7500000000000000e-01;
        
        weight_[113] = -3.5025386136497250e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[113])[0] = 1.2500000000000000e-01;
        (*gp_[113])[1] = 6.2500000000000000e-01;
        (*gp_[113])[2] = 1.2500000000000000e-01;
        
        weight_[114] = -3.5025386136497250e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[114])[0] = 3.7500000000000000e-01;
        (*gp_[114])[1] = 1.2500000000000000e-01;
        (*gp_[114])[2] = 3.7500000000000000e-01;
        
        weight_[115] = -3.5025386136497250e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[115])[0] = 3.7500000000000000e-01;
        (*gp_[115])[1] = 3.7500000000000000e-01;
        (*gp_[115])[2] = 1.2500000000000000e-01;
        
        weight_[116] = -3.5025386136497250e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[116])[0] = 6.2500000000000000e-01;
        (*gp_[116])[1] = 1.2500000000000000e-01;
        (*gp_[116])[2] = 1.2500000000000000e-01;
        
        weight_[117] = -3.5025386136497250e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[117])[0] = 1.2500000000000000e-01;
        (*gp_[117])[1] = 1.2500000000000000e-01;
        (*gp_[117])[2] = 3.7500000000000000e-01;
        
        weight_[118] = -3.5025386136497250e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[118])[0] = 1.2500000000000000e-01;
        (*gp_[118])[1] = 3.7500000000000000e-01;
        (*gp_[118])[2] = 1.2500000000000000e-01;
        
        weight_[119] = -3.5025386136497250e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[119])[0] = 3.7500000000000000e-01;
        (*gp_[119])[1] = 1.2500000000000000e-01;
        (*gp_[119])[2] = 1.2500000000000000e-01;
        
        weight_[120] = -3.5025386136497250e-02;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[120])[0] = 1.2500000000000000e-01;
        (*gp_[120])[1] = 1.2500000000000000e-01;
        (*gp_[120])[2] = 1.2500000000000000e-01;
        
        weight_[121] = 4.0680803571428751e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[121])[0] = 1.6666666666666652e-01;
        (*gp_[121])[1] = 1.6666666666666652e-01;
        (*gp_[121])[2] = 5.0000000000000000e-01;
        
        weight_[122] = 4.0680803571428751e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[122])[0] = 1.6666666666666652e-01;
        (*gp_[122])[1] = 5.0000000000000000e-01;
        (*gp_[122])[2] = 1.6666666666666652e-01;
        
        weight_[123] = 4.0680803571428751e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[123])[0] = 5.0000000000000000e-01;
        (*gp_[123])[1] = 1.6666666666666652e-01;
        (*gp_[123])[2] = 1.6666666666666652e-01;
        
        weight_[124] = 4.0680803571428751e-03;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[124])[0] = 1.6666666666666652e-01;
        (*gp_[124])[1] = 1.6666666666666652e-01;
        (*gp_[124])[2] = 1.6666666666666652e-01;
        
        weight_[125] = -9.4062316284499996e-05;
        gp_.push_back(new Geometry::PointReference(3));
        (*gp_[125])[0] = 2.5000000000000000e-01;
        (*gp_[125])[1] = 2.5000000000000000e-01;
        (*gp_[125])[2] = 2.5000000000000000e-01;
        
    }
//---------------------------------------------------------------------------
}// close namespace IntegrationRules

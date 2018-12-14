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
#include "Integration/QuadratureRules/GaussQuadratureRulesForTriangle.h"
#include "Geometry/ReferenceTriangle.h"
#include "Geometry/PointReference.h"
//---------------------------------------------------------------------------
namespace QuadratureRules
{
//---------------------------------------------------------------------------
    std::string Tn2_1_1::getName() const
    {
        return name_;
    }
    
    std::size_t Tn2_1_1::order() const
    {
        return 1;
    }
    
    std::size_t Tn2_1_1::dimension() const
    {
        return 2;
    }
    
    std::size_t Tn2_1_1::getNumberOfPoints() const
    {
        return 1;
    }
    
    double Tn2_1_1::weight(std::size_t i) const
    {
        logger.assert_debug(i < getNumberOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReferenceBase&
    Tn2_1_1::getPoint(std::size_t i) const
    {
        logger.assert_debug(i < getNumberOfPoints(), "%::getPoint - wrong index!", name_);
        return gp_[i];
    }
    
    Geometry::ReferenceGeometry*
    Tn2_1_1::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    Tn2_1_1::Tn2_1_1()
            : name_("Tn2_1_1"), refGeoPtr_(&Geometry::ReferenceTriangle::Instance()), gp_(0)
    {
        weight_[0] = 5.0000000000000000e-01;
        gp_.push_back({3.3333333333333348e-01, 3.3333333333333348e-01});
        
    }
    
//---------------------------------------------------------------------------
    std::string Tn2_2_3::getName() const
    {
        return name_;
    }
    
    std::size_t Tn2_2_3::order() const
    {
        return 2;
    }
    
    std::size_t Tn2_2_3::dimension() const
    {
        return 2;
    }
    
    std::size_t Tn2_2_3::getNumberOfPoints() const
    {
        return 3;
    }
    
    double Tn2_2_3::weight(std::size_t i) const
    {
        logger.assert_debug(i < getNumberOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReferenceBase&
    Tn2_2_3::getPoint(std::size_t i) const
    {
        logger.assert_debug(i < getNumberOfPoints(), "%::getPoint - wrong index!", name_);
        return gp_[i];
    }
    
    Geometry::ReferenceGeometry*
    Tn2_2_3::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    Tn2_2_3::Tn2_2_3()
            : name_("Tn2_2_1"), refGeoPtr_(&Geometry::ReferenceTriangle::Instance()), gp_(0)
    {
        weight_[0] = 1.6666666666666674e-01;
        gp_.push_back({1.6666666666666652e-01, 1.6666666666666652e-01});
        
        weight_[1] = 1.6666666666666674e-01;
        gp_.push_back({1.6666666666666652e-01, 6.6666666666666652e-01});
        
        weight_[2] = 1.6666666666666674e-01;
        gp_.push_back({6.6666666666666652e-01, 1.6666666666666652e-01});
        
    }
    
//---------------------------------------------------------------------------
    std::string Tn2_3_4::getName() const
    {
        return name_;
    }
    
    std::size_t Tn2_3_4::order() const
    {
        return 3;
    }
    
    std::size_t Tn2_3_4::dimension() const
    {
        return 2;
    }
    
    std::size_t Tn2_3_4::getNumberOfPoints() const
    {
        return 4;
    }
    
    double Tn2_3_4::weight(std::size_t i) const
    {
        logger.assert_debug(i < getNumberOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReferenceBase&
    Tn2_3_4::getPoint(std::size_t i) const
    {
        logger.assert_debug(i < getNumberOfPoints(), "%::getPoint - wrong index!", name_);
        return gp_[i];
    }
    
    Geometry::ReferenceGeometry*
    Tn2_3_4::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    Tn2_3_4::Tn2_3_4()
            : name_("Tn2_3_1"), refGeoPtr_(&Geometry::ReferenceTriangle::Instance()), gp_(0)
    {
        weight_[0] = -2.8125000000000000e-01;
        gp_.push_back({3.3333333333333348e-01, 3.3333333333333348e-01});

        weight_[1] = 2.6041666666666674e-01;
        gp_.push_back({2.0000000000000001e-01, 2.0000000000000001e-01});
        
        weight_[2] = 2.6041666666666674e-01;
        gp_.push_back({2.0000000000000001e-01, 5.9999999999999998e-01});
        
        weight_[3] = 2.6041666666666674e-01;
        gp_.push_back({5.9999999999999998e-01, 2.0000000000000001e-01});
        
    }
    
//---------------------------------------------------------------------------
    std::string Tn2_4_6::getName() const
    {
        return name_;
    }
    
    std::size_t Tn2_4_6::order() const
    {
        return 4;
    }
    
    std::size_t Tn2_4_6::dimension() const
    {
        return 2;
    }
    
    std::size_t Tn2_4_6::getNumberOfPoints() const
    {
        return 6;
    }
    
    double Tn2_4_6::weight(std::size_t i) const
    {
        logger.assert_debug(i < getNumberOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReferenceBase&
    Tn2_4_6::getPoint(std::size_t i) const
    {
        logger.assert_debug(i < getNumberOfPoints(), "%::getPoint - wrong index!", name_);
        return gp_[i];
    }
    
    Geometry::ReferenceGeometry*
    Tn2_4_6::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    Tn2_4_6::Tn2_4_6()
            : name_("Tn2_4_1"), refGeoPtr_(&Geometry::ReferenceTriangle::Instance()), gp_(0)
    {
        weight_[0] = 1.1169079483900550e-01;
        gp_.push_back({4.4594849091596500e-01, 4.4594849091596500e-01});
        
        weight_[1] = 1.1169079483900550e-01;
        gp_.push_back({4.4594849091596500e-01, 1.0810301816807000e-01});
        
        weight_[2] = 1.1169079483900550e-01;
        gp_.push_back({1.0810301816807000e-01, 4.4594849091596500e-01});
        
        weight_[3] = 5.4975871827660998e-02;
        gp_.push_back({9.1576213509771021e-02, 9.1576213509771021e-02});
        
        weight_[4] = 5.4975871827660998e-02;
        gp_.push_back({9.1576213509771021e-02, 8.1684757298045896e-01});
        
        weight_[5] = 5.4975871827660998e-02;
        gp_.push_back({8.1684757298045896e-01, 9.1576213509771021e-02});
        
    }
    
//---------------------------------------------------------------------------
    std::string T2_5_7::getName() const
    {
        return name_;
    }
    
    std::size_t T2_5_7::order() const
    {
        return 5;
    }
    
    std::size_t T2_5_7::dimension() const
    {
        return 2;
    }
    
    std::size_t T2_5_7::getNumberOfPoints() const
    {
        return 7;
    }
    
    double T2_5_7::weight(std::size_t i) const
    {
        logger.assert_debug(i < getNumberOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReferenceBase&
    T2_5_7::getPoint(std::size_t i) const
    {
        logger.assert_debug(i < getNumberOfPoints(), "%::getPoint - wrong index!", name_);
        return gp_[i];
    }
    
    Geometry::ReferenceGeometry*
    T2_5_7::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    T2_5_7::T2_5_7()
            : name_("T2_5_1"), refGeoPtr_(&Geometry::ReferenceTriangle::Instance()), gp_(0)
    {
        weight_[0] = 1.1250000000000000e-01;
        gp_.push_back({3.3333333333333348e-01, 3.3333333333333348e-01});
        
        weight_[1] = 6.6197076394252999e-02;
        gp_.push_back({4.7014206410511500e-01, 4.7014206410511500e-01});
        
        weight_[2] = 6.6197076394252999e-02;
        gp_.push_back({4.7014206410511500e-01, 5.9715871789770003e-02});
        
        weight_[3] = 6.6197076394252999e-02;
        gp_.push_back({5.9715871789770003e-02, 4.7014206410511500e-01});
        
        weight_[4] = 6.2969590272413500e-02;
        gp_.push_back({1.0128650732345601e-01, 1.0128650732345601e-01});
        
        weight_[5] = 6.2969590272413500e-02;
        gp_.push_back({1.0128650732345601e-01, 7.9742698535308698e-01});
        
        weight_[6] = 6.2969590272413500e-02;
        gp_.push_back({7.9742698535308698e-01, 1.0128650732345601e-01});
        
    }
    
//---------------------------------------------------------------------------
    std::string T2_6_12::getName() const
    {
        return name_;
    }
    
    std::size_t T2_6_12::order() const
    {
        return 6;
    }
    
    std::size_t T2_6_12::dimension() const
    {
        return 2;
    }
    
    std::size_t T2_6_12::getNumberOfPoints() const
    {
        return 12;
    }
    
    double T2_6_12::weight(std::size_t i) const
    {
        logger.assert_debug(i < getNumberOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReferenceBase&
    T2_6_12::getPoint(std::size_t i) const
    {
        logger.assert_debug(i < getNumberOfPoints(), "%::getPoint - wrong index!", name_);
        return gp_[i];
    }
    
    Geometry::ReferenceGeometry*
    T2_6_12::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    T2_6_12::T2_6_12()
            : name_("T2_6_1"), refGeoPtr_(&Geometry::ReferenceTriangle::Instance()), gp_(0)
    {
        weight_[0] = 5.8393137863189497e-02;
        gp_.push_back({2.4928674517090998e-01, 2.4928674517090998e-01});
        
        weight_[1] = 5.8393137863189497e-02;
        gp_.push_back({2.4928674517090998e-01, 5.0142650965817903e-01});
        
        weight_[2] = 5.8393137863189497e-02;
        gp_.push_back({5.0142650965817903e-01, 2.4928674517090998e-01});
        
        weight_[3] = 2.5422453185103500e-02;
        gp_.push_back({6.3089014491502005e-02, 6.3089014491502005e-02});
        
        weight_[4] = 2.5422453185103500e-02;
        gp_.push_back({6.3089014491502005e-02, 8.7382197101699599e-01});
        
        weight_[5] = 2.5422453185103500e-02;
        gp_.push_back({8.7382197101699599e-01, 6.3089014491502005e-02});
        
        weight_[6] = 4.1425537809187001e-02;
        gp_.push_back({3.1035245103378400e-01, 6.3650249912139900e-01});
        
        weight_[7] = 4.1425537809187001e-02;
        gp_.push_back({6.3650249912139900e-01, 5.3145049844816994e-02});
        
        weight_[8] = 4.1425537809187001e-02;
        gp_.push_back({5.3145049844816994e-02, 3.1035245103378400e-01});
        
        weight_[9] = 4.1425537809187001e-02;
        gp_.push_back({3.1035245103378400e-01, 5.3145049844816994e-02});
        
        weight_[10] = 4.1425537809187001e-02;
        gp_.push_back({6.3650249912139900e-01, 3.1035245103378400e-01});
        
        weight_[11] = 4.1425537809187001e-02;
        gp_.push_back({5.3145049844816994e-02, 6.3650249912139900e-01});
        
    }
    
//---------------------------------------------------------------------------
    std::string T2_7_13::getName() const
    {
        return name_;
    }
    
    std::size_t T2_7_13::order() const
    {
        return 7;
    }
    
    std::size_t T2_7_13::dimension() const
    {
        return 2;
    }
    
    std::size_t T2_7_13::getNumberOfPoints() const
    {
        return 13;
    }
    
    double T2_7_13::weight(std::size_t i) const
    {
        logger.assert_debug(i < getNumberOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReferenceBase&
    T2_7_13::getPoint(std::size_t i) const
    {
        logger.assert_debug(i < getNumberOfPoints(), "%::getPoint - wrong index!", name_);
        return gp_[i];
    }
    
    Geometry::ReferenceGeometry*
    T2_7_13::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    T2_7_13::T2_7_13()
            : name_("T2_7_1"), refGeoPtr_(&Geometry::ReferenceTriangle::Instance()), gp_(0)
    {
        weight_[0] = -7.4785022233840995e-02;
        gp_.push_back({3.3333333333333348e-01, 3.3333333333333348e-01});
        
        weight_[1] = 8.7807628716603997e-02;
        gp_.push_back({2.6034596607904004e-01, 2.6034596607904004e-01});
        
        weight_[2] = 8.7807628716603997e-02;
        gp_.push_back({2.6034596607904004e-01, 4.7930806784191998e-01});
        
        weight_[3] = 8.7807628716603997e-02;
        gp_.push_back({4.7930806784191998e-01, 2.6034596607904004e-01});
        
        weight_[4] = 2.6673617804419000e-02;
        gp_.push_back({6.5130102902215992e-02, 6.5130102902215992e-02});
        
        weight_[5] = 2.6673617804419000e-02;
        gp_.push_back({6.5130102902215992e-02, 8.6973979419556802e-01});
        
        weight_[6] = 2.6673617804419000e-02;
        gp_.push_back({8.6973979419556802e-01, 6.5130102902215992e-02});
        
        weight_[7] = 3.8556880445128498e-02;
        gp_.push_back({3.1286549600487401e-01, 6.3844418856981000e-01});
        
        weight_[8] = 3.8556880445128498e-02;
        gp_.push_back({6.3844418856981000e-01, 4.8690315425315989e-02});
        
        weight_[9] = 3.8556880445128498e-02;
        gp_.push_back({4.8690315425315989e-02, 3.1286549600487401e-01});
        
        weight_[10] = 3.8556880445128498e-02;
        gp_.push_back({3.1286549600487401e-01, 4.8690315425315989e-02});
        
        weight_[11] = 3.8556880445128498e-02;
        gp_.push_back({6.3844418856981000e-01, 3.1286549600487401e-01});
        
        weight_[12] = 3.8556880445128498e-02;
        gp_.push_back({4.8690315425315989e-02, 6.3844418856981000e-01});
        
    }
    
//---------------------------------------------------------------------------
    std::string T2_8_16::getName() const
    {
        return name_;
    }
    
    std::size_t T2_8_16::order() const
    {
        return 8;
    }
    
    std::size_t T2_8_16::dimension() const
    {
        return 2;
    }
    
    std::size_t T2_8_16::getNumberOfPoints() const
    {
        return 16;
    }
    
    double T2_8_16::weight(std::size_t i) const
    {
        logger.assert_debug(i < getNumberOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReferenceBase&
    T2_8_16::getPoint(std::size_t i) const
    {
        logger.assert_debug(i < getNumberOfPoints(), "%::getPoint - wrong index!", name_);
        return gp_[i];
    }
    
    Geometry::ReferenceGeometry*
    T2_8_16::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    T2_8_16::T2_8_16()
            : name_("T2_8_1"), refGeoPtr_(&Geometry::ReferenceTriangle::Instance()), gp_(0)
    {
        weight_[0] = 7.2157803838893503e-02;
        gp_.push_back({3.3333333333333348e-01, 3.3333333333333348e-01});
        
        weight_[1] = 4.7545817133642497e-02;
        gp_.push_back({4.5929258829272301e-01, 4.5929258829272301e-01});
        
        weight_[2] = 4.7545817133642497e-02;
        gp_.push_back({4.5929258829272301e-01, 8.1414823414554027e-02});
        
        weight_[3] = 4.7545817133642497e-02;
        gp_.push_back({8.1414823414554027e-02, 4.5929258829272301e-01});
        
        weight_[4] = 5.1608685267358997e-02;
        gp_.push_back({1.7056930775175999e-01, 1.7056930775175999e-01});
        
        weight_[5] = 5.1608685267358997e-02;
        gp_.push_back({1.7056930775175999e-01, 6.5886138449648002e-01});
        
        weight_[6] = 5.1608685267358997e-02;
        gp_.push_back({6.5886138449648002e-01, 1.7056930775175999e-01});
        
        weight_[7] = 1.6229248811599001e-02;
        gp_.push_back({5.0547228317031012e-02, 5.0547228317031012e-02});
        
        weight_[8] = 1.6229248811598752e-02;
        gp_.push_back({5.0547228317031012e-02, 8.9890554336593798e-01});
        
        weight_[9] = 1.6229248811599001e-02;
        gp_.push_back({8.9890554336593798e-01, 5.0547228317031012e-02});
        
        weight_[10] = 1.3615157087217500e-02;
        gp_.push_back({2.6311282963463800e-01, 7.2849239295540402e-01});
        
        weight_[11] = 1.3615157087217500e-02;
        gp_.push_back({7.2849239295540402e-01, 8.3947774099579764e-03});
        
        weight_[12] = 1.3615157087217500e-02;
        gp_.push_back({8.3947774099579764e-03, 2.6311282963463800e-01});
        
        weight_[13] = 1.3615157087217500e-02;
        gp_.push_back({2.6311282963463800e-01, 8.3947774099579764e-03});
        
        weight_[14] = 1.3615157087217500e-02;
        gp_.push_back({7.2849239295540402e-01, 2.6311282963463800e-01});
        
        weight_[15] = 1.3615157087217500e-02;
        gp_.push_back({8.3947774099579764e-03, 7.2849239295540402e-01});
        
    }
    
//---------------------------------------------------------------------------
    std::string T2_9_19::getName() const
    {
        return name_;
    }
    
    std::size_t T2_9_19::order() const
    {
        return 9;
    }
    
    std::size_t T2_9_19::dimension() const
    {
        return 2;
    }
    
    std::size_t T2_9_19::getNumberOfPoints() const
    {
        return 19;
    }
    
    double T2_9_19::weight(std::size_t i) const
    {
        logger.assert_debug(i < getNumberOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReferenceBase&
    T2_9_19::getPoint(std::size_t i) const
    {
        logger.assert_debug(i < getNumberOfPoints(), "%::getPoint - wrong index!", name_);
        return gp_[i];
    }
    
    Geometry::ReferenceGeometry*
    T2_9_19::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    T2_9_19::T2_9_19()
            : name_("T2_9_1"), refGeoPtr_(&Geometry::ReferenceTriangle::Instance()), gp_(0)
    {
        weight_[0] = 4.8567898141399501e-02;
        gp_.push_back({3.3333333333333348e-01, 3.3333333333333348e-01});
        
        weight_[1] = 1.5667350113569501e-02;
        gp_.push_back({4.8968251919873801e-01, 4.8968251919873801e-01});
        
        weight_[2] = 1.5667350113569501e-02;
        gp_.push_back({4.8968251919873801e-01, 2.0634961602524982e-02});
        
        weight_[3] = 1.5667350113569501e-02;
        gp_.push_back({2.0634961602524982e-02, 4.8968251919873801e-01});
        
        weight_[4] = 3.8913770502387000e-02;
        gp_.push_back({4.3708959149293702e-01, 4.3708959149293702e-01});
        
        weight_[5] = 3.8913770502387000e-02;
        gp_.push_back({4.3708959149293702e-01, 1.2582081701412701e-01});
        
        weight_[6] = 3.8913770502387000e-02;
        gp_.push_back({1.2582081701412701e-01, 4.3708959149293702e-01});
        
        weight_[7] = 3.9823869463604999e-02;
        gp_.push_back({1.8820353561903302e-01, 1.8820353561903302e-01});
        
        weight_[8] = 3.9823869463604999e-02;
        gp_.push_back({1.8820353561903302e-01, 6.2359292876193506e-01});
        
        weight_[9] = 3.9823869463604999e-02;
        gp_.push_back({6.2359292876193506e-01, 1.8820353561903302e-01});
        
        weight_[10] = 1.2788837829349000e-02;
        gp_.push_back({4.4729513394453024e-02, 4.4729513394453024e-02});
        
        weight_[11] = 1.2788837829349000e-02;
        gp_.push_back({4.4729513394453024e-02, 9.1054097321109495e-01});
        
        weight_[12] = 1.2788837829349000e-02;
        gp_.push_back({9.1054097321109495e-01, 4.4729513394453024e-02});
        
        weight_[13] = 2.1641769688644501e-02;
        gp_.push_back({2.2196298916076601e-01, 7.4119859878449801e-01});
        
        weight_[14] = 2.1641769688644501e-02;
        gp_.push_back({7.4119859878449801e-01, 3.6838412054735981e-02});
        
        weight_[15] = 2.1641769688644501e-02;
        gp_.push_back({3.6838412054735981e-02, 2.2196298916076601e-01});
        
        weight_[16] = 2.1641769688644501e-02;
        gp_.push_back({2.2196298916076601e-01, 3.6838412054735981e-02});
        
        weight_[17] = 2.1641769688644501e-02;
        gp_.push_back({7.4119859878449801e-01, 2.2196298916076601e-01});
        
        weight_[18] = 2.1641769688644501e-02;
        gp_.push_back({3.6838412054735981e-02, 7.4119859878449801e-01});
        
    }
    
//---------------------------------------------------------------------------
    std::string T2_10_25::getName() const
    {
        return name_;
    }
    
    std::size_t T2_10_25::order() const
    {
        return 10;
    }
    
    std::size_t T2_10_25::dimension() const
    {
        return 2;
    }
    
    std::size_t T2_10_25::getNumberOfPoints() const
    {
        return 25;
    }
    
    double T2_10_25::weight(std::size_t i) const
    {
        logger.assert_debug(i < getNumberOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReferenceBase&
    T2_10_25::getPoint(std::size_t i) const
    {
        logger.assert_debug(i < getNumberOfPoints(), "%::getPoint - wrong index!", name_);
        return gp_[i];
    }
    
    Geometry::ReferenceGeometry*
    T2_10_25::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    T2_10_25::T2_10_25()
            : name_("T2_10_1"), refGeoPtr_(&Geometry::ReferenceTriangle::Instance()), gp_(0)
    {
        weight_[0] = 4.5408995191376998e-02;
        gp_.push_back({3.3333333333333348e-01, 3.3333333333333348e-01});
        
        weight_[1] = 1.8362978878233498e-02;
        gp_.push_back({4.8557763338365700e-01, 4.8557763338365700e-01});
        
        weight_[2] = 1.8362978878233498e-02;
        gp_.push_back({4.8557763338365700e-01, 2.8844733232684994e-02});
        
        weight_[3] = 1.8362978878233498e-02;
        gp_.push_back({2.8844733232684994e-02, 4.8557763338365700e-01});
        
        weight_[4] = 2.2660529717764000e-02;
        gp_.push_back({1.0948157548503701e-01, 1.0948157548503701e-01});
        
        weight_[5] = 2.2660529717764000e-02;
        gp_.push_back({1.0948157548503701e-01, 7.8103684902992598e-01});
        
        weight_[6] = 2.2660529717764000e-02;
        gp_.push_back({7.8103684902992598e-01, 1.0948157548503701e-01});
        
        weight_[7] = 3.6378958422710002e-02;
        gp_.push_back({3.0793983876412101e-01, 5.5035294182099903e-01});
        
        weight_[8] = 3.6378958422710002e-02;
        gp_.push_back({5.5035294182099903e-01, 1.4170721941488001e-01});
        
        weight_[9] = 3.6378958422710002e-02;
        gp_.push_back({1.4170721941488001e-01, 3.0793983876412101e-01});
        
        weight_[10] = 3.6378958422710002e-02;
        gp_.push_back({3.0793983876412101e-01, 1.4170721941488001e-01});
        
        weight_[11] = 3.6378958422710002e-02;
        gp_.push_back({5.5035294182099903e-01, 3.0793983876412101e-01});
        
        weight_[12] = 3.6378958422710002e-02;
        gp_.push_back({1.4170721941488001e-01, 5.5035294182099903e-01});
        
        weight_[13] = 1.4163621265528500e-02;
        gp_.push_back({2.4667256063990300e-01, 7.2832390459741103e-01});
        
        weight_[14] = 1.4163621265528500e-02;
        gp_.push_back({7.2832390459741103e-01, 2.5003534762686019e-02});
        
        weight_[15] = 1.4163621265528500e-02;
        gp_.push_back({2.5003534762686019e-02, 2.4667256063990300e-01});
        
        weight_[16] = 1.4163621265528500e-02;
        gp_.push_back({2.4667256063990300e-01, 2.5003534762686019e-02});
        
        weight_[17] = 1.4163621265528500e-02;
        gp_.push_back({7.2832390459741103e-01, 2.4667256063990300e-01});
        
        weight_[18] = 1.4163621265528500e-02;
        gp_.push_back({2.5003534762686019e-02, 7.2832390459741103e-01});
        
        weight_[19] = 4.7108334818665000e-03;
        gp_.push_back({6.6803251012200027e-02, 9.2365593358749998e-01});
        
        weight_[20] = 4.7108334818665000e-03;
        gp_.push_back({9.2365593358749998e-01, 9.5408154002989964e-03});
        
        weight_[21] = 4.7108334818665000e-03;
        gp_.push_back({9.5408154002989964e-03, 6.6803251012200027e-02});
        
        weight_[22] = 4.7108334818665000e-03;
        gp_.push_back({6.6803251012200027e-02, 9.5408154002989964e-03});
        
        weight_[23] = 4.7108334818665000e-03;
        gp_.push_back({9.2365593358749998e-01, 6.6803251012200027e-02});
        
        weight_[24] = 4.7108334818665000e-03;
        gp_.push_back({9.5408154002989964e-03, 9.2365593358749998e-01});
        
    }
    
//---------------------------------------------------------------------------
    std::string T2_11_28::getName() const
    {
        return name_;
    }
    
    std::size_t T2_11_28::order() const
    {
        return 11;
    }
    
    std::size_t T2_11_28::dimension() const
    {
        return 2;
    }
    
    std::size_t T2_11_28::getNumberOfPoints() const
    {
        return 28;
    }
    
    double T2_11_28::weight(std::size_t i) const
    {
        logger.assert_debug(i < getNumberOfPoints(), "%::weight - wrong index!", name_);
        return weight_[i];
    }
    
    const Geometry::PointReferenceBase&
    T2_11_28::getPoint(std::size_t i) const
    {
        logger.assert_debug(i < getNumberOfPoints(), "%::getPoint - wrong index!", name_);
        return gp_[i];
    }
    
    Geometry::ReferenceGeometry*
    T2_11_28::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    T2_11_28::T2_11_28()
            : //warning: points and weights found vary wildly between runs of the quadrature rule generator
            name_("T2_11_1"), refGeoPtr_(&Geometry::ReferenceTriangle::Instance()), gp_(0)
    {
        weight_[0] = 8.005121009414541471823946e-02 / 2.;
        gp_.push_back({3.3333333333333333e-01, 3.3333333333333333e-01});
        
        weight_[1] = 6.71206859217794936920366e-02 / 2.;
        gp_.push_back({2.154717529698554e-01, 2.154717529698554e-01});
        
        weight_[2] = 6.71206859217794936920366e-02 / 2.;
        gp_.push_back({2.154717529698554e-01, 1. - 2. * 2.154717529698554e-01});
        
        weight_[3] = 6.71206859217794936920366e-02 / 2.;
        gp_.push_back({1. - 2. * 2.154717529698554e-01, 2.154717529698554e-01});
        
        weight_[4] = 1.260585310278218701321e-02 / 2.;
        gp_.push_back({3.12795522110139299208957e-02, 3.12795522110139299208957e-02});
        
        weight_[5] = 1.260585310278218701321e-02 / 2.;
        gp_.push_back({3.12795522110139299208957e-02, 1. - 2. * 3.12795522110139299208957e-02});
        
        weight_[6] = 1.260585310278218701321e-02 / 2.;
        gp_.push_back({1. - 2. * 3.12795522110139299208957e-02, 3.12795522110139299208957e-02});
        
        weight_[7] = 4.05824609257477118939e-02 / 2.;
        gp_.push_back({1.155035357737921489271e-01, 1.155035357737921489271e-01});
        
        weight_[8] = 4.05824609257477118939e-02 / 2.;
        gp_.push_back({1.155035357737921489271e-01, 1. - 2. * 1.155035357737921489271e-01});
        
        weight_[9] = 4.05824609257477118939e-02 / 2.;
        gp_.push_back({1. - 2. * 1.155035357737921489271e-01, 1.155035357737921489271e-01});
        
        weight_[10] = 6.17123685629664404606197e-02 / 2.;
        gp_.push_back({4.357774152823301051418e-01, 4.357774152823301051418e-01});
        
        weight_[11] = 6.17123685629664404606197e-02 / 2.;
        gp_.push_back({4.357774152823301051418e-01, 1. - 2. * 4.357774152823301051418e-01});
        
        weight_[12] = 6.17123685629664404606197e-02 / 2.;
        gp_.push_back({1. - 2. * 4.357774152823301051418e-01, 4.357774152823301051418e-01});
        
        weight_[13] = 1.12084059305228051993e-02 / 2.;
        gp_.push_back({4.9988403867731539816248e-01, 4.9988403867731539816248e-01});
        
        weight_[14] = 1.12084059305228051993e-02 / 2.;
        gp_.push_back({4.9988403867731539816248e-01, 1. - 2. * 4.9988403867731539816248e-01});
        
        weight_[15] = 1.12084059305228051993e-02 / 2.;
        gp_.push_back({1. - 2. * 4.9988403867731539816248e-01, 4.9988403867731539816248e-01});
        
        weight_[16] = 4.115295964681454289321412e-02 / 2.;
        gp_.push_back({3.164262022308527e-01, 4.8055796920866524676254359e-02});
        
        weight_[17] = 4.115295964681454289321412e-02 / 2.;
        gp_.push_back({4.8055796920866524676254359e-02, 3.164262022308527e-01});
        
        weight_[18] = 4.115295964681454289321412e-02 / 2.;
        gp_.push_back({3.164262022308527e-01, 1. - 4.8055796920866524676254359e-02 - 3.164262022308527e-01});
        
        weight_[19] = 4.115295964681454289321412e-02 / 2.;
        gp_.push_back({1. - 4.8055796920866524676254359e-02 - 3.164262022308527e-01, 3.164262022308527e-01});
        
        weight_[20] = 4.115295964681454289321412e-02 / 2.;
        gp_.push_back({4.8055796920866524676254359e-02, 1. - 4.8055796920866524676254359e-02 - 3.164262022308527e-01});
        
        weight_[21] = 4.115295964681454289321412e-02 / 2.;
        gp_.push_back({1. - 4.8055796920866524676254359e-02 - 3.164262022308527e-01, 4.8055796920866524676254359e-02});
        
        weight_[22] = 1.55569514489285688575e-02 / 2.;
        gp_.push_back({1.614188489108937e-01, 1.56544537223683037837e-02});
        
        weight_[23] = 1.55569514489285688575e-02 / 2.;
        gp_.push_back({1.56544537223683037837e-02, 1.614188489108937e-01});
        
        weight_[24] = 1.55569514489285688575e-02 / 2.;
        gp_.push_back({1.614188489108937e-01, 1. - 1.614188489108937e-01 - 1.56544537223683037837e-02});
        
        weight_[25] = 1.55569514489285688575e-02 / 2.;
        gp_.push_back({1. - 1.614188489108937e-01 - 1.56544537223683037837e-02, 1.614188489108937e-01});
        
        weight_[26] = 1.55569514489285688575e-02 / 2.;
        gp_.push_back({1.56544537223683037837e-02, 1. - 1.614188489108937e-01 - 1.56544537223683037837e-02});
        
        weight_[27] = 1.55569514489285688575e-02 / 2.;
        gp_.push_back({1. - 1.614188489108937e-01 - 1.56544537223683037837e-02, 1.56544537223683037837e-02});
        
    }

//---------------------------------------------------------------------------
}// close namespace IntegrationRules

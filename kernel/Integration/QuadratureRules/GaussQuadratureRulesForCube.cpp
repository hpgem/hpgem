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
// File: GaussQuadratureRulesForCube.cpp 
// Implementation of Gauss quadrature rules for reference cube.
// Lars Pesch, Fri Mar  3 12:59:11 CET 2006
//----
// Modified from original file allGaussQuadratureRules.h
// by M.T. Julianto, Wed Feb 23 10:45:06 UTC 2013
//---------------------------------------------------------------------------
// System includes and names imported from them:
#include <cmath>
//---------------------------------------------------------------------------
// Package includes:
#include "Integration/QuadratureRules/GaussQuadratureRulesForCube.h"
#include "Geometry/ReferenceCube.h"
#include "GaussQuadratureRulesForLine.h"
#include "Geometry/PointReference.h"
using Geometry::ReferenceCube;

//---------------------------------------------------------------------------
namespace QuadratureRules
{
//---------------------------------------------------------------------------
    std::string Cn3_1_1::getName() const
    {
        return name_;
    }
    
    std::size_t Cn3_1_1::order() const
    {
        return 1;
    }
    
    std::size_t Cn3_1_1::dimension() const
    {
        return 3;
    }
    
    std::size_t Cn3_1_1::nrOfPoints() const
    {
        return 1;
    }
    
    double Cn3_1_1::weight(std::size_t i) const
    {
        if (i < 1)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }
    
    const Geometry::PointReference&
    Cn3_1_1::getPoint(std::size_t i) const
    {
        if (i < 1)
            return gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }
    
    Cn3_1_1::ReferenceGeometryT*
    Cn3_1_1::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    Cn3_1_1::Cn3_1_1()
            : name_("Cn3_1_1"), refGeoPtr_(&ReferenceCube::Instance()), gp_(1, 3)
    {
        weight_[0] = ((2.0) * (2.0)) * (2.0);
        gp_[0][0] = 0.0;
        gp_[0][1] = 0.0;
        gp_[0][2] = 0.0;
        
    }
    
    Cn3_1_1::~Cn3_1_1()
    {
    }
    
//---------------------------------------------------------------------------
    std::string Cn3_3_4::getName() const
    {
        return name_;
    }
    
    std::size_t Cn3_3_4::order() const
    {
        return 3;
    }
    
    std::size_t Cn3_3_4::dimension() const
    {
        return 3;
    }
    
    std::size_t Cn3_3_4::nrOfPoints() const
    {
        return 8;
    }
    
    double Cn3_3_4::weight(std::size_t i) const
    {
        if (i < 8)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }
    
    const Geometry::PointReference&
    Cn3_3_4::getPoint(std::size_t i) const
    {
        if (i < 8)
            return gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }
    
    Cn3_3_4::ReferenceGeometryT*
    Cn3_3_4::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    Cn3_3_4::Cn3_3_4()
            : name_("Cn3_3_4"), refGeoPtr_(&ReferenceCube::Instance()), gp_(8, 3)
    {
        weight_[0] = ((1.0) * (1.0)) * (1.0);
        gp_[0][0] = -std::sqrt(3.0) / 3.0;
        gp_[0][1] = -std::sqrt(3.0) / 3.0;
        gp_[0][2] = -std::sqrt(3.0) / 3.0;
        
        weight_[1] = ((1.0) * (1.0)) * (1.0);
        gp_[1][0] = +std::sqrt(3.0) / 3.0;
        gp_[1][1] = -std::sqrt(3.0) / 3.0;
        gp_[1][2] = -std::sqrt(3.0) / 3.0;
        
        weight_[2] = ((1.0) * (1.0)) * (1.0);
        gp_[2][0] = -std::sqrt(3.0) / 3.0;
        gp_[2][1] = +std::sqrt(3.0) / 3.0;
        gp_[2][2] = -std::sqrt(3.0) / 3.0;
        
        weight_[3] = ((1.0) * (1.0)) * (1.0);
        gp_[3][0] = +std::sqrt(3.0) / 3.0;
        gp_[3][1] = +std::sqrt(3.0) / 3.0;
        gp_[3][2] = -std::sqrt(3.0) / 3.0;
        
        weight_[4] = ((1.0) * (1.0)) * (1.0);
        gp_[4][0] = -std::sqrt(3.0) / 3.0;
        gp_[4][1] = -std::sqrt(3.0) / 3.0;
        gp_[4][2] = +std::sqrt(3.0) / 3.0;
        
        weight_[5] = ((1.0) * (1.0)) * (1.0);
        gp_[5][0] = +std::sqrt(3.0) / 3.0;
        gp_[5][1] = -std::sqrt(3.0) / 3.0;
        gp_[5][2] = +std::sqrt(3.0) / 3.0;
        
        weight_[6] = ((1.0) * (1.0)) * (1.0);
        gp_[6][0] = -std::sqrt(3.0) / 3.0;
        gp_[6][1] = +std::sqrt(3.0) / 3.0;
        gp_[6][2] = +std::sqrt(3.0) / 3.0;
        
        weight_[7] = ((1.0) * (1.0)) * (1.0);
        gp_[7][0] = +std::sqrt(3.0) / 3.0;
        gp_[7][1] = +std::sqrt(3.0) / 3.0;
        gp_[7][2] = +std::sqrt(3.0) / 3.0;
        
    }
    
    Cn3_3_4::~Cn3_3_4()
    {
        
    }
    
//---------------------------------------------------------------------------
    std::string Cn3_5_9::getName() const
    {
        return name_;
    }
    
    std::size_t Cn3_5_9::order() const
    {
        return 5;
    }
    
    std::size_t Cn3_5_9::dimension() const
    {
        return 3;
    }
    
    std::size_t Cn3_5_9::nrOfPoints() const
    {
        return 27;
    }
    
    double Cn3_5_9::weight(std::size_t i) const
    {
        if (i < 27)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }
    
    const Geometry::PointReference&
    Cn3_5_9::getPoint(std::size_t i) const
    {
        if (i < 27)
            return gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }
    
    Cn3_5_9::ReferenceGeometryT*
    Cn3_5_9::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    Cn3_5_9::Cn3_5_9()
            : name_("Cn3_5_9"), refGeoPtr_(&ReferenceCube::Instance()), gp_(27, 3)
    {
        weight_[0] = ((5. / 9.) * (5. / 9.)) * (5. / 9.);
        gp_[0][0] = -std::sqrt(3.0 / 5.0);
        gp_[0][1] = -std::sqrt(3.0 / 5.0);
        gp_[0][2] = -std::sqrt(3.0 / 5.0);
        
        weight_[1] = ((8. / 9.) * (5. / 9.)) * (5. / 9.);
        gp_[1][0] = 0.0;
        gp_[1][1] = -std::sqrt(3.0 / 5.0);
        gp_[1][2] = -std::sqrt(3.0 / 5.0);
        
        weight_[2] = ((5. / 9.) * (5. / 9.)) * (5. / 9.);
        gp_[2][0] = +std::sqrt(3.0 / 5.0);
        gp_[2][1] = -std::sqrt(3.0 / 5.0);
        gp_[2][2] = -std::sqrt(3.0 / 5.0);
        
        weight_[3] = ((5. / 9.) * (8. / 9.)) * (5. / 9.);
        gp_[3][0] = -std::sqrt(3.0 / 5.0);
        gp_[3][1] = 0.0;
        gp_[3][2] = -std::sqrt(3.0 / 5.0);
        
        weight_[4] = ((8. / 9.) * (8. / 9.)) * (5. / 9.);
        gp_[4][0] = 0.0;
        gp_[4][1] = 0.0;
        gp_[4][2] = -std::sqrt(3.0 / 5.0);
        
        weight_[5] = ((5. / 9.) * (8. / 9.)) * (5. / 9.);
        gp_[5][0] = +std::sqrt(3.0 / 5.0);
        gp_[5][1] = 0.0;
        gp_[5][2] = -std::sqrt(3.0 / 5.0);
        
        weight_[6] = ((5. / 9.) * (5. / 9.)) * (5. / 9.);
        gp_[6][0] = -std::sqrt(3.0 / 5.0);
        gp_[6][1] = +std::sqrt(3.0 / 5.0);
        gp_[6][2] = -std::sqrt(3.0 / 5.0);
        
        weight_[7] = ((8. / 9.) * (5. / 9.)) * (5. / 9.);
        gp_[7][0] = 0.0;
        gp_[7][1] = +std::sqrt(3.0 / 5.0);
        gp_[7][2] = -std::sqrt(3.0 / 5.0);
        
        weight_[8] = ((5. / 9.) * (5. / 9.)) * (5. / 9.);
        gp_[8][0] = +std::sqrt(3.0 / 5.0);
        gp_[8][1] = +std::sqrt(3.0 / 5.0);
        gp_[8][2] = -std::sqrt(3.0 / 5.0);
        
        weight_[9] = ((5. / 9.) * (5. / 9.)) * (8. / 9.);
        gp_[9][0] = -std::sqrt(3.0 / 5.0);
        gp_[9][1] = -std::sqrt(3.0 / 5.0);
        gp_[9][2] = 0.0;
        
        weight_[10] = ((8. / 9.) * (5. / 9.)) * (8. / 9.);
        gp_[10][0] = 0.0;
        gp_[10][1] = -std::sqrt(3.0 / 5.0);
        gp_[10][2] = 0.0;
        
        weight_[11] = ((5. / 9.) * (5. / 9.)) * (8. / 9.);
        gp_[11][0] = +std::sqrt(3.0 / 5.0);
        gp_[11][1] = -std::sqrt(3.0 / 5.0);
        gp_[11][2] = 0.0;
        
        weight_[12] = ((5. / 9.) * (8. / 9.)) * (8. / 9.);
        gp_[12][0] = -std::sqrt(3.0 / 5.0);
        gp_[12][1] = 0.0;
        gp_[12][2] = 0.0;
        
        weight_[13] = ((8. / 9.) * (8. / 9.)) * (8. / 9.);
        gp_[13][0] = 0.0;
        gp_[13][1] = 0.0;
        gp_[13][2] = 0.0;
        
        weight_[14] = ((5. / 9.) * (8. / 9.)) * (8. / 9.);
        gp_[14][0] = +std::sqrt(3.0 / 5.0);
        gp_[14][1] = 0.0;
        gp_[14][2] = 0.0;
        
        weight_[15] = ((5. / 9.) * (5. / 9.)) * (8. / 9.);
        gp_[15][0] = -std::sqrt(3.0 / 5.0);
        gp_[15][1] = +std::sqrt(3.0 / 5.0);
        gp_[15][2] = 0.0;
        
        weight_[16] = ((8. / 9.) * (5. / 9.)) * (8. / 9.);
        gp_[16][0] = 0.0;
        gp_[16][1] = +std::sqrt(3.0 / 5.0);
        gp_[16][2] = 0.0;
        
        weight_[17] = ((5. / 9.) * (5. / 9.)) * (8. / 9.);
        gp_[17][0] = +std::sqrt(3.0 / 5.0);
        gp_[17][1] = +std::sqrt(3.0 / 5.0);
        gp_[17][2] = 0.0;
        
        weight_[18] = ((5. / 9.) * (5. / 9.)) * (5. / 9.);
        gp_[18][0] = -std::sqrt(3.0 / 5.0);
        gp_[18][1] = -std::sqrt(3.0 / 5.0);
        gp_[18][2] = +std::sqrt(3.0 / 5.0);
        
        weight_[19] = ((8. / 9.) * (5. / 9.)) * (5. / 9.);
        gp_[19][0] = 0.0;
        gp_[19][1] = -std::sqrt(3.0 / 5.0);
        gp_[19][2] = +std::sqrt(3.0 / 5.0);
        
        weight_[20] = ((5. / 9.) * (5. / 9.)) * (5. / 9.);
        gp_[20][0] = +std::sqrt(3.0 / 5.0);
        gp_[20][1] = -std::sqrt(3.0 / 5.0);
        gp_[20][2] = +std::sqrt(3.0 / 5.0);
        
        weight_[21] = ((5. / 9.) * (8. / 9.)) * (5. / 9.);
        gp_[21][0] = -std::sqrt(3.0 / 5.0);
        gp_[21][1] = 0.0;
        gp_[21][2] = +std::sqrt(3.0 / 5.0);
        
        weight_[22] = ((8. / 9.) * (8. / 9.)) * (5. / 9.);
        gp_[22][0] = 0.0;
        gp_[22][1] = 0.0;
        gp_[22][2] = +std::sqrt(3.0 / 5.0);
        
        weight_[23] = ((5. / 9.) * (8. / 9.)) * (5. / 9.);
        gp_[23][0] = +std::sqrt(3.0 / 5.0);
        gp_[23][1] = 0.0;
        gp_[23][2] = +std::sqrt(3.0 / 5.0);
        
        weight_[24] = ((5. / 9.) * (5. / 9.)) * (5. / 9.);
        gp_[24][0] = -std::sqrt(3.0 / 5.0);
        gp_[24][1] = +std::sqrt(3.0 / 5.0);
        gp_[24][2] = +std::sqrt(3.0 / 5.0);
        
        weight_[25] = ((8. / 9.) * (5. / 9.)) * (5. / 9.);
        gp_[25][0] = 0.0;
        gp_[25][1] = +std::sqrt(3.0 / 5.0);
        gp_[25][2] = +std::sqrt(3.0 / 5.0);
        
        weight_[26] = ((5. / 9.) * (5. / 9.)) * (5. / 9.);
        gp_[26][0] = +std::sqrt(3.0 / 5.0);
        gp_[26][1] = +std::sqrt(3.0 / 5.0);
        gp_[26][2] = +std::sqrt(3.0 / 5.0);
        
    }
    
    Cn3_5_9::~Cn3_5_9()
    {
    }
    
//---------------------------------------------------------------------------
    std::string C3_7_2::getName() const
    {
        return name_;
    }
    
    std::size_t C3_7_2::order() const
    {
        return 7;
    }
    
    std::size_t C3_7_2::dimension() const
    {
        return 3;
    }
    
    std::size_t C3_7_2::nrOfPoints() const
    {
        return 64;
    }
    
    double C3_7_2::weight(std::size_t i) const
    {
        if (i < 64)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }
    
    const Geometry::PointReference&
    C3_7_2::getPoint(std::size_t i) const
    {
        if (i < 64)
            return gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }
    
    C3_7_2::ReferenceGeometryT*
    C3_7_2::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    C3_7_2::C3_7_2()
            : name_("C3_7_2"), refGeoPtr_(&ReferenceCube::Instance()), gp_(64, 3)
    {
        int position(0);
        C1_7_x& ruleForLine = C1_7_x::Instance();
        Geometry::PointReference point1D(1);
        for (int i = 0; i < ruleForLine.nrOfPoints(); ++i)
        {
            for (int j = 0; j < ruleForLine.nrOfPoints(); ++j)
            {
                for (int k = 0; k < ruleForLine.nrOfPoints(); ++k)
                {
                    weight_[position] = ruleForLine.weight(i) * ruleForLine.weight(j) * ruleForLine.weight(k);
                    point1D = ruleForLine.getPoint(i);
                    gp_[position][0] = point1D[0];
                    point1D = ruleForLine.getPoint(j);
                    gp_[position][1] = point1D[0];
                    point1D = ruleForLine.getPoint(k);
                    gp_[position][2] = point1D[0];
                    ++position;
                }
            }
        }
        /*weight_[0] = 1078. / 3645.;//rule breaks for (1-z)(1-x^2)(1-y^2)(5x^2-1)(5y^2-1)
         gp_[0][0] = +std::sqrt((6. / 7.));//and some others
         gp_[0][1] = 0.;
         gp_[0][2] = 0.;

         weight_[1] = 1078. / 3645.;
         gp_[1][0] = -std::sqrt((6. / 7.));
         gp_[1][1] = 0.;
         gp_[1][2] = 0.;

         weight_[2] = 1078. / 3645.;
         gp_[2][0] = 0.;
         gp_[2][1] = +std::sqrt((6. / 7.));
         gp_[2][2] = 0.;

         weight_[3] = 1078. / 3645.;
         gp_[3][0] = 0.;
         gp_[3][1] = -std::sqrt((6. / 7.));
         gp_[3][2] = 0.;

         weight_[4] = 1078. / 3645.;
         gp_[4][0] = 0.;
         gp_[4][1] = 0.;
         gp_[4][2] = +std::sqrt((6. / 7.));

         weight_[5] = 1078. / 3645.;
         gp_[5][0] = 0.;
         gp_[5][1] = 0.;
         gp_[5][2] = -std::sqrt((6. / 7.));

         weight_[6] = 343. / 3645.;
         gp_[6][0] = +std::sqrt((6. / 7.));
         gp_[6][1] = +std::sqrt((6. / 7.));
         gp_[6][2] = 0.;

         weight_[7] = 343. / 3645.;
         gp_[7][0] = -std::sqrt((6. / 7.));
         gp_[7][1] = +std::sqrt((6. / 7.));
         gp_[7][2] = 0.;

         weight_[8] = 343. / 3645.;
         gp_[8][0] = +std::sqrt((6. / 7.));
         gp_[8][1] = -std::sqrt((6. / 7.));
         gp_[8][2] = 0.;

         weight_[9] = 343. / 3645.;
         gp_[9][0] = -std::sqrt((6. / 7.));
         gp_[9][1] = -std::sqrt((6. / 7.));
         gp_[9][2] = 0.;

         weight_[10] = 343. / 3645.;
         gp_[10][0] = +std::sqrt((6. / 7.));
         gp_[10][1] = 0.;
         gp_[10][2] = +std::sqrt((6. / 7.));

         weight_[11] = 343. / 3645.;
         gp_[11][0] = -std::sqrt((6. / 7.));
         gp_[11][1] = 0.;
         gp_[11][2] = +std::sqrt((6. / 7.));

         weight_[12] = 343. / 3645.;
         gp_[12][0] = +std::sqrt((6. / 7.));
         gp_[12][1] = 0.;
         gp_[12][2] = -std::sqrt((6. / 7.));

         weight_[13] = 343. / 3645.;
         gp_[13][0] = -std::sqrt((6. / 7.));
         gp_[13][1] = 0.;
         gp_[13][2] = -std::sqrt((6. / 7.));

         weight_[14] = 343. / 3645.;
         gp_[14][0] = 0.;
         gp_[14][1] = +std::sqrt((6. / 7.));
         gp_[14][2] = +std::sqrt((6. / 7.));

         weight_[15] = 343. / 3645.;
         gp_[15][0] = 0.;
         gp_[15][1] = -std::sqrt((6. / 7.));
         gp_[15][2] = +std::sqrt((6. / 7.));

         weight_[16] = 343. / 3645.;
         gp_[16][0] = 0.;
         gp_[16][1] = +std::sqrt((6. / 7.));
         gp_[16][2] = -std::sqrt((6. / 7.));

         weight_[17] = 343. / 3645.;
         gp_[17][0] = 0.;
         gp_[17][1] = -std::sqrt((6. / 7.));
         gp_[17][2] = -std::sqrt((6. / 7.));

         weight_[18] = (774. * ((960. + 3. * std::sqrt(28798.)) / 2726.) - 230.) / (1215. * (((960. + 3. * std::sqrt(28798.)) / 2726.) - ((960. - 3. * std::sqrt(28798.)) / 2726.)));
         gp_[18][0] = +std::sqrt(((960. - 3. * std::sqrt(28798.)) / 2726.));
         gp_[18][1] = +std::sqrt(((960. - 3. * std::sqrt(28798.)) / 2726.));
         gp_[18][2] = +std::sqrt(((960. - 3. * std::sqrt(28798.)) / 2726.));

         weight_[19] = (774. * ((960. + 3. * std::sqrt(28798.)) / 2726.) - 230.) / (1215. * (((960. + 3. * std::sqrt(28798.)) / 2726.) - ((960. - 3. * std::sqrt(28798.)) / 2726.)));
         gp_[19][0] = -std::sqrt(((960. - 3. * std::sqrt(28798.)) / 2726.));
         gp_[19][1] = +std::sqrt(((960. - 3. * std::sqrt(28798.)) / 2726.));
         gp_[19][2] = +std::sqrt(((960. - 3. * std::sqrt(28798.)) / 2726.));

         weight_[20] = (774. * ((960. + 3. * std::sqrt(28798.)) / 2726.) - 230.) / (1215. * (((960. + 3. * std::sqrt(28798.)) / 2726.) - ((960. - 3. * std::sqrt(28798.)) / 2726.)));
         gp_[20][0] = +std::sqrt(((960. - 3. * std::sqrt(28798.)) / 2726.));
         gp_[20][1] = -std::sqrt(((960. - 3. * std::sqrt(28798.)) / 2726.));
         gp_[20][2] = +std::sqrt(((960. - 3. * std::sqrt(28798.)) / 2726.));

         weight_[21] = (774. * ((960. + 3. * std::sqrt(28798.)) / 2726.) - 230.) / (1215. * (((960. + 3. * std::sqrt(28798.)) / 2726.) - ((960. - 3. * std::sqrt(28798.)) / 2726.)));
         gp_[21][0] = -std::sqrt(((960. - 3. * std::sqrt(28798.)) / 2726.));
         gp_[21][1] = -std::sqrt(((960. - 3. * std::sqrt(28798.)) / 2726.));
         gp_[21][2] = +std::sqrt(((960. - 3. * std::sqrt(28798.)) / 2726.));

         weight_[22] = (774. * ((960. + 3. * std::sqrt(28798.)) / 2726.) - 230.) / (1215. * (((960. + 3. * std::sqrt(28798.)) / 2726.) - ((960. - 3. * std::sqrt(28798.)) / 2726.)));
         gp_[22][0] = +std::sqrt(((960. - 3. * std::sqrt(28798.)) / 2726.));
         gp_[22][1] = +std::sqrt(((960. - 3. * std::sqrt(28798.)) / 2726.));
         gp_[22][2] = -std::sqrt(((960. - 3. * std::sqrt(28798.)) / 2726.));

         weight_[23] = (774. * ((960. + 3. * std::sqrt(28798.)) / 2726.) - 230.) / (1215. * (((960. + 3. * std::sqrt(28798.)) / 2726.) - ((960. - 3. * std::sqrt(28798.)) / 2726.)));
         gp_[23][0] = -std::sqrt(((960. - 3. * std::sqrt(28798.)) / 2726.));
         gp_[23][1] = +std::sqrt(((960. - 3. * std::sqrt(28798.)) / 2726.));
         gp_[23][2] = -std::sqrt(((960. - 3. * std::sqrt(28798.)) / 2726.));

         weight_[24] = (774. * ((960. + 3. * std::sqrt(28798.)) / 2726.) - 230.) / (1215. * (((960. + 3. * std::sqrt(28798.)) / 2726.) - ((960. - 3. * std::sqrt(28798.)) / 2726.)));
         gp_[24][0] = +std::sqrt(((960. - 3. * std::sqrt(28798.)) / 2726.));
         gp_[24][1] = -std::sqrt(((960. - 3. * std::sqrt(28798.)) / 2726.));
         gp_[24][2] = -std::sqrt(((960. - 3. * std::sqrt(28798.)) / 2726.));

         weight_[25] = (774. * ((960. + 3. * std::sqrt(28798.)) / 2726.) - 230.) / (1215. * (((960. + 3. * std::sqrt(28798.)) / 2726.) - ((960. - 3. * std::sqrt(28798.)) / 2726.)));
         gp_[25][0] = -std::sqrt(((960. - 3. * std::sqrt(28798.)) / 2726.));
         gp_[25][1] = -std::sqrt(((960. - 3. * std::sqrt(28798.)) / 2726.));
         gp_[25][2] = -std::sqrt(((960. - 3. * std::sqrt(28798.)) / 2726.));

         weight_[26] = (230. - 774. * ((960. - 3. * std::sqrt(28798.)) / 2726.)) / (1215. * (((960. + 3. * std::sqrt(28798.)) / 2726.) - ((960. - 3. * std::sqrt(28798.)) / 2726.)));
         gp_[26][0] = +std::sqrt(((960. + 3. * std::sqrt(28798.)) / 2726.));
         gp_[26][1] = +std::sqrt(((960. + 3. * std::sqrt(28798.)) / 2726.));
         gp_[26][2] = +std::sqrt(((960. + 3. * std::sqrt(28798.)) / 2726.));

         weight_[27] = (230. - 774. * ((960. - 3. * std::sqrt(28798.)) / 2726.)) / (1215. * (((960. + 3. * std::sqrt(28798.)) / 2726.) - ((960. - 3. * std::sqrt(28798.)) / 2726.)));
         gp_[27][0] = -std::sqrt(((960. + 3. * std::sqrt(28798.)) / 2726.));
         gp_[27][1] = +std::sqrt(((960. + 3. * std::sqrt(28798.)) / 2726.));
         gp_[27][2] = +std::sqrt(((960. + 3. * std::sqrt(28798.)) / 2726.));

         weight_[28] = (230. - 774. * ((960. - 3. * std::sqrt(28798.)) / 2726.)) / (1215. * (((960. + 3. * std::sqrt(28798.)) / 2726.) - ((960. - 3. * std::sqrt(28798.)) / 2726.)));
         gp_[28][0] = +std::sqrt(((960. + 3. * std::sqrt(28798.)) / 2726.));
         gp_[28][1] = -std::sqrt(((960. + 3. * std::sqrt(28798.)) / 2726.));
         gp_[28][2] = +std::sqrt(((960. + 3. * std::sqrt(28798.)) / 2726.));

         weight_[29] = (230. - 774. * ((960. - 3. * std::sqrt(28798.)) / 2726.)) / (1215. * (((960. + 3. * std::sqrt(28798.)) / 2726.) - ((960. - 3. * std::sqrt(28798.)) / 2726.)));
         gp_[29][0] = -std::sqrt(((960. + 3. * std::sqrt(28798.)) / 2726.));
         gp_[29][1] = -std::sqrt(((960. + 3. * std::sqrt(28798.)) / 2726.));
         gp_[29][2] = +std::sqrt(((960. + 3. * std::sqrt(28798.)) / 2726.));

         weight_[30] = (230. - 774. * ((960. - 3. * std::sqrt(28798.)) / 2726.)) / (1215. * (((960. + 3. * std::sqrt(28798.)) / 2726.) - ((960. - 3. * std::sqrt(28798.)) / 2726.)));
         gp_[30][0] = +std::sqrt(((960. + 3. * std::sqrt(28798.)) / 2726.));
         gp_[30][1] = +std::sqrt(((960. + 3. * std::sqrt(28798.)) / 2726.));
         gp_[30][2] = -std::sqrt(((960. + 3. * std::sqrt(28798.)) / 2726.));

         weight_[31] = (230. - 774. * ((960. - 3. * std::sqrt(28798.)) / 2726.)) / (1215. * (((960. + 3. * std::sqrt(28798.)) / 2726.) - ((960. - 3. * std::sqrt(28798.)) / 2726.)));
         gp_[31][0] = -std::sqrt(((960. + 3. * std::sqrt(28798.)) / 2726.));
         gp_[31][1] = +std::sqrt(((960. + 3. * std::sqrt(28798.)) / 2726.));
         gp_[31][2] = -std::sqrt(((960. + 3. * std::sqrt(28798.)) / 2726.));

         weight_[32] = (230. - 774. * ((960. - 3. * std::sqrt(28798.)) / 2726.)) / (1215. * (((960. + 3. * std::sqrt(28798.)) / 2726.) - ((960. - 3. * std::sqrt(28798.)) / 2726.)));
         gp_[32][0] = +std::sqrt(((960. + 3. * std::sqrt(28798.)) / 2726.));
         gp_[32][1] = -std::sqrt(((960. + 3. * std::sqrt(28798.)) / 2726.));
         gp_[32][2] = -std::sqrt(((960. + 3. * std::sqrt(28798.)) / 2726.));

         weight_[33] = (230. - 774. * ((960. - 3. * std::sqrt(28798.)) / 2726.)) / (1215. * (((960. + 3. * std::sqrt(28798.)) / 2726.) - ((960. - 3. * std::sqrt(28798.)) / 2726.)));
         gp_[33][0] = -std::sqrt(((960. + 3. * std::sqrt(28798.)) / 2726.));
         gp_[33][1] = -std::sqrt(((960. + 3. * std::sqrt(28798.)) / 2726.));
         gp_[33][2] = -std::sqrt(((960. + 3. * std::sqrt(28798.)) / 2726.));
         */
    }
    
    C3_7_2::~C3_7_2()
    {
    }
    
//---------------------------------------------------------------------------
    std::string C3_9_2::getName() const
    {
        return name_;
    }
    
    std::size_t C3_9_2::order() const
    {
        return 9;
    }
    
    std::size_t C3_9_2::dimension() const
    {
        return 3;
    }
    
    std::size_t C3_9_2::nrOfPoints() const
    {
        return 125;
    }
    
    double C3_9_2::weight(std::size_t i) const
    {
        if (i < 125)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }
    
    const Geometry::PointReference&
    C3_9_2::getPoint(std::size_t i) const
    {
        if (i < 125)
            return gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }
    
    C3_9_2::ReferenceGeometryT*
    C3_9_2::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    C3_9_2::C3_9_2()
            : name_("C3_9_2"), refGeoPtr_(&ReferenceCube::Instance()), gp_(125, 3)
    {
        int position(0);
        C1_9_25& ruleForLine = C1_9_25::Instance();
        Geometry::PointReference point1D(1);
        for (int i = 0; i < ruleForLine.nrOfPoints(); ++i)
        {
            for (int j = 0; j < ruleForLine.nrOfPoints(); ++j)
            {
                for (int k = 0; k < ruleForLine.nrOfPoints(); ++k)
                {
                    weight_[position] = ruleForLine.weight(i) * ruleForLine.weight(j) * ruleForLine.weight(k);
                    point1D = ruleForLine.getPoint(i);
                    gp_[position][0] = point1D[0];
                    point1D = ruleForLine.getPoint(j);
                    gp_[position][1] = point1D[0];
                    point1D = ruleForLine.getPoint(k);
                    gp_[position][2] = point1D[0];
                    ++position;
                }
            }
        }
    }
    
    C3_9_2::~C3_9_2()
    {
    }
    
//---------------------------------------------------------------------------
    std::string C3_11_2::getName() const
    {
        return name_;
    }
    
    std::size_t C3_11_2::order() const
    {
        return 11;
    }
    
    std::size_t C3_11_2::dimension() const
    {
        return 3;
    }
    
    std::size_t C3_11_2::nrOfPoints() const
    {
        return 216;
    }
    
    double C3_11_2::weight(std::size_t i) const
    {
        if (i < 216)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }
    
    const Geometry::PointReference&
    C3_11_2::getPoint(std::size_t i) const
    {
        if (i < 216)
            return gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }
    
    C3_11_2::ReferenceGeometryT*
    C3_11_2::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }
    
    C3_11_2::C3_11_2()
            : name_("C3_11_2"), refGeoPtr_(&ReferenceCube::Instance()), gp_(216, 3)
    {
        int position(0);
        C1_11_36& ruleForLine = C1_11_36::Instance();
        Geometry::PointReference point1D(1);
        for (int i = 0; i < ruleForLine.nrOfPoints(); ++i)
        {
            for (int j = 0; j < ruleForLine.nrOfPoints(); ++j)
            {
                for (int k = 0; k < ruleForLine.nrOfPoints(); ++k)
                {
                    weight_[position] = ruleForLine.weight(i) * ruleForLine.weight(j) * ruleForLine.weight(k);
                    point1D = ruleForLine.getPoint(i);
                    gp_[position][0] = point1D[0];
                    point1D = ruleForLine.getPoint(j);
                    gp_[position][1] = point1D[0];
                    point1D = ruleForLine.getPoint(k);
                    gp_[position][2] = point1D[0];
                    ++position;
                }
            }
        }
    }
    
    C3_11_2::~C3_11_2()
    {
    }

//---------------------------------------------------------------------------
}// close namespace QuadratureRules

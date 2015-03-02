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
#include "Integration/GlobalNamespaceIntegration.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRulesForSquare.hpp"
#include "Geometry/ReferenceSquare.hpp"
#include "GaussQuadratureRulesForLine.hpp"
#include "Geometry/PointReference.hpp"
using Geometry::ReferenceSquare;

//---------------------------------------------------------------------------
namespace QuadratureRules
{
//---------------------------------------------------------------------------
    std::string
    Cn2_1_1::getName() const
    {
        return name_;
    }

    unsigned int
    Cn2_1_1::order() const
    {
        return 1;
    }

    unsigned int
    Cn2_1_1::dimension() const
    {
        return 2;
    }

    unsigned int
    Cn2_1_1::nrOfPoints() const
    {
        return 1;
    }

    double
    Cn2_1_1::weight(unsigned int i) const
    {
        if (i < 1)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }

    void
    Cn2_1_1::getPoint(unsigned int i, PointReferenceT& p) const
    {
        if (i < 1)
            p=gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }

    Cn2_1_1::ReferenceGeometryT*
    Cn2_1_1::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    Cn2_1_1::Cn2_1_1():
        name_("Cn2_1_1"),
        refGeoPtr_(&ReferenceSquare::Instance()),gp_(1,2)
    {
        weight_[0] = ( 2.0 ) * ( 2.0 );
        gp_[0][0] = 0.0;
        gp_[0][1] = 0.0;

    }

    Cn2_1_1::~Cn2_1_1()
    {
    }


//---------------------------------------------------------------------------
    std::string
    Cn2_3_4::getName() const
    {
        return name_;
    }

    unsigned int
    Cn2_3_4::order() const
    {
        return 3;
    }

    unsigned int
    Cn2_3_4::dimension() const
    {
        return 2;
    }

    unsigned int
    Cn2_3_4::nrOfPoints() const
    {
        return 4;
    }

    double
    Cn2_3_4::weight(unsigned int i) const
    {
        if (i < 4)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }

    void
    Cn2_3_4::getPoint(unsigned int i, PointReferenceT& p) const
    {
        if (i < 4)
            p=gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }

    Cn2_3_4::ReferenceGeometryT*
    Cn2_3_4::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    Cn2_3_4::Cn2_3_4():
        name_("Cn2_3_4"),
        refGeoPtr_(&ReferenceSquare::Instance()),gp_(4,2)
    {
        weight_[0] = ( 1.0 ) * ( 1.0 );
        gp_[0][0] = -std::sqrt(3.0) / 3.0;
        gp_[0][1] = -std::sqrt(3.0) / 3.0;

        weight_[1] = ( 1.0 ) * ( 1.0 );
        gp_[1][0] = +std::sqrt(3.0) / 3.0;
        gp_[1][1] = -std::sqrt(3.0) / 3.0;

        weight_[2] = ( 1.0 ) * ( 1.0 );
        gp_[2][0] = -std::sqrt(3.0) / 3.0;
        gp_[2][1] = +std::sqrt(3.0) / 3.0;

        weight_[3] = ( 1.0 ) * ( 1.0 );
        gp_[3][0] = +std::sqrt(3.0) / 3.0;
        gp_[3][1] = +std::sqrt(3.0) / 3.0;

    }

    Cn2_3_4::~Cn2_3_4()
    {}


//---------------------------------------------------------------------------
    std::string
    Cn2_5_9::getName() const
    {
        return name_;
    }

    unsigned int
    Cn2_5_9::order() const
    {
        return 5;
    }

    unsigned int
    Cn2_5_9::dimension() const
    {
        return 2;
    }

    unsigned int
    Cn2_5_9::nrOfPoints() const
    {
        return 9;
    }

    double
    Cn2_5_9::weight(unsigned int i) const
    {
        if (i < 9)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }

    void
    Cn2_5_9::getPoint(unsigned int i, PointReferenceT& p) const
    {
        if (i < 9)
            p=gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }

    Cn2_5_9::ReferenceGeometryT*
    Cn2_5_9::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    Cn2_5_9::Cn2_5_9():
        name_("Cn2_5_9"),
        refGeoPtr_(&ReferenceSquare::Instance()),gp_(9,2)
    {
        weight_[0] = ( 5. / 9. ) * ( 5. / 9. );
        gp_[0][0] = -std::sqrt(3.0 / 5.0);
        gp_[0][1] = -std::sqrt(3.0 / 5.0);

        weight_[1] = ( 8. / 9. ) * ( 5. / 9. );
        gp_[1][0] = 0.0;
        gp_[1][1] = -std::sqrt(3.0 / 5.0);

        weight_[2] = ( 5. / 9. ) * ( 5. / 9. );
        gp_[2][0] = +std::sqrt(3.0 / 5.0);
        gp_[2][1] = -std::sqrt(3.0 / 5.0);

        weight_[3] = ( 5. / 9. ) * ( 8. / 9. );
        gp_[3][0] = -std::sqrt(3.0 / 5.0);
        gp_[3][1] = 0.0;

        weight_[4] = ( 8. / 9. ) * ( 8. / 9. );
        gp_[4][0] = 0.0;
        gp_[4][1] = 0.0;

        weight_[5] = ( 5. / 9. ) * ( 8. / 9. );
        gp_[5][0] = +std::sqrt(3.0 / 5.0);
        gp_[5][1] = 0.0;

        weight_[6] = ( 5. / 9. ) * ( 5. / 9. );
        gp_[6][0] = -std::sqrt(3.0 / 5.0);
        gp_[6][1] = +std::sqrt(3.0 / 5.0);

        weight_[7] = ( 8. / 9. ) * ( 5. / 9. );
        gp_[7][0] = 0.0;
        gp_[7][1] = +std::sqrt(3.0 / 5.0);

        weight_[8] = ( 5. / 9. ) * ( 5. / 9. );
        gp_[8][0] = +std::sqrt(3.0 / 5.0);
        gp_[8][1] = +std::sqrt(3.0 / 5.0);

    }

    Cn2_5_9::~Cn2_5_9()
    {
    
    }


//---------------------------------------------------------------------------
    std::string
    C2_7_4::getName() const
    {
        return name_;
    }

    unsigned int
    C2_7_4::order() const
    {
        return 7;
    }

    unsigned int
    C2_7_4::dimension() const
    {
        return 2;
    }

    unsigned int
    C2_7_4::nrOfPoints() const
    {
        return 16;
    }

    double
    C2_7_4::weight(unsigned int i) const
    {
        if (i < 16)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }

    void
    C2_7_4::getPoint(unsigned int i, PointReferenceT& p) const
    {
        if (i < 16)
            p=gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }

    C2_7_4::ReferenceGeometryT*
    C2_7_4::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    C2_7_4::C2_7_4():
        name_("C2_7_4"),
        refGeoPtr_(&ReferenceSquare::Instance()),gp_(16,2)
    {
        weight_[0] = (59. + 6. * std::sqrt(30.)) / 216.;
        gp_[0][0] = +std::sqrt((15. - 2. * std::sqrt(30.)) / 35.);
        gp_[0][1] = +std::sqrt((15. - 2. * std::sqrt(30.)) / 35.);

        weight_[1] = (59. + 6. * std::sqrt(30.)) / 216.;
        gp_[1][0] = -std::sqrt((15. - 2. * std::sqrt(30.)) / 35.);
        gp_[1][1] = +std::sqrt((15. - 2. * std::sqrt(30.)) / 35.);

        weight_[2] = (59. + 6. * std::sqrt(30.)) / 216.;
        gp_[2][0] = +std::sqrt((15. - 2. * std::sqrt(30.)) / 35.);
        gp_[2][1] = -std::sqrt((15. - 2. * std::sqrt(30.)) / 35.);

        weight_[3] = (59. + 6. * std::sqrt(30.)) / 216.;
        gp_[3][0] = -std::sqrt((15. - 2. * std::sqrt(30.)) / 35.);
        gp_[3][1] = -std::sqrt((15. - 2. * std::sqrt(30.)) / 35.);

        weight_[4] = (59. - 6. * std::sqrt(30.)) / 216.;
        gp_[4][0] = +std::sqrt((15. + 2. * std::sqrt(30.)) / 35.);
        gp_[4][1] = +std::sqrt((15. + 2. * std::sqrt(30.)) / 35.);

        weight_[5] = (59. - 6. * std::sqrt(30.)) / 216.;
        gp_[5][0] = -std::sqrt((15. + 2. * std::sqrt(30.)) / 35.);
        gp_[5][1] = +std::sqrt((15. + 2. * std::sqrt(30.)) / 35.);

        weight_[6] = (59. - 6. * std::sqrt(30.)) / 216.;
        gp_[6][0] = +std::sqrt((15. + 2. * std::sqrt(30.)) / 35.);
        gp_[6][1] = -std::sqrt((15. + 2. * std::sqrt(30.)) / 35.);

        weight_[7] = (59. - 6. * std::sqrt(30.)) / 216.;
        gp_[7][0] = -std::sqrt((15. + 2. * std::sqrt(30.)) / 35.);
        gp_[7][1] = -std::sqrt((15. + 2. * std::sqrt(30.)) / 35.);

        weight_[8] = 49. / 216.;
        gp_[8][0] = +std::sqrt((15. - 2. * std::sqrt(30.)) / 35.);
        gp_[8][1] = +std::sqrt((15. + 2. * std::sqrt(30.)) / 35.);

        weight_[9] = 49. / 216.;
        gp_[9][0] = -std::sqrt((15. - 2. * std::sqrt(30.)) / 35.);
        gp_[9][1] = +std::sqrt((15. + 2. * std::sqrt(30.)) / 35.);

        weight_[10] = 49. / 216.;
        gp_[10][0] = +std::sqrt((15. - 2. * std::sqrt(30.)) / 35.);
        gp_[10][1] = -std::sqrt((15. + 2. * std::sqrt(30.)) / 35.);

        weight_[11] = 49. / 216.;
        gp_[11][0] = -std::sqrt((15. - 2. * std::sqrt(30.)) / 35.);
        gp_[11][1] = -std::sqrt((15. + 2. * std::sqrt(30.)) / 35.);

        weight_[12] = 49. / 216.;
        gp_[12][0] = +std::sqrt((15. + 2. * std::sqrt(30.)) / 35.);
        gp_[12][1] = +std::sqrt((15. - 2. * std::sqrt(30.)) / 35.);

        weight_[13] = 49. / 216.;
        gp_[13][0] = -std::sqrt((15. + 2. * std::sqrt(30.)) / 35.);
        gp_[13][1] = +std::sqrt((15. - 2. * std::sqrt(30.)) / 35.);

        weight_[14] = 49. / 216.;
        gp_[14][0] = +std::sqrt((15. + 2. * std::sqrt(30.)) / 35.);
        gp_[14][1] = -std::sqrt((15. - 2. * std::sqrt(30.)) / 35.);

        weight_[15] = 49. / 216.;
        gp_[15][0] = -std::sqrt((15. + 2. * std::sqrt(30.)) / 35.);
        gp_[15][1] = -std::sqrt((15. - 2. * std::sqrt(30.)) / 35.);

    }

    C2_7_4::~C2_7_4()
    {
    }


//---------------------------------------------------------------------------
    std::string
    C2_9_5::getName() const
    {
        return name_;
    }

    unsigned int
    C2_9_5::order() const
    {
        return 9;
    }

    unsigned int
    C2_9_5::dimension() const
    {
        return 2;
    }

    unsigned int
    C2_9_5::nrOfPoints() const
    {
        return 25;
    }

    double
    C2_9_5::weight(unsigned int i) const
    {
        if (i < 25)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }

    void
    C2_9_5::getPoint(unsigned int i, PointReferenceT& p) const
    {
        if (i < 25)
            p=gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }

    C2_9_5::ReferenceGeometryT*
    C2_9_5::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    C2_9_5::C2_9_5():
        name_("C2_9_5"),
        refGeoPtr_(&ReferenceSquare::Instance()),gp_(25,2)
    {
    	int position(0);
        C1_9_25& ruleForLine = C1_9_25::Instance();
        Geometry::PointReference point1D(1);
        for(int i=0;i<ruleForLine.nrOfPoints();++i){
        	for(int j=0;j<ruleForLine.nrOfPoints();++j){
        		weight_[position]=ruleForLine.weight(i)*ruleForLine.weight(j);
        		ruleForLine.getPoint(i,point1D);
        		gp_[position][0]=point1D[0];
        		ruleForLine.getPoint(j,point1D);
        		gp_[position][1]=point1D[0];
        		++position;
        	}
        }

    }

    C2_9_5::~C2_9_5()
    {
    }


//---------------------------------------------------------------------------
    std::string
    C2_11_6::getName() const
    {
        return name_;
    }

    unsigned int
    C2_11_6::order() const
    {
        return 11;
    }

    unsigned int
    C2_11_6::dimension() const
    {
        return 2;
    }

    unsigned int
    C2_11_6::nrOfPoints() const
    {
        return 36;
    }

    double
    C2_11_6::weight(unsigned int i) const
    {
        if (i < 36)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }

    void
    C2_11_6::getPoint(unsigned int i, PointReferenceT& p) const
    {
        if (i < 36)
            p=gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }

    C2_11_6::ReferenceGeometryT*
    C2_11_6::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    C2_11_6::C2_11_6():
        name_("C2_11_6"),
        refGeoPtr_(&ReferenceSquare::Instance()),gp_(36,2)
    {
    	int position(0);
        C1_11_36& ruleForLine = C1_11_36::Instance();
        Geometry::PointReference point1D(1);
        for(int i=0;i<ruleForLine.nrOfPoints();++i){
        	for(int j=0;j<ruleForLine.nrOfPoints();++j){
        		weight_[position]=ruleForLine.weight(i)*ruleForLine.weight(j);
        		ruleForLine.getPoint(i,point1D);
        		gp_[position][0]=point1D[0];
        		ruleForLine.getPoint(j,point1D);
        		gp_[position][1]=point1D[0];
        		++position;
        	}
        }

    }

    C2_11_6::~C2_11_6()
    {
    }


//---------------------------------------------------------------------------
} // close namespace IntegrationRules

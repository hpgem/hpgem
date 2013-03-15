//------------------------------------------------------------------------------
// File: GaussQuadratureRulesForSquare.cpp 
// Implementation of Gauss quadrature rules for reference square.
// Lars Pesch, Fri Mar  3 12:59:11 CET 2006
//----
// Modified from original file allGaussQuadratureRules.hpp
// by M.T. Julianto, Wed Feb 23 10:45:06 UTC 2013
//---------------------------------------------------------------------------
// System includes and names imported from them:
#include <cmath>
//---------------------------------------------------------------------------
// Package includes:
#include "Integration/GlobalNamespaceIntegration.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRulesForSquare.hpp"

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
    Cn2_1_1::getPoint(unsigned int i, PointReference<2>& p) const
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
        refGeoPtr_(&ReferenceSquare::Instance())
    {
        weight_[0] = ( 2.0 ) * ( 2.0 );
        gp_[0][0] = 0.0;
        gp_[0][1] = 0.0;

        refGeoPtr_->addGaussQuadratureRule(this);
    }

    Cn2_1_1::~Cn2_1_1()
    {
    }

    namespace
    {
        const Cn2_1_1& instantiateCn2_1_1 = Cn2_1_1::Instance();
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
    Cn2_3_4::getPoint(unsigned int i, PointReference<2>& p) const
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
        refGeoPtr_(&ReferenceSquare::Instance())
    {
        weight_[0] = ( 1.0 ) * ( 1.0 );
        gp_[0][0] = -sqrt(3.0) / 3.0;
        gp_[0][1] = -sqrt(3.0) / 3.0;

        weight_[1] = ( 1.0 ) * ( 1.0 );
        gp_[1][0] = +sqrt(3.0) / 3.0;
        gp_[1][1] = -sqrt(3.0) / 3.0;

        weight_[2] = ( 1.0 ) * ( 1.0 );
        gp_[2][0] = -sqrt(3.0) / 3.0;
        gp_[2][1] = +sqrt(3.0) / 3.0;

        weight_[3] = ( 1.0 ) * ( 1.0 );
        gp_[3][0] = +sqrt(3.0) / 3.0;
        gp_[3][1] = +sqrt(3.0) / 3.0;

        refGeoPtr_->addGaussQuadratureRule(this);
    }

    Cn2_3_4::~Cn2_3_4()
    {}

    namespace
    {
        const Cn2_3_4& instantiateCn2_3_4 = Cn2_3_4::Instance();
    }
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
    Cn2_5_9::getPoint(unsigned int i, PointReference<2>& p) const
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
        refGeoPtr_(&ReferenceSquare::Instance())
    {
        weight_[0] = ( 5. / 9. ) * ( 5. / 9. );
        gp_[0][0] = -sqrt(3.0 / 5.0);
        gp_[0][1] = -sqrt(3.0 / 5.0);

        weight_[1] = ( 8. / 9. ) * ( 5. / 9. );
        gp_[1][0] = 0.0;
        gp_[1][1] = -sqrt(3.0 / 5.0);

        weight_[2] = ( 5. / 9. ) * ( 5. / 9. );
        gp_[2][0] = +sqrt(3.0 / 5.0);
        gp_[2][1] = -sqrt(3.0 / 5.0);

        weight_[3] = ( 5. / 9. ) * ( 8. / 9. );
        gp_[3][0] = -sqrt(3.0 / 5.0);
        gp_[3][1] = 0.0;

        weight_[4] = ( 8. / 9. ) * ( 8. / 9. );
        gp_[4][0] = 0.0;
        gp_[4][1] = 0.0;

        weight_[5] = ( 5. / 9. ) * ( 8. / 9. );
        gp_[5][0] = +sqrt(3.0 / 5.0);
        gp_[5][1] = 0.0;

        weight_[6] = ( 5. / 9. ) * ( 5. / 9. );
        gp_[6][0] = -sqrt(3.0 / 5.0);
        gp_[6][1] = +sqrt(3.0 / 5.0);

        weight_[7] = ( 8. / 9. ) * ( 5. / 9. );
        gp_[7][0] = 0.0;
        gp_[7][1] = +sqrt(3.0 / 5.0);

        weight_[8] = ( 5. / 9. ) * ( 5. / 9. );
        gp_[8][0] = +sqrt(3.0 / 5.0);
        gp_[8][1] = +sqrt(3.0 / 5.0);

        refGeoPtr_->addGaussQuadratureRule(this);
    }

    Cn2_5_9::~Cn2_5_9()
    {
    
    }

    namespace
    {
        const Cn2_5_9& instantiateCn2_5_9 = Cn2_5_9::Instance();
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
    C2_7_4::getPoint(unsigned int i, PointReference<2>& p) const
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
        refGeoPtr_(&ReferenceSquare::Instance())
    {
        weight_[0] = (59. + 6. * sqrt(30.)) / 216.;
        gp_[0][0] = +sqrt((15. - 2. * sqrt(30.)) / 35.);
        gp_[0][1] = +sqrt((15. - 2. * sqrt(30.)) / 35.);

        weight_[1] = (59. + 6. * sqrt(30.)) / 216.;
        gp_[1][0] = -sqrt((15. - 2. * sqrt(30.)) / 35.);
        gp_[1][1] = +sqrt((15. - 2. * sqrt(30.)) / 35.);

        weight_[2] = (59. + 6. * sqrt(30.)) / 216.;
        gp_[2][0] = +sqrt((15. - 2. * sqrt(30.)) / 35.);
        gp_[2][1] = -sqrt((15. - 2. * sqrt(30.)) / 35.);

        weight_[3] = (59. + 6. * sqrt(30.)) / 216.;
        gp_[3][0] = -sqrt((15. - 2. * sqrt(30.)) / 35.);
        gp_[3][1] = -sqrt((15. - 2. * sqrt(30.)) / 35.);

        weight_[4] = (59. - 6. * sqrt(30.)) / 216.;
        gp_[4][0] = +sqrt((15. + 2. * sqrt(30.)) / 35.);
        gp_[4][1] = +sqrt((15. + 2. * sqrt(30.)) / 35.);

        weight_[5] = (59. - 6. * sqrt(30.)) / 216.;
        gp_[5][0] = -sqrt((15. + 2. * sqrt(30.)) / 35.);
        gp_[5][1] = +sqrt((15. + 2. * sqrt(30.)) / 35.);

        weight_[6] = (59. - 6. * sqrt(30.)) / 216.;
        gp_[6][0] = +sqrt((15. + 2. * sqrt(30.)) / 35.);
        gp_[6][1] = -sqrt((15. + 2. * sqrt(30.)) / 35.);

        weight_[7] = (59. - 6. * sqrt(30.)) / 216.;
        gp_[7][0] = -sqrt((15. + 2. * sqrt(30.)) / 35.);
        gp_[7][1] = -sqrt((15. + 2. * sqrt(30.)) / 35.);

        weight_[8] = 49. / 216.;
        gp_[8][0] = +sqrt((15. - 2. * sqrt(30.)) / 35.);
        gp_[8][1] = +sqrt((15. + 2. * sqrt(30.)) / 35.);

        weight_[9] = 49. / 216.;
        gp_[9][0] = -sqrt((15. - 2. * sqrt(30.)) / 35.);
        gp_[9][1] = +sqrt((15. + 2. * sqrt(30.)) / 35.);

        weight_[10] = 49. / 216.;
        gp_[10][0] = +sqrt((15. - 2. * sqrt(30.)) / 35.);
        gp_[10][1] = -sqrt((15. + 2. * sqrt(30.)) / 35.);

        weight_[11] = 49. / 216.;
        gp_[11][0] = -sqrt((15. - 2. * sqrt(30.)) / 35.);
        gp_[11][1] = -sqrt((15. + 2. * sqrt(30.)) / 35.);

        weight_[12] = 49. / 216.;
        gp_[12][0] = +sqrt((15. + 2. * sqrt(30.)) / 35.);
        gp_[12][1] = +sqrt((15. - 2. * sqrt(30.)) / 35.);

        weight_[13] = 49. / 216.;
        gp_[13][0] = -sqrt((15. + 2. * sqrt(30.)) / 35.);
        gp_[13][1] = +sqrt((15. - 2. * sqrt(30.)) / 35.);

        weight_[14] = 49. / 216.;
        gp_[14][0] = +sqrt((15. + 2. * sqrt(30.)) / 35.);
        gp_[14][1] = -sqrt((15. - 2. * sqrt(30.)) / 35.);

        weight_[15] = 49. / 216.;
        gp_[15][0] = -sqrt((15. + 2. * sqrt(30.)) / 35.);
        gp_[15][1] = -sqrt((15. - 2. * sqrt(30.)) / 35.);

        refGeoPtr_->addGaussQuadratureRule(this);
    }

    C2_7_4::~C2_7_4()
    {
    }

    namespace
    {
        const C2_7_4& instantiateC2_7_4 = C2_7_4::Instance();
    }
//---------------------------------------------------------------------------
} // close namespace IntegrationRules

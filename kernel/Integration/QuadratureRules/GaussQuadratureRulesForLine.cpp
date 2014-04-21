//------------------------------------------------------------------------------
// File: GaussQuadratureRulesForLine.cpp 
// Implementation of Gauss quadrature rules for reference line.
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
#include "Integration/QuadratureRules/GaussQuadratureRulesForLine.hpp"

#include "Geometry/ReferenceLine.hpp"
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

    unsigned int
    Cn1_1_1::order() const
    {
        return 1;
    }

    unsigned int
    Cn1_1_1::dimension() const
    {
        return 1;
    }

    unsigned int
    Cn1_1_1::nrOfPoints() const
    {
        return 1;
    }

    double
    Cn1_1_1::weight(unsigned int i) const
    {
        if (i < 1)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }

    void
    Cn1_1_1::getPoint(unsigned int i, PointReferenceT& p) const
    {
        if (i < 1)
            p=gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }

    Cn1_1_1::ReferenceGeometryT* Cn1_1_1::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    Cn1_1_1::Cn1_1_1():
        name_("Cn1_1_1"),
        refGeoPtr_(&ReferenceLine::Instance()),gp_(1,1)
    {
        weight_[0] = 2.0;
        gp_[0][0] = 0.0;

        refGeoPtr_->addGaussQuadratureRule(this);
    }

    Cn1_1_1::~Cn1_1_1()
    {
    }

//---------------------------------------------------------------------------
    std::string Cn1_3_4::getName() const
    {
        return name_;
    }

    unsigned int
    Cn1_3_4::order() const
    {
        return 3;
    }

    unsigned int
    Cn1_3_4::dimension() const
    {
        return 1;
    }

    unsigned int
    Cn1_3_4::nrOfPoints() const
    {
        return 2;
    }

    double
    Cn1_3_4::weight(unsigned int i) const
    {
        if (i < 2)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }

    void
    Cn1_3_4::getPoint(unsigned int i, PointReferenceT& p) const
    {
        if (i < 2)
            p=gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }

    Cn1_3_4::ReferenceGeometryT*
    Cn1_3_4::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    Cn1_3_4::Cn1_3_4():
        name_("Cn1_3_4"),
        refGeoPtr_(&ReferenceLine::Instance()),gp_(2,1)
    {
        weight_[0] = 1.0;
        gp_[0][0] = -sqrt(3.0) / 3.0;

        weight_[1] = 1.0;
        gp_[1][0] = +sqrt(3.0) / 3.0;

        refGeoPtr_->addGaussQuadratureRule(this);
    }

    Cn1_3_4::~Cn1_3_4()
    {
    }

//---------------------------------------------------------------------------
    std::string
    Cn1_5_9::getName() const
    {
        return name_;
    }

    unsigned int
    Cn1_5_9::order() const
    {
        return 5;
    }

    unsigned int
    Cn1_5_9::dimension() const
    {
        return 1;
    }

    unsigned int
    Cn1_5_9::nrOfPoints() const
    {
        return 3;
    }

    double
    Cn1_5_9::weight(unsigned int i) const
    {
        if (i < 3)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }

    void
    Cn1_5_9::getPoint(unsigned int i, PointReferenceT& p) const
    {
        if (i < 3)
            p=gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }

    Cn1_5_9::ReferenceGeometryT*
    Cn1_5_9::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    Cn1_5_9::Cn1_5_9():
        name_("Cn1_5_9"),
        refGeoPtr_(&ReferenceLine::Instance()),gp_(3,1)
    {
        weight_[0] = 5. / 9.;
        gp_[0][0] = -sqrt(3.0 / 5.0);

        weight_[1] = 8. / 9.;
        gp_[1][0] = 0.0;

        weight_[2] = 5. / 9.;
        gp_[2][0] = +sqrt(3.0 / 5.0);

        refGeoPtr_->addGaussQuadratureRule(this);
    }

    Cn1_5_9::~Cn1_5_9()
    {
    }


//---------------------------------------------------------------------------
    std::string
    C1_7_x::getName() const
    {
        return name_;
    }

    unsigned int
    C1_7_x::order() const
    {
        return 7;
    }

    unsigned int
    C1_7_x::dimension() const
    {
        return 1;
    }

    unsigned int
    C1_7_x::nrOfPoints() const
    {
        return 4;
    }

    double
    C1_7_x::weight(unsigned int i) const
    {
        if (i < 4)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }

    void
    C1_7_x::getPoint(unsigned int i, PointReferenceT& p) const
    {
        if (i < 4)
            p=gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }

    C1_7_x::ReferenceGeometryT*
    C1_7_x::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    C1_7_x::C1_7_x():
        name_("C1_7_x"),
        refGeoPtr_(&ReferenceLine::Instance()),gp_(4,1)
    {
        weight_[0] = (0.347854845137453);
        gp_[0][0] = (-0.861136311594053);

        weight_[1] = (0.652145154862546);
        gp_[1][0] = (-0.339981043584856);

        weight_[2] = (0.652145154862546);
        gp_[2][0] = (+0.339981043584856);

        weight_[3] = (0.347854845137453);
        gp_[3][0] = (+0.861136311594053);

        refGeoPtr_->addGaussQuadratureRule(this);
    }

    C1_7_x::~C1_7_x()
    {
    }


//---------------------------------------------------------------------------
    std::string
    C1_9_25::getName() const
    {
        return name_;
    }

    unsigned int
    C1_9_25::order() const
    {
        return 9;
    }

    unsigned int
    C1_9_25::dimension() const
    {
        return 1;
    }

    unsigned int
    C1_9_25::nrOfPoints() const
    {
        return 5;
    }

    double
    C1_9_25::weight(unsigned int i) const
    {
        if (i < 5)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }

    void
    C1_9_25::getPoint(unsigned int i, PointReferenceT& p) const
    {
        if (i < 5)
            p=gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }

    C1_9_25::ReferenceGeometryT*
    C1_9_25::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    C1_9_25::C1_9_25():
        name_("C1_9_25"),
        refGeoPtr_(&ReferenceLine::Instance()),gp_(5,1)
    {
        weight_[0] = ( 0.236926885056189);
        gp_[0][0] = (-0.906179845938663);

        weight_[1] = ( 0.478628670499366);
        gp_[1][0] = (-0.538469310105683);

        weight_[2]=(0.56888888888888888888888888);
        gp_[2][0]=0.0;

        weight_[3] = ( 0.478628670499366);
        gp_[3][0] = (0.538469310105683);

        weight_[4] = ( 0.236926885056189);
        gp_[4][0] = (0.906179845938663);


        refGeoPtr_->addGaussQuadratureRule(this);
    }

    C1_9_25::~C1_9_25()
    {
    }


//---------------------------------------------------------------------------
    std::string
    C1_11_36::getName() const
    {
        return name_;
    }

    unsigned int
    C1_11_36::order() const
    {
        return 11;
    }

    unsigned int
    C1_11_36::dimension() const
    {
        return 1;
    }

    unsigned int
    C1_11_36::nrOfPoints() const
    {
        return 6;
    }

    double
    C1_11_36::weight(unsigned int i) const
    {
        if (i < 6)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }

    void
    C1_11_36::getPoint(unsigned int i, PointReferenceT& p) const
    {
        if (i < 6)
            p=gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }

    C1_11_36::ReferenceGeometryT*
    C1_11_36::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    C1_11_36::C1_11_36():
        name_("C1_11_36"),
        refGeoPtr_(&ReferenceLine::Instance()),gp_(6,1)
    {
        weight_[0] = ( 0.171324492379170);
        gp_[0][0] = (-0.932469514203152);

        weight_[1] = ( 0.360761573048138);
        gp_[1][0] = (-0.661209386466264);

        weight_[2]=( 0.467913934572691);
        gp_[2][0]=(-0.238619186083196);

        weight_[3]=( 0.467913934572691);
        gp_[3][0]=(0.238619186083196);

        weight_[4] = ( 0.360761573048138);
        gp_[4][0] = (0.661209386466264);

        weight_[5] = ( 0.171324492379170);
        gp_[5][0] = (0.932469514203152);



        refGeoPtr_->addGaussQuadratureRule(this);
    }

    C1_11_36::~C1_11_36()
    {
    }


//---------------------------------------------------------------------------
} // close namespace IntegrationRules

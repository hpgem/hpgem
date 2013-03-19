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
        refGeoPtr_(&ReferenceLine::Instance())
    {
        weight_[0] = 2.0;
        gp_[0][0] = 0.0;

        refGeoPtr_->addGaussQuadratureRule(this);
    }

    Cn1_1_1::~Cn1_1_1()
    {
    }

    namespace
    {
        const Cn1_1_1& instantiateCn1_1_1 = Cn1_1_1::Instance();
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
        refGeoPtr_(&ReferenceLine::Instance())
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

    namespace
    {
        const Cn1_3_4& instantiateCn1_3_4 = Cn1_3_4::Instance();
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
        refGeoPtr_(&ReferenceLine::Instance())
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

    namespace
    {
        const Cn1_5_9& instantiateCn1_5_9 = Cn1_5_9::Instance();
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
        refGeoPtr_(&ReferenceLine::Instance())
    {
        weight_[0] = (2. * 0.1739274226);
        gp_[0][0] = (-0.861136312);

        weight_[1] = (2. * 0.3260725774);
        gp_[1][0] = (-0.3399810436);

        weight_[2] = (2. * 0.3260725774);
        gp_[2][0] = (+0.3399810436);

        weight_[3] = (2. * 0.1739274226);
        gp_[3][0] = (+0.861136312);

        refGeoPtr_->addGaussQuadratureRule(this);
    }

    C1_7_x::~C1_7_x()
    {
    }

    namespace
    {
        const C1_7_x& instantiateC1_7_x = C1_7_x::Instance();
    }
//---------------------------------------------------------------------------
} // close namespace IntegrationRules

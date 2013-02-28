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

//---------------------------------------------------------------------------
namespace QuadratureRules
{
//---------------------------------------------------------------------------
    std::string Cn1_1_1::getName() const
        {
            return _name;
        }

    DimType  Cn1_1_1::order() const
        {
            return 1;
        }

    DimType Cn1_1_1::dimension() const
        {
            return 1;
        }

    DimType Cn1_1_1::nrOfPoints() const
        {
            return 1;
        }

    NumType Cn1_1_1::weight(DimType i) const
        {
            if (i < 1)
                return _weight[i];
            else
                throw _name + "::weight - wrong index!";
        }

    void Cn1_1_1::getPoint(DimType i, PointReference<1>& p) const
        {
            if (i < 1)
                p=_gp[i];
            else
                throw _name + "::getPoint -  wrong index!";
        }

    Geometry::ReferenceGeometry<1>* Cn1_1_1::forReferenceGeometry() const
        {
            return _refGeoPtr;
        }

    Cn1_1_1::Cn1_1_1()
     : _name("Cn1_1_1"),
       _refGeoPtr(&ReferenceLine::Instance())
            {
                _weight[0] = 2.0;
                _gp[0][0] = 0.0;

                _refGeoPtr->addGaussQuadratureRule(this);
            }

    Cn1_1_1::~Cn1_1_1()
          {}

    namespace
    {
        const Cn1_1_1& instantiateCn1_1_1 = Cn1_1_1::Instance();
    }

//---------------------------------------------------------------------------
    std::string Cn1_3_4::getName() const
        {
            return _name;
        }

    DimType  Cn1_3_4::order() const
        {
            return 3;
        }

    DimType Cn1_3_4::dimension() const
        {
            return 1;
        }

    DimType Cn1_3_4::nrOfPoints() const
        {
            return 2;
        }

    NumType Cn1_3_4::weight(DimType i) const
        {
            if (i < 2)
                return _weight[i];
            else
                throw _name + "::weight - wrong index!";
        }

    void Cn1_3_4::getPoint(DimType i, PointReference<1>& p) const
        {
            if (i < 2)
                p=_gp[i];
            else
                throw _name + "::getPoint -  wrong index!";
        }

    Geometry::ReferenceGeometry<1>* Cn1_3_4::forReferenceGeometry() const
        {
            return _refGeoPtr;
        }

    Cn1_3_4::Cn1_3_4()
     : _name("Cn1_3_4"),
       _refGeoPtr(&ReferenceLine::Instance())
            {
                _weight[0] = 1.0;
                _gp[0][0] = -sqrt(3.0) / 3.0;

                _weight[1] = 1.0;
                _gp[1][0] = +sqrt(3.0) / 3.0;

                _refGeoPtr->addGaussQuadratureRule(this);
            }

    Cn1_3_4::~Cn1_3_4()
          {}

    namespace
    {
        const Cn1_3_4& instantiateCn1_3_4 = Cn1_3_4::Instance();
    }
//---------------------------------------------------------------------------
    std::string Cn1_5_9::getName() const
        {
            return _name;
        }

    DimType  Cn1_5_9::order() const
        {
            return 5;
        }

    DimType Cn1_5_9::dimension() const
        {
            return 1;
        }

    DimType Cn1_5_9::nrOfPoints() const
        {
            return 3;
        }

    NumType Cn1_5_9::weight(DimType i) const
        {
            if (i < 3)
                return _weight[i];
            else
                throw _name + "::weight - wrong index!";
        }

    void Cn1_5_9::getPoint(DimType i, PointReference<1>& p) const
        {
            if (i < 3)
                p=_gp[i];
            else
                throw _name + "::getPoint -  wrong index!";
        }

    Geometry::ReferenceGeometry<1>* Cn1_5_9::forReferenceGeometry() const
        {
            return _refGeoPtr;
        }

    Cn1_5_9::Cn1_5_9()
     : _name("Cn1_5_9"),
       _refGeoPtr(&ReferenceLine::Instance())
            {
                _weight[0] = 5. / 9.;
                _gp[0][0] = -sqrt(3.0 / 5.0);

                _weight[1] = 8. / 9.;
                _gp[1][0] = 0.0;

                _weight[2] = 5. / 9.;
                _gp[2][0] = +sqrt(3.0 / 5.0);

                _refGeoPtr->addGaussQuadratureRule(this);
            }

    Cn1_5_9::~Cn1_5_9()
          {}

    namespace
    {
        const Cn1_5_9& instantiateCn1_5_9 = Cn1_5_9::Instance();
    }
//---------------------------------------------------------------------------
    std::string C1_7_x::getName() const
        {
            return _name;
        }

    DimType  C1_7_x::order() const
        {
            return 7;
        }

    DimType C1_7_x::dimension() const
        {
            return 1;
        }

    DimType C1_7_x::nrOfPoints() const
        {
            return 4;
        }

    NumType C1_7_x::weight(DimType i) const
        {
            if (i < 4)
                return _weight[i];
            else
                throw _name + "::weight - wrong index!";
        }

    void C1_7_x::getPoint(DimType i, PointReference<1>& p) const
        {
            if (i < 4)
                p=_gp[i];
            else
                throw _name + "::getPoint -  wrong index!";
        }

    Geometry::ReferenceGeometry<1>* C1_7_x::forReferenceGeometry() const
        {
            return _refGeoPtr;
        }

    C1_7_x::C1_7_x()
     : _name("C1_7_x"),
       _refGeoPtr(&ReferenceLine::Instance())
            {
                _weight[0] = (2. * 0.1739274226);
                _gp[0][0] = (-0.861136312);

                _weight[1] = (2. * 0.3260725774);
                _gp[1][0] = (-0.3399810436);

                _weight[2] = (2. * 0.3260725774);
                _gp[2][0] = (+0.3399810436);

                _weight[3] = (2. * 0.1739274226);
                _gp[3][0] = (+0.861136312);

                _refGeoPtr->addGaussQuadratureRule(this);
            }

    C1_7_x::~C1_7_x()
          {}

    namespace
    {
        const C1_7_x& instantiateC1_7_x = C1_7_x::Instance();
    }
//---------------------------------------------------------------------------
} // close namespace IntegrationRules

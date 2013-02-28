//------------------------------------------------------------------------------
// File: GaussQuadratureRulesForHypercube.cpp 
// Implementation of Gauss quadrature rules for reference hypercube.
// Lars Pesch, Fri Mar  3 12:59:11 CET 2006
//----
// Modified from original file allGaussQuadratureRules.hpp
// by M.T. Julianto, Wed Feb 25 10:45:06 UTC 2013
//---------------------------------------------------------------------------
// System includes and names imported from them:
#include <cmath>
//---------------------------------------------------------------------------
// Package includes:
#include "Integration/GlobalNamespaceIntegration.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRulesForHypercube.hpp"

//---------------------------------------------------------------------------
namespace QuadratureRules
{
//---------------------------------------------------------------------------
    std::string Cn4_1_1::getName() const
        {
            return _name;
        }

    DimType  Cn4_1_1::order() const
        {
            return 1;
        }

    DimType Cn4_1_1::dimension() const
        {
            return 4;
        }

    DimType Cn4_1_1::nrOfPoints() const
        {
            return 1;
        }

    NumType Cn4_1_1::weight(DimType i) const
        {
            if (i < 1)
                return _weight[i];
            else
                throw _name + "::weight - wrong index!";
        }

    void Cn4_1_1::getPoint(DimType i, PointReference<4>& p) const
        {
            if (i < 1)
                p=_gp[i];
            else
                throw _name + "::getPoint -  wrong index!";
        }

    Geometry::ReferenceGeometry<4>* Cn4_1_1::forReferenceGeometry() const
        {
            return _refGeoPtr;
        }

    Cn4_1_1::Cn4_1_1()
     : _name("Cn4_1_1"),
       _refGeoPtr(&ReferenceHypercube::Instance())
            {
                _weight[0] = 16.0;
                _gp[0][0] = 0.0;
                _gp[0][1] = 0.0;
                _gp[0][2] = 0.0;
                _gp[0][3] = 0.0;

                _refGeoPtr->addGaussQuadratureRule(this);
            }

    Cn4_1_1::~Cn4_1_1()
          {}

    namespace
    {
        const Cn4_1_1& instantiateCn4_1_1 = Cn4_1_1::Instance();
    }
//---------------------------------------------------------------------------
    std::string Cn4_3_4::getName() const
        {
            return _name;
        }

    DimType  Cn4_3_4::order() const
        {
            return 3;
        }

    DimType Cn4_3_4::dimension() const
        {
            return 4;
        }

    DimType Cn4_3_4::nrOfPoints() const
        {
            return 16;
        }

    NumType Cn4_3_4::weight(DimType i) const
        {
            if (i < 16)
                return _weight[i];
            else
                throw _name + "::weight - wrong index!";
        }

    void Cn4_3_4::getPoint(DimType i, PointReference<4>& p) const
        {
            if (i < 16)
                p=_gp[i];
            else
                throw _name + "::getPoint -  wrong index!";
        }

    Geometry::ReferenceGeometry<4>* Cn4_3_4::forReferenceGeometry() const
        {
            return _refGeoPtr;
        }

    Cn4_3_4::Cn4_3_4()
     : _name("Cn4_3_4"),
       _refGeoPtr(&ReferenceHypercube::Instance())
            {
                _weight[0] = 1.0;
                _gp[0][0] = -sqrt(3.0) / 3.0;
                _gp[0][1] = -sqrt(3.0) / 3.0;
                _gp[0][2] = -sqrt(3.0) / 3.0;
                _gp[0][3] = -sqrt(3.0) / 3.0;

                _weight[1] = 1.0;
                _gp[1][0] = +sqrt(3.0) / 3.0;
                _gp[1][1] = -sqrt(3.0) / 3.0;
                _gp[1][2] = -sqrt(3.0) / 3.0;
                _gp[1][3] = -sqrt(3.0) / 3.0;

                _weight[2] = 1.0;
                _gp[2][0] = -sqrt(3.0) / 3.0;
                _gp[2][1] = +sqrt(3.0) / 3.0;
                _gp[2][2] = -sqrt(3.0) / 3.0;
                _gp[2][3] = -sqrt(3.0) / 3.0;

                _weight[3] = 1.0;
                _gp[3][0] = +sqrt(3.0) / 3.0;
                _gp[3][1] = +sqrt(3.0) / 3.0;
                _gp[3][2] = -sqrt(3.0) / 3.0;
                _gp[3][3] = -sqrt(3.0) / 3.0;

                _weight[4] = 1.0;
                _gp[4][0] = -sqrt(3.0) / 3.0;
                _gp[4][1] = -sqrt(3.0) / 3.0;
                _gp[4][2] = +sqrt(3.0) / 3.0;
                _gp[4][3] = -sqrt(3.0) / 3.0;

                _weight[5] = 1.0;
                _gp[5][0] = +sqrt(3.0) / 3.0;
                _gp[5][1] = -sqrt(3.0) / 3.0;
                _gp[5][2] = +sqrt(3.0) / 3.0;
                _gp[5][3] = -sqrt(3.0) / 3.0;

                _weight[6] = 1.0;
                _gp[6][0] = -sqrt(3.0) / 3.0;
                _gp[6][1] = +sqrt(3.0) / 3.0;
                _gp[6][2] = +sqrt(3.0) / 3.0;
                _gp[6][3] = -sqrt(3.0) / 3.0;

                _weight[7] = 1.0;
                _gp[7][0] = +sqrt(3.0) / 3.0;
                _gp[7][1] = +sqrt(3.0) / 3.0;
                _gp[7][2] = +sqrt(3.0) / 3.0;
                _gp[7][3] = -sqrt(3.0) / 3.0;

                _weight[8] = 1.0;
                _gp[8][0] = -sqrt(3.0) / 3.0;
                _gp[8][1] = -sqrt(3.0) / 3.0;
                _gp[8][2] = -sqrt(3.0) / 3.0;
                _gp[8][3] = +sqrt(3.0) / 3.0;

                _weight[9] = 1.0;
                _gp[9][0] = +sqrt(3.0) / 3.0;
                _gp[9][1] = -sqrt(3.0) / 3.0;
                _gp[9][2] = -sqrt(3.0) / 3.0;
                _gp[9][3] = +sqrt(3.0) / 3.0;

                _weight[10] = 1.0;
                _gp[10][0] = -sqrt(3.0) / 3.0;
                _gp[10][1] = +sqrt(3.0) / 3.0;
                _gp[10][2] = -sqrt(3.0) / 3.0;
                _gp[10][3] = +sqrt(3.0) / 3.0;

                _weight[11] = 1.0;
                _gp[11][0] = +sqrt(3.0) / 3.0;
                _gp[11][1] = +sqrt(3.0) / 3.0;
                _gp[11][2] = -sqrt(3.0) / 3.0;
                _gp[11][3] = +sqrt(3.0) / 3.0;

                _weight[12] = 1.0;
                _gp[12][0] = -sqrt(3.0) / 3.0;
                _gp[12][1] = -sqrt(3.0) / 3.0;
                _gp[12][2] = +sqrt(3.0) / 3.0;
                _gp[12][3] = +sqrt(3.0) / 3.0;

                _weight[13] = 1.0;
                _gp[13][0] = +sqrt(3.0) / 3.0;
                _gp[13][1] = -sqrt(3.0) / 3.0;
                _gp[13][2] = +sqrt(3.0) / 3.0;
                _gp[13][3] = +sqrt(3.0) / 3.0;

                _weight[14] = 1.0;
                _gp[14][0] = -sqrt(3.0) / 3.0;
                _gp[14][1] = +sqrt(3.0) / 3.0;
                _gp[14][2] = +sqrt(3.0) / 3.0;
                _gp[14][3] = +sqrt(3.0) / 3.0;

                _weight[15] = 1.0;
                _gp[15][0] = +sqrt(3.0) / 3.0;
                _gp[15][1] = +sqrt(3.0) / 3.0;
                _gp[15][2] = +sqrt(3.0) / 3.0;
                _gp[15][3] = +sqrt(3.0) / 3.0;

                _refGeoPtr->addGaussQuadratureRule(this);
            }

    Cn4_3_4::~Cn4_3_4()
          {}

    namespace
    {
        const Cn4_3_4& instantiateCn4_3_4 = Cn4_3_4::Instance();
    }

//---------------------------------------------------------------------------
} // close namespace IntegrationRules

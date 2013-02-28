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
    std::string Cn2_1_1::getName() const
        {
            return _name;
        }

    DimType  Cn2_1_1::order() const
        {
            return 1;
        }

    DimType Cn2_1_1::dimension() const
        {
            return 2;
        }

    DimType Cn2_1_1::nrOfPoints() const
        {
            return 1;
        }

    NumType Cn2_1_1::weight(DimType i) const
        {
            if (i < 1)
                return _weight[i];
            else
                throw _name + "::weight - wrong index!";
        }

    void Cn2_1_1::getPoint(DimType i, PointReference<2>& p) const
        {
            if (i < 1)
                p=_gp[i];
            else
                throw _name + "::getPoint -  wrong index!";
        }

    Geometry::ReferenceGeometry<2>* Cn2_1_1::forReferenceGeometry() const
        {
            return _refGeoPtr;
        }

    Cn2_1_1::Cn2_1_1()
     : _name("Cn2_1_1"),
       _refGeoPtr(&ReferenceSquare::Instance())
            {
                _weight[0] = ( 2.0 ) * ( 2.0 );
                _gp[0][0] = 0.0;
                _gp[0][1] = 0.0;

                _refGeoPtr->addGaussQuadratureRule(this);
           }

    Cn2_1_1::~Cn2_1_1()
          {}

    namespace
    {
        const Cn2_1_1& instantiateCn2_1_1 = Cn2_1_1::Instance();
    }
//---------------------------------------------------------------------------
    std::string Cn2_3_4::getName() const
        {
            return _name;
        }

    DimType  Cn2_3_4::order() const
        {
            return 3;
        }

    DimType Cn2_3_4::dimension() const
        {
            return 2;
        }

    DimType Cn2_3_4::nrOfPoints() const
        {
            return 4;
        }

    NumType Cn2_3_4::weight(DimType i) const
        {
            if (i < 4)
                return _weight[i];
            else
                throw _name + "::weight - wrong index!";
        }

    void Cn2_3_4::getPoint(DimType i, PointReference<2>& p) const
        {
            if (i < 4)
                p=_gp[i];
            else
                throw _name + "::getPoint -  wrong index!";
        }

    Geometry::ReferenceGeometry<2>* Cn2_3_4::forReferenceGeometry() const
        {
            return _refGeoPtr;
        }

    Cn2_3_4::Cn2_3_4()
     : _name("Cn2_3_4"),
       _refGeoPtr(&ReferenceSquare::Instance())
            {
                _weight[0] = ( 1.0 ) * ( 1.0 );
                _gp[0][0] = -sqrt(3.0) / 3.0;
                _gp[0][1] = -sqrt(3.0) / 3.0;

                _weight[1] = ( 1.0 ) * ( 1.0 );
                _gp[1][0] = +sqrt(3.0) / 3.0;
                _gp[1][1] = -sqrt(3.0) / 3.0;

                _weight[2] = ( 1.0 ) * ( 1.0 );
                _gp[2][0] = -sqrt(3.0) / 3.0;
                _gp[2][1] = +sqrt(3.0) / 3.0;

                _weight[3] = ( 1.0 ) * ( 1.0 );
                _gp[3][0] = +sqrt(3.0) / 3.0;
                _gp[3][1] = +sqrt(3.0) / 3.0;

                _refGeoPtr->addGaussQuadratureRule(this);
            }

    Cn2_3_4::~Cn2_3_4()
          {}

    namespace
    {
        const Cn2_3_4& instantiateCn2_3_4 = Cn2_3_4::Instance();
    }
//---------------------------------------------------------------------------
    std::string Cn2_5_9::getName() const
        {
            return _name;
        }

    DimType  Cn2_5_9::order() const
        {
            return 5;
        }

    DimType Cn2_5_9::dimension() const
        {
            return 2;
        }

    DimType Cn2_5_9::nrOfPoints() const
        {
            return 9;
        }

    NumType Cn2_5_9::weight(DimType i) const
        {
            if (i < 9)
                return _weight[i];
            else
                throw _name + "::weight - wrong index!";
        }

    void Cn2_5_9::getPoint(DimType i, PointReference<2>& p) const
        {
            if (i < 9)
                p=_gp[i];
            else
                throw _name + "::getPoint -  wrong index!";
        }

    Geometry::ReferenceGeometry<2>* Cn2_5_9::forReferenceGeometry() const
        {
            return _refGeoPtr;
        }

    Cn2_5_9::Cn2_5_9()
     : _name("Cn2_5_9"),
       _refGeoPtr(&ReferenceSquare::Instance())
            {
                _weight[0] = ( 5. / 9. ) * ( 5. / 9. );
                _gp[0][0] = -sqrt(3.0 / 5.0);
                _gp[0][1] = -sqrt(3.0 / 5.0);

                _weight[1] = ( 8. / 9. ) * ( 5. / 9. );
                _gp[1][0] = 0.0;
                _gp[1][1] = -sqrt(3.0 / 5.0);

                _weight[2] = ( 5. / 9. ) * ( 5. / 9. );
                _gp[2][0] = +sqrt(3.0 / 5.0);
                _gp[2][1] = -sqrt(3.0 / 5.0);

                _weight[3] = ( 5. / 9. ) * ( 8. / 9. );
                _gp[3][0] = -sqrt(3.0 / 5.0);
                _gp[3][1] = 0.0;

                _weight[4] = ( 8. / 9. ) * ( 8. / 9. );
                _gp[4][0] = 0.0;
                _gp[4][1] = 0.0;

                _weight[5] = ( 5. / 9. ) * ( 8. / 9. );
                _gp[5][0] = +sqrt(3.0 / 5.0);
                _gp[5][1] = 0.0;

                _weight[6] = ( 5. / 9. ) * ( 5. / 9. );
                _gp[6][0] = -sqrt(3.0 / 5.0);
                _gp[6][1] = +sqrt(3.0 / 5.0);

                _weight[7] = ( 8. / 9. ) * ( 5. / 9. );
                _gp[7][0] = 0.0;
                _gp[7][1] = +sqrt(3.0 / 5.0);

                _weight[8] = ( 5. / 9. ) * ( 5. / 9. );
                _gp[8][0] = +sqrt(3.0 / 5.0);
                _gp[8][1] = +sqrt(3.0 / 5.0);

                _refGeoPtr->addGaussQuadratureRule(this);
            }

    Cn2_5_9::~Cn2_5_9()
          {}

    namespace
    {
        const Cn2_5_9& instantiateCn2_5_9 = Cn2_5_9::Instance();
    }
//---------------------------------------------------------------------------
    std::string C2_7_4::getName() const
        {
            return _name;
        }

    DimType  C2_7_4::order() const
        {
            return 7;
        }

    DimType C2_7_4::dimension() const
        {
            return 2;
        }

    DimType C2_7_4::nrOfPoints() const
        {
            return 16;
        }

    NumType C2_7_4::weight(DimType i) const
        {
            if (i < 16)
                return _weight[i];
            else
                throw _name + "::weight - wrong index!";
        }

    void C2_7_4::getPoint(DimType i, PointReference<2>& p) const
        {
            if (i < 16)
                p=_gp[i];
            else
                throw _name + "::getPoint -  wrong index!";
        }

    Geometry::ReferenceGeometry<2>* C2_7_4::forReferenceGeometry() const
        {
            return _refGeoPtr;
        }

    C2_7_4::C2_7_4()
     : _name("C2_7_4"),
       _refGeoPtr(&ReferenceSquare::Instance())
            {
                _weight[0] = (59. + 6. * sqrt(30.)) / 216.;
                _gp[0][0] = +sqrt((15. - 2. * sqrt(30.)) / 35.);
                _gp[0][1] = +sqrt((15. - 2. * sqrt(30.)) / 35.);

                _weight[1] = (59. + 6. * sqrt(30.)) / 216.;
                _gp[1][0] = -sqrt((15. - 2. * sqrt(30.)) / 35.);
                _gp[1][1] = +sqrt((15. - 2. * sqrt(30.)) / 35.);

                _weight[2] = (59. + 6. * sqrt(30.)) / 216.;
                _gp[2][0] = +sqrt((15. - 2. * sqrt(30.)) / 35.);
                _gp[2][1] = -sqrt((15. - 2. * sqrt(30.)) / 35.);

                _weight[3] = (59. + 6. * sqrt(30.)) / 216.;
                _gp[3][0] = -sqrt((15. - 2. * sqrt(30.)) / 35.);
                _gp[3][1] = -sqrt((15. - 2. * sqrt(30.)) / 35.);

                _weight[4] = (59. - 6. * sqrt(30.)) / 216.;
                _gp[4][0] = +sqrt((15. + 2. * sqrt(30.)) / 35.);
                _gp[4][1] = +sqrt((15. + 2. * sqrt(30.)) / 35.);

                _weight[5] = (59. - 6. * sqrt(30.)) / 216.;
                _gp[5][0] = -sqrt((15. + 2. * sqrt(30.)) / 35.);
                _gp[5][1] = +sqrt((15. + 2. * sqrt(30.)) / 35.);

                _weight[6] = (59. - 6. * sqrt(30.)) / 216.;
                _gp[6][0] = +sqrt((15. + 2. * sqrt(30.)) / 35.);
                _gp[6][1] = -sqrt((15. + 2. * sqrt(30.)) / 35.);

                _weight[7] = (59. - 6. * sqrt(30.)) / 216.;
                _gp[7][0] = -sqrt((15. + 2. * sqrt(30.)) / 35.);
                _gp[7][1] = -sqrt((15. + 2. * sqrt(30.)) / 35.);

                _weight[8] = 49. / 216.;
                _gp[8][0] = +sqrt((15. - 2. * sqrt(30.)) / 35.);
                _gp[8][1] = +sqrt((15. + 2. * sqrt(30.)) / 35.);

                _weight[9] = 49. / 216.;
                _gp[9][0] = -sqrt((15. - 2. * sqrt(30.)) / 35.);
                _gp[9][1] = +sqrt((15. + 2. * sqrt(30.)) / 35.);

                _weight[10] = 49. / 216.;
                _gp[10][0] = +sqrt((15. - 2. * sqrt(30.)) / 35.);
                _gp[10][1] = -sqrt((15. + 2. * sqrt(30.)) / 35.);

                _weight[11] = 49. / 216.;
                _gp[11][0] = -sqrt((15. - 2. * sqrt(30.)) / 35.);
                _gp[11][1] = -sqrt((15. + 2. * sqrt(30.)) / 35.);

                _weight[12] = 49. / 216.;
                _gp[12][0] = +sqrt((15. + 2. * sqrt(30.)) / 35.);
                _gp[12][1] = +sqrt((15. - 2. * sqrt(30.)) / 35.);

                _weight[13] = 49. / 216.;
                _gp[13][0] = -sqrt((15. + 2. * sqrt(30.)) / 35.);
                _gp[13][1] = +sqrt((15. - 2. * sqrt(30.)) / 35.);

                _weight[14] = 49. / 216.;
                _gp[14][0] = +sqrt((15. + 2. * sqrt(30.)) / 35.);
                _gp[14][1] = -sqrt((15. - 2. * sqrt(30.)) / 35.);

                _weight[15] = 49. / 216.;
                _gp[15][0] = -sqrt((15. + 2. * sqrt(30.)) / 35.);
                _gp[15][1] = -sqrt((15. - 2. * sqrt(30.)) / 35.);

                _refGeoPtr->addGaussQuadratureRule(this);
            }

    C2_7_4::~C2_7_4()
          {}

    namespace
    {
        const C2_7_4& instantiateC2_7_4 = C2_7_4::Instance();
    }
//---------------------------------------------------------------------------
} // close namespace IntegrationRules

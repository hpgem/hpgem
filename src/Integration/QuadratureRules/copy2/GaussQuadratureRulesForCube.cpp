//------------------------------------------------------------------------------
// File: GaussQuadratureRulesForCube.cpp 
// Implementation of Gauss quadrature rules for reference cube.
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
#include "Integration/QuadratureRules/GaussQuadratureRulesForCube.hpp"

//---------------------------------------------------------------------------
namespace QuadratureRules
{
//---------------------------------------------------------------------------
    std::string Cn3_1_1::getName() const
        {
            return _name;
        }

    unsigned int  Cn3_1_1::order() const
        {
            return 1;
        }

    unsigned int Cn3_1_1::dimension() const
        {
            return 3;
        }

    unsigned int Cn3_1_1::nrOfPoints() const
        {
            return 1;
        }

    NumType Cn3_1_1::weight(unsigned int i) const
        {
            if (i < 1)
                return _weight[i];
            else
                throw _name + "::weight - wrong index!";
        }

    void Cn3_1_1::getPoint(unsigned int i, PointReference<3>& p) const
        {
            if (i < 1)
                p=_gp[i];
            else
                throw _name + "::getPoint -  wrong index!";
        }

    Geometry::ReferenceGeometry<3>* Cn3_1_1::forReferenceGeometry() const
        {
            return _refGeoPtr;
        }

    Cn3_1_1::Cn3_1_1()
     : _name("Cn3_1_1"),
       _refGeoPtr(&ReferenceCube::Instance())
            {
                _weight[0] = ( ( 2.0 ) * ( 2.0 ) ) * ( 2.0 );
                _gp[0][0] = 0.0;
                _gp[0][1] = 0.0;
                _gp[0][2] = 0.0;

                _refGeoPtr->addGaussQuadratureRule(this);
            }

    Cn3_1_1::~Cn3_1_1()
          {}

    namespace
    {
        const Cn3_1_1& instantiateCn3_1_1 = Cn3_1_1::Instance();
    }
//---------------------------------------------------------------------------
    std::string Cn3_3_4::getName() const
        {
            return _name;
        }

    unsigned int  Cn3_3_4::order() const
        {
            return 3;
        }

    unsigned int Cn3_3_4::dimension() const
        {
            return 3;
        }

    unsigned int Cn3_3_4::nrOfPoints() const
        {
            return 8;
        }

    NumType Cn3_3_4::weight(unsigned int i) const
        {
            if (i < 8)
                return _weight[i];
            else
                throw _name + "::weight - wrong index!";
        }

    void Cn3_3_4::getPoint(unsigned int i, PointReference<3>& p) const
        {
            if (i < 8)
                p=_gp[i];
            else
                throw _name + "::getPoint -  wrong index!";
        }

    Geometry::ReferenceGeometry<3>* Cn3_3_4::forReferenceGeometry() const
        {
            return _refGeoPtr;
        }

    Cn3_3_4::Cn3_3_4()
     : _name("Cn3_3_4"),
       _refGeoPtr(&ReferenceCube::Instance())
            {
                _weight[0] = ( ( 1.0 ) * ( 1.0 ) ) * ( 1.0 );
                _gp[0][0] = -sqrt(3.0) / 3.0;
                _gp[0][1] = -sqrt(3.0) / 3.0;
                _gp[0][2] = -sqrt(3.0) / 3.0;

                _weight[1] = ( ( 1.0 ) * ( 1.0 ) ) * ( 1.0 );
                _gp[1][0] = +sqrt(3.0) / 3.0;
                _gp[1][1] = -sqrt(3.0) / 3.0;
                _gp[1][2] = -sqrt(3.0) / 3.0;

                _weight[2] = ( ( 1.0 ) * ( 1.0 ) ) * ( 1.0 );
                _gp[2][0] = -sqrt(3.0) / 3.0;
                _gp[2][1] = +sqrt(3.0) / 3.0;
                _gp[2][2] = -sqrt(3.0) / 3.0;

                _weight[3] = ( ( 1.0 ) * ( 1.0 ) ) * ( 1.0 );
                _gp[3][0] = +sqrt(3.0) / 3.0;
                _gp[3][1] = +sqrt(3.0) / 3.0;
                _gp[3][2] = -sqrt(3.0) / 3.0;

                _weight[4] = ( ( 1.0 ) * ( 1.0 ) ) * ( 1.0 );
                _gp[4][0] = -sqrt(3.0) / 3.0;
                _gp[4][1] = -sqrt(3.0) / 3.0;
                _gp[4][2] = +sqrt(3.0) / 3.0;

                _weight[5] = ( ( 1.0 ) * ( 1.0 ) ) * ( 1.0 );
                _gp[5][0] = +sqrt(3.0) / 3.0;
                _gp[5][1] = -sqrt(3.0) / 3.0;
                _gp[5][2] = +sqrt(3.0) / 3.0;

                _weight[6] = ( ( 1.0 ) * ( 1.0 ) ) * ( 1.0 );
                _gp[6][0] = -sqrt(3.0) / 3.0;
                _gp[6][1] = +sqrt(3.0) / 3.0;
                _gp[6][2] = +sqrt(3.0) / 3.0;

                _weight[7] = ( ( 1.0 ) * ( 1.0 ) ) * ( 1.0 );
                _gp[7][0] = +sqrt(3.0) / 3.0;
                _gp[7][1] = +sqrt(3.0) / 3.0;
                _gp[7][2] = +sqrt(3.0) / 3.0;

                _refGeoPtr->addGaussQuadratureRule(this);
            }

    Cn3_3_4::~Cn3_3_4()
          {}

    namespace
    {
        const Cn3_3_4& instantiateCn3_3_4 = Cn3_3_4::Instance();
    }
//---------------------------------------------------------------------------
    std::string Cn3_5_9::getName() const
        {
            return _name;
        }

    unsigned int  Cn3_5_9::order() const
        {
            return 5;
        }

    unsigned int Cn3_5_9::dimension() const
        {
            return 3;
        }

    unsigned int Cn3_5_9::nrOfPoints() const
        {
            return 27;
        }

    NumType Cn3_5_9::weight(unsigned int i) const
        {
            if (i < 27)
                return _weight[i];
            else
                throw _name + "::weight - wrong index!";
        }

    void Cn3_5_9::getPoint(unsigned int i, PointReference<3>& p) const
        {
            if (i < 27)
                p=_gp[i];
            else
                throw _name + "::getPoint -  wrong index!";
        }

    Geometry::ReferenceGeometry<3>* Cn3_5_9::forReferenceGeometry() const
        {
            return _refGeoPtr;
        }

    Cn3_5_9::Cn3_5_9()
     : _name("Cn3_5_9"),
       _refGeoPtr(&ReferenceCube::Instance())
            {
                _weight[0] = ( ( 5. / 9. ) * ( 5. / 9. ) ) * ( 5. / 9. );
                _gp[0][0] = -sqrt(3.0 / 5.0);
                _gp[0][1] = -sqrt(3.0 / 5.0);
                _gp[0][2] = -sqrt(3.0 / 5.0);

                _weight[1] = ( ( 8. / 9. ) * ( 5. / 9. ) ) * ( 5. / 9. );
                _gp[1][0] = 0.0;
                _gp[1][1] = -sqrt(3.0 / 5.0);
                _gp[1][2] = -sqrt(3.0 / 5.0);

                _weight[2] = ( ( 5. / 9. ) * ( 5. / 9. ) ) * ( 5. / 9. );
                _gp[2][0] = +sqrt(3.0 / 5.0);
                _gp[2][1] = -sqrt(3.0 / 5.0);
                _gp[2][2] = -sqrt(3.0 / 5.0);

                _weight[3] = ( ( 5. / 9. ) * ( 8. / 9. ) ) * ( 5. / 9. );
                _gp[3][0] = -sqrt(3.0 / 5.0);
                _gp[3][1] = 0.0;
                _gp[3][2] = -sqrt(3.0 / 5.0);

                _weight[4] = ( ( 8. / 9. ) * ( 8. / 9. ) ) * ( 5. / 9. );
                _gp[4][0] = 0.0;
                _gp[4][1] = 0.0;
                _gp[4][2] = -sqrt(3.0 / 5.0);

                _weight[5] = ( ( 5. / 9. ) * ( 8. / 9. ) ) * ( 5. / 9. );
                _gp[5][0] = +sqrt(3.0 / 5.0);
                _gp[5][1] = 0.0;
                _gp[5][2] = -sqrt(3.0 / 5.0);

                _weight[6] = ( ( 5. / 9. ) * ( 5. / 9. ) ) * ( 5. / 9. );
                _gp[6][0] = -sqrt(3.0 / 5.0);
                _gp[6][1] = +sqrt(3.0 / 5.0);
                _gp[6][2] = -sqrt(3.0 / 5.0);

                _weight[7] = ( ( 8. / 9. ) * ( 5. / 9. ) ) * ( 5. / 9. );
                _gp[7][0] = 0.0;
                _gp[7][1] = +sqrt(3.0 / 5.0);
                _gp[7][2] = -sqrt(3.0 / 5.0);

                _weight[8] = ( ( 5. / 9. ) * ( 5. / 9. ) ) * ( 5. / 9. );
                _gp[8][0] = +sqrt(3.0 / 5.0);
                _gp[8][1] = +sqrt(3.0 / 5.0);
                _gp[8][2] = -sqrt(3.0 / 5.0);

                _weight[9] = ( ( 5. / 9. ) * ( 5. / 9. ) ) * ( 8. / 9. );
                _gp[9][0] = -sqrt(3.0 / 5.0);
                _gp[9][1] = -sqrt(3.0 / 5.0);
                _gp[9][2] = 0.0;

                _weight[10] = ( ( 8. / 9. ) * ( 5. / 9. ) ) * ( 8. / 9. );
                _gp[10][0] = 0.0;
                _gp[10][1] = -sqrt(3.0 / 5.0);
                _gp[10][2] = 0.0;

                _weight[11] = ( ( 5. / 9. ) * ( 5. / 9. ) ) * ( 8. / 9. );
                _gp[11][0] = +sqrt(3.0 / 5.0);
                _gp[11][1] = -sqrt(3.0 / 5.0);
                _gp[11][2] = 0.0;

                _weight[12] = ( ( 5. / 9. ) * ( 8. / 9. ) ) * ( 8. / 9. );
                _gp[12][0] = -sqrt(3.0 / 5.0);
                _gp[12][1] = 0.0;
                _gp[12][2] = 0.0;

                _weight[13] = ( ( 8. / 9. ) * ( 8. / 9. ) ) * ( 8. / 9. );
                _gp[13][0] = 0.0;
                _gp[13][1] = 0.0;
                _gp[13][2] = 0.0;

                _weight[14] = ( ( 5. / 9. ) * ( 8. / 9. ) ) * ( 8. / 9. );
                _gp[14][0] = +sqrt(3.0 / 5.0);
                _gp[14][1] = 0.0;
                _gp[14][2] = 0.0;

                _weight[15] = ( ( 5. / 9. ) * ( 5. / 9. ) ) * ( 8. / 9. );
                _gp[15][0] = -sqrt(3.0 / 5.0);
                _gp[15][1] = +sqrt(3.0 / 5.0);
                _gp[15][2] = 0.0;

                _weight[16] = ( ( 8. / 9. ) * ( 5. / 9. ) ) * ( 8. / 9. );
                _gp[16][0] = 0.0;
                _gp[16][1] = +sqrt(3.0 / 5.0);
                _gp[16][2] = 0.0;

                _weight[17] = ( ( 5. / 9. ) * ( 5. / 9. ) ) * ( 8. / 9. );
                _gp[17][0] = +sqrt(3.0 / 5.0);
                _gp[17][1] = +sqrt(3.0 / 5.0);
                _gp[17][2] = 0.0;

                _weight[18] = ( ( 5. / 9. ) * ( 5. / 9. ) ) * ( 5. / 9. );
                _gp[18][0] = -sqrt(3.0 / 5.0);
                _gp[18][1] = -sqrt(3.0 / 5.0);
                _gp[18][2] = +sqrt(3.0 / 5.0);

                _weight[19] = ( ( 8. / 9. ) * ( 5. / 9. ) ) * ( 5. / 9. );
                _gp[19][0] = 0.0;
                _gp[19][1] = -sqrt(3.0 / 5.0);
                _gp[19][2] = +sqrt(3.0 / 5.0);

                _weight[20] = ( ( 5. / 9. ) * ( 5. / 9. ) ) * ( 5. / 9. );
                _gp[20][0] = +sqrt(3.0 / 5.0);
                _gp[20][1] = -sqrt(3.0 / 5.0);
                _gp[20][2] = +sqrt(3.0 / 5.0);

                _weight[21] = ( ( 5. / 9. ) * ( 8. / 9. ) ) * ( 5. / 9. );
                _gp[21][0] = -sqrt(3.0 / 5.0);
                _gp[21][1] = 0.0;
                _gp[21][2] = +sqrt(3.0 / 5.0);

                _weight[22] = ( ( 8. / 9. ) * ( 8. / 9. ) ) * ( 5. / 9. );
                _gp[22][0] = 0.0;
                _gp[22][1] = 0.0;
                _gp[22][2] = +sqrt(3.0 / 5.0);

                _weight[23] = ( ( 5. / 9. ) * ( 8. / 9. ) ) * ( 5. / 9. );
                _gp[23][0] = +sqrt(3.0 / 5.0);
                _gp[23][1] = 0.0;
                _gp[23][2] = +sqrt(3.0 / 5.0);

                _weight[24] = ( ( 5. / 9. ) * ( 5. / 9. ) ) * ( 5. / 9. );
                _gp[24][0] = -sqrt(3.0 / 5.0);
                _gp[24][1] = +sqrt(3.0 / 5.0);
                _gp[24][2] = +sqrt(3.0 / 5.0);

                _weight[25] = ( ( 8. / 9. ) * ( 5. / 9. ) ) * ( 5. / 9. );
                _gp[25][0] = 0.0;
                _gp[25][1] = +sqrt(3.0 / 5.0);
                _gp[25][2] = +sqrt(3.0 / 5.0);

                _weight[26] = ( ( 5. / 9. ) * ( 5. / 9. ) ) * ( 5. / 9. );
                _gp[26][0] = +sqrt(3.0 / 5.0);
                _gp[26][1] = +sqrt(3.0 / 5.0);
                _gp[26][2] = +sqrt(3.0 / 5.0);

                _refGeoPtr->addGaussQuadratureRule(this);
            }

    Cn3_5_9::~Cn3_5_9()
          {}

    namespace
    {
        const Cn3_5_9& instantiateCn3_5_9 = Cn3_5_9::Instance();
    }
//---------------------------------------------------------------------------
    std::string C3_7_2::getName() const
        {
            return _name;
        }

    unsigned int  C3_7_2::order() const
        {
            return 7;
        }

    unsigned int C3_7_2::dimension() const
        {
            return 3;
        }

    unsigned int C3_7_2::nrOfPoints() const
        {
            return 34;
        }

    NumType C3_7_2::weight(unsigned int i) const
        {
            if (i < 34)
                return _weight[i];
            else
                throw _name + "::weight - wrong index!";
        }

    void C3_7_2::getPoint(unsigned int i, PointReference<3>& p) const
        {
            if (i < 34)
                p=_gp[i];
            else
                throw _name + "::getPoint -  wrong index!";
        }

    Geometry::ReferenceGeometry<3>* C3_7_2::forReferenceGeometry() const
        {
            return _refGeoPtr;
        }

    C3_7_2::C3_7_2()
     : _name("C3_7_2"),
       _refGeoPtr(&ReferenceCube::Instance())
            {
                _weight[0] = 1078. / 3645.;
                _gp[0][0] = +sqrt((6. / 7.));
                _gp[0][1] = 0.;
                _gp[0][2] = 0.;

                _weight[1] = 1078. / 3645.;
                _gp[1][0] = -sqrt((6. / 7.));
                _gp[1][1] = 0.;
                _gp[1][2] = 0.;

                _weight[2] = 1078. / 3645.;
                _gp[2][0] = 0.;
                _gp[2][1] = +sqrt((6. / 7.));
                _gp[2][2] = 0.;

                _weight[3] = 1078. / 3645.;
                _gp[3][0] = 0.;
                _gp[3][1] = -sqrt((6. / 7.));
                _gp[3][2] = 0.;

                _weight[4] = 1078. / 3645.;
                _gp[4][0] = 0.;
                _gp[4][1] = 0.;
                _gp[4][2] = +sqrt((6. / 7.));

                _weight[5] = 1078. / 3645.;
                _gp[5][0] = 0.;
                _gp[5][1] = 0.;
                _gp[5][2] = -sqrt((6. / 7.));

                _weight[6] = 343. / 3645.;
                _gp[6][0] = +sqrt((6. / 7.));
                _gp[6][1] = +sqrt((6. / 7.));
                _gp[6][2] = 0.;

                _weight[7] = 343. / 3645.;
                _gp[7][0] = -sqrt((6. / 7.));
                _gp[7][1] = +sqrt((6. / 7.));
                _gp[7][2] = 0.;

                _weight[8] = 343. / 3645.;
                _gp[8][0] = +sqrt((6. / 7.));
                _gp[8][1] = -sqrt((6. / 7.));
                _gp[8][2] = 0.;

                _weight[9] = 343. / 3645.;
                _gp[9][0] = -sqrt((6. / 7.));
                _gp[9][1] = -sqrt((6. / 7.));
                _gp[9][2] = 0.;

                _weight[10] = 343. / 3645.;
                _gp[10][0] = +sqrt((6. / 7.));
                _gp[10][1] = 0.;
                _gp[10][2] = +sqrt((6. / 7.));

                _weight[11] = 343. / 3645.;
                _gp[11][0] = -sqrt((6. / 7.));
                _gp[11][1] = 0.;
                _gp[11][2] = +sqrt((6. / 7.));

                _weight[12] = 343. / 3645.;
                _gp[12][0] = +sqrt((6. / 7.));
                _gp[12][1] = 0.;
                _gp[12][2] = -sqrt((6. / 7.));

                _weight[13] = 343. / 3645.;
                _gp[13][0] = -sqrt((6. / 7.));
                _gp[13][1] = 0.;
                _gp[13][2] = -sqrt((6. / 7.));

                _weight[14] = 343. / 3645.;
                _gp[14][0] = 0.;
                _gp[14][1] = +sqrt((6. / 7.));
                _gp[14][2] = +sqrt((6. / 7.));

                _weight[15] = 343. / 3645.;
                _gp[15][0] = 0.;
                _gp[15][1] = -sqrt((6. / 7.));
                _gp[15][2] = +sqrt((6. / 7.));

                _weight[16] = 343. / 3645.;
                _gp[16][0] = 0.;
                _gp[16][1] = +sqrt((6. / 7.));
                _gp[16][2] = -sqrt((6. / 7.));

                _weight[17] = 343. / 3645.;
                _gp[17][0] = 0.;
                _gp[17][1] = -sqrt((6. / 7.));
                _gp[17][2] = -sqrt((6. / 7.));

                _weight[18] = (774. * ((960. + 3. * sqrt(28798.)) / 2726.) - 230.) / (1215. * (((960. + 3. * sqrt(28798.)) / 2726.) - ((960. - 3. * sqrt(28798.)) / 2726.)));
                _gp[18][0] = +sqrt(((960. - 3. * sqrt(28798.)) / 2726.));
                _gp[18][1] = +sqrt(((960. - 3. * sqrt(28798.)) / 2726.));
                _gp[18][2] = +sqrt(((960. - 3. * sqrt(28798.)) / 2726.));

                _weight[19] = (774. * ((960. + 3. * sqrt(28798.)) / 2726.) - 230.) / (1215. * (((960. + 3. * sqrt(28798.)) / 2726.) - ((960. - 3. * sqrt(28798.)) / 2726.)));
                _gp[19][0] = -sqrt(((960. - 3. * sqrt(28798.)) / 2726.));
                _gp[19][1] = +sqrt(((960. - 3. * sqrt(28798.)) / 2726.));
                _gp[19][2] = +sqrt(((960. - 3. * sqrt(28798.)) / 2726.));

                _weight[20] = (774. * ((960. + 3. * sqrt(28798.)) / 2726.) - 230.) / (1215. * (((960. + 3. * sqrt(28798.)) / 2726.) - ((960. - 3. * sqrt(28798.)) / 2726.)));
                _gp[20][0] = +sqrt(((960. - 3. * sqrt(28798.)) / 2726.));
                _gp[20][1] = -sqrt(((960. - 3. * sqrt(28798.)) / 2726.));
                _gp[20][2] = +sqrt(((960. - 3. * sqrt(28798.)) / 2726.));

                _weight[21] = (774. * ((960. + 3. * sqrt(28798.)) / 2726.) - 230.) / (1215. * (((960. + 3. * sqrt(28798.)) / 2726.) - ((960. - 3. * sqrt(28798.)) / 2726.)));
                _gp[21][0] = -sqrt(((960. - 3. * sqrt(28798.)) / 2726.));
                _gp[21][1] = -sqrt(((960. - 3. * sqrt(28798.)) / 2726.));
                _gp[21][2] = +sqrt(((960. - 3. * sqrt(28798.)) / 2726.));

                _weight[22] = (774. * ((960. + 3. * sqrt(28798.)) / 2726.) - 230.) / (1215. * (((960. + 3. * sqrt(28798.)) / 2726.) - ((960. - 3. * sqrt(28798.)) / 2726.)));
                _gp[22][0] = +sqrt(((960. - 3. * sqrt(28798.)) / 2726.));
                _gp[22][1] = +sqrt(((960. - 3. * sqrt(28798.)) / 2726.));
                _gp[22][2] = -sqrt(((960. - 3. * sqrt(28798.)) / 2726.));

                _weight[23] = (774. * ((960. + 3. * sqrt(28798.)) / 2726.) - 230.) / (1215. * (((960. + 3. * sqrt(28798.)) / 2726.) - ((960. - 3. * sqrt(28798.)) / 2726.)));
                _gp[23][0] = -sqrt(((960. - 3. * sqrt(28798.)) / 2726.));
                _gp[23][1] = +sqrt(((960. - 3. * sqrt(28798.)) / 2726.));
                _gp[23][2] = -sqrt(((960. - 3. * sqrt(28798.)) / 2726.));

                _weight[24] = (774. * ((960. + 3. * sqrt(28798.)) / 2726.) - 230.) / (1215. * (((960. + 3. * sqrt(28798.)) / 2726.) - ((960. - 3. * sqrt(28798.)) / 2726.)));
                _gp[24][0] = +sqrt(((960. - 3. * sqrt(28798.)) / 2726.));
                _gp[24][1] = -sqrt(((960. - 3. * sqrt(28798.)) / 2726.));
                _gp[24][2] = -sqrt(((960. - 3. * sqrt(28798.)) / 2726.));

                _weight[25] = (774. * ((960. + 3. * sqrt(28798.)) / 2726.) - 230.) / (1215. * (((960. + 3. * sqrt(28798.)) / 2726.) - ((960. - 3. * sqrt(28798.)) / 2726.)));
                _gp[25][0] = -sqrt(((960. - 3. * sqrt(28798.)) / 2726.));
                _gp[25][1] = -sqrt(((960. - 3. * sqrt(28798.)) / 2726.));
                _gp[25][2] = -sqrt(((960. - 3. * sqrt(28798.)) / 2726.));

                _weight[26] = (230. - 774. * ((960. - 3. * sqrt(28798.)) / 2726.)) / (1215. * (((960. + 3. * sqrt(28798.)) / 2726.) - ((960. - 3. * sqrt(28798.)) / 2726.)));
                _gp[26][0] = +sqrt(((960. + 3. * sqrt(28798.)) / 2726.));
                _gp[26][1] = +sqrt(((960. + 3. * sqrt(28798.)) / 2726.));
                _gp[26][2] = +sqrt(((960. + 3. * sqrt(28798.)) / 2726.));

                _weight[27] = (230. - 774. * ((960. - 3. * sqrt(28798.)) / 2726.)) / (1215. * (((960. + 3. * sqrt(28798.)) / 2726.) - ((960. - 3. * sqrt(28798.)) / 2726.)));
                _gp[27][0] = -sqrt(((960. + 3. * sqrt(28798.)) / 2726.));
                _gp[27][1] = +sqrt(((960. + 3. * sqrt(28798.)) / 2726.));
                _gp[27][2] = +sqrt(((960. + 3. * sqrt(28798.)) / 2726.));

                _weight[28] = (230. - 774. * ((960. - 3. * sqrt(28798.)) / 2726.)) / (1215. * (((960. + 3. * sqrt(28798.)) / 2726.) - ((960. - 3. * sqrt(28798.)) / 2726.)));
                _gp[28][0] = +sqrt(((960. + 3. * sqrt(28798.)) / 2726.));
                _gp[28][1] = -sqrt(((960. + 3. * sqrt(28798.)) / 2726.));
                _gp[28][2] = +sqrt(((960. + 3. * sqrt(28798.)) / 2726.));

                _weight[29] = (230. - 774. * ((960. - 3. * sqrt(28798.)) / 2726.)) / (1215. * (((960. + 3. * sqrt(28798.)) / 2726.) - ((960. - 3. * sqrt(28798.)) / 2726.)));
                _gp[29][0] = -sqrt(((960. + 3. * sqrt(28798.)) / 2726.));
                _gp[29][1] = -sqrt(((960. + 3. * sqrt(28798.)) / 2726.));
                _gp[29][2] = +sqrt(((960. + 3. * sqrt(28798.)) / 2726.));

                _weight[30] = (230. - 774. * ((960. - 3. * sqrt(28798.)) / 2726.)) / (1215. * (((960. + 3. * sqrt(28798.)) / 2726.) - ((960. - 3. * sqrt(28798.)) / 2726.)));
                _gp[30][0] = +sqrt(((960. + 3. * sqrt(28798.)) / 2726.));
                _gp[30][1] = +sqrt(((960. + 3. * sqrt(28798.)) / 2726.));
                _gp[30][2] = -sqrt(((960. + 3. * sqrt(28798.)) / 2726.));

                _weight[31] = (230. - 774. * ((960. - 3. * sqrt(28798.)) / 2726.)) / (1215. * (((960. + 3. * sqrt(28798.)) / 2726.) - ((960. - 3. * sqrt(28798.)) / 2726.)));
                _gp[31][0] = -sqrt(((960. + 3. * sqrt(28798.)) / 2726.));
                _gp[31][1] = +sqrt(((960. + 3. * sqrt(28798.)) / 2726.));
                _gp[31][2] = -sqrt(((960. + 3. * sqrt(28798.)) / 2726.));

                _weight[32] = (230. - 774. * ((960. - 3. * sqrt(28798.)) / 2726.)) / (1215. * (((960. + 3. * sqrt(28798.)) / 2726.) - ((960. - 3. * sqrt(28798.)) / 2726.)));
                _gp[32][0] = +sqrt(((960. + 3. * sqrt(28798.)) / 2726.));
                _gp[32][1] = -sqrt(((960. + 3. * sqrt(28798.)) / 2726.));
                _gp[32][2] = -sqrt(((960. + 3. * sqrt(28798.)) / 2726.));

                _weight[33] = (230. - 774. * ((960. - 3. * sqrt(28798.)) / 2726.)) / (1215. * (((960. + 3. * sqrt(28798.)) / 2726.) - ((960. - 3. * sqrt(28798.)) / 2726.)));
                _gp[33][0] = -sqrt(((960. + 3. * sqrt(28798.)) / 2726.));
                _gp[33][1] = -sqrt(((960. + 3. * sqrt(28798.)) / 2726.));
                _gp[33][2] = -sqrt(((960. + 3. * sqrt(28798.)) / 2726.));

                _refGeoPtr->addGaussQuadratureRule(this);
            }

    C3_7_2::~C3_7_2()
          {}

    namespace
    {
        const C3_7_2& instantiateC3_7_2 = C3_7_2::Instance();
    }
//---------------------------------------------------------------------------
} // close namespace QuadratureRules

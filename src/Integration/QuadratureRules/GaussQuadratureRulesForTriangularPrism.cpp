//------------------------------------------------------------------------------
// File: GaussQuadratureRulesForTriangularPrism.cpp 
// Implementation of Gauss quadrature rules for reference triangular-prism.
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
#include "Integration/QuadratureRules/GaussQuadratureRulesForTriangularPrism.hpp"
#include "Geometry/ReferenceTriangularPrism.hpp"
using Geometry::ReferenceTriangularPrism;

//---------------------------------------------------------------------------
namespace QuadratureRules
{
//---------------------------------------------------------------------------
    std::string
    TriPrism_1_1::getName() const
    {
        return name_;
    }

    unsigned int
    TriPrism_1_1::order() const
    {
        return 1;
    }

    unsigned int
    TriPrism_1_1::dimension() const
    {
        return 3;
    }

    unsigned int
    TriPrism_1_1::nrOfPoints() const
    {
        return 1;
    }

    double
    TriPrism_1_1::weight(unsigned int i) const
    {
        if (i < 1)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }

    void
    TriPrism_1_1::getPoint(unsigned int i, PointReference<3>& p) const
    {
        if (i < 1)
            p=gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }

    TriPrism_1_1::ReferenceGeometryT*
    TriPrism_1_1::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    TriPrism_1_1::TriPrism_1_1():
        name_("TriPrism_1_1"),
        refGeoPtr_(&ReferenceTriangularPrism::Instance())
    {
        weight_[0] = ( 0.5 ) * ( 2.0 );
        gp_[0][0] = 1.0 / 3.0;
        gp_[0][1] = 1.0 / 3.0;
        gp_[0][2] = 0.0;

        refGeoPtr_->addGaussQuadratureRule(this);
    }

    TriPrism_1_1::~TriPrism_1_1()
    {
    }

    namespace
    {
        const TriPrism_1_1& instantiateTriPrism_1_1 = TriPrism_1_1::Instance();
    }
//---------------------------------------------------------------------------
    std::string 
    TriPrism_3_1::getName() const
    {
        return name_;
    }

    unsigned int
    TriPrism_3_1::order() const
    {
        return 3;
    }

    unsigned int
    TriPrism_3_1::dimension() const
    {
        return 3;
    }

    unsigned int
    TriPrism_3_1::nrOfPoints() const
    {
        return 8;
    }

    double
    TriPrism_3_1::weight(unsigned int i) const
    {
        if (i < 8)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }

    void
    TriPrism_3_1::getPoint(unsigned int i, PointReferenceT& p) const
    {
        if (i < 8)
            p=gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }

    TriPrism_3_1::ReferenceGeometryT* 
    TriPrism_3_1::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    TriPrism_3_1::TriPrism_3_1():
        name_("TriPrism_3_1"),
        refGeoPtr_(&ReferenceTriangularPrism::Instance())
    {
        weight_[0] = ( -9. / 32. ) * ( 1.0 );
        gp_[0][0] = 1. / 3.;
        gp_[0][1] = 1. / 3.;
        gp_[0][2] = -sqrt(3.0) / 3.0;

        weight_[1] = ( 25. / 96. ) * ( 1.0 );
        gp_[1][0] = 1. / 5.;
        gp_[1][1] = 1. / 5.;
        gp_[1][2] = -sqrt(3.0) / 3.0;

        weight_[2] = ( 25. / 96. ) * ( 1.0 );
        gp_[2][0] = 1. / 5.;
        gp_[2][1] = 3. / 5.;
        gp_[2][2] = -sqrt(3.0) / 3.0;

        weight_[3] = ( 25. / 96. ) * ( 1.0 );
        gp_[3][0] = 3. / 5.;
        gp_[3][1] = 1. / 5.;
        gp_[3][2] = -sqrt(3.0) / 3.0;

        weight_[4] = ( -9. / 32. ) * ( 1.0 );
        gp_[4][0] = 1. / 3.;
        gp_[4][1] = 1. / 3.;
        gp_[4][2] = +sqrt(3.0) / 3.0;

        weight_[5] = ( 25. / 96. ) * ( 1.0 );
        gp_[5][0] = 1. / 5.;
        gp_[5][1] = 1. / 5.;
        gp_[5][2] = +sqrt(3.0) / 3.0;

        weight_[6] = ( 25. / 96. ) * ( 1.0 );
        gp_[6][0] = 1. / 5.;
        gp_[6][1] = 3. / 5.;
        gp_[6][2] = +sqrt(3.0) / 3.0;

        weight_[7] = ( 25. / 96. ) * ( 1.0 );
        gp_[7][0] = 3. / 5.;
        gp_[7][1] = 1. / 5.;
        gp_[7][2] = +sqrt(3.0) / 3.0;

        refGeoPtr_->addGaussQuadratureRule(this);
    }

    TriPrism_3_1::~TriPrism_3_1()
    {
    }

    namespace
    {
        const TriPrism_3_1& instantiateTriPrism_3_1 = TriPrism_3_1::Instance();
    }
//---------------------------------------------------------------------------
    std::string
    TriPrism_5_1::getName() const
    {
        return name_;
    }

    unsigned int
    TriPrism_5_1::order() const
    {
        return 5;
    }

    unsigned int
    TriPrism_5_1::dimension() const
    {
        return 3;
    }

    unsigned int
    TriPrism_5_1::nrOfPoints() const
    {
        return 21;
    }

    double
    TriPrism_5_1::weight(unsigned int i) const
    {
        if (i < 21)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }

    void
    TriPrism_5_1::getPoint(unsigned int i, PointReferenceT& p) const
    {
        if (i < 21)
            p=gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }

    TriPrism_5_1::ReferenceGeometryT* 
    TriPrism_5_1::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    TriPrism_5_1::TriPrism_5_1():
        name_("TriPrism_5_1"),
        refGeoPtr_(&ReferenceTriangularPrism::Instance())
    {
        weight_[0] = ( 9. / 80. ) * ( 5. / 9. );
        gp_[0][0] = 1. / 3.;
        gp_[0][1] = 1. / 3.;
        gp_[0][2] = -sqrt(3.0 / 5.0);

        weight_[1] = ( (155. - sqrt(15.)) / 2400. ) * ( 5. / 9. );
        gp_[1][0] = (6. - sqrt(15.)) / 21.;
        gp_[1][1] = (6. - sqrt(15.)) / 21.;
        gp_[1][2] = -sqrt(3.0 / 5.0);

        weight_[2] = ( (155. - sqrt(15.)) / 2400. ) * ( 5. / 9. );
        gp_[2][0] = (9. + 2. * sqrt(15.)) / 21.;
        gp_[2][1] = (6. - sqrt(15.)) / 21.;
        gp_[2][2] = -sqrt(3.0 / 5.0);

        weight_[3] = ( (155. - sqrt(15.)) / 2400. ) * ( 5. / 9. );
        gp_[3][0] = (6. - sqrt(15.)) / 21.;
        gp_[3][1] = (9. + 2. * sqrt(15.)) / 21.;
        gp_[3][2] = -sqrt(3.0 / 5.0);

        weight_[4] = ( (155. + sqrt(15.)) / 2400. ) * ( 5. / 9. );
        gp_[4][0] = (6. + sqrt(15.)) / 21.;
        gp_[4][1] = (6. + sqrt(15.)) / 21.;
        gp_[4][2] = -sqrt(3.0 / 5.0);

        weight_[5] = ( (155. + sqrt(15.)) / 2400. ) * ( 5. / 9. );
        gp_[5][0] = (9. - 2. * sqrt(15.)) / 21.;
        gp_[5][1] = (6. + sqrt(15.)) / 21.;
        gp_[5][2] = -sqrt(3.0 / 5.0);

        weight_[6] = ( (155. + sqrt(15.)) / 2400. ) * ( 5. / 9. );
        gp_[6][0] = (6. + sqrt(15.)) / 21.;
        gp_[6][1] = (9. - 2. * sqrt(15.)) / 21.;
        gp_[6][2] = -sqrt(3.0 / 5.0);

        weight_[7] = ( 9. / 80. ) * ( 8. / 9. );
        gp_[7][0] = 1. / 3.;
        gp_[7][1] = 1. / 3.;
        gp_[7][2] = 0.0;

        weight_[8] = ( (155. - sqrt(15.)) / 2400. ) * ( 8. / 9. );
        gp_[8][0] = (6. - sqrt(15.)) / 21.;
        gp_[8][1] = (6. - sqrt(15.)) / 21.;
        gp_[8][2] = 0.0;

        weight_[9] = ( (155. - sqrt(15.)) / 2400. ) * ( 8. / 9. );
        gp_[9][0] = (9. + 2. * sqrt(15.)) / 21.;
        gp_[9][1] = (6. - sqrt(15.)) / 21.;
        gp_[9][2] = 0.0;

        weight_[10] = ( (155. - sqrt(15.)) / 2400. ) * ( 8. / 9. );
        gp_[10][0] = (6. - sqrt(15.)) / 21.;
        gp_[10][1] = (9. + 2. * sqrt(15.)) / 21.;
        gp_[10][2] = 0.0;

        weight_[11] = ( (155. + sqrt(15.)) / 2400. ) * ( 8. / 9. );
        gp_[11][0] = (6. + sqrt(15.)) / 21.;
        gp_[11][1] = (6. + sqrt(15.)) / 21.;
        gp_[11][2] = 0.0;

        weight_[12] = ( (155. + sqrt(15.)) / 2400. ) * ( 8. / 9. );
        gp_[12][0] = (9. - 2. * sqrt(15.)) / 21.;
        gp_[12][1] = (6. + sqrt(15.)) / 21.;
        gp_[12][2] = 0.0;

        weight_[13] = ( (155. + sqrt(15.)) / 2400. ) * ( 8. / 9. );
        gp_[13][0] = (6. + sqrt(15.)) / 21.;
        gp_[13][1] = (9. - 2. * sqrt(15.)) / 21.;
        gp_[13][2] = 0.0;

        weight_[14] = ( 9. / 80. ) * ( 5. / 9. );
        gp_[14][0] = 1. / 3.;
        gp_[14][1] = 1. / 3.;
        gp_[14][2] = +sqrt(3.0 / 5.0);

        weight_[15] = ( (155. - sqrt(15.)) / 2400. ) * ( 5. / 9. );
        gp_[15][0] = (6. - sqrt(15.)) / 21.;
        gp_[15][1] = (6. - sqrt(15.)) / 21.;
        gp_[15][2] = +sqrt(3.0 / 5.0);

        weight_[16] = ( (155. - sqrt(15.)) / 2400. ) * ( 5. / 9. );
        gp_[16][0] = (9. + 2. * sqrt(15.)) / 21.;
        gp_[16][1] = (6. - sqrt(15.)) / 21.;
        gp_[16][2] = +sqrt(3.0 / 5.0);

        weight_[17] = ( (155. - sqrt(15.)) / 2400. ) * ( 5. / 9. );
        gp_[17][0] = (6. - sqrt(15.)) / 21.;
        gp_[17][1] = (9. + 2. * sqrt(15.)) / 21.;
        gp_[17][2] = +sqrt(3.0 / 5.0);

        weight_[18] = ( (155. + sqrt(15.)) / 2400. ) * ( 5. / 9. );
        gp_[18][0] = (6. + sqrt(15.)) / 21.;
        gp_[18][1] = (6. + sqrt(15.)) / 21.;
        gp_[18][2] = +sqrt(3.0 / 5.0);

        weight_[19] = ( (155. + sqrt(15.)) / 2400. ) * ( 5. / 9. );
        gp_[19][0] = (9. - 2. * sqrt(15.)) / 21.;
        gp_[19][1] = (6. + sqrt(15.)) / 21.;
        gp_[19][2] = +sqrt(3.0 / 5.0);

        weight_[20] = ( (155. + sqrt(15.)) / 2400. ) * ( 5. / 9. );
        gp_[20][0] = (6. + sqrt(15.)) / 21.;
        gp_[20][1] = (9. - 2. * sqrt(15.)) / 21.;
        gp_[20][2] = +sqrt(3.0 / 5.0);

        refGeoPtr_->addGaussQuadratureRule(this);
    }

    TriPrism_5_1::~TriPrism_5_1()
    {}

    namespace
    {
        const TriPrism_5_1& instantiateTriPrism_5_1 = TriPrism_5_1::Instance();
    }
//---------------------------------------------------------------------------
    std::string
    TriPrism_7_1::getName() const
    {
        return name_;
    }

    unsigned int
    TriPrism_7_1::order() const
    {
        return 7;
    }

    unsigned int 
    TriPrism_7_1::dimension() const
    {
        return 3;
    }

    unsigned int
    TriPrism_7_1::nrOfPoints() const
    {
        return 64;
    }

    double
    TriPrism_7_1::weight(unsigned int i) const
    {
        if (i < 64)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }

    void
    TriPrism_7_1::getPoint(unsigned int i, PointReference<3>& p) const
    {
        if (i < 64)
            p=gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }

    TriPrism_7_1::ReferenceGeometryT*
    TriPrism_7_1::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    TriPrism_7_1::TriPrism_7_1():
        name_("TriPrism_7_1"),
        refGeoPtr_(&ReferenceTriangularPrism::Instance())
    {
        weight_[0] = ( (0.1739274226) * (0.1355069134) ) * ( (2. * 0.1739274226) );
        gp_[0][0] = ((-0.861136312 + 1.) / 2.) * (1. - (0.0571041961));
        gp_[0][1] = 0.0571041961;
        gp_[0][2] = (-0.861136312);

        weight_[1] = ( (0.3260725774) * (0.1355069134) ) * ( (2. * 0.1739274226) );
        gp_[1][0] = ((-0.3399810436 + 1.)/ 2.) * (1. - (0.0571041961));
        gp_[1][1] = 0.0571041961;
        gp_[1][2] = (-0.861136312);

        weight_[2] = ( (0.3260725774) * (0.1355069134) ) * ( (2. * 0.1739274226) );
        gp_[2][0] = ((+0.3399810436 + 1.)/ 2.) * (1. - (0.0571041961));
        gp_[2][1] = 0.0571041961;
        gp_[2][2] = (-0.861136312);

        weight_[3] = ( (0.1739274226) * (0.1355069134) ) * ( (2. * 0.1739274226) );
        gp_[3][0] = ((+0.861136312 + 1.) / 2.) * (1. - (0.0571041961));
        gp_[3][1] = 0.0571041961;
        gp_[3][2] = (-0.861136312);

        weight_[4] = ( (0.1739274226) * (0.2034645680) ) * ( (2. * 0.1739274226) );
        gp_[4][0] = ((-0.861136312 + 1.) / 2.) * (1. - (0.2768430136));
        gp_[4][1] = 0.2768430136;
        gp_[4][2] = (-0.861136312);

        weight_[5] = ( (0.3260725774) * (0.2034645680) ) * ( (2. * 0.1739274226) );
        gp_[5][0] = ((-0.3399810436 + 1.)/ 2.) * (1. - (0.2768430136));
        gp_[5][1] = 0.2768430136;
        gp_[5][2] = (-0.861136312);

        weight_[6] = ( (0.3260725774) * (0.2034645680) ) * ( (2. * 0.1739274226) );
        gp_[6][0] = ((+0.3399810436 + 1.)/ 2.) * (1. - (0.2768430136));
        gp_[6][1] = 0.2768430136;
        gp_[6][2] = (-0.861136312);

        weight_[7] = ( (0.1739274226) * (0.2034645680) ) * ( (2. * 0.1739274226) );
        gp_[7][0] = ((+0.861136312 + 1.) / 2.) * (1. - (0.2768430136));
        gp_[7][1] = 0.2768430136;
        gp_[7][2] = (-0.861136312);

        weight_[8] = ( (0.1739274226) * (0.1298475476) ) * ( (2. * 0.1739274226) );
        gp_[8][0] = ((-0.861136312 + 1.) / 2.) * (1. - (0.5835904324));
        gp_[8][1] = 0.5835904324;
        gp_[8][2] = (-0.861136312);

        weight_[9] = ( (0.3260725774) * (0.1298475476) ) * ( (2. * 0.1739274226) );
        gp_[9][0] = ((-0.3399810436 + 1.)/ 2.) * (1. - (0.5835904324));
        gp_[9][1] = 0.5835904324;
        gp_[9][2] = (-0.861136312);

        weight_[10] = ( (0.3260725774) * (0.1298475476) ) * ( (2. * 0.1739274226) );
        gp_[10][0] = ((+0.3399810436 + 1.)/ 2.) * (1. - (0.5835904324));
        gp_[10][1] = 0.5835904324;
        gp_[10][2] = (-0.861136312);

        weight_[11] = ( (0.1739274226) * (0.1298475476) ) * ( (2. * 0.1739274226) );
        gp_[11][0] = ((+0.861136312 + 1.) / 2.) * (1. - (0.5835904324));
        gp_[11][1] = 0.5835904324;
        gp_[11][2] = (-0.861136312);

        weight_[12] = ( (0.1739274226) * (0.0311809709) ) * ( (2. * 0.1739274226) );
        gp_[12][0] = ((-0.861136312 + 1.) / 2.) * (1. - (0.8602401357));
        gp_[12][1] = 0.8602401357;
        gp_[12][2] = (-0.861136312);

        weight_[13] = ( (0.3260725774) * (0.0311809709) ) * ( (2. * 0.1739274226) );
        gp_[13][0] = ((-0.3399810436 + 1.)/ 2.) * (1. - (0.8602401357));
        gp_[13][1] = 0.8602401357;
        gp_[13][2] = (-0.861136312);

        weight_[14] = ( (0.3260725774) * (0.0311809709) ) * ( (2. * 0.1739274226) );
        gp_[14][0] = ((+0.3399810436 + 1.)/ 2.) * (1. - (0.8602401357));
        gp_[14][1] = 0.8602401357;
        gp_[14][2] = (-0.861136312);

        weight_[15] = ( (0.1739274226) * (0.0311809709) ) * ( (2. * 0.1739274226) );
        gp_[15][0] = ((+0.861136312 + 1.) / 2.) * (1. - (0.8602401357));
        gp_[15][1] = 0.8602401357;
        gp_[15][2] = (-0.861136312);

        weight_[16] = ( (0.1739274226) * (0.1355069134) ) * ( (2. * 0.3260725774) );
        gp_[16][0] = ((-0.861136312 + 1.) / 2.) * (1. - (0.0571041961));
        gp_[16][1] = 0.0571041961;
        gp_[16][2] = (-0.3399810436);

        weight_[17] = ( (0.3260725774) * (0.1355069134) ) * ( (2. * 0.3260725774) );
        gp_[17][0] = ((-0.3399810436 + 1.)/ 2.) * (1. - (0.0571041961));
        gp_[17][1] = 0.0571041961;
        gp_[17][2] = (-0.3399810436);

        weight_[18] = ( (0.3260725774) * (0.1355069134) ) * ( (2. * 0.3260725774) );
        gp_[18][0] = ((+0.3399810436 + 1.)/ 2.) * (1. - (0.0571041961));
        gp_[18][1] = 0.0571041961;
        gp_[18][2] = (-0.3399810436);

        weight_[19] = ( (0.1739274226) * (0.1355069134) ) * ( (2. * 0.3260725774) );
        gp_[19][0] = ((+0.861136312 + 1.) / 2.) * (1. - (0.0571041961));
        gp_[19][1] = 0.0571041961;
        gp_[19][2] = (-0.3399810436);

        weight_[20] = ( (0.1739274226) * (0.2034645680) ) * ( (2. * 0.3260725774) );
        gp_[20][0] = ((-0.861136312 + 1.) / 2.) * (1. - (0.2768430136));
        gp_[20][1] = 0.2768430136;
        gp_[20][2] = (-0.3399810436);

        weight_[21] = ( (0.3260725774) * (0.2034645680) ) * ( (2. * 0.3260725774) );
        gp_[21][0] = ((-0.3399810436 + 1.)/ 2.) * (1. - (0.2768430136));
        gp_[21][1] = 0.2768430136;
        gp_[21][2] = (-0.3399810436);

        weight_[22] = ( (0.3260725774) * (0.2034645680) ) * ( (2. * 0.3260725774) );
        gp_[22][0] = ((+0.3399810436 + 1.)/ 2.) * (1. - (0.2768430136));
        gp_[22][1] = 0.2768430136;
        gp_[22][2] = (-0.3399810436);

        weight_[23] = ( (0.1739274226) * (0.2034645680) ) * ( (2. * 0.3260725774) );
        gp_[23][0] = ((+0.861136312 + 1.) / 2.) * (1. - (0.2768430136));
        gp_[23][1] = 0.2768430136;
        gp_[23][2] = (-0.3399810436);

        weight_[24] = ( (0.1739274226) * (0.1298475476) ) * ( (2. * 0.3260725774) );
        gp_[24][0] = ((-0.861136312 + 1.) / 2.) * (1. - (0.5835904324));
        gp_[24][1] = 0.5835904324;
        gp_[24][2] = (-0.3399810436);

        weight_[25] = ( (0.3260725774) * (0.1298475476) ) * ( (2. * 0.3260725774) );
        gp_[25][0] = ((-0.3399810436 + 1.)/ 2.) * (1. - (0.5835904324));
        gp_[25][1] = 0.5835904324;
        gp_[25][2] = (-0.3399810436);

        weight_[26] = ( (0.3260725774) * (0.1298475476) ) * ( (2. * 0.3260725774) );
        gp_[26][0] = ((+0.3399810436 + 1.)/ 2.) * (1. - (0.5835904324));
        gp_[26][1] = 0.5835904324;
        gp_[26][2] = (-0.3399810436);

        weight_[27] = ( (0.1739274226) * (0.1298475476) ) * ( (2. * 0.3260725774) );
        gp_[27][0] = ((+0.861136312 + 1.) / 2.) * (1. - (0.5835904324));
        gp_[27][1] = 0.5835904324;
        gp_[27][2] = (-0.3399810436);

        weight_[28] = ( (0.1739274226) * (0.0311809709) ) * ( (2. * 0.3260725774) );
        gp_[28][0] = ((-0.861136312 + 1.) / 2.) * (1. - (0.8602401357));
        gp_[28][1] = 0.8602401357;
        gp_[28][2] = (-0.3399810436);

        weight_[29] = ( (0.3260725774) * (0.0311809709) ) * ( (2. * 0.3260725774) );
        gp_[29][0] = ((-0.3399810436 + 1.)/ 2.) * (1. - (0.8602401357));
        gp_[29][1] = 0.8602401357;
        gp_[29][2] = (-0.3399810436);

        weight_[30] = ( (0.3260725774) * (0.0311809709) ) * ( (2. * 0.3260725774) );
        gp_[30][0] = ((+0.3399810436 + 1.)/ 2.) * (1. - (0.8602401357));
        gp_[30][1] = 0.8602401357;
        gp_[30][2] = (-0.3399810436);

        weight_[31] = ( (0.1739274226) * (0.0311809709) ) * ( (2. * 0.3260725774) );
        gp_[31][0] = ((+0.861136312 + 1.) / 2.) * (1. - (0.8602401357));
        gp_[31][1] = 0.8602401357;
        gp_[31][2] = (-0.3399810436);

        weight_[32] = ( (0.1739274226) * (0.1355069134) ) * ( (2. * 0.3260725774) );
        gp_[32][0] = ((-0.861136312 + 1.) / 2.) * (1. - (0.0571041961));
        gp_[32][1] = 0.0571041961;
        gp_[32][2] = (+0.3399810436);

        weight_[33] = ( (0.3260725774) * (0.1355069134) ) * ( (2. * 0.3260725774) );
        gp_[33][0] = ((-0.3399810436 + 1.)/ 2.) * (1. - (0.0571041961));
        gp_[33][1] = 0.0571041961;
        gp_[33][2] = (+0.3399810436);

        weight_[34] = ( (0.3260725774) * (0.1355069134) ) * ( (2. * 0.3260725774) );
        gp_[34][0] = ((+0.3399810436 + 1.)/ 2.) * (1. - (0.0571041961));
        gp_[34][1] = 0.0571041961;
        gp_[34][2] = (+0.3399810436);

        weight_[35] = ( (0.1739274226) * (0.1355069134) ) * ( (2. * 0.3260725774) );
        gp_[35][0] = ((+0.861136312 + 1.) / 2.) * (1. - (0.0571041961));
        gp_[35][1] = 0.0571041961;
        gp_[35][2] = (+0.3399810436);

        weight_[36] = ( (0.1739274226) * (0.2034645680) ) * ( (2. * 0.3260725774) );
        gp_[36][0] = ((-0.861136312 + 1.) / 2.) * (1. - (0.2768430136));
        gp_[36][1] = 0.2768430136;
        gp_[36][2] = (+0.3399810436);

        weight_[37] = ( (0.3260725774) * (0.2034645680) ) * ( (2. * 0.3260725774) );
        gp_[37][0] = ((-0.3399810436 + 1.)/ 2.) * (1. - (0.2768430136));
        gp_[37][1] = 0.2768430136;
        gp_[37][2] = (+0.3399810436);

        weight_[38] = ( (0.3260725774) * (0.2034645680) ) * ( (2. * 0.3260725774) );
        gp_[38][0] = ((+0.3399810436 + 1.)/ 2.) * (1. - (0.2768430136));
        gp_[38][1] = 0.2768430136;
        gp_[38][2] = (+0.3399810436);

        weight_[39] = ( (0.1739274226) * (0.2034645680) ) * ( (2. * 0.3260725774) );
        gp_[39][0] = ((+0.861136312 + 1.) / 2.) * (1. - (0.2768430136));
        gp_[39][1] = 0.2768430136;
        gp_[39][2] = (+0.3399810436);

        weight_[40] = ( (0.1739274226) * (0.1298475476) ) * ( (2. * 0.3260725774) );
        gp_[40][0] = ((-0.861136312 + 1.) / 2.) * (1. - (0.5835904324));
        gp_[40][1] = 0.5835904324;
        gp_[40][2] = (+0.3399810436);

        weight_[41] = ( (0.3260725774) * (0.1298475476) ) * ( (2. * 0.3260725774) );
        gp_[41][0] = ((-0.3399810436 + 1.)/ 2.) * (1. - (0.5835904324));
        gp_[41][1] = 0.5835904324;
        gp_[41][2] = (+0.3399810436);

        weight_[42] = ( (0.3260725774) * (0.1298475476) ) * ( (2. * 0.3260725774) );
        gp_[42][0] = ((+0.3399810436 + 1.)/ 2.) * (1. - (0.5835904324));
        gp_[42][1] = 0.5835904324;
        gp_[42][2] = (+0.3399810436);

        weight_[43] = ( (0.1739274226) * (0.1298475476) ) * ( (2. * 0.3260725774) );
        gp_[43][0] = ((+0.861136312 + 1.) / 2.) * (1. - (0.5835904324));
        gp_[43][1] = 0.5835904324;
        gp_[43][2] = (+0.3399810436);

        weight_[44] = ( (0.1739274226) * (0.0311809709) ) * ( (2. * 0.3260725774) );
        gp_[44][0] = ((-0.861136312 + 1.) / 2.) * (1. - (0.8602401357));
        gp_[44][1] = 0.8602401357;
        gp_[44][2] = (+0.3399810436);

        weight_[45] = ( (0.3260725774) * (0.0311809709) ) * ( (2. * 0.3260725774) );
        gp_[45][0] = ((-0.3399810436 + 1.)/ 2.) * (1. - (0.8602401357));
        gp_[45][1] = 0.8602401357;
        gp_[45][2] = (+0.3399810436);

        weight_[46] = ( (0.3260725774) * (0.0311809709) ) * ( (2. * 0.3260725774) );
        gp_[46][0] = ((+0.3399810436 + 1.)/ 2.) * (1. - (0.8602401357));
        gp_[46][1] = 0.8602401357;
        gp_[46][2] = (+0.3399810436);

        weight_[47] = ( (0.1739274226) * (0.0311809709) ) * ( (2. * 0.3260725774) );
        gp_[47][0] = ((+0.861136312 + 1.) / 2.) * (1. - (0.8602401357));
        gp_[47][1] = 0.8602401357;
        gp_[47][2] = (+0.3399810436);

        weight_[48] = ( (0.1739274226) * (0.1355069134) ) * ( (2. * 0.1739274226) );
        gp_[48][0] = ((-0.861136312 + 1.) / 2.) * (1. - (0.0571041961));
        gp_[48][1] = 0.0571041961;
        gp_[48][2] = (+0.861136312);

        weight_[49] = ( (0.3260725774) * (0.1355069134) ) * ( (2. * 0.1739274226) );
        gp_[49][0] = ((-0.3399810436 + 1.)/ 2.) * (1. - (0.0571041961));
        gp_[49][1] = 0.0571041961;
        gp_[49][2] = (+0.861136312);

        weight_[50] = ( (0.3260725774) * (0.1355069134) ) * ( (2. * 0.1739274226) );
        gp_[50][0] = ((+0.3399810436 + 1.)/ 2.) * (1. - (0.0571041961));
        gp_[50][1] = 0.0571041961;
        gp_[50][2] = (+0.861136312);

        weight_[51] = ( (0.1739274226) * (0.1355069134) ) * ( (2. * 0.1739274226) );
        gp_[51][0] = ((+0.861136312 + 1.) / 2.) * (1. - (0.0571041961));
        gp_[51][1] = 0.0571041961;
        gp_[51][2] = (+0.861136312);

        weight_[52] = ( (0.1739274226) * (0.2034645680) ) * ( (2. * 0.1739274226) );
        gp_[52][0] = ((-0.861136312 + 1.) / 2.) * (1. - (0.2768430136));
        gp_[52][1] = 0.2768430136;
        gp_[52][2] = (+0.861136312);

        weight_[53] = ( (0.3260725774) * (0.2034645680) ) * ( (2. * 0.1739274226) );
        gp_[53][0] = ((-0.3399810436 + 1.)/ 2.) * (1. - (0.2768430136));
        gp_[53][1] = 0.2768430136;
        gp_[53][2] = (+0.861136312);

        weight_[54] = ( (0.3260725774) * (0.2034645680) ) * ( (2. * 0.1739274226) );
        gp_[54][0] = ((+0.3399810436 + 1.)/ 2.) * (1. - (0.2768430136));
        gp_[54][1] = 0.2768430136;
        gp_[54][2] = (+0.861136312);

        weight_[55] = ( (0.1739274226) * (0.2034645680) ) * ( (2. * 0.1739274226) );
        gp_[55][0] = ((+0.861136312 + 1.) / 2.) * (1. - (0.2768430136));
        gp_[55][1] = 0.2768430136;
        gp_[55][2] = (+0.861136312);

        weight_[56] = ( (0.1739274226) * (0.1298475476) ) * ( (2. * 0.1739274226) );
        gp_[56][0] = ((-0.861136312 + 1.) / 2.) * (1. - (0.5835904324));
        gp_[56][1] = 0.5835904324;
        gp_[56][2] = (+0.861136312);

        weight_[57] = ( (0.3260725774) * (0.1298475476) ) * ( (2. * 0.1739274226) );
        gp_[57][0] = ((-0.3399810436 + 1.)/ 2.) * (1. - (0.5835904324));
        gp_[57][1] = 0.5835904324;
        gp_[57][2] = (+0.861136312);

        weight_[58] = ( (0.3260725774) * (0.1298475476) ) * ( (2. * 0.1739274226) );
        gp_[58][0] = ((+0.3399810436 + 1.)/ 2.) * (1. - (0.5835904324));
        gp_[58][1] = 0.5835904324;
        gp_[58][2] = (+0.861136312);

        weight_[59] = ( (0.1739274226) * (0.1298475476) ) * ( (2. * 0.1739274226) );
        gp_[59][0] = ((+0.861136312 + 1.) / 2.) * (1. - (0.5835904324));
        gp_[59][1] = 0.5835904324;
        gp_[59][2] = (+0.861136312);

        weight_[60] = ( (0.1739274226) * (0.0311809709) ) * ( (2. * 0.1739274226) );
        gp_[60][0] = ((-0.861136312 + 1.) / 2.) * (1. - (0.8602401357));
        gp_[60][1] = 0.8602401357;
        gp_[60][2] = (+0.861136312);

        weight_[61] = ( (0.3260725774) * (0.0311809709) ) * ( (2. * 0.1739274226) );
        gp_[61][0] = ((-0.3399810436 + 1.)/ 2.) * (1. - (0.8602401357));
        gp_[61][1] = 0.8602401357;
        gp_[61][2] = (+0.861136312);

        weight_[62] = ( (0.3260725774) * (0.0311809709) ) * ( (2. * 0.1739274226) );
        gp_[62][0] = ((+0.3399810436 + 1.)/ 2.) * (1. - (0.8602401357));
        gp_[62][1] = 0.8602401357;
        gp_[62][2] = (+0.861136312);

        weight_[63] = ( (0.1739274226) * (0.0311809709) ) * ( (2. * 0.1739274226) );
        gp_[63][0] = ((+0.861136312 + 1.) / 2.) * (1. - (0.8602401357));
        gp_[63][1] = 0.8602401357;
        gp_[63][2] = (+0.861136312);

        refGeoPtr_->addGaussQuadratureRule(this);
    }

    TriPrism_7_1::~TriPrism_7_1()
    {
    }

    namespace
    {
        const TriPrism_7_1& instantiateTriPrism_7_1 = TriPrism_7_1::Instance();
    }
//---------------------------------------------------------------------------
} // close namespace IntegrationRules

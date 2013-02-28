//------------------------------------------------------------------------------
// File: GaussQuadratureRulesForTriangle.cpp 
// Implementation of Gauss quadrature rules for reference triangle.
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
#include "Integration/QuadratureRules/GaussQuadratureRulesForTriangle.hpp"

//---------------------------------------------------------------------------
namespace QuadratureRules
{
//---------------------------------------------------------------------------
    std::string Tn2_1_1::getName() const
        {
            return _name;
        }

    DimType  Tn2_1_1::order() const
        {
            return 1;
        }

    DimType Tn2_1_1::dimension() const
        {
            return 2;
        }

    DimType Tn2_1_1::nrOfPoints() const
        {
            return 1;
        }

    NumType Tn2_1_1::weight(DimType i) const
        {
            if (i < 1)
                return _weight[i];
            else
                throw _name + "::weight - wrong index!";
        }

    void Tn2_1_1::getPoint(DimType i, PointReference<2>& p) const
        {
            if (i < 1)
                p=_gp[i];
            else
                throw _name + "::getPoint -  wrong index!";
        }

    Geometry::ReferenceGeometry<2>* Tn2_1_1::forReferenceGeometry() const
        {
            return _refGeoPtr;
        }

    Tn2_1_1::Tn2_1_1()
     : _name("Tn2_1_1"),
       _refGeoPtr(&ReferenceTriangle::Instance())
            {
                _weight[0] = 5.0000000000000000e-01;
                _gp[0][0] = 3.3333333333333348e-01;
                _gp[0][1] = 3.3333333333333348e-01;

                _refGeoPtr->addGaussQuadratureRule(this);
            }

    Tn2_1_1::~Tn2_1_1()
          {}

    namespace
    {
        const Tn2_1_1& instantiateTn2_1_1 = Tn2_1_1::Instance();
    }
//---------------------------------------------------------------------------
    std::string Tn2_2_1::getName() const
        {
            return _name;
        }

    DimType  Tn2_2_1::order() const
        {
            return 2;
        }

    DimType Tn2_2_1::dimension() const
        {
            return 2;
        }

    DimType Tn2_2_1::nrOfPoints() const
        {
            return 3;
        }

    NumType Tn2_2_1::weight(DimType i) const
        {
            if (i < 3)
                return _weight[i];
            else
                throw _name + "::weight - wrong index!";
        }

    void Tn2_2_1::getPoint(DimType i, PointReference<2>& p) const
        {
            if (i < 3)
                p=_gp[i];
            else
                throw _name + "::getPoint -  wrong index!";
        }

    Geometry::ReferenceGeometry<2>* Tn2_2_1::forReferenceGeometry() const
        {
            return _refGeoPtr;
        }

    Tn2_2_1::Tn2_2_1()
     : _name("Tn2_2_1"),
       _refGeoPtr(&ReferenceTriangle::Instance())
            {
                _weight[0] = 1.6666666666666674e-01;
                _gp[0][0] = 1.6666666666666652e-01;
                _gp[0][1] = 1.6666666666666652e-01;
                
                _weight[1] = 1.6666666666666674e-01;
                _gp[1][0] = 1.6666666666666652e-01;
                _gp[1][1] = 6.6666666666666652e-01;
                
                _weight[2] = 1.6666666666666674e-01;
                _gp[2][0] = 6.6666666666666652e-01;
                _gp[2][1] = 1.6666666666666652e-01;

                _refGeoPtr->addGaussQuadratureRule(this);
            }

    Tn2_2_1::~Tn2_2_1()
          {}

    namespace
    {
        const Tn2_2_1& instantiateTn2_2_1 = Tn2_2_1::Instance();
    }
//---------------------------------------------------------------------------
    std::string Tn2_3_1::getName() const
        {
            return _name;
        }

    DimType  Tn2_3_1::order() const
        {
            return 3;
        }

    DimType Tn2_3_1::dimension() const
        {
            return 2;
        }

    DimType Tn2_3_1::nrOfPoints() const
        {
            return 4;
        }

    NumType Tn2_3_1::weight(DimType i) const
        {
            if (i < 4)
                return _weight[i];
            else
                throw _name + "::weight - wrong index!";
        }

    void Tn2_3_1::getPoint(DimType i, PointReference<2>& p) const
        {
            if (i < 4)
                p=_gp[i];
            else
                throw _name + "::getPoint -  wrong index!";
        }

    Geometry::ReferenceGeometry<2>* Tn2_3_1::forReferenceGeometry() const
        {
            return _refGeoPtr;
        }

    Tn2_3_1::Tn2_3_1()
     : _name("Tn2_3_1"),
       _refGeoPtr(&ReferenceTriangle::Instance())
            {
                _weight[0] = -2.8125000000000000e-01;
                _gp[0][0] = 3.3333333333333348e-01;
                _gp[0][1] = 3.3333333333333348e-01;

                _weight[1] = 2.6041666666666674e-01;
                _gp[1][0] = 2.0000000000000001e-01;
                _gp[1][1] = 2.0000000000000001e-01;
                
                _weight[2] = 2.6041666666666674e-01;
                _gp[2][0] = 2.0000000000000001e-01;
                _gp[2][1] = 5.9999999999999998e-01;

                _weight[3] = 2.6041666666666674e-01;
                _gp[3][0] = 5.9999999999999998e-01;
                _gp[3][1] = 2.0000000000000001e-01;

                _refGeoPtr->addGaussQuadratureRule(this);
            }

    Tn2_3_1::~Tn2_3_1()
          {}

    namespace
    {
        const Tn2_3_1& instantiateTn2_3_1 = Tn2_3_1::Instance();
    }
//---------------------------------------------------------------------------
    std::string Tn2_4_1::getName() const
        {
            return _name;
        }

    DimType  Tn2_4_1::order() const
        {
            return 4;
        }

    DimType Tn2_4_1::dimension() const
        {
            return 2;
        }

    DimType Tn2_4_1::nrOfPoints() const
        {
            return 6;
        }

    NumType Tn2_4_1::weight(DimType i) const
        {
            if (i < 6)
                return _weight[i];
            else
                throw _name + "::weight - wrong index!";
        }

    void Tn2_4_1::getPoint(DimType i, PointReference<2>& p) const
        {
            if (i < 6)
                p=_gp[i];
            else
                throw _name + "::getPoint -  wrong index!";
        }

    Geometry::ReferenceGeometry<2>* Tn2_4_1::forReferenceGeometry() const
        {
            return _refGeoPtr;
        }

    Tn2_4_1::Tn2_4_1()
     : _name("Tn2_4_1"),
       _refGeoPtr(&ReferenceTriangle::Instance())
            {
                _weight[0] = 1.1169079483900550e-01;
                _gp[0][0] = 4.4594849091596500e-01;
                _gp[0][1] = 4.4594849091596500e-01;

                _weight[1] = 1.1169079483900550e-01;
                _gp[1][0] = 4.4594849091596500e-01;
                _gp[1][1] = 1.0810301816807000e-01;

                _weight[2] = 1.1169079483900550e-01;
                _gp[2][0] = 1.0810301816807000e-01;
                _gp[2][1] = 4.4594849091596500e-01;

                _weight[3] = 5.4975871827660998e-02;
                _gp[3][0] = 9.1576213509771021e-02;
                _gp[3][1] = 9.1576213509771021e-02;

                _weight[4] = 5.4975871827660998e-02;
                _gp[4][0] = 9.1576213509771021e-02;
                _gp[4][1] = 8.1684757298045896e-01;

                _weight[5] = 5.4975871827660998e-02;
                _gp[5][0] = 8.1684757298045896e-01;
                _gp[5][1] = 9.1576213509771021e-02;

                _refGeoPtr->addGaussQuadratureRule(this);
            }

    Tn2_4_1::~Tn2_4_1()
          {}

    namespace
    {
        const Tn2_4_1& instantiateTn2_4_1 = Tn2_4_1::Instance();
    }
//---------------------------------------------------------------------------
    std::string T2_5_1::getName() const
        {
            return _name;
        }

    DimType  T2_5_1::order() const
        {
            return 5;
        }

    DimType T2_5_1::dimension() const
        {
            return 2;
        }

    DimType T2_5_1::nrOfPoints() const
        {
            return 7;
        }

    NumType T2_5_1::weight(DimType i) const
        {
            if (i < 7)
                return _weight[i];
            else
                throw _name + "::weight - wrong index!";
        }

    void T2_5_1::getPoint(DimType i, PointReference<2>& p) const
        {
            if (i < 7)
                p=_gp[i];
            else
                throw _name + "::getPoint -  wrong index!";
        }

    Geometry::ReferenceGeometry<2>* T2_5_1::forReferenceGeometry() const
        {
            return _refGeoPtr;
        }

    T2_5_1::T2_5_1()
     : _name("T2_5_1"),
       _refGeoPtr(&ReferenceTriangle::Instance())
            {
                _weight[0] = 1.1250000000000000e-01;
                _gp[0][0] = 3.3333333333333348e-01;
                _gp[0][1] = 3.3333333333333348e-01;

                _weight[1] = 6.6197076394252999e-02;
                _gp[1][0] = 4.7014206410511500e-01;
                _gp[1][1] = 4.7014206410511500e-01;
                
                _weight[2] = 6.6197076394252999e-02;
                _gp[2][0] = 4.7014206410511500e-01;
                _gp[2][1] = 5.9715871789770003e-02;

                _weight[3] = 6.6197076394252999e-02;
                _gp[3][0] = 5.9715871789770003e-02;
                _gp[3][1] = 4.7014206410511500e-01;

                _weight[4] = 6.2969590272413500e-02;
                _gp[4][0] = 1.0128650732345601e-01;
                _gp[4][1] = 1.0128650732345601e-01;

                _weight[5] = 6.2969590272413500e-02;
                _gp[5][0] = 1.0128650732345601e-01;
                _gp[5][1] = 7.9742698535308698e-01;

                _weight[6] = 6.2969590272413500e-02;
                _gp[6][0] = 7.9742698535308698e-01;
                _gp[6][1] = 1.0128650732345601e-01;

                _refGeoPtr->addGaussQuadratureRule(this);
            }

    T2_5_1::~T2_5_1()
          {}

    namespace
    {
        const T2_5_1& instantiateT2_5_1 = T2_5_1::Instance();
    }
//---------------------------------------------------------------------------
    std::string T2_6_1::getName() const
        {
            return _name;
        }

    DimType  T2_6_1::order() const
        {
            return 6;
        }

    DimType T2_6_1::dimension() const
        {
            return 2;
        }

    DimType T2_6_1::nrOfPoints() const
        {
            return 12;
        }

    NumType T2_6_1::weight(DimType i) const
        {
            if (i < 12)
                return _weight[i];
            else
                throw _name + "::weight - wrong index!";
        }

    void T2_6_1::getPoint(DimType i, PointReference<2>& p) const
        {
            if (i < 12)
                p=_gp[i];
            else
                throw _name + "::getPoint -  wrong index!";
        }

    Geometry::ReferenceGeometry<2>* T2_6_1::forReferenceGeometry() const
        {
            return _refGeoPtr;
        }

    T2_6_1::T2_6_1()
     : _name("T2_6_1"),
       _refGeoPtr(&ReferenceTriangle::Instance())
            {
                _weight[0] = 5.8393137863189497e-02;
                _gp[0][0] = 2.4928674517090998e-01;
                _gp[0][1] = 2.4928674517090998e-01;

                _weight[1] = 5.8393137863189497e-02;
                _gp[1][0] = 2.4928674517090998e-01;
                _gp[1][1] = 5.0142650965817903e-01;

                _weight[2] = 5.8393137863189497e-02;
                _gp[2][0] = 5.0142650965817903e-01;
                _gp[2][1] = 2.4928674517090998e-01;

                _weight[3] = 2.5422453185103500e-02;
                _gp[3][0] = 6.3089014491502005e-02;
                _gp[3][1] = 6.3089014491502005e-02;

                _weight[4] = 2.5422453185103500e-02;
                _gp[4][0] = 6.3089014491502005e-02;
                _gp[4][1] = 8.7382197101699599e-01;

                _weight[5] = 2.5422453185103500e-02;
                _gp[5][0] = 8.7382197101699599e-01;
                _gp[5][1] = 6.3089014491502005e-02;

                _weight[6] = 4.1425537809187001e-02;
                _gp[6][0] = 3.1035245103378400e-01;
                _gp[6][1] = 6.3650249912139900e-01;

                _weight[7] = 4.1425537809187001e-02;
                _gp[7][0] = 6.3650249912139900e-01;
                _gp[7][1] = 5.3145049844816994e-02;

                _weight[8] = 4.1425537809187001e-02;
                _gp[8][0] = 5.3145049844816994e-02;
                _gp[8][1] = 3.1035245103378400e-01;

                _weight[9] = 4.1425537809187001e-02;
                _gp[9][0] = 3.1035245103378400e-01;
                _gp[9][1] = 5.3145049844816994e-02;

                _weight[10] = 4.1425537809187001e-02;
                _gp[10][0] = 6.3650249912139900e-01;
                _gp[10][1] = 3.1035245103378400e-01;

                _weight[11] = 4.1425537809187001e-02;
                _gp[11][0] = 5.3145049844816994e-02;
                _gp[11][1] = 6.3650249912139900e-01;

                _refGeoPtr->addGaussQuadratureRule(this);
            }

    T2_6_1::~T2_6_1()
          {}

    namespace
    {
        const T2_6_1& instantiateT2_6_1 = T2_6_1::Instance();
    }
//---------------------------------------------------------------------------
    std::string T2_7_1::getName() const
        {
            return _name;
        }

    DimType  T2_7_1::order() const
        {
            return 7;
        }

    DimType T2_7_1::dimension() const
        {
            return 2;
        }

    DimType T2_7_1::nrOfPoints() const
        {
            return 13;
        }

    NumType T2_7_1::weight(DimType i) const
        {
            if (i < 13)
                return _weight[i];
            else
                throw _name + "::weight - wrong index!";
        }

    void T2_7_1::getPoint(DimType i, PointReference<2>& p) const
        {
            if (i < 13)
                p=_gp[i];
            else
                throw _name + "::getPoint -  wrong index!";
        }

    Geometry::ReferenceGeometry<2>* T2_7_1::forReferenceGeometry() const
        {
            return _refGeoPtr;
        }

    T2_7_1::T2_7_1()
     : _name("T2_7_1"),
       _refGeoPtr(&ReferenceTriangle::Instance())
            {
                _weight[0] = -7.4785022233840995e-02;
                _gp[0][0] = 3.3333333333333348e-01;
                _gp[0][1] = 3.3333333333333348e-01;

                _weight[1] = 8.7807628716603997e-02;
                _gp[1][0] = 2.6034596607904004e-01;
                _gp[1][1] = 2.6034596607904004e-01;

                _weight[2] = 8.7807628716603997e-02;
                _gp[2][0] = 2.6034596607904004e-01;
                _gp[2][1] = 4.7930806784191998e-01;

                _weight[3] = 8.7807628716603997e-02;
                _gp[3][0] = 4.7930806784191998e-01;
                _gp[3][1] = 2.6034596607904004e-01;

                _weight[4] = 2.6673617804419000e-02;
                _gp[4][0] = 6.5130102902215992e-02;
                _gp[4][1] = 6.5130102902215992e-02;

                _weight[5] = 2.6673617804419000e-02;
                _gp[5][0] = 6.5130102902215992e-02;
                _gp[5][1] = 8.6973979419556802e-01;

                _weight[6] = 2.6673617804419000e-02;
                _gp[6][0] = 8.6973979419556802e-01;
                _gp[6][1] = 6.5130102902215992e-02;

                _weight[7] = 3.8556880445128498e-02;
                _gp[7][0] = 3.1286549600487401e-01;
                _gp[7][1] = 6.3844418856981000e-01;

                _weight[8] = 3.8556880445128498e-02;
                _gp[8][0] = 6.3844418856981000e-01;
                _gp[8][1] = 4.8690315425315989e-02;

                _weight[9] = 3.8556880445128498e-02;
                _gp[9][0] = 4.8690315425315989e-02;
                _gp[9][1] = 3.1286549600487401e-01;

                _weight[10] = 3.8556880445128498e-02;
                _gp[10][0] = 3.1286549600487401e-01;
                _gp[10][1] = 4.8690315425315989e-02;

                _weight[11] = 3.8556880445128498e-02;
                _gp[11][0] = 6.3844418856981000e-01;
                _gp[11][1] = 3.1286549600487401e-01;
 
                _weight[12] = 3.8556880445128498e-02;
                _gp[12][0] = 4.8690315425315989e-02;
                _gp[12][1] = 6.3844418856981000e-01;

                _refGeoPtr->addGaussQuadratureRule(this);
            }

    T2_7_1::~T2_7_1()
          {}

    namespace
    {
        const T2_7_1& instantiateT2_7_1 = T2_7_1::Instance();
    }
//---------------------------------------------------------------------------
    std::string T2_8_1::getName() const
        {
            return _name;
        }

    DimType  T2_8_1::order() const
        {
            return 8;
        }

    DimType T2_8_1::dimension() const
        {
            return 2;
        }

    DimType T2_8_1::nrOfPoints() const
        {
            return 16;
        }

    NumType T2_8_1::weight(DimType i) const
        {
            if (i < 16)
                return _weight[i];
            else
                throw _name + "::weight - wrong index!";
        }

    void T2_8_1::getPoint(DimType i, PointReference<2>& p) const
        {
            if (i < 16)
                p=_gp[i];
            else
                throw _name + "::getPoint -  wrong index!";
        }

    Geometry::ReferenceGeometry<2>* T2_8_1::forReferenceGeometry() const
        {
            return _refGeoPtr;
        }

    T2_8_1::T2_8_1()
     : _name("T2_8_1"),
       _refGeoPtr(&ReferenceTriangle::Instance())
            {
                _weight[0] = 7.2157803838893503e-02;
                _gp[0][0] = 3.3333333333333348e-01;
                _gp[0][1] = 3.3333333333333348e-01;
                
                _weight[1] = 4.7545817133642497e-02;
                _gp[1][0] = 4.5929258829272301e-01;
                _gp[1][1] = 4.5929258829272301e-01;

                _weight[2] = 4.7545817133642497e-02;
                _gp[2][0] = 4.5929258829272301e-01;
                _gp[2][1] = 8.1414823414554027e-02;
                
                _weight[3] = 4.7545817133642497e-02;
                _gp[3][0] = 8.1414823414554027e-02;
                _gp[3][1] = 4.5929258829272301e-01;
                
                _weight[4] = 5.1608685267358997e-02;
                _gp[4][0] = 1.7056930775175999e-01;
                _gp[4][1] = 1.7056930775175999e-01;

                _weight[5] = 5.1608685267358997e-02;
                _gp[5][0] = 1.7056930775175999e-01;
                _gp[5][1] = 6.5886138449648002e-01;

                _weight[6] = 5.1608685267358997e-02;
                _gp[6][0] = 6.5886138449648002e-01;
                _gp[6][1] = 1.7056930775175999e-01;

                _weight[7] = 1.6229248811599001e-02;
                _gp[7][0] = 5.0547228317031012e-02;
                _gp[7][1] = 5.0547228317031012e-02;

                _weight[8] = 1.6229248811598752e-02;
                _gp[8][0] = 5.0547228317031012e-02;
                _gp[8][1] = 8.9890554336593798e-01;

                _weight[9] = 1.6229248811599001e-02;
                _gp[9][0] = 8.9890554336593798e-01;
                _gp[9][1] = 5.0547228317031012e-02;

                _weight[10] = 1.3615157087217500e-02;
                _gp[10][0] = 2.6311282963463800e-01;
                _gp[10][1] = 7.2849239295540402e-01;

                _weight[11] = 1.3615157087217500e-02;
                _gp[11][0] = 7.2849239295540402e-01;
                _gp[11][1] = 8.3947774099579764e-03;

                _weight[12] = 1.3615157087217500e-02;
                _gp[12][0] = 8.3947774099579764e-03;
                _gp[12][1] = 2.6311282963463800e-01;

                _weight[13] = 1.3615157087217500e-02;
                _gp[13][0] = 2.6311282963463800e-01;
                _gp[13][1] = 8.3947774099579764e-03;

                _weight[14] = 1.3615157087217500e-02;
                _gp[14][0] = 7.2849239295540402e-01;
                _gp[14][1] = 2.6311282963463800e-01;

                _weight[15] = 1.3615157087217500e-02;
                _gp[15][0] = 8.3947774099579764e-03;
                _gp[15][1] = 7.2849239295540402e-01;

                _refGeoPtr->addGaussQuadratureRule(this);
            }

    T2_8_1::~T2_8_1()
          {}

    namespace
    {
        const T2_8_1& instantiateT2_8_1 = T2_8_1::Instance();
    }
//---------------------------------------------------------------------------
    std::string T2_9_1::getName() const
        {
            return _name;
        }

    DimType  T2_9_1::order() const
        {
            return 9;
        }

    DimType T2_9_1::dimension() const
        {
            return 2;
        }

    DimType T2_9_1::nrOfPoints() const
        {
            return 19;
        }

    NumType T2_9_1::weight(DimType i) const
        {
            if (i < 19)
                return _weight[i];
            else
                throw _name + "::weight - wrong index!";
        }

    void T2_9_1::getPoint(DimType i, PointReference<2>& p) const
        {
            if (i < 19)
                p=_gp[i];
            else
                throw _name + "::getPoint -  wrong index!";
        }

    Geometry::ReferenceGeometry<2>* T2_9_1::forReferenceGeometry() const
        {
            return _refGeoPtr;
        }

    T2_9_1::T2_9_1()
     : _name("T2_9_1"),
       _refGeoPtr(&ReferenceTriangle::Instance())
            {
                _weight[0] = 4.8567898141399501e-02;
                _gp[0][0] = 3.3333333333333348e-01;
                _gp[0][1] = 3.3333333333333348e-01;
                
                _weight[1] = 1.5667350113569501e-02;
                _gp[1][0] = 4.8968251919873801e-01;
                _gp[1][1] = 4.8968251919873801e-01;
                
                _weight[2] = 1.5667350113569501e-02;
                _gp[2][0] = 4.8968251919873801e-01;
                _gp[2][1] = 2.0634961602524982e-02;
                
                _weight[3] = 1.5667350113569501e-02;
                _gp[3][0] = 2.0634961602524982e-02;
                _gp[3][1] = 4.8968251919873801e-01;
                
                _weight[4] = 3.8913770502387000e-02;
                _gp[4][0] = 4.3708959149293702e-01;
                _gp[4][1] = 4.3708959149293702e-01;

                _weight[5] = 3.8913770502387000e-02;
                _gp[5][0] = 4.3708959149293702e-01;
                _gp[5][1] = 1.2582081701412701e-01;

                _weight[6] = 3.8913770502387000e-02;
                _gp[6][0] = 1.2582081701412701e-01;
                _gp[6][1] = 4.3708959149293702e-01;

                _weight[7] = 3.9823869463604999e-02;
                _gp[7][0] = 1.8820353561903302e-01;
                _gp[7][1] = 1.8820353561903302e-01;

                _weight[8] = 3.9823869463604999e-02;
                _gp[8][0] = 1.8820353561903302e-01;
                _gp[8][1] = 6.2359292876193506e-01;
                
                _weight[9] = 3.9823869463604999e-02;
                _gp[9][0] = 6.2359292876193506e-01;
                _gp[9][1] = 1.8820353561903302e-01;
                
                _weight[10] = 1.2788837829349000e-02;
                _gp[10][0] = 4.4729513394453024e-02;
                _gp[10][1] = 4.4729513394453024e-02;

                _weight[11] = 1.2788837829349000e-02;
                _gp[11][0] = 4.4729513394453024e-02;
                _gp[11][1] = 9.1054097321109495e-01;

                _weight[12] = 1.2788837829349000e-02;
                _gp[12][0] = 9.1054097321109495e-01;
                _gp[12][1] = 4.4729513394453024e-02;

                _weight[13] = 2.1641769688644501e-02;
                _gp[13][0] = 2.2196298916076601e-01;
                _gp[13][1] = 7.4119859878449801e-01;
                
                _weight[14] = 2.1641769688644501e-02;
                _gp[14][0] = 7.4119859878449801e-01;
                _gp[14][1] = 3.6838412054735981e-02;

                _weight[15] = 2.1641769688644501e-02;
                _gp[15][0] = 3.6838412054735981e-02;
                _gp[15][1] = 2.2196298916076601e-01;
                
                _weight[16] = 2.1641769688644501e-02;
                _gp[16][0] = 2.2196298916076601e-01;
                _gp[16][1] = 3.6838412054735981e-02;
                
                _weight[17] = 2.1641769688644501e-02;
                _gp[17][0] = 7.4119859878449801e-01;
                _gp[17][1] = 2.2196298916076601e-01;
                
                _weight[18] = 2.1641769688644501e-02;
                _gp[18][0] = 3.6838412054735981e-02;
                _gp[18][1] = 7.4119859878449801e-01;

                _refGeoPtr->addGaussQuadratureRule(this);
            }

    T2_9_1::~T2_9_1()
          {}

    namespace
    {
        const T2_9_1& instantiateT2_9_1 = T2_9_1::Instance();
    }
//---------------------------------------------------------------------------
    std::string T2_10_1::getName() const
        {
            return _name;
        }

    DimType  T2_10_1::order() const
        {
            return 10;
        }

    DimType T2_10_1::dimension() const
        {
            return 2;
        }

    DimType T2_10_1::nrOfPoints() const
        {
            return 25;
        }

    NumType T2_10_1::weight(DimType i) const
        {
            if (i < 25)
                return _weight[i];
            else
                throw _name + "::weight - wrong index!";
        }

    void T2_10_1::getPoint(DimType i, PointReference<2>& p) const
        {
            if (i < 25)
                p=_gp[i];
            else
                throw _name + "::getPoint -  wrong index!";
        }

    Geometry::ReferenceGeometry<2>* T2_10_1::forReferenceGeometry() const
        {
            return _refGeoPtr;
        }

    T2_10_1::T2_10_1()
     : _name("T2_10_1"),
       _refGeoPtr(&ReferenceTriangle::Instance())
            {
                _weight[0] = 4.5408995191376998e-02;
                _gp[0][0] = 3.3333333333333348e-01;
                _gp[0][1] = 3.3333333333333348e-01;

                _weight[1] = 1.8362978878233498e-02;
                _gp[1][0] = 4.8557763338365700e-01;
                _gp[1][1] = 4.8557763338365700e-01;
                
                _weight[2] = 1.8362978878233498e-02;
                _gp[2][0] = 4.8557763338365700e-01;
                _gp[2][1] = 2.8844733232684994e-02;
                
                _weight[3] = 1.8362978878233498e-02;
                _gp[3][0] = 2.8844733232684994e-02;
                _gp[3][1] = 4.8557763338365700e-01;
                
                _weight[4] = 2.2660529717764000e-02;
                _gp[4][0] = 1.0948157548503701e-01;
                _gp[4][1] = 1.0948157548503701e-01;

                _weight[5] = 2.2660529717764000e-02;
                _gp[5][0] = 1.0948157548503701e-01;
                _gp[5][1] = 7.8103684902992598e-01;

                _weight[6] = 2.2660529717764000e-02;
                _gp[6][0] = 7.8103684902992598e-01;
                _gp[6][1] = 1.0948157548503701e-01;

                _weight[7] = 3.6378958422710002e-02;
                _gp[7][0] = 3.0793983876412101e-01;
                _gp[7][1] = 5.5035294182099903e-01;

                _weight[8] = 3.6378958422710002e-02;
                _gp[8][0] = 5.5035294182099903e-01;
                _gp[8][1] = 1.4170721941488001e-01;

                _weight[9] = 3.6378958422710002e-02;
                _gp[9][0] = 1.4170721941488001e-01;
                _gp[9][1] = 3.0793983876412101e-01;

                _weight[10] = 3.6378958422710002e-02;
                _gp[10][0] = 3.0793983876412101e-01;
                _gp[10][1] = 1.4170721941488001e-01;
                
                _weight[11] = 3.6378958422710002e-02;
                _gp[11][0] = 5.5035294182099903e-01;
                _gp[11][1] = 3.0793983876412101e-01;

                _weight[12] = 3.6378958422710002e-02;
                _gp[12][0] = 1.4170721941488001e-01;
                _gp[12][1] = 5.5035294182099903e-01;

                _weight[13] = 1.4163621265528500e-02;
                _gp[13][0] = 2.4667256063990300e-01;
                _gp[13][1] = 7.2832390459741103e-01;

                _weight[14] = 1.4163621265528500e-02;
                _gp[14][0] = 7.2832390459741103e-01;
                _gp[14][1] = 2.5003534762686019e-02;
 
                _weight[15] = 1.4163621265528500e-02;
                _gp[15][0] = 2.5003534762686019e-02;
                _gp[15][1] = 2.4667256063990300e-01;

                _weight[16] = 1.4163621265528500e-02;
                _gp[16][0] = 2.4667256063990300e-01;
                _gp[16][1] = 2.5003534762686019e-02;

                _weight[17] = 1.4163621265528500e-02;
                _gp[17][0] = 7.2832390459741103e-01;
                _gp[17][1] = 2.4667256063990300e-01;

                _weight[18] = 1.4163621265528500e-02;
                _gp[18][0] = 2.5003534762686019e-02;
                _gp[18][1] = 7.2832390459741103e-01;

                _weight[19] = 4.7108334818665000e-03;
                _gp[19][0] = 6.6803251012200027e-02;
                _gp[19][1] = 9.2365593358749998e-01;

                _weight[20] = 4.7108334818665000e-03;
                _gp[20][0] = 9.2365593358749998e-01;
                _gp[20][1] = 9.5408154002989964e-03;

                _weight[21] = 4.7108334818665000e-03;
                _gp[21][0] = 9.5408154002989964e-03;
                _gp[21][1] = 6.6803251012200027e-02;

                _weight[22] = 4.7108334818665000e-03;
                _gp[22][0] = 6.6803251012200027e-02;
                _gp[22][1] = 9.5408154002989964e-03;

                _weight[23] = 4.7108334818665000e-03;
                _gp[23][0] = 9.2365593358749998e-01;
                _gp[23][1] = 6.6803251012200027e-02;

                _weight[24] = 4.7108334818665000e-03;
                _gp[24][0] = 9.5408154002989964e-03;
                _gp[24][1] = 9.2365593358749998e-01;

                _refGeoPtr->addGaussQuadratureRule(this);
            }

    T2_10_1::~T2_10_1()
          {}

    namespace
    {
        const T2_10_1& instantiateT2_10_1 = T2_10_1::Instance();
    }

//---------------------------------------------------------------------------
} // close namespace IntegrationRules

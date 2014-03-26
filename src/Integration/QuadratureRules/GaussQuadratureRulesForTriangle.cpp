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
#include "Geometry/ReferenceTriangle.hpp"
using Geometry::ReferenceTriangle;
//---------------------------------------------------------------------------
namespace QuadratureRules
{
//---------------------------------------------------------------------------
    std::string
    Tn2_1_1::getName() const
    {
        return name_;
    }

    unsigned int
    Tn2_1_1::order() const
    {
        return 1;
    }

    unsigned int
    Tn2_1_1::dimension() const
    {
        return 2;
    }

    unsigned int
    Tn2_1_1::nrOfPoints() const
    {
        return 1;
    }

    double
    Tn2_1_1::weight(unsigned int i) const
    {
        if (i < 1)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }

    void
    Tn2_1_1::getPoint(unsigned int i, PointReferenceT& p) const
    {
        if (i < 1)
            p=gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }

    Tn2_1_1::ReferenceGeometryT*
    Tn2_1_1::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    Tn2_1_1::Tn2_1_1():
        name_("Tn2_1_1"),
        refGeoPtr_(&ReferenceTriangle::Instance()),gp_(1,2)
    {
        weight_[0] = 5.0000000000000000e-01;
        gp_[0][0] = 3.3333333333333348e-01;
        gp_[0][1] = 3.3333333333333348e-01;

        refGeoPtr_->addGaussQuadratureRule(this);
    }

    Tn2_1_1::~Tn2_1_1()
    {
    }

 
//---------------------------------------------------------------------------
    std::string
    Tn2_2_1::getName() const
    {
        return name_;
    }

    unsigned int
    Tn2_2_1::order() const
    {
        return 2;
    }

    unsigned int
    Tn2_2_1::dimension() const
    {
        return 2;
    }

    unsigned int
    Tn2_2_1::nrOfPoints() const
    {
        return 3;
    }

    double
    Tn2_2_1::weight(unsigned int i) const
    {
        if (i < 3)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }

    void
    Tn2_2_1::getPoint(unsigned int i, PointReferenceT& p) const
    {
        if (i < 3)
            p=gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }

    Tn2_2_1::ReferenceGeometryT*
    Tn2_2_1::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    Tn2_2_1::Tn2_2_1():
        name_("Tn2_2_1"),
        refGeoPtr_(&ReferenceTriangle::Instance()),gp_(3,2)
    {
        weight_[0] = 1.6666666666666674e-01;
        gp_[0][0] = 1.6666666666666652e-01;
        gp_[0][1] = 1.6666666666666652e-01;
        
        weight_[1] = 1.6666666666666674e-01;
        gp_[1][0] = 1.6666666666666652e-01;
        gp_[1][1] = 6.6666666666666652e-01;
        
        weight_[2] = 1.6666666666666674e-01;
        gp_[2][0] = 6.6666666666666652e-01;
        gp_[2][1] = 1.6666666666666652e-01;

        refGeoPtr_->addGaussQuadratureRule(this);
     }

    Tn2_2_1::~Tn2_2_1()
    {
    }


//---------------------------------------------------------------------------
    std::string
    Tn2_3_1::getName() const
    {
        return name_;
    }

    unsigned int
    Tn2_3_1::order() const
    {
        return 3;
    }

    unsigned int
    Tn2_3_1::dimension() const
    {
        return 2;
    }

    unsigned int
    Tn2_3_1::nrOfPoints() const
    {
        return 4;
    }

    double
    Tn2_3_1::weight(unsigned int i) const
    {
        if (i < 4)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }

    void
    Tn2_3_1::getPoint(unsigned int i, PointReferenceT& p) const
    {
        if (i < 4)
            p=gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }

    Tn2_3_1::ReferenceGeometryT*
    Tn2_3_1::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    Tn2_3_1::Tn2_3_1():
        name_("Tn2_3_1"),
        refGeoPtr_(&ReferenceTriangle::Instance()),gp_(4,2)
    {
        weight_[0] = -2.8125000000000000e-01;
        gp_[0][0] = 3.3333333333333348e-01;
        gp_[0][1] = 3.3333333333333348e-01;

        weight_[1] = 2.6041666666666674e-01;
        gp_[1][0] = 2.0000000000000001e-01;
        gp_[1][1] = 2.0000000000000001e-01;
        
        weight_[2] = 2.6041666666666674e-01;
        gp_[2][0] = 2.0000000000000001e-01;
        gp_[2][1] = 5.9999999999999998e-01;

        weight_[3] = 2.6041666666666674e-01;
        gp_[3][0] = 5.9999999999999998e-01;
        gp_[3][1] = 2.0000000000000001e-01;

        refGeoPtr_->addGaussQuadratureRule(this);
    }

    Tn2_3_1::~Tn2_3_1()
    {
    }


//---------------------------------------------------------------------------
    std::string
    Tn2_4_1::getName() const
    {
        return name_;
    }

    unsigned int
    Tn2_4_1::order() const
    {
        return 4;
    }

    unsigned int
    Tn2_4_1::dimension() const
    {
        return 2;
    }

    unsigned int
    Tn2_4_1::nrOfPoints() const
    {
        return 6;
    }

    double
    Tn2_4_1::weight(unsigned int i) const
    {
        if (i < 6)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }

    void
    Tn2_4_1::getPoint(unsigned int i, PointReferenceT& p) const
    {
        if (i < 6)
            p=gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }

    Tn2_4_1::ReferenceGeometryT*
    Tn2_4_1::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    Tn2_4_1::Tn2_4_1():
        name_("Tn2_4_1"),
        refGeoPtr_(&ReferenceTriangle::Instance()),gp_(6,2)
    {
        weight_[0] = 1.1169079483900550e-01;
        gp_[0][0] = 4.4594849091596500e-01;
        gp_[0][1] = 4.4594849091596500e-01;

        weight_[1] = 1.1169079483900550e-01;
        gp_[1][0] = 4.4594849091596500e-01;
        gp_[1][1] = 1.0810301816807000e-01;

        weight_[2] = 1.1169079483900550e-01;
        gp_[2][0] = 1.0810301816807000e-01;
        gp_[2][1] = 4.4594849091596500e-01;

        weight_[3] = 5.4975871827660998e-02;
        gp_[3][0] = 9.1576213509771021e-02;
        gp_[3][1] = 9.1576213509771021e-02;

        weight_[4] = 5.4975871827660998e-02;
        gp_[4][0] = 9.1576213509771021e-02;
        gp_[4][1] = 8.1684757298045896e-01;

        weight_[5] = 5.4975871827660998e-02;
        gp_[5][0] = 8.1684757298045896e-01;
        gp_[5][1] = 9.1576213509771021e-02;

        refGeoPtr_->addGaussQuadratureRule(this);
    }

    Tn2_4_1::~Tn2_4_1()
    {
    }

//---------------------------------------------------------------------------
    std::string
    T2_5_1::getName() const
    {
        return name_;
    }

    unsigned int
    T2_5_1::order() const
    {
        return 5;
    }

    unsigned int
    T2_5_1::dimension() const
    {
        return 2;
    }

    unsigned int
    T2_5_1::nrOfPoints() const
    {
        return 7;
    }

    double
    T2_5_1::weight(unsigned int i) const
    {
        if (i < 7)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }

    void
    T2_5_1::getPoint(unsigned int i, PointReferenceT& p) const
    {
        if (i < 7)
            p=gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }

    T2_5_1::ReferenceGeometryT*
    T2_5_1::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    T2_5_1::T2_5_1():
        name_("T2_5_1"),
        refGeoPtr_(&ReferenceTriangle::Instance()),gp_(7,2)
    {
        weight_[0] = 1.1250000000000000e-01;
        gp_[0][0] = 3.3333333333333348e-01;
        gp_[0][1] = 3.3333333333333348e-01;

        weight_[1] = 6.6197076394252999e-02;
        gp_[1][0] = 4.7014206410511500e-01;
        gp_[1][1] = 4.7014206410511500e-01;
        
        weight_[2] = 6.6197076394252999e-02;
        gp_[2][0] = 4.7014206410511500e-01;
        gp_[2][1] = 5.9715871789770003e-02;

        weight_[3] = 6.6197076394252999e-02;
        gp_[3][0] = 5.9715871789770003e-02;
        gp_[3][1] = 4.7014206410511500e-01;

        weight_[4] = 6.2969590272413500e-02;
        gp_[4][0] = 1.0128650732345601e-01;
        gp_[4][1] = 1.0128650732345601e-01;

        weight_[5] = 6.2969590272413500e-02;
        gp_[5][0] = 1.0128650732345601e-01;
        gp_[5][1] = 7.9742698535308698e-01;

        weight_[6] = 6.2969590272413500e-02;
        gp_[6][0] = 7.9742698535308698e-01;
        gp_[6][1] = 1.0128650732345601e-01;

        refGeoPtr_->addGaussQuadratureRule(this);
    }

    T2_5_1::~T2_5_1()
    {
    }

//---------------------------------------------------------------------------
    std::string
    T2_6_1::getName() const
    {
        return name_;
    }

    unsigned int
    T2_6_1::order() const
    {
        return 6;
    }

    unsigned int
    T2_6_1::dimension() const
    {
        return 2;
    }

    unsigned int
    T2_6_1::nrOfPoints() const
    {
        return 12;
    }

    double
    T2_6_1::weight(unsigned int i) const
    {
        if (i < 12)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }

    void
    T2_6_1::getPoint(unsigned int i, PointReferenceT& p) const
    {
        if (i < 12)
            p=gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }

    T2_6_1::ReferenceGeometryT*
    T2_6_1::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    T2_6_1::T2_6_1():
        name_("T2_6_1"),
        refGeoPtr_(&ReferenceTriangle::Instance()),gp_(12,2)
    {
        weight_[0] = 5.8393137863189497e-02;
        gp_[0][0] = 2.4928674517090998e-01;
        gp_[0][1] = 2.4928674517090998e-01;

        weight_[1] = 5.8393137863189497e-02;
        gp_[1][0] = 2.4928674517090998e-01;
        gp_[1][1] = 5.0142650965817903e-01;

        weight_[2] = 5.8393137863189497e-02;
        gp_[2][0] = 5.0142650965817903e-01;
        gp_[2][1] = 2.4928674517090998e-01;

        weight_[3] = 2.5422453185103500e-02;
        gp_[3][0] = 6.3089014491502005e-02;
        gp_[3][1] = 6.3089014491502005e-02;

        weight_[4] = 2.5422453185103500e-02;
        gp_[4][0] = 6.3089014491502005e-02;
        gp_[4][1] = 8.7382197101699599e-01;

        weight_[5] = 2.5422453185103500e-02;
        gp_[5][0] = 8.7382197101699599e-01;
        gp_[5][1] = 6.3089014491502005e-02;

        weight_[6] = 4.1425537809187001e-02;
        gp_[6][0] = 3.1035245103378400e-01;
        gp_[6][1] = 6.3650249912139900e-01;

        weight_[7] = 4.1425537809187001e-02;
        gp_[7][0] = 6.3650249912139900e-01;
        gp_[7][1] = 5.3145049844816994e-02;

        weight_[8] = 4.1425537809187001e-02;
        gp_[8][0] = 5.3145049844816994e-02;
        gp_[8][1] = 3.1035245103378400e-01;

        weight_[9] = 4.1425537809187001e-02;
        gp_[9][0] = 3.1035245103378400e-01;
        gp_[9][1] = 5.3145049844816994e-02;

        weight_[10] = 4.1425537809187001e-02;
        gp_[10][0] = 6.3650249912139900e-01;
        gp_[10][1] = 3.1035245103378400e-01;

        weight_[11] = 4.1425537809187001e-02;
        gp_[11][0] = 5.3145049844816994e-02;
        gp_[11][1] = 6.3650249912139900e-01;

        refGeoPtr_->addGaussQuadratureRule(this);
    }

    T2_6_1::~T2_6_1()
    {
    }


//---------------------------------------------------------------------------
    std::string
    T2_7_1::getName() const
    {
        return name_;
    }

    unsigned int
    T2_7_1::order() const
    {
        return 7;
    }

    unsigned int
    T2_7_1::dimension() const
    {
        return 2;
    }

    unsigned int
    T2_7_1::nrOfPoints() const
    {
        return 13;
    }

    double
    T2_7_1::weight(unsigned int i) const
    {
        if (i < 13)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }

    void
    T2_7_1::getPoint(unsigned int i, PointReferenceT& p) const
    {
        if (i < 13)
            p=gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }

    T2_7_1::ReferenceGeometryT*
    T2_7_1::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    T2_7_1::T2_7_1():
        name_("T2_7_1"),
        refGeoPtr_(&ReferenceTriangle::Instance()),gp_(13,2)
    {
        weight_[0] = -7.4785022233840995e-02;
        gp_[0][0] = 3.3333333333333348e-01;
        gp_[0][1] = 3.3333333333333348e-01;

        weight_[1] = 8.7807628716603997e-02;
        gp_[1][0] = 2.6034596607904004e-01;
        gp_[1][1] = 2.6034596607904004e-01;

        weight_[2] = 8.7807628716603997e-02;
        gp_[2][0] = 2.6034596607904004e-01;
        gp_[2][1] = 4.7930806784191998e-01;

        weight_[3] = 8.7807628716603997e-02;
        gp_[3][0] = 4.7930806784191998e-01;
        gp_[3][1] = 2.6034596607904004e-01;

        weight_[4] = 2.6673617804419000e-02;
        gp_[4][0] = 6.5130102902215992e-02;
        gp_[4][1] = 6.5130102902215992e-02;

        weight_[5] = 2.6673617804419000e-02;
        gp_[5][0] = 6.5130102902215992e-02;
        gp_[5][1] = 8.6973979419556802e-01;

        weight_[6] = 2.6673617804419000e-02;
        gp_[6][0] = 8.6973979419556802e-01;
        gp_[6][1] = 6.5130102902215992e-02;

        weight_[7] = 3.8556880445128498e-02;
        gp_[7][0] = 3.1286549600487401e-01;
        gp_[7][1] = 6.3844418856981000e-01;

        weight_[8] = 3.8556880445128498e-02;
        gp_[8][0] = 6.3844418856981000e-01;
        gp_[8][1] = 4.8690315425315989e-02;

        weight_[9] = 3.8556880445128498e-02;
        gp_[9][0] = 4.8690315425315989e-02;
        gp_[9][1] = 3.1286549600487401e-01;

        weight_[10] = 3.8556880445128498e-02;
        gp_[10][0] = 3.1286549600487401e-01;
        gp_[10][1] = 4.8690315425315989e-02;

        weight_[11] = 3.8556880445128498e-02;
        gp_[11][0] = 6.3844418856981000e-01;
        gp_[11][1] = 3.1286549600487401e-01;

        weight_[12] = 3.8556880445128498e-02;
        gp_[12][0] = 4.8690315425315989e-02;
        gp_[12][1] = 6.3844418856981000e-01;

        refGeoPtr_->addGaussQuadratureRule(this);
    }

    T2_7_1::~T2_7_1()
    {
    }


//---------------------------------------------------------------------------
    std::string
    T2_8_1::getName() const
    {
        return name_;
    }

    unsigned int
    T2_8_1::order() const
    {
        return 8;
    }

    unsigned int
    T2_8_1::dimension() const
    {
        return 2;
    }

    unsigned int
    T2_8_1::nrOfPoints() const
    {
        return 16;
    }

    double
    T2_8_1::weight(unsigned int i) const
    {
        if (i < 16)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }

    void
    T2_8_1::getPoint(unsigned int i, PointReferenceT& p) const
    {
        if (i < 16)
            p=gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }

    T2_8_1::ReferenceGeometryT*
    T2_8_1::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    T2_8_1::T2_8_1():
        name_("T2_8_1"),
        refGeoPtr_(&ReferenceTriangle::Instance()),gp_(16,2)
    {
        weight_[0] = 7.2157803838893503e-02;
        gp_[0][0] = 3.3333333333333348e-01;
        gp_[0][1] = 3.3333333333333348e-01;
        
        weight_[1] = 4.7545817133642497e-02;
        gp_[1][0] = 4.5929258829272301e-01;
        gp_[1][1] = 4.5929258829272301e-01;

        weight_[2] = 4.7545817133642497e-02;
        gp_[2][0] = 4.5929258829272301e-01;
        gp_[2][1] = 8.1414823414554027e-02;
        
        weight_[3] = 4.7545817133642497e-02;
        gp_[3][0] = 8.1414823414554027e-02;
        gp_[3][1] = 4.5929258829272301e-01;
        
        weight_[4] = 5.1608685267358997e-02;
        gp_[4][0] = 1.7056930775175999e-01;
        gp_[4][1] = 1.7056930775175999e-01;

        weight_[5] = 5.1608685267358997e-02;
        gp_[5][0] = 1.7056930775175999e-01;
        gp_[5][1] = 6.5886138449648002e-01;

        weight_[6] = 5.1608685267358997e-02;
        gp_[6][0] = 6.5886138449648002e-01;
        gp_[6][1] = 1.7056930775175999e-01;

        weight_[7] = 1.6229248811599001e-02;
        gp_[7][0] = 5.0547228317031012e-02;
        gp_[7][1] = 5.0547228317031012e-02;

        weight_[8] = 1.6229248811598752e-02;
        gp_[8][0] = 5.0547228317031012e-02;
        gp_[8][1] = 8.9890554336593798e-01;

        weight_[9] = 1.6229248811599001e-02;
        gp_[9][0] = 8.9890554336593798e-01;
        gp_[9][1] = 5.0547228317031012e-02;

        weight_[10] = 1.3615157087217500e-02;
        gp_[10][0] = 2.6311282963463800e-01;
        gp_[10][1] = 7.2849239295540402e-01;

        weight_[11] = 1.3615157087217500e-02;
        gp_[11][0] = 7.2849239295540402e-01;
        gp_[11][1] = 8.3947774099579764e-03;

        weight_[12] = 1.3615157087217500e-02;
        gp_[12][0] = 8.3947774099579764e-03;
        gp_[12][1] = 2.6311282963463800e-01;

        weight_[13] = 1.3615157087217500e-02;
        gp_[13][0] = 2.6311282963463800e-01;
        gp_[13][1] = 8.3947774099579764e-03;

        weight_[14] = 1.3615157087217500e-02;
        gp_[14][0] = 7.2849239295540402e-01;
        gp_[14][1] = 2.6311282963463800e-01;

        weight_[15] = 1.3615157087217500e-02;
        gp_[15][0] = 8.3947774099579764e-03;
        gp_[15][1] = 7.2849239295540402e-01;

        refGeoPtr_->addGaussQuadratureRule(this);
    }

    T2_8_1::~T2_8_1()
    {
    }

//---------------------------------------------------------------------------
    std::string
    T2_9_1::getName() const
    {
        return name_;
    }

    unsigned int
    T2_9_1::order() const
    {
        return 9;
    }

    unsigned int
    T2_9_1::dimension() const
    {
        return 2;
    }

    unsigned int
    T2_9_1::nrOfPoints() const
    {
        return 19;
    }

    double
    T2_9_1::weight(unsigned int i) const
    {
        if (i < 19)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }

    void
    T2_9_1::getPoint(unsigned int i, PointReferenceT& p) const
    {
        if (i < 19)
            p=gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }

    T2_9_1::ReferenceGeometryT*
    T2_9_1::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    T2_9_1::T2_9_1():
        name_("T2_9_1"),
        refGeoPtr_(&ReferenceTriangle::Instance()),gp_(19,2)
    {
        weight_[0] = 4.8567898141399501e-02;
        gp_[0][0] = 3.3333333333333348e-01;
        gp_[0][1] = 3.3333333333333348e-01;
        
        weight_[1] = 1.5667350113569501e-02;
        gp_[1][0] = 4.8968251919873801e-01;
        gp_[1][1] = 4.8968251919873801e-01;
        
        weight_[2] = 1.5667350113569501e-02;
        gp_[2][0] = 4.8968251919873801e-01;
        gp_[2][1] = 2.0634961602524982e-02;
        
        weight_[3] = 1.5667350113569501e-02;
        gp_[3][0] = 2.0634961602524982e-02;
        gp_[3][1] = 4.8968251919873801e-01;
        
        weight_[4] = 3.8913770502387000e-02;
        gp_[4][0] = 4.3708959149293702e-01;
        gp_[4][1] = 4.3708959149293702e-01;

        weight_[5] = 3.8913770502387000e-02;
        gp_[5][0] = 4.3708959149293702e-01;
        gp_[5][1] = 1.2582081701412701e-01;

        weight_[6] = 3.8913770502387000e-02;
        gp_[6][0] = 1.2582081701412701e-01;
        gp_[6][1] = 4.3708959149293702e-01;

        weight_[7] = 3.9823869463604999e-02;
        gp_[7][0] = 1.8820353561903302e-01;
        gp_[7][1] = 1.8820353561903302e-01;

        weight_[8] = 3.9823869463604999e-02;
        gp_[8][0] = 1.8820353561903302e-01;
        gp_[8][1] = 6.2359292876193506e-01;
        
        weight_[9] = 3.9823869463604999e-02;
        gp_[9][0] = 6.2359292876193506e-01;
        gp_[9][1] = 1.8820353561903302e-01;
        
        weight_[10] = 1.2788837829349000e-02;
        gp_[10][0] = 4.4729513394453024e-02;
        gp_[10][1] = 4.4729513394453024e-02;

        weight_[11] = 1.2788837829349000e-02;
        gp_[11][0] = 4.4729513394453024e-02;
        gp_[11][1] = 9.1054097321109495e-01;

        weight_[12] = 1.2788837829349000e-02;
        gp_[12][0] = 9.1054097321109495e-01;
        gp_[12][1] = 4.4729513394453024e-02;

        weight_[13] = 2.1641769688644501e-02;
        gp_[13][0] = 2.2196298916076601e-01;
        gp_[13][1] = 7.4119859878449801e-01;
        
        weight_[14] = 2.1641769688644501e-02;
        gp_[14][0] = 7.4119859878449801e-01;
        gp_[14][1] = 3.6838412054735981e-02;

        weight_[15] = 2.1641769688644501e-02;
        gp_[15][0] = 3.6838412054735981e-02;
        gp_[15][1] = 2.2196298916076601e-01;
        
        weight_[16] = 2.1641769688644501e-02;
        gp_[16][0] = 2.2196298916076601e-01;
        gp_[16][1] = 3.6838412054735981e-02;
        
        weight_[17] = 2.1641769688644501e-02;
        gp_[17][0] = 7.4119859878449801e-01;
        gp_[17][1] = 2.2196298916076601e-01;
        
        weight_[18] = 2.1641769688644501e-02;
        gp_[18][0] = 3.6838412054735981e-02;
        gp_[18][1] = 7.4119859878449801e-01;

        refGeoPtr_->addGaussQuadratureRule(this);
    }

    T2_9_1::~T2_9_1()
    {
    }


//---------------------------------------------------------------------------
    std::string
    T2_10_1::getName() const
    {
        return name_;
    }

    unsigned int
    T2_10_1::order() const
    {
        return 10;
    }

    unsigned int 
    T2_10_1::dimension() const
    {
        return 2;
    }

    unsigned int 
    T2_10_1::nrOfPoints() const
    {
        return 25;
    }

    double
    T2_10_1::weight(unsigned int i) const
    {
        if (i < 25)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }

    void
    T2_10_1::getPoint(unsigned int i, PointReferenceT& p) const
    {
        if (i < 25)
            p=gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }

    T2_10_1::ReferenceGeometryT*
    T2_10_1::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    T2_10_1::T2_10_1():
        name_("T2_10_1"),
        refGeoPtr_(&ReferenceTriangle::Instance()),gp_(25,2)
    {
        weight_[0] = 4.5408995191376998e-02;
        gp_[0][0] = 3.3333333333333348e-01;
        gp_[0][1] = 3.3333333333333348e-01;

        weight_[1] = 1.8362978878233498e-02;
        gp_[1][0] = 4.8557763338365700e-01;
        gp_[1][1] = 4.8557763338365700e-01;
        
        weight_[2] = 1.8362978878233498e-02;
        gp_[2][0] = 4.8557763338365700e-01;
        gp_[2][1] = 2.8844733232684994e-02;
        
        weight_[3] = 1.8362978878233498e-02;
        gp_[3][0] = 2.8844733232684994e-02;
        gp_[3][1] = 4.8557763338365700e-01;
        
        weight_[4] = 2.2660529717764000e-02;
        gp_[4][0] = 1.0948157548503701e-01;
        gp_[4][1] = 1.0948157548503701e-01;

        weight_[5] = 2.2660529717764000e-02;
        gp_[5][0] = 1.0948157548503701e-01;
        gp_[5][1] = 7.8103684902992598e-01;

        weight_[6] = 2.2660529717764000e-02;
        gp_[6][0] = 7.8103684902992598e-01;
        gp_[6][1] = 1.0948157548503701e-01;

        weight_[7] = 3.6378958422710002e-02;
        gp_[7][0] = 3.0793983876412101e-01;
        gp_[7][1] = 5.5035294182099903e-01;

        weight_[8] = 3.6378958422710002e-02;
        gp_[8][0] = 5.5035294182099903e-01;
        gp_[8][1] = 1.4170721941488001e-01;

        weight_[9] = 3.6378958422710002e-02;
        gp_[9][0] = 1.4170721941488001e-01;
        gp_[9][1] = 3.0793983876412101e-01;

        weight_[10] = 3.6378958422710002e-02;
        gp_[10][0] = 3.0793983876412101e-01;
        gp_[10][1] = 1.4170721941488001e-01;
        
        weight_[11] = 3.6378958422710002e-02;
        gp_[11][0] = 5.5035294182099903e-01;
        gp_[11][1] = 3.0793983876412101e-01;

        weight_[12] = 3.6378958422710002e-02;
        gp_[12][0] = 1.4170721941488001e-01;
        gp_[12][1] = 5.5035294182099903e-01;

        weight_[13] = 1.4163621265528500e-02;
        gp_[13][0] = 2.4667256063990300e-01;
        gp_[13][1] = 7.2832390459741103e-01;

        weight_[14] = 1.4163621265528500e-02;
        gp_[14][0] = 7.2832390459741103e-01;
        gp_[14][1] = 2.5003534762686019e-02;

        weight_[15] = 1.4163621265528500e-02;
        gp_[15][0] = 2.5003534762686019e-02;
        gp_[15][1] = 2.4667256063990300e-01;

        weight_[16] = 1.4163621265528500e-02;
        gp_[16][0] = 2.4667256063990300e-01;
        gp_[16][1] = 2.5003534762686019e-02;

        weight_[17] = 1.4163621265528500e-02;
        gp_[17][0] = 7.2832390459741103e-01;
        gp_[17][1] = 2.4667256063990300e-01;

        weight_[18] = 1.4163621265528500e-02;
        gp_[18][0] = 2.5003534762686019e-02;
        gp_[18][1] = 7.2832390459741103e-01;

        weight_[19] = 4.7108334818665000e-03;
        gp_[19][0] = 6.6803251012200027e-02;
        gp_[19][1] = 9.2365593358749998e-01;

        weight_[20] = 4.7108334818665000e-03;
        gp_[20][0] = 9.2365593358749998e-01;
        gp_[20][1] = 9.5408154002989964e-03;

        weight_[21] = 4.7108334818665000e-03;
        gp_[21][0] = 9.5408154002989964e-03;
        gp_[21][1] = 6.6803251012200027e-02;

        weight_[22] = 4.7108334818665000e-03;
        gp_[22][0] = 6.6803251012200027e-02;
        gp_[22][1] = 9.5408154002989964e-03;

        weight_[23] = 4.7108334818665000e-03;
        gp_[23][0] = 9.2365593358749998e-01;
        gp_[23][1] = 6.6803251012200027e-02;

        weight_[24] = 4.7108334818665000e-03;
        gp_[24][0] = 9.5408154002989964e-03;
        gp_[24][1] = 9.2365593358749998e-01;

        refGeoPtr_->addGaussQuadratureRule(this);
    }

    T2_10_1::~T2_10_1()
    {
    }

//---------------------------------------------------------------------------
    std::string
    T2_11_1::getName() const
    {
        return name_;
    }

    unsigned int
    T2_11_1::order() const
    {
        return 11;
    }

    unsigned int
    T2_11_1::dimension() const
    {
        return 2;
    }

    unsigned int
    T2_11_1::nrOfPoints() const
    {
        return 28;
    }

    double
    T2_11_1::weight(unsigned int i) const
    {
        if (i < 28)
            return weight_[i];
        else
            throw name_ + "::weight - wrong index!";
    }

    void
    T2_11_1::getPoint(unsigned int i, PointReferenceT& p) const
    {
        if (i < 28)
            p=gp_[i];
        else
            throw name_ + "::getPoint -  wrong index!";
    }

    T2_11_1::ReferenceGeometryT*
    T2_11_1::forReferenceGeometry() const
    {
        return refGeoPtr_;
    }

    T2_11_1::T2_11_1()://warning: points and weights found vary wildly between runs of the quadrature rule generator
        name_("T2_11_1"),
        refGeoPtr_(&ReferenceTriangle::Instance()),gp_(28,2)
    {
        weight_[0] = 8.005121009414541471823946e-02/2.;
        gp_[0][0] = 3.3333333333333333e-01;
        gp_[0][1] = 3.3333333333333333e-01;

        weight_[1] = 6.71206859217794936920366e-02/2.;
        gp_[1][0] = 2.154717529698554e-01;
        gp_[1][1] = 2.154717529698554e-01;

        weight_[2] = 6.71206859217794936920366e-02/2.;
        gp_[2][0] = 2.154717529698554e-01;
        gp_[2][1] = 1.-2.*2.154717529698554e-01;

        weight_[3] = 6.71206859217794936920366e-02/2.;
        gp_[3][0] = 1.-2.*2.154717529698554e-01;
        gp_[3][1] = 2.154717529698554e-01;

        weight_[4] = 1.260585310278218701321e-02/2.;
        gp_[4][0] = 3.12795522110139299208957e-02;
        gp_[4][1] = 3.12795522110139299208957e-02;

        weight_[5] = 1.260585310278218701321e-02/2.;
        gp_[5][0] = 3.12795522110139299208957e-02;
        gp_[5][1] = 1.-2.*3.12795522110139299208957e-02;

        weight_[6] = 1.260585310278218701321e-02/2.;
        gp_[6][0] = 1.-2.*3.12795522110139299208957e-02;
        gp_[6][1] = 3.12795522110139299208957e-02;

        weight_[7] = 4.05824609257477118939e-02/2.;
        gp_[7][0] = 1.155035357737921489271e-01;
        gp_[7][1] = 1.155035357737921489271e-01;

        weight_[8] = 4.05824609257477118939e-02/2.;
        gp_[8][0] = 1.155035357737921489271e-01;
        gp_[8][1] = 1.-2.*1.155035357737921489271e-01;

        weight_[9] = 4.05824609257477118939e-02/2.;
        gp_[9][0] = 1.-2.*1.155035357737921489271e-01;
        gp_[9][1] = 1.155035357737921489271e-01;

        weight_[10] = 6.17123685629664404606197e-02/2.;
        gp_[10][0] = 4.357774152823301051418e-01;
        gp_[10][1] = 4.357774152823301051418e-01;

        weight_[11] = 6.17123685629664404606197e-02/2.;
        gp_[11][0] = 4.357774152823301051418e-01;
        gp_[11][1] = 1.-2.*4.357774152823301051418e-01;

        weight_[12] = 6.17123685629664404606197e-02/2.;
        gp_[12][0] = 1.-2.*4.357774152823301051418e-01;
        gp_[12][1] = 4.357774152823301051418e-01;

        weight_[13] = 1.12084059305228051993e-02/2.;
        gp_[13][0] = 4.9988403867731539816248e-01;
        gp_[13][1] = 4.9988403867731539816248e-01;

        weight_[14] = 1.12084059305228051993e-02/2.;
        gp_[14][0] = 4.9988403867731539816248e-01;
        gp_[14][1] = 1.-2.*4.9988403867731539816248e-01;

        weight_[15] = 1.12084059305228051993e-02/2.;
        gp_[15][0] = 1.-2.*4.9988403867731539816248e-01;
        gp_[15][1] = 4.9988403867731539816248e-01;

        weight_[16] = 4.115295964681454289321412e-02/2.;
        gp_[16][0] = 3.164262022308527e-01;
        gp_[16][1] = 4.8055796920866524676254359e-02;

        weight_[17] = 4.115295964681454289321412e-02/2.;
        gp_[17][0] = 4.8055796920866524676254359e-02;
        gp_[17][1] = 3.164262022308527e-01;

        weight_[18] = 4.115295964681454289321412e-02/2.;
        gp_[18][0] = 3.164262022308527e-01;
        gp_[18][1] = 1.-4.8055796920866524676254359e-02-3.164262022308527e-01;

        weight_[19] = 4.115295964681454289321412e-02/2.;
        gp_[19][0] = 1.-4.8055796920866524676254359e-02-3.164262022308527e-01;
        gp_[19][1] = 3.164262022308527e-01;

        weight_[20] = 4.115295964681454289321412e-02/2.;
        gp_[20][0] = 4.8055796920866524676254359e-02;
        gp_[20][1] = 1.-4.8055796920866524676254359e-02-3.164262022308527e-01;

        weight_[21] = 4.115295964681454289321412e-02/2.;
        gp_[21][0] = 1.-4.8055796920866524676254359e-02-3.164262022308527e-01;
        gp_[21][1] = 4.8055796920866524676254359e-02;

        weight_[22] = 1.55569514489285688575e-02/2.;
        gp_[22][0] = 1.614188489108937e-01;
        gp_[22][1] = 1.56544537223683037837e-02;

        weight_[23] = 1.55569514489285688575e-02/2.;
        gp_[23][0] = 1.56544537223683037837e-02;
        gp_[23][1] = 1.614188489108937e-01;

        weight_[24] = 1.55569514489285688575e-02/2.;
        gp_[24][0] = 1.614188489108937e-01;
        gp_[24][1] = 1.-1.614188489108937e-01-1.56544537223683037837e-02;

        weight_[25] = 1.55569514489285688575e-02/2.;
        gp_[25][0] = 1.-1.614188489108937e-01-1.56544537223683037837e-02;
        gp_[25][1] = 1.614188489108937e-01;

        weight_[26] = 1.55569514489285688575e-02/2.;
        gp_[26][0] = 1.56544537223683037837e-02;
        gp_[26][1] = 1.-1.614188489108937e-01-1.56544537223683037837e-02;

        weight_[27] = 1.55569514489285688575e-02/2.;
        gp_[27][0] = 1.-1.614188489108937e-01-1.56544537223683037837e-02;
        gp_[27][1] = 1.56544537223683037837e-02;

        refGeoPtr_->addGaussQuadratureRule(this);
    }

    T2_11_1::~T2_11_1()
    {
    }

//---------------------------------------------------------------------------
} // close namespace IntegrationRules

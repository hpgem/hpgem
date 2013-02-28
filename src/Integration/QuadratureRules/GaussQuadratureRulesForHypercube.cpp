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
            return name_;
        }

    unsigned int  Cn4_1_1::order() const
        {
            return 1;
        }

    unsigned int Cn4_1_1::dimension() const
        {
            return 4;
        }

    unsigned int Cn4_1_1::nrOfPoints() const
        {
            return 1;
        }

    NumType Cn4_1_1::weight(unsigned int i) const
        {
            if (i < 1)
                return weight_[i];
            else
                throw name_ + "::weight - wrong index!";
        }

    void Cn4_1_1::getPoint(unsigned int i, PointReference<4>& p) const
        {
            if (i < 1)
                p=gp_[i];
            else
                throw name_ + "::getPoint -  wrong index!";
        }

    Geometry::ReferenceGeometry<4>* Cn4_1_1::forReferenceGeometry() const
        {
            return refGeoPtr_;
        }

    Cn4_1_1::Cn4_1_1()
     : name_("Cn4_1_1"),
       refGeoPtr_(&ReferenceHypercube::Instance())
            {
                weight_[0] = 16.0;
                gp_[0][0] = 0.0;
                gp_[0][1] = 0.0;
                gp_[0][2] = 0.0;
                gp_[0][3] = 0.0;

                refGeoPtr_->addGaussQuadratureRule(this);
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
            return name_;
        }

    unsigned int  Cn4_3_4::order() const
        {
            return 3;
        }

    unsigned int Cn4_3_4::dimension() const
        {
            return 4;
        }

    unsigned int Cn4_3_4::nrOfPoints() const
        {
            return 16;
        }

    NumType Cn4_3_4::weight(unsigned int i) const
        {
            if (i < 16)
                return weight_[i];
            else
                throw name_ + "::weight - wrong index!";
        }

    void Cn4_3_4::getPoint(unsigned int i, PointReference<4>& p) const
        {
            if (i < 16)
                p=gp_[i];
            else
                throw name_ + "::getPoint -  wrong index!";
        }

    Geometry::ReferenceGeometry<4>* Cn4_3_4::forReferenceGeometry() const
        {
            return refGeoPtr_;
        }

    Cn4_3_4::Cn4_3_4()
     : name_("Cn4_3_4"),
       refGeoPtr_(&ReferenceHypercube::Instance())
            {
                weight_[0] = 1.0;
                gp_[0][0] = -sqrt(3.0) / 3.0;
                gp_[0][1] = -sqrt(3.0) / 3.0;
                gp_[0][2] = -sqrt(3.0) / 3.0;
                gp_[0][3] = -sqrt(3.0) / 3.0;

                weight_[1] = 1.0;
                gp_[1][0] = +sqrt(3.0) / 3.0;
                gp_[1][1] = -sqrt(3.0) / 3.0;
                gp_[1][2] = -sqrt(3.0) / 3.0;
                gp_[1][3] = -sqrt(3.0) / 3.0;

                weight_[2] = 1.0;
                gp_[2][0] = -sqrt(3.0) / 3.0;
                gp_[2][1] = +sqrt(3.0) / 3.0;
                gp_[2][2] = -sqrt(3.0) / 3.0;
                gp_[2][3] = -sqrt(3.0) / 3.0;

                weight_[3] = 1.0;
                gp_[3][0] = +sqrt(3.0) / 3.0;
                gp_[3][1] = +sqrt(3.0) / 3.0;
                gp_[3][2] = -sqrt(3.0) / 3.0;
                gp_[3][3] = -sqrt(3.0) / 3.0;

                weight_[4] = 1.0;
                gp_[4][0] = -sqrt(3.0) / 3.0;
                gp_[4][1] = -sqrt(3.0) / 3.0;
                gp_[4][2] = +sqrt(3.0) / 3.0;
                gp_[4][3] = -sqrt(3.0) / 3.0;

                weight_[5] = 1.0;
                gp_[5][0] = +sqrt(3.0) / 3.0;
                gp_[5][1] = -sqrt(3.0) / 3.0;
                gp_[5][2] = +sqrt(3.0) / 3.0;
                gp_[5][3] = -sqrt(3.0) / 3.0;

                weight_[6] = 1.0;
                gp_[6][0] = -sqrt(3.0) / 3.0;
                gp_[6][1] = +sqrt(3.0) / 3.0;
                gp_[6][2] = +sqrt(3.0) / 3.0;
                gp_[6][3] = -sqrt(3.0) / 3.0;

                weight_[7] = 1.0;
                gp_[7][0] = +sqrt(3.0) / 3.0;
                gp_[7][1] = +sqrt(3.0) / 3.0;
                gp_[7][2] = +sqrt(3.0) / 3.0;
                gp_[7][3] = -sqrt(3.0) / 3.0;

                weight_[8] = 1.0;
                gp_[8][0] = -sqrt(3.0) / 3.0;
                gp_[8][1] = -sqrt(3.0) / 3.0;
                gp_[8][2] = -sqrt(3.0) / 3.0;
                gp_[8][3] = +sqrt(3.0) / 3.0;

                weight_[9] = 1.0;
                gp_[9][0] = +sqrt(3.0) / 3.0;
                gp_[9][1] = -sqrt(3.0) / 3.0;
                gp_[9][2] = -sqrt(3.0) / 3.0;
                gp_[9][3] = +sqrt(3.0) / 3.0;

                weight_[10] = 1.0;
                gp_[10][0] = -sqrt(3.0) / 3.0;
                gp_[10][1] = +sqrt(3.0) / 3.0;
                gp_[10][2] = -sqrt(3.0) / 3.0;
                gp_[10][3] = +sqrt(3.0) / 3.0;

                weight_[11] = 1.0;
                gp_[11][0] = +sqrt(3.0) / 3.0;
                gp_[11][1] = +sqrt(3.0) / 3.0;
                gp_[11][2] = -sqrt(3.0) / 3.0;
                gp_[11][3] = +sqrt(3.0) / 3.0;

                weight_[12] = 1.0;
                gp_[12][0] = -sqrt(3.0) / 3.0;
                gp_[12][1] = -sqrt(3.0) / 3.0;
                gp_[12][2] = +sqrt(3.0) / 3.0;
                gp_[12][3] = +sqrt(3.0) / 3.0;

                weight_[13] = 1.0;
                gp_[13][0] = +sqrt(3.0) / 3.0;
                gp_[13][1] = -sqrt(3.0) / 3.0;
                gp_[13][2] = +sqrt(3.0) / 3.0;
                gp_[13][3] = +sqrt(3.0) / 3.0;

                weight_[14] = 1.0;
                gp_[14][0] = -sqrt(3.0) / 3.0;
                gp_[14][1] = +sqrt(3.0) / 3.0;
                gp_[14][2] = +sqrt(3.0) / 3.0;
                gp_[14][3] = +sqrt(3.0) / 3.0;

                weight_[15] = 1.0;
                gp_[15][0] = +sqrt(3.0) / 3.0;
                gp_[15][1] = +sqrt(3.0) / 3.0;
                gp_[15][2] = +sqrt(3.0) / 3.0;
                gp_[15][3] = +sqrt(3.0) / 3.0;

                refGeoPtr_->addGaussQuadratureRule(this);
            }

    Cn4_3_4::~Cn4_3_4()
          {}

    namespace
    {
        const Cn4_3_4& instantiateCn4_3_4 = Cn4_3_4::Instance();
    }

//---------------------------------------------------------------------------
} // close namespace IntegrationRules

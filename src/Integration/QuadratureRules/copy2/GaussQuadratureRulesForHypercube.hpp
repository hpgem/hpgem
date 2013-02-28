//------------------------------------------------------------------------------
// File: GaussQuadratureRulesForHypercube.hpp 
// Headers of Gauss quadrature rules for reference hypercube.
// Lars Pesch, Fri Mar  3 12:59:11 CET 2006
//----
// Modified from original file allGaussQuadratureRules.hpp
// by M.T. Julianto, Wed Feb 25 10:45:06 UTC 2013
//------------------------------------------------------------------------------
#ifndef GaussQuadratureRulesForHypercube_hpp
#define GaussQuadratureRulesForHypercube_hpp
//---------------------------------------------------------------------------
#include "Geometry/PointReference.hpp"
#include "Geometry/ReferenceGeometry.hpp"
#include "Geometry/ReferenceHypercube.hpp"
#include "Integration/GlobalNamespaceIntegration.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"

//---------------------------------------------------------------------------
namespace QuadratureRules
{
    using Integration::NumType;
    using Geometry::PointReference;
    using Geometry::ReferenceGeometry;
    using Geometry::ReferenceHypercube;

//---------------------------------------------------------------------------
    class Cn4_1_1
        : public GaussQuadratureRule<4>
    {
    public:
        static Cn4_1_1& Instance()
            {
                static Cn4_1_1 theInstance;
                return theInstance;
            }

        virtual std::string getName() const;
        virtual unsigned int order() const;
        virtual unsigned int dimension() const;
        virtual unsigned int nrOfPoints() const;
        virtual NumType weight(unsigned int i) const;
        virtual void getPoint(unsigned int i, PointReference<4>& p) const;
        virtual ReferenceGeometry<4>* forReferenceGeometry() const;

    private:
        Cn4_1_1();
        Cn4_1_1(const Cn4_1_1&);
        virtual ~Cn4_1_1();

        const std::string _name;
        NumType _weight[1];
        ReferenceGeometry<4>* const _refGeoPtr;
        PointReference<4> _gp[1];
    };
//---------------------------------------------------------------------------
    class Cn4_3_4
        : public GaussQuadratureRule<4>
    {
    public:
        static Cn4_3_4& Instance()
            {
                static Cn4_3_4 theInstance;
                return theInstance;
            }

        virtual std::string getName() const;
        virtual unsigned int order() const;
        virtual unsigned int dimension() const;
        virtual unsigned int nrOfPoints() const;
        virtual NumType weight(unsigned int i) const;
        virtual void getPoint(unsigned int i, PointReference<4>& p) const;
        virtual ReferenceGeometry<4>* forReferenceGeometry() const;

    private:
        Cn4_3_4();
        Cn4_3_4(const Cn4_3_4&);
        virtual ~Cn4_3_4();

        const std::string _name;
        NumType _weight[16];
        ReferenceGeometry<4>* const _refGeoPtr;
        PointReference<4> _gp[16];
    };

//---------------------------------------------------------------------------
} // close namespace QuadratureRules
#endif

//------------------------------------------------------------------------------
// File: GaussQuadratureRulesForSquare.hpp 
// Headers of Gauss quadrature rules for reference square.
// Lars Pesch, Fri Mar  3 12:59:11 CET 2006
//----
// Modified from original file allGaussQuadratureRules.hpp
// by M.T. Julianto, Wed Feb 25 10:45:06 UTC 2013
//------------------------------------------------------------------------------
#ifndef GaussQuadratureRulesForSquare_hpp
#define GaussQuadratureRulesForSquare_hpp
//---------------------------------------------------------------------------
#include "Geometry/PointReference.hpp"
#include "Geometry/ReferenceGeometry.hpp"
#include "Geometry/ReferenceSquare.hpp"
#include "Integration/GlobalNamespaceIntegration.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"

//---------------------------------------------------------------------------
namespace QuadratureRules
{
    using Integration::NumType;
    using Integration::unsigned int;
    using Geometry::PointReference;
    using Geometry::ReferenceGeometry;
    using Geometry::ReferenceSquare;

//---------------------------------------------------------------------------
    class Cn2_1_1
        : public GaussQuadratureRule<2>
    {
    public:
        static Cn2_1_1& Instance()
            {
                static Cn2_1_1 theInstance;
                return theInstance;
            }

        virtual std::string getName() const;
        virtual unsigned int order() const;
        virtual unsigned int dimension() const;
        virtual unsigned int nrOfPoints() const;
        virtual NumType weight(unsigned int i) const;
        virtual void getPoint(unsigned int i, PointReference<2>& p) const;
        virtual ReferenceGeometry<2>* forReferenceGeometry() const;

    private:
        Cn2_1_1();
        Cn2_1_1(const Cn2_1_1&);
        virtual ~Cn2_1_1();

        const std::string _name;
        NumType _weight[1];
        ReferenceGeometry<2>* const _refGeoPtr;
        PointReference<2> _gp[1];
    };

//---------------------------------------------------------------------------
    class Cn2_3_4
        : public GaussQuadratureRule<2>
    {
    public:
        static Cn2_3_4& Instance()
            {
                static Cn2_3_4 theInstance;
                return theInstance;
            }

        virtual std::string getName() const;
        virtual unsigned int order() const;
        virtual unsigned int dimension() const;
        virtual unsigned int nrOfPoints() const;
        virtual NumType weight(unsigned int i) const;
        virtual void getPoint(unsigned int i, PointReference<2>& p) const;
        virtual ReferenceGeometry<2>* forReferenceGeometry() const;

    private:
        Cn2_3_4();
        Cn2_3_4(const Cn2_3_4&);
        virtual ~Cn2_3_4();

        const std::string _name;
        NumType _weight[4];
        ReferenceGeometry<2>* const _refGeoPtr;
        PointReference<2> _gp[4];
    };

//---------------------------------------------------------------------------
    class Cn2_5_9
        : public GaussQuadratureRule<2>
    {
    public:
        static Cn2_5_9& Instance()
            {
                static Cn2_5_9 theInstance;
                return theInstance;
            }

        virtual std::string getName() const;
        virtual unsigned int order() const;
        virtual unsigned int dimension() const;
        virtual unsigned int nrOfPoints() const;
        virtual NumType weight(unsigned int i) const;
        virtual void getPoint(unsigned int i, PointReference<2>& p) const;
        virtual ReferenceGeometry<2>* forReferenceGeometry() const;

    private:
        Cn2_5_9();
        Cn2_5_9(const Cn2_5_9&);
        virtual ~Cn2_5_9();

        const std::string _name;
        NumType _weight[9];
        ReferenceGeometry<2>* const _refGeoPtr;
        PointReference<2> _gp[9];
    };

//---------------------------------------------------------------------------
    class C2_7_4
        : public GaussQuadratureRule<2>
    {
    public:
        static C2_7_4& Instance()
            {
                static C2_7_4 theInstance;
                return theInstance;
            }

        virtual std::string getName() const;
        virtual unsigned int order() const;
        virtual unsigned int dimension() const;
        virtual unsigned int nrOfPoints() const;
        virtual NumType weight(unsigned int i) const;
        virtual void getPoint(unsigned int i, PointReference<2>& p) const;
        virtual ReferenceGeometry<2>* forReferenceGeometry() const;

    private:
        C2_7_4();
        C2_7_4(const C2_7_4&);
        virtual ~C2_7_4();

        const std::string _name;
        NumType _weight[16];
        ReferenceGeometry<2>* const _refGeoPtr;
        PointReference<2> _gp[16];
    };

//---------------------------------------------------------------------------
} // close namespace QuadratureRules
#endif

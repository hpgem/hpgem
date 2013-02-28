//------------------------------------------------------------------------------
// File: GaussQuadratureRulesForLine.hpp 
// Headers of Gauss quadrature rules for reference line.
// Lars Pesch, Fri Mar  3 12:59:11 CET 2006
//----
// Modified from original file allGaussQuadratureRules.hpp
// by M.T. Julianto, Wed Feb 25 10:45:06 UTC 2013
//------------------------------------------------------------------------------
#ifndef GaussQuadratureRulesForLine_hpp
#define GaussQuadratureRulesForLine_hpp
//---------------------------------------------------------------------------
#include "Geometry/PointReference.hpp"
#include "Geometry/ReferenceGeometry.hpp"
#include "Geometry/ReferenceLine.hpp"
#include "Integration/GlobalNamespaceIntegration.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"

//---------------------------------------------------------------------------
namespace QuadratureRules
{
    using Integration::NumType;
    using Geometry::PointReference;
    using Geometry::ReferenceGeometry;
    using Geometry::ReferenceLine;

//---------------------------------------------------------------------------
    class Cn1_1_1
        : public GaussQuadratureRule<1>
    {
    public:
        static Cn1_1_1& Instance()
            {
                static Cn1_1_1 theInstance;
                return theInstance;
            }

        virtual std::string getName() const;
        virtual unsigned int order() const;
        virtual unsigned int dimension() const;
        virtual unsigned int nrOfPoints() const;
        virtual NumType weight(unsigned int i) const;
        virtual void getPoint(unsigned int i, PointReference<1>& p) const;
        virtual ReferenceGeometry<1>* forReferenceGeometry() const;

    private:
        Cn1_1_1();
        Cn1_1_1(const Cn1_1_1&);
        virtual ~Cn1_1_1();

        const std::string _name;
        NumType _weight[1];
        ReferenceGeometry<1>* const _refGeoPtr;
        PointReference<1> _gp[1];
    };

//---------------------------------------------------------------------------
    class Cn1_3_4
        : public GaussQuadratureRule<1>
    {
    public:
        static Cn1_3_4& Instance()
            {
                static Cn1_3_4 theInstance;
                return theInstance;
            }

        virtual std::string getName() const;
        virtual unsigned int order() const;
        virtual unsigned int dimension() const;
        virtual unsigned int nrOfPoints() const;
        virtual NumType weight(unsigned int i) const;
        virtual void getPoint(unsigned int i, PointReference<1>& p) const;
        virtual ReferenceGeometry<1>* forReferenceGeometry() const;

    private:
        Cn1_3_4();
        Cn1_3_4(const Cn1_3_4&);
        virtual ~Cn1_3_4();

        const std::string _name;
        NumType _weight[2];
        ReferenceGeometry<1>* const _refGeoPtr;
        PointReference<1> _gp[2];
    };

//---------------------------------------------------------------------------
    class Cn1_5_9
        : public GaussQuadratureRule<1>
    {
    public:
        static Cn1_5_9& Instance()
            {
                static Cn1_5_9 theInstance;
                return theInstance;
            }

        virtual std::string getName() const;
        virtual unsigned int order() const;
        virtual unsigned int dimension() const;
        virtual unsigned int nrOfPoints() const;
        virtual NumType weight(unsigned int i) const;
        virtual void getPoint(unsigned int i, PointReference<1>& p) const;
        virtual ReferenceGeometry<1>* forReferenceGeometry() const;

    private:
        Cn1_5_9();
        Cn1_5_9(const Cn1_5_9&);
        virtual ~Cn1_5_9();

        const std::string _name;
        NumType _weight[3];
        ReferenceGeometry<1>* const _refGeoPtr;
        PointReference<1> _gp[3];
    };

//---------------------------------------------------------------------------
    class C1_7_x
        : public GaussQuadratureRule<1>
    {
    public:
        static C1_7_x& Instance()
            {
                static C1_7_x theInstance;
                return theInstance;
            }

        virtual std::string getName() const;
        virtual unsigned int order() const;
        virtual unsigned int dimension() const;
        virtual unsigned int nrOfPoints() const;
        virtual NumType weight(unsigned int i) const;
        virtual void getPoint(unsigned int i, PointReference<1>& p) const;
        virtual ReferenceGeometry<1>* forReferenceGeometry() const;

    private:
        C1_7_x();
        C1_7_x(const C1_7_x&);
        virtual ~C1_7_x();

        const std::string _name;
        NumType _weight[4];
        ReferenceGeometry<1>* const _refGeoPtr;
        PointReference<1> _gp[4];
    };

//---------------------------------------------------------------------------
} // close namespace QuadratureRules
#endif

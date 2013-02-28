//------------------------------------------------------------------------------
// File: GaussQuadratureRulesForTetrahedron.hpp 
// Headers of Gauss quadrature rules for reference tetrahedron.
// Lars Pesch, Fri Mar  3 12:59:11 CET 2006
//----
// Modified from original file allGaussQuadratureRules.hpp
// by M.T. Julianto, Wed Feb 25 10:45:06 UTC 2013
//------------------------------------------------------------------------------
#ifndef GaussQuadratureRulesForTetrahedron_hpp
#define GaussQuadratureRulesForTetrahedron_hpp
//---------------------------------------------------------------------------
#include "Geometry/PointReference.hpp"
#include "Geometry/ReferenceGeometry.hpp"
#include "Geometry/ReferenceTetrahedron.hpp"
#include "Integration/GlobalNamespaceIntegration.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"

//---------------------------------------------------------------------------
namespace QuadratureRules
{
    using Integration::NumType;
    using Integration::unsigned int;
    using Geometry::PointReference;
    using Geometry::ReferenceGeometry;
    using Geometry::ReferenceTetrahedron;

//---------------------------------------------------------------------------
    class Tn3_1_1
        : public GaussQuadratureRule<3>
    {
    public:
        static Tn3_1_1& Instance()
            {
                static Tn3_1_1 theInstance;
                return theInstance;
            }

        virtual std::string getName() const;
        virtual unsigned int order() const;
        virtual unsigned int dimension() const;
        virtual unsigned int nrOfPoints() const;
        virtual NumType weight(unsigned int i) const;
        virtual void getPoint(unsigned int i, PointReference<3>& p) const;
        virtual ReferenceGeometry<3>* forReferenceGeometry() const;

    private:
        Tn3_1_1();
        Tn3_1_1(const Tn3_1_1&);
        virtual ~Tn3_1_1();

        const std::string _name;
        NumType _weight[1];
        ReferenceGeometry<3>* const _refGeoPtr;
        PointReference<3> _gp[1];
    };

//---------------------------------------------------------------------------
    class Tn3_2_1
        : public GaussQuadratureRule<3>
    {
    public:
        static Tn3_2_1& Instance()
            {
                static Tn3_2_1 theInstance;
                return theInstance;
            }

        virtual std::string getName() const;
        virtual unsigned int order() const;
        virtual unsigned int dimension() const;
        virtual unsigned int nrOfPoints() const;
        virtual NumType weight(unsigned int i) const;
        virtual void getPoint(unsigned int i, PointReference<3>& p) const;
        virtual ReferenceGeometry<3>* forReferenceGeometry() const;

    private:
        Tn3_2_1();
        Tn3_2_1(const Tn3_2_1&);
        virtual ~Tn3_2_1();

        const std::string _name;
        NumType _weight[4];
        ReferenceGeometry<3>* const _refGeoPtr;
        PointReference<3> _gp[4];
    };

//---------------------------------------------------------------------------
    class Tn3_3_1
        : public GaussQuadratureRule<3>
    {
    public:
        static Tn3_3_1& Instance()
            {
                static Tn3_3_1 theInstance;
                return theInstance;
            }

        virtual std::string getName() const;
        virtual unsigned int order() const;
        virtual unsigned int dimension() const;
        virtual unsigned int nrOfPoints() const;
        virtual NumType weight(unsigned int i) const;
        virtual void getPoint(unsigned int i, PointReference<3>& p) const;
        virtual ReferenceGeometry<3>* forReferenceGeometry() const;

    private:
        Tn3_3_1();
        Tn3_3_1(const Tn3_3_1&);
        virtual ~Tn3_3_1();

        const std::string _name;
        NumType _weight[5];
        ReferenceGeometry<3>* const _refGeoPtr;
        PointReference<3> _gp[5];
    };

//---------------------------------------------------------------------------
    class Tn3_4_1
        : public GaussQuadratureRule<3>
    {
    public:
        static Tn3_4_1& Instance()
            {
                static Tn3_4_1 theInstance;
                return theInstance;
            }

        virtual std::string getName() const;
        virtual unsigned int order() const;
        virtual unsigned int dimension() const;
        virtual unsigned int nrOfPoints() const;
        virtual NumType weight(unsigned int i) const;
        virtual void getPoint(unsigned int i, PointReference<3>& p) const;
        virtual ReferenceGeometry<3>* forReferenceGeometry() const;

    private:
        Tn3_4_1();
        Tn3_4_1(const Tn3_4_1&);
        virtual ~Tn3_4_1();

        const std::string _name;
        NumType _weight[11];
        ReferenceGeometry<3>* const _refGeoPtr;
        PointReference<3> _gp[11];
    };

//---------------------------------------------------------------------------
    class T3_5_1
        : public GaussQuadratureRule<3>
    {
    public:
        static T3_5_1& Instance()
            {
                static T3_5_1 theInstance;
                return theInstance;
            }

        virtual std::string getName() const;
        virtual unsigned int order() const;
        virtual unsigned int dimension() const;
        virtual unsigned int nrOfPoints() const;
        virtual NumType weight(unsigned int i) const;
        virtual void getPoint(unsigned int i, PointReference<3>& p) const;
        virtual ReferenceGeometry<3>* forReferenceGeometry() const;

    private:
        T3_5_1();
        T3_5_1(const T3_5_1&);
        virtual ~T3_5_1();

        const std::string _name;
        NumType _weight[14];
        ReferenceGeometry<3>* const _refGeoPtr;
        PointReference<3> _gp[14];
    };

//---------------------------------------------------------------------------
    class T3_6_1
        : public GaussQuadratureRule<3>
    {
    public:
        static T3_6_1& Instance()
            {
                static T3_6_1 theInstance;
                return theInstance;
            }

        virtual std::string getName() const;
        virtual unsigned int order() const;
        virtual unsigned int dimension() const;
        virtual unsigned int nrOfPoints() const;
        virtual NumType weight(unsigned int i) const;
        virtual void getPoint(unsigned int i, PointReference<3>& p) const;
        virtual ReferenceGeometry<3>* forReferenceGeometry() const;

    private:
        T3_6_1();
        T3_6_1(const T3_6_1&);
        virtual ~T3_6_1();

        const std::string _name;
        NumType _weight[24];
        ReferenceGeometry<3>* const _refGeoPtr;
        PointReference<3> _gp[24];
    };

//---------------------------------------------------------------------------
    class T3_7_1
        : public GaussQuadratureRule<3>
    {
    public:
        static T3_7_1& Instance()
            {
                static T3_7_1 theInstance;
                return theInstance;
            }

        virtual std::string getName() const;
        virtual unsigned int order() const;
        virtual unsigned int dimension() const;
        virtual unsigned int nrOfPoints() const;
        virtual NumType weight(unsigned int i) const;
        virtual void getPoint(unsigned int i, PointReference<3>& p) const;
        virtual ReferenceGeometry<3>* forReferenceGeometry() const;

    private:
        T3_7_1();
        T3_7_1(const T3_7_1&);
        virtual ~T3_7_1();

        const std::string _name;
        NumType _weight[31];
        ReferenceGeometry<3>* const _refGeoPtr;
        PointReference<3> _gp[31];
    };

//---------------------------------------------------------------------------
    class T3_8_1
        : public GaussQuadratureRule<3>
    {
    public:
        static T3_8_1& Instance()
            {
                static T3_8_1 theInstance;
                return theInstance;
            }

        virtual std::string getName() const;
        virtual unsigned int order() const;
        virtual unsigned int dimension() const;
        virtual unsigned int nrOfPoints() const;
        virtual NumType weight(unsigned int i) const;
        virtual void getPoint(unsigned int i, PointReference<3>& p) const;
        virtual ReferenceGeometry<3>* forReferenceGeometry() const;

    private:
        T3_8_1();
        T3_8_1(const T3_8_1&);
        virtual ~T3_8_1();

        const std::string _name;
        NumType _weight[43];
        ReferenceGeometry<3>* const _refGeoPtr;
        PointReference<3> _gp[43];
    };

//---------------------------------------------------------------------------
    class T3_9_1
        : public GaussQuadratureRule<3>
    {
    public:
        static T3_9_1& Instance()
            {
                static T3_9_1 theInstance;
                return theInstance;
            }

        virtual std::string getName() const;
        virtual unsigned int order() const;
        virtual unsigned int dimension() const;
        virtual unsigned int nrOfPoints() const;
        virtual NumType weight(unsigned int i) const;
        virtual void getPoint(unsigned int i, PointReference<3>& p) const;
        virtual ReferenceGeometry<3>* forReferenceGeometry() const;

    private:
        T3_9_1();
        T3_9_1(const T3_9_1&);
        virtual ~T3_9_1();

        const std::string _name;
        NumType _weight[53];
        ReferenceGeometry<3>* const _refGeoPtr;
        PointReference<3> _gp[53];
    };

//---------------------------------------------------------------------------
    class T3_10_1
        : public GaussQuadratureRule<3>
    {
    public:
        static T3_10_1& Instance()
            {
                static T3_10_1 theInstance;
                return theInstance;
            }

        virtual std::string getName() const;
        virtual unsigned int order() const;
        virtual unsigned int dimension() const;
        virtual unsigned int nrOfPoints() const;
        virtual NumType weight(unsigned int i) const;
        virtual void getPoint(unsigned int i, PointReference<3>& p) const;
        virtual ReferenceGeometry<3>* forReferenceGeometry() const;

    private:
        T3_10_1();
        T3_10_1(const T3_10_1&);
        virtual ~T3_10_1();

        const std::string _name;
        NumType _weight[126];
        ReferenceGeometry<3>* const _refGeoPtr;
        PointReference<3> _gp[126];
    };

//---------------------------------------------------------------------------
} // close namespace QuadratureRules
#endif

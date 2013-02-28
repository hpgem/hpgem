//------------------------------------------------------------------------------
// File: GaussQuadratureRulesForTriangle.hpp 
// Headers of Gauss quadrature rules for reference triangle.
// Lars Pesch, Fri Mar  3 12:59:11 CET 2006
//----
// Modified from original file allGaussQuadratureRules.hpp
// by M.T. Julianto, Wed Feb 25 10:45:06 UTC 2013
//------------------------------------------------------------------------------
#ifndef GaussQuadratureRulesForTriangle_hpp
#define GaussQuadratureRulesForTriangle_hpp
//---------------------------------------------------------------------------
#include "Geometry/PointReference.hpp"
#include "Geometry/ReferenceGeometry.hpp"
#include "Geometry/ReferenceTriangle.hpp"
#include "Integration/GlobalNamespaceIntegration.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"

//---------------------------------------------------------------------------
namespace QuadratureRules
{
    using Integration::NumType;
    using Geometry::PointReference;
    using Geometry::ReferenceGeometry;
    using Geometry::ReferenceTriangle;

//---------------------------------------------------------------------------
    class Tn2_1_1
        : public GaussQuadratureRule<2>
    {
    public:
        static Tn2_1_1& Instance()
            {
                static Tn2_1_1 theInstance;
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
        Tn2_1_1();
        Tn2_1_1(const Tn2_1_1&);
        virtual ~Tn2_1_1();

        const std::string name_;
        NumType weight_[1];
        ReferenceGeometry<2>* const refGeoPtr_;
        PointReference<2> gp_[1];
    };

//---------------------------------------------------------------------------
    class Tn2_2_1
        : public GaussQuadratureRule<2>
    {
    public:
        static Tn2_2_1& Instance()
            {
                static Tn2_2_1 theInstance;
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
        Tn2_2_1();
        Tn2_2_1(const Tn2_2_1&);
        virtual ~Tn2_2_1();

        const std::string name_;
        NumType weight_[3];
        ReferenceGeometry<2>* const refGeoPtr_;
        PointReference<2> gp_[3];
    };

//---------------------------------------------------------------------------
    class Tn2_3_1
        : public GaussQuadratureRule<2>
    {
    public:
        static Tn2_3_1& Instance()
            {
                static Tn2_3_1 theInstance;
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
        Tn2_3_1();
        Tn2_3_1(const Tn2_3_1&);
        virtual ~Tn2_3_1();

        const std::string name_;
        NumType weight_[4];
        ReferenceGeometry<2>* const refGeoPtr_;
        PointReference<2> gp_[4];
    };

//---------------------------------------------------------------------------
    class Tn2_4_1
        : public GaussQuadratureRule<2>
    {
    public:
        static Tn2_4_1& Instance()
            {
                static Tn2_4_1 theInstance;
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
        Tn2_4_1();
        Tn2_4_1(const Tn2_4_1&);
        virtual ~Tn2_4_1();

        const std::string name_;
        NumType weight_[6];
        ReferenceGeometry<2>* const refGeoPtr_;
        PointReference<2> gp_[6];
    };

//---------------------------------------------------------------------------
    class T2_5_1
        : public GaussQuadratureRule<2>
    {
    public:
        static T2_5_1& Instance()
            {
                static T2_5_1 theInstance;
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
        T2_5_1();
        T2_5_1(const T2_5_1&);
        virtual ~T2_5_1();

        const std::string name_;
        NumType weight_[7];
        ReferenceGeometry<2>* const refGeoPtr_;
        PointReference<2> gp_[7];
    };

//---------------------------------------------------------------------------
    class T2_6_1
        : public GaussQuadratureRule<2>
    {
    public:
        static T2_6_1& Instance()
            {
                static T2_6_1 theInstance;
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
        T2_6_1();
        T2_6_1(const T2_6_1&);
        virtual ~T2_6_1();

        const std::string name_;
        NumType weight_[12];
        ReferenceGeometry<2>* const refGeoPtr_;
        PointReference<2> gp_[12];
    };

//---------------------------------------------------------------------------
    class T2_7_1
        : public GaussQuadratureRule<2>
    {
    public:
        static T2_7_1& Instance()
            {
                static T2_7_1 theInstance;
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
        T2_7_1();
        T2_7_1(const T2_7_1&);
        virtual ~T2_7_1();

        const std::string name_;
        NumType weight_[13];
        ReferenceGeometry<2>* const refGeoPtr_;
        PointReference<2> gp_[13];
    };

//---------------------------------------------------------------------------
    class T2_8_1
        : public GaussQuadratureRule<2>
    {
    public:
        static T2_8_1& Instance()
            {
                static T2_8_1 theInstance;
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
        T2_8_1();
        T2_8_1(const T2_8_1&);
        virtual ~T2_8_1();

        const std::string name_;
        NumType weight_[16];
        ReferenceGeometry<2>* const refGeoPtr_;
        PointReference<2> gp_[16];
    };

//---------------------------------------------------------------------------
    class T2_9_1
        : public GaussQuadratureRule<2>
    {
    public:
        static T2_9_1& Instance()
            {
                static T2_9_1 theInstance;
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
        T2_9_1();
        T2_9_1(const T2_9_1&);
        virtual ~T2_9_1();

        const std::string name_;
        NumType weight_[19];
        ReferenceGeometry<2>* const refGeoPtr_;
        PointReference<2> gp_[19];
    };

//---------------------------------------------------------------------------
    class T2_10_1
        : public GaussQuadratureRule<2>
    {
    public:
        static T2_10_1& Instance()
            {
                static T2_10_1 theInstance;
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
        T2_10_1();
        T2_10_1(const T2_10_1&);
        virtual ~T2_10_1();

        const std::string name_;
        NumType weight_[25];
        ReferenceGeometry<2>* const refGeoPtr_;
        PointReference<2> gp_[25];
    };

//---------------------------------------------------------------------------
} // close namespace QuadratureRules
#endif

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
#include "Integration/GlobalNamespaceIntegration.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"

//---------------------------------------------------------------------------
namespace QuadratureRules
{
    using Geometry::PointReference;
    using Geometry::ReferenceGeometry;

//---------------------------------------------------------------------------
    class Tn2_1_1: public GaussQuadratureRule
    {
    public:
        typedef ReferenceGeometry    ReferenceGeometryT;
        typedef PointReference      PointReferenceT;
    public:
        static Tn2_1_1& Instance()
            {
                static Tn2_1_1 theInstance;
                return theInstance;
            }

        virtual std::string             getName() const;
        virtual unsigned int            order() const;
        virtual unsigned int            dimension() const;
        virtual unsigned int            nrOfPoints() const;
        virtual double                  weight(unsigned int i) const;
        virtual void                    getPoint(unsigned int i, PointReferenceT& p) const;
        virtual ReferenceGeometryT*     forReferenceGeometry() const;

    private:
        Tn2_1_1();
        Tn2_1_1(const Tn2_1_1&);
        virtual ~Tn2_1_1();
    private:
        const std::string               name_;
        double                          weight_[1];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                 gp_;
    };

//---------------------------------------------------------------------------
    class Tn2_2_1: public GaussQuadratureRule
    {
    public:
        typedef ReferenceGeometry    ReferenceGeometryT;
        typedef PointReference       PointReferenceT;
    public:
        static Tn2_2_1& Instance()
            {
                static Tn2_2_1 theInstance;
                return theInstance;
            }

        virtual std::string             getName() const;
        virtual unsigned int            order() const;
        virtual unsigned int            dimension() const;
        virtual unsigned int            nrOfPoints() const;
        virtual double                  weight(unsigned int i) const;
        virtual void                    getPoint(unsigned int i, PointReferenceT& p) const;
        virtual ReferenceGeometryT*     forReferenceGeometry() const;

    private:
        Tn2_2_1();
        Tn2_2_1(const Tn2_2_1&);
        virtual ~Tn2_2_1();
    private:
        const std::string               name_;
        double                          weight_[3];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };

//---------------------------------------------------------------------------
    class Tn2_3_1: public GaussQuadratureRule
    {
    public:
        typedef ReferenceGeometry    ReferenceGeometryT;
        typedef PointReference       PointReferenceT;
    public:
        static Tn2_3_1& Instance()
            {
                static Tn2_3_1 theInstance;
                return theInstance;
            }

        virtual std::string             getName() const;
        virtual unsigned int            order() const;
        virtual unsigned int            dimension() const;
        virtual unsigned int            nrOfPoints() const;
        virtual double                  weight(unsigned int i) const;
        virtual void                    getPoint(unsigned int i, PointReferenceT& p) const;
        virtual ReferenceGeometryT*     forReferenceGeometry() const;

    private:
        Tn2_3_1();
        Tn2_3_1(const Tn2_3_1&);
        virtual ~Tn2_3_1();
    private:
        const std::string               name_;
        double                          weight_[4];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };

//---------------------------------------------------------------------------
    class Tn2_4_1: public GaussQuadratureRule
    {
    public:
        typedef ReferenceGeometry    ReferenceGeometryT;
        typedef PointReference       PointReferenceT;
    public:
        static Tn2_4_1& Instance()
            {
                static Tn2_4_1 theInstance;
                return theInstance;
            }

        virtual std::string             getName() const;
        virtual unsigned int            order() const;
        virtual unsigned int            dimension() const;
        virtual unsigned int            nrOfPoints() const;
        virtual double                  weight(unsigned int i) const;
        virtual void                    getPoint(unsigned int i, PointReferenceT& p) const;
        virtual ReferenceGeometryT*     forReferenceGeometry() const;

    private:
        Tn2_4_1();
        Tn2_4_1(const Tn2_4_1&);
        virtual ~Tn2_4_1();
    private:
        const std::string               name_;
        double                          weight_[6];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };

//---------------------------------------------------------------------------
    class T2_5_1: public GaussQuadratureRule
    {
    public:
        typedef ReferenceGeometry    ReferenceGeometryT;
        typedef PointReference       PointReferenceT;
    public:
        static T2_5_1& Instance()
            {
                static T2_5_1 theInstance;
                return theInstance;
            }

        virtual std::string             getName() const;
        virtual unsigned int            order() const;
        virtual unsigned int            dimension() const;
        virtual unsigned int            nrOfPoints() const;
        virtual double                  weight(unsigned int i) const;
        virtual void                    getPoint(unsigned int i, PointReferenceT& p) const;
        virtual ReferenceGeometryT*     forReferenceGeometry() const;

    private:
        T2_5_1();
        T2_5_1(const T2_5_1&);
        virtual ~T2_5_1();

        const std::string               name_;
        double                          weight_[7];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };

//---------------------------------------------------------------------------
    class T2_6_1: public GaussQuadratureRule
    {
    public:
        typedef ReferenceGeometry    ReferenceGeometryT;
        typedef PointReference       PointReferenceT;
    public:
        static T2_6_1& Instance()
            {
                static T2_6_1 theInstance;
                return theInstance;
            }

        virtual std::string             getName() const;
        virtual unsigned int            order() const;
        virtual unsigned int            dimension() const;
        virtual unsigned int            nrOfPoints() const;
        virtual double                  weight(unsigned int i) const;
        virtual void                    getPoint(unsigned int i, PointReferenceT& p) const;
        virtual ReferenceGeometryT*     forReferenceGeometry() const;

    private:
        T2_6_1();
        T2_6_1(const T2_6_1&);
        virtual ~T2_6_1();
    private:
        const std::string               name_;
        double                          weight_[12];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };

//---------------------------------------------------------------------------
    class T2_7_1: public GaussQuadratureRule
    {
    public:
        typedef ReferenceGeometry    ReferenceGeometryT;
        typedef PointReference       PointReferenceT;
    public:
        static T2_7_1& Instance()
            {
                static T2_7_1 theInstance;
                return theInstance;
            }

        virtual std::string             getName() const;
        virtual unsigned int            order() const;
        virtual unsigned int            dimension() const;
        virtual unsigned int            nrOfPoints() const;
        virtual double                  weight(unsigned int i) const;
        virtual void                    getPoint(unsigned int i, PointReferenceT& p) const;
        virtual ReferenceGeometryT*     forReferenceGeometry() const;

    private:
        T2_7_1();
        T2_7_1(const T2_7_1&);
        virtual ~T2_7_1();
    private:
        const std::string               name_;
        double                          weight_[13];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };

//---------------------------------------------------------------------------
    class T2_8_1: public GaussQuadratureRule
    {
    public:
        typedef ReferenceGeometry    ReferenceGeometryT;
        typedef PointReference       PointReferenceT;
    public:
        static T2_8_1& Instance()
            {
                static T2_8_1 theInstance;
                return theInstance;
            }

        virtual std::string             getName() const;
        virtual unsigned int            order() const;
        virtual unsigned int            dimension() const;
        virtual unsigned int            nrOfPoints() const;
        virtual double                  weight(unsigned int i) const;
        virtual void                    getPoint(unsigned int i, PointReferenceT& p) const;
        virtual ReferenceGeometryT*     forReferenceGeometry() const;

    private:
        T2_8_1();
        T2_8_1(const T2_8_1&);
        virtual ~T2_8_1();
    private:
        const std::string               name_;
        double                          weight_[16];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };

//---------------------------------------------------------------------------
    class T2_9_1: public GaussQuadratureRule
    {
    public:
        typedef ReferenceGeometry    ReferenceGeometryT;
        typedef PointReference       PointReferenceT;
    public:
        static T2_9_1& Instance()
            {
                static T2_9_1 theInstance;
                return theInstance;
            }

        virtual std::string             getName() const;
        virtual unsigned int            order() const;
        virtual unsigned int            dimension() const;
        virtual unsigned int            nrOfPoints() const;
        virtual double                  weight(unsigned int i) const;
        virtual void                    getPoint(unsigned int i, PointReferenceT& p) const;
        virtual ReferenceGeometryT*     forReferenceGeometry() const;

    private:
        T2_9_1();
        T2_9_1(const T2_9_1&);
        virtual ~T2_9_1();
    private:
        const std::string               name_;
        double                          weight_[19];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };

//---------------------------------------------------------------------------
    class T2_10_1: public GaussQuadratureRule
    {
    public:
        typedef ReferenceGeometry    ReferenceGeometryT;
        typedef PointReference       PointReferenceT;
    public:
        static T2_10_1& Instance()
            {
                static T2_10_1 theInstance;
                return theInstance;
            }

        virtual std::string             getName() const;
        virtual unsigned int            order() const;
        virtual unsigned int            dimension() const;
        virtual unsigned int            nrOfPoints() const;
        virtual double                  weight(unsigned int i) const;
        virtual void                    getPoint(unsigned int i, PointReferenceT& p) const;
        virtual ReferenceGeometryT*     forReferenceGeometry() const;

    private:
        T2_10_1();
        T2_10_1(const T2_10_1&);
        virtual ~T2_10_1();
    private:
        const std::string               name_;
        double                          weight_[25];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };

//---------------------------------------------------------------------------
    class T2_11_1: public GaussQuadratureRule
    {
    public:
        typedef ReferenceGeometry    ReferenceGeometryT;
        typedef PointReference       PointReferenceT;
    public:
        static T2_11_1& Instance()
            {
                static T2_11_1 theInstance;
                return theInstance;
            }

        virtual std::string             getName() const;
        virtual unsigned int            order() const;
        virtual unsigned int            dimension() const;
        virtual unsigned int            nrOfPoints() const;
        virtual double                  weight(unsigned int i) const;
        virtual void                    getPoint(unsigned int i, PointReferenceT& p) const;
        virtual ReferenceGeometryT*     forReferenceGeometry() const;

    private:
        T2_11_1();
        T2_11_1(const T2_11_1&);
        virtual ~T2_11_1();
    private:
        const std::string               name_;
        double                          weight_[28];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };

//---------------------------------------------------------------------------
} // close namespace QuadratureRules
#endif

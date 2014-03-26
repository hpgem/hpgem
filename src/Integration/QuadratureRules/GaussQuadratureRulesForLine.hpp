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
#include "Integration/GlobalNamespaceIntegration.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"

//---------------------------------------------------------------------------
namespace QuadratureRules
{
    using Geometry::PointReference;
    using Geometry::ReferenceGeometry;
  
//---------------------------------------------------------------------------
    class Cn1_1_1: public GaussQuadratureRule
    {
    public:
        typedef PointReference       PointReferenceT;
        typedef ReferenceGeometry    ReferenceGeometryT;
    public:
        static Cn1_1_1& Instance()
        {
            static Cn1_1_1 theInstance;
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
        Cn1_1_1();
        Cn1_1_1(const Cn1_1_1&);
        virtual ~Cn1_1_1();
    private:
        const std::string           name_;
        double                      weight_[1];
        ReferenceGeometryT* const   refGeoPtr_;
        std::vector<PointReferenceT>             gp_;
    };

//---------------------------------------------------------------------------
    class Cn1_3_4: public GaussQuadratureRule
    {
    public:
        typedef PointReference       PointReferenceT;
        typedef ReferenceGeometry    ReferenceGeometryT;
    public:
        static Cn1_3_4& Instance()
        {
            static Cn1_3_4 theInstance;
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
        Cn1_3_4();
        Cn1_3_4(const Cn1_3_4&);
        virtual ~Cn1_3_4();
    private:
        const std::string               name_;
        double                          weight_[2];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                 gp_;
    };

//---------------------------------------------------------------------------
    class Cn1_5_9: public GaussQuadratureRule
    {
    public:
        typedef PointReference       PointReferenceT;
        typedef ReferenceGeometry    ReferenceGeometryT;
    public:
        static Cn1_5_9& Instance()
        {
            static Cn1_5_9 theInstance;
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
        Cn1_5_9();
        Cn1_5_9(const Cn1_5_9&);
        virtual ~Cn1_5_9();

    private:
        const std::string           name_;
        double                      weight_[3];
        ReferenceGeometryT* const   refGeoPtr_;
        std::vector<PointReferenceT>             gp_;
    };

//---------------------------------------------------------------------------
    class C1_7_x: public GaussQuadratureRule
    {
    public:
        typedef PointReference       PointReferenceT;
        typedef ReferenceGeometry    ReferenceGeometryT;
    public:
        static C1_7_x& Instance()
        {
            static C1_7_x theInstance;
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
        C1_7_x();
        C1_7_x(const C1_7_x&);
        virtual ~C1_7_x();
    private:
        const std::string               name_;
        double                          weight_[4];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                 gp_;
    };

//---------------------------------------------------------------------------
    class C1_9_25: public GaussQuadratureRule///What is the magic number at the end?? -FB
    {
    public:
        typedef PointReference       PointReferenceT;
        typedef ReferenceGeometry    ReferenceGeometryT;
    public:
        static C1_9_25& Instance()
        {
            static C1_9_25 theInstance;
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
        C1_9_25();
        C1_9_25(const C1_9_25&);
        virtual ~C1_9_25();
    private:
        const std::string               name_;
        double                          weight_[5];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                 gp_;
    };

//---------------------------------------------------------------------------
    class C1_11_36: public GaussQuadratureRule///What is the magic number at the end?? -FB
    {
    public:
        typedef PointReference       PointReferenceT;
        typedef ReferenceGeometry    ReferenceGeometryT;
    public:
        static C1_11_36& Instance()
        {
            static C1_11_36 theInstance;
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
        C1_11_36();
        C1_11_36(const C1_11_36&);
        virtual ~C1_11_36();
    private:
        const std::string               name_;
        double                          weight_[6];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                 gp_;
    };

//---------------------------------------------------------------------------
} // close namespace QuadratureRules
#endif

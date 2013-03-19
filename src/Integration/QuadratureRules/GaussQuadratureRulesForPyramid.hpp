//------------------------------------------------------------------------------
// File: GaussQuadratureRulesForPyramid.hpp 
// Headers of Gauss quadrature rules for reference Pyramid.
// Lars Pesch, Fri Mar  3 12:59:11 CET 2006
//----
// Modified from original file allGaussQuadratureRules.hpp
// by M.T. Julianto, Wed Feb 25 10:45:06 UTC 2013
//------------------------------------------------------------------------------
#ifndef GaussQuadratureRulesForPyramid_hpp
#define GaussQuadratureRulesForPyramid_hpp
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
    class Pyramid_1_1: public GaussQuadratureRule<3>
    {
    public:
        typedef PointReference<3>       PointReferenceT;
        typedef ReferenceGeometry<3>    ReferenceGeometryT;
    public:
        static Pyramid_1_1& Instance()
            {
                static Pyramid_1_1 theInstance;
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
        Pyramid_1_1();
        Pyramid_1_1(const Pyramid_1_1&);
        virtual ~Pyramid_1_1();
    private:

        const std::string               name_;
        double                          weight_[4];
        ReferenceGeometryT* const       refGeoPtr_;
        PointReferenceT                 gp_[4];
    };

//---------------------------------------------------------------------------
    class Pyramid_3_1: public GaussQuadratureRule<3>
    {
    public:
        typedef PointReference<3>       PointReferenceT;
        typedef ReferenceGeometry<3>    ReferenceGeometryT;
    public:
        static Pyramid_3_1& Instance()
            {
                static Pyramid_3_1 theInstance;
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
        Pyramid_3_1();
        Pyramid_3_1(const Pyramid_3_1&);
        virtual ~Pyramid_3_1();
    private:
        const std::string               name_;
        double                          weight_[16];
        ReferenceGeometryT* const       refGeoPtr_;
        PointReferenceT                 gp_[16];
    };

//---------------------------------------------------------------------------
    class Pyramid_5_1: public GaussQuadratureRule<3>
    {
     public:
        typedef PointReference<3>       PointReferenceT;
        typedef ReferenceGeometry<3>    ReferenceGeometryT;
    public:
        static Pyramid_5_1& Instance()
            {
                static Pyramid_5_1 theInstance;
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
        Pyramid_5_1();
        Pyramid_5_1(const Pyramid_5_1&);
        virtual ~Pyramid_5_1();
    private:
        const std::string               name_;
        double                          weight_[36];
        ReferenceGeometryT* const       refGeoPtr_;
        PointReferenceT                 gp_[36];
    };

//---------------------------------------------------------------------------
    class Pyramid_7_1: public GaussQuadratureRule<3>
    {
     public:
        typedef PointReference<3>       PointReferenceT;
        typedef ReferenceGeometry<3>    ReferenceGeometryT;
    public:
        static Pyramid_7_1& Instance()
            {
                static Pyramid_7_1 theInstance;
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
        Pyramid_7_1();
        Pyramid_7_1(const Pyramid_7_1&);
        virtual ~Pyramid_7_1();
    private:
        const std::string               name_;
        double                          weight_[48];
        ReferenceGeometryT* const       refGeoPtr_;
        PointReferenceT                 gp_[48];
    };

//---------------------------------------------------------------------------
} // close namespace QuadratureRules
#endif

//------------------------------------------------------------------------------
// File: GaussQuadratureRulesForCube.hpp 
// Headers of Gauss quadrature rules for reference cube.
// Lars Pesch, Fri Mar  3 12:59:11 CET 2006
//----
// Modified from original file allGaussQuadratureRules.hpp
// by M.T. Julianto, Wed Feb 25 10:45:06 UTC 2013
//------------------------------------------------------------------------------
#ifndef GaussQuadratureRulesForCube_hpp
#define GaussQuadratureRulesForCube_hpp
//---------------------------------------------------------------------------
#include "Geometry/PointReference.hpp"
#include "Geometry/ReferenceGeometry.hpp"
#include "Geometry/ReferenceCube.hpp"
#include "Integration/GlobalNamespaceIntegration.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"

//---------------------------------------------------------------------------
namespace QuadratureRules
{
    using Geometry::PointReference;
    using Geometry::ReferenceGeometry;
    using Geometry::ReferenceCube;

//---------------------------------------------------------------------------
    class Cn3_1_1
        : public GaussQuadratureRule<3>
    {
    public:
        typedef ReferenceGeometry<3>    ReferenceGeometryT;
        typedef PointReference<3>       PointReferenceT;
    public:
        static Cn3_1_1& Instance()
            {
                static Cn3_1_1 theInstance;
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
        Cn3_1_1();
        Cn3_1_1(const Cn3_1_1&);
        virtual ~Cn3_1_1();
    private:
        const std::string           name_;
        double                      weight_[1];
        ReferenceGeometryT* const   refGeoPtr_;
        PointReferenceT             gp_[1];
    };

//---------------------------------------------------------------------------
    class Cn3_3_4: public GaussQuadratureRule<3>
    {
    public:
        typedef ReferenceGeometry<3>    ReferenceGeometryT;
        typedef PointReference<3>       PointReferenceT;
    public:
        static Cn3_3_4& Instance()
        {
            static Cn3_3_4 theInstance;
            return theInstance;
        }

        virtual std::string                 getName() const;
        virtual unsigned int                order() const;
        virtual unsigned int                dimension() const;
        virtual unsigned int                nrOfPoints() const;
        virtual double                      weight(unsigned int i) const;
        virtual void                        getPoint(unsigned int i, PointReferenceT& p) const;
        virtual ReferenceGeometryT*         forReferenceGeometry() const;

    private:
        Cn3_3_4();
        Cn3_3_4(const Cn3_3_4&);
        virtual ~Cn3_3_4();
    private:
        const std::string               name_;
        double                          weight_[8];
        ReferenceGeometryT* const       refGeoPtr_;
        PointReferenceT                 gp_[8];
    };

//---------------------------------------------------------------------------
    class Cn3_5_9: public GaussQuadratureRule<3>
    {
    public:
        typedef ReferenceGeometry<3>    ReferenceGeometryT;
        typedef PointReference<3>       PointReferenceT;
    public:
        static Cn3_5_9& Instance()
        {
            static Cn3_5_9 theInstance;
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
        Cn3_5_9();
        Cn3_5_9(const Cn3_5_9&);
        virtual ~Cn3_5_9();
    private:
        const std::string               name_;
        double                          weight_[27];
        ReferenceGeometryT* const       refGeoPtr_;
        PointReferenceT                 gp_[27];
    };

//---------------------------------------------------------------------------
    class C3_7_2: public GaussQuadratureRule<3>
    {
    public:
        typedef ReferenceGeometry<3>    ReferenceGeometryT;
        typedef PointReference<3>       PointReferenceT;
    public:
        static C3_7_2& Instance()
        {
            static C3_7_2 theInstance;
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
        C3_7_2();
        C3_7_2(const C3_7_2&);
        virtual ~C3_7_2();

    private:
        const std::string               name_;
        double                          weight_[34];
        ReferenceGeometryT* const       refGeoPtr_;
        PointReferenceT                 gp_[34];
    };

//---------------------------------------------------------------------------
} // close namespace QuadratureRules
#endif

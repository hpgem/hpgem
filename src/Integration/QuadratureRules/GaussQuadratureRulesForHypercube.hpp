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
#include "Integration/GlobalNamespaceIntegration.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"

//---------------------------------------------------------------------------
namespace QuadratureRules
{
    using Geometry::PointReference;
    using Geometry::ReferenceGeometry;
//---------------------------------------------------------------------------
    class Cn4_1_1: public GaussQuadratureRule
    {
    public:    
        typedef PointReference       PointReferenceT;
        typedef ReferenceGeometry    ReferenceGeometryT;
    public:
        static Cn4_1_1& Instance()
        {
            static Cn4_1_1 theInstance;
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
        Cn4_1_1();
        Cn4_1_1(const Cn4_1_1&);
        virtual ~Cn4_1_1();
    private:
        const std::string               name_;
        double                          weight_[1];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };
//---------------------------------------------------------------------------
    class Cn4_3_4: public GaussQuadratureRule
    {
    public:    
        typedef PointReference       PointReferenceT;
        typedef ReferenceGeometry    ReferenceGeometryT;
    public:
        static Cn4_3_4& Instance()
        {
            static Cn4_3_4 theInstance;
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
        Cn4_3_4();
        Cn4_3_4(const Cn4_3_4&);
        virtual ~Cn4_3_4();
    private:
        const std::string               name_;
        double                          weight_[16];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };

//---------------------------------------------------------------------------
} // close namespace QuadratureRules
#endif

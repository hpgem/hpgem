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
    using Geometry::PointReference;
    using Geometry::ReferenceGeometry;
    using Geometry::ReferenceSquare;

//---------------------------------------------------------------------------
    class Cn2_1_1:public GaussQuadratureRule<2>
    {
    public:
        typedef ReferenceGeometry<2>    ReferenceGeometryT;
        typedef PointReference<2>       PointReferenceT;
        
    public:
        static Cn2_1_1& Instance()
        {
            static Cn2_1_1 theInstance;
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
        Cn2_1_1();
        Cn2_1_1(const Cn2_1_1&);
        virtual ~Cn2_1_1();

    private:
        const std::string               name_;
        double                          weight_[1];
        ReferenceGeometryT* const       refGeoPtr_;
        PointReferenceT                 gp_[1];
    };

//---------------------------------------------------------------------------
    class Cn2_3_4: public GaussQuadratureRule<2>
    {
    public:
        typedef ReferenceGeometry<2>    ReferenceGeometryT;
        typedef PointReference<2>       PointReferenceT;
    public:
        static Cn2_3_4& Instance()
        {
            static Cn2_3_4 theInstance;
            return theInstance;
        }

        virtual std::string         getName() const;
        
        virtual unsigned int        order() const;
        
        virtual unsigned int        dimension() const;
        
        virtual unsigned int        nrOfPoints() const;
        
        virtual double              weight(unsigned int i) const;
        
        virtual void                getPoint(unsigned int i, PointReferenceT& p) const;
        
        virtual ReferenceGeometryT* forReferenceGeometry() const;

    private:
        Cn2_3_4();
        Cn2_3_4(const Cn2_3_4&);
        virtual ~Cn2_3_4();
    private:
        const std::string               name_;
        double                          weight_[4];
        ReferenceGeometryT* const       refGeoPtr_;
        PointReferenceT                 gp_[4];
    };

//---------------------------------------------------------------------------
    class Cn2_5_9: public GaussQuadratureRule<2>
    {
    public:
        typedef ReferenceGeometry<2>    ReferenceGeometryT;
        typedef PointReference<2>       PointReferenceT;
    public:
        static Cn2_5_9& Instance()
            {
                static Cn2_5_9 theInstance;
                return theInstance;
            }

        virtual std::string         getName() const;
        
        virtual unsigned int        order() const;
        
        virtual unsigned int        dimension() const;
        
        virtual unsigned int        nrOfPoints() const;
        
        virtual double              weight(unsigned int i) const;
        
        virtual void                getPoint(unsigned int i, PointReferenceT& p) const;
        
        virtual ReferenceGeometryT* forReferenceGeometry() const;

    private:
        Cn2_5_9();
        Cn2_5_9(const Cn2_5_9&);
        virtual ~Cn2_5_9();
    private:
        const std::string               name_;
        double                          weight_[9];
        ReferenceGeometryT* const       refGeoPtr_;
        PointReferenceT                 gp_[9];
    };

//---------------------------------------------------------------------------
    class C2_7_4:public GaussQuadratureRule<2>
    {
    public:
        typedef ReferenceGeometry<2>    ReferenceGeometryT;
        typedef PointReference<2>       PointReferenceT;
    public:
        static C2_7_4& Instance()
        {
            static C2_7_4 theInstance;
            return theInstance;
        }

        virtual std::string getName() const;
        
        virtual unsigned int        order() const;
        
        virtual unsigned int        dimension() const;
        
        virtual unsigned int        nrOfPoints() const;
        
        virtual double              weight(unsigned int i) const;
        
        virtual void                getPoint(unsigned int i, PointReferenceT& p) const;
        
        virtual ReferenceGeometryT* forReferenceGeometry() const;

    private:
        C2_7_4();
        C2_7_4(const C2_7_4&);
        virtual ~C2_7_4();
    private:
        const std::string               name_;
        double                          weight_[16];
        ReferenceGeometryT* const       refGeoPtr_;
        PointReferenceT                 gp_[16];
    };
//---------------------------------------------------------------------------
} // close namespace QuadratureRules
#endif

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
#include "Integration/GlobalNamespaceIntegration.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"

//---------------------------------------------------------------------------
namespace QuadratureRules
{
    using Geometry::PointReference;
    using Geometry::ReferenceGeometry;
    
//---------------------------------------------------------------------------
    class Cn2_1_1:public GaussQuadratureRule
    {
    public:
        typedef ReferenceGeometry    ReferenceGeometryT;
        typedef PointReference       PointReferenceT;
        
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
        std::vector<PointReferenceT>                   gp_;
    };

//---------------------------------------------------------------------------
    class Cn2_3_4: public GaussQuadratureRule
    {
    public:
        typedef ReferenceGeometry    ReferenceGeometryT;
        typedef PointReference       PointReferenceT;
    public:
        static Cn2_3_4& Instance()
        {
            static Cn2_3_4 theInstance;
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
        Cn2_3_4();
        Cn2_3_4(const Cn2_3_4&);
        virtual ~Cn2_3_4();
    private:
        const std::string               name_;
        double                          weight_[4];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };

//---------------------------------------------------------------------------
    class Cn2_5_9: public GaussQuadratureRule
    {
    public:
        typedef ReferenceGeometry    ReferenceGeometryT;
        typedef PointReference       PointReferenceT;
    public:
        static Cn2_5_9& Instance()
            {
                static Cn2_5_9 theInstance;
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
        Cn2_5_9();
        Cn2_5_9(const Cn2_5_9&);
        virtual ~Cn2_5_9();
    private:
        const std::string               name_;
        double                          weight_[9];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };

//---------------------------------------------------------------------------
    class C2_7_4:public GaussQuadratureRule
    {
    public:
        typedef ReferenceGeometry    ReferenceGeometryT;
        typedef PointReference       PointReferenceT;
    public:
        static C2_7_4& Instance()
        {
            static C2_7_4 theInstance;
            return theInstance;
        }

        virtual std::string getName() const;
        
        virtual unsigned int            order() const;
        
        virtual unsigned int            dimension() const;
        
        virtual unsigned int            nrOfPoints() const;
        
        virtual double                  weight(unsigned int i) const;
        
        virtual void                    getPoint(unsigned int i, PointReferenceT& p) const;
        
        virtual ReferenceGeometryT*     forReferenceGeometry() const;

    private:
        C2_7_4();
        C2_7_4(const C2_7_4&);
        virtual ~C2_7_4();
    private:
        const std::string               name_;
        double                          weight_[16];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };
//---------------------------------------------------------------------------
    class C2_9_5:public GaussQuadratureRule
    {
    public:
        typedef ReferenceGeometry    ReferenceGeometryT;
        typedef PointReference       PointReferenceT;
    public:
        static C2_9_5& Instance()
        {
            static C2_9_5 theInstance;
            return theInstance;
        }

        virtual std::string getName() const;

        virtual unsigned int            order() const;

        virtual unsigned int            dimension() const;

        virtual unsigned int            nrOfPoints() const;

        virtual double                  weight(unsigned int i) const;

        virtual void                    getPoint(unsigned int i, PointReferenceT& p) const;

        virtual ReferenceGeometryT*     forReferenceGeometry() const;

    private:
        C2_9_5();
        C2_9_5(const C2_9_5&);
        virtual ~C2_9_5();
    private:
        const std::string               name_;
        double                          weight_[25];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };
//---------------------------------------------------------------------------
    class C2_11_6:public GaussQuadratureRule
    {
    public:
        typedef ReferenceGeometry    ReferenceGeometryT;
        typedef PointReference       PointReferenceT;
    public:
        static C2_11_6& Instance()
        {
            static C2_11_6 theInstance;
            return theInstance;
        }

        virtual std::string getName() const;

        virtual unsigned int            order() const;

        virtual unsigned int            dimension() const;

        virtual unsigned int            nrOfPoints() const;

        virtual double                  weight(unsigned int i) const;

        virtual void                    getPoint(unsigned int i, PointReferenceT& p) const;

        virtual ReferenceGeometryT*     forReferenceGeometry() const;

    private:
        C2_11_6();
        C2_11_6(const C2_11_6&);
        virtual ~C2_11_6();
    private:
        const std::string               name_;
        double                          weight_[36];
        ReferenceGeometryT* const       refGeoPtr_;
        std::vector<PointReferenceT>                   gp_;
    };
//---------------------------------------------------------------------------
} // close namespace QuadratureRules
#endif

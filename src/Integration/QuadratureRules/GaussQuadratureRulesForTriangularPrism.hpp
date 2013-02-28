//------------------------------------------------------------------------------
// File: GaussQuadratureRulesForTriangularPrism.hpp 
// Headers of Gauss quadrature rules for reference triangular-prism.
// Lars Pesch, Fri Mar  3 12:59:11 CET 2006
//----
// Modified from original file allGaussQuadratureRules.hpp
// by M.T. Julianto, Wed Feb 25 10:45:06 UTC 2013
//------------------------------------------------------------------------------
#ifndef GaussQuadratureRulesForTriangularPrism_hpp
#define GaussQuadratureRulesForTriangularPrism_hpp
//---------------------------------------------------------------------------
#include "Geometry/PointReference.hpp"
#include "Geometry/ReferenceGeometry.hpp"
#include "Geometry/ReferenceTriangularPrism.hpp"
#include "Integration/GlobalNamespaceIntegration.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"

//---------------------------------------------------------------------------
namespace QuadratureRules
{
    using Integration::NumType;
    using Geometry::PointReference;
    using Geometry::ReferenceGeometry;
    using Geometry::ReferenceTriangularPrism;

//---------------------------------------------------------------------------
    class TriPrism_1_1
        : public GaussQuadratureRule<3>
    {
    public:
        static TriPrism_1_1& Instance()
            {
                static TriPrism_1_1 theInstance;
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
        TriPrism_1_1();
        TriPrism_1_1(const TriPrism_1_1&);
        virtual ~TriPrism_1_1();

        const std::string name_;
        NumType weight_[1];
        ReferenceGeometry<3>* const refGeoPtr_;
        PointReference<3> gp_[1];
    };

//---------------------------------------------------------------------------
    class TriPrism_3_1
        : public GaussQuadratureRule<3>
    {
    public:
        static TriPrism_3_1& Instance()
            {
                static TriPrism_3_1 theInstance;
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
        TriPrism_3_1();
        TriPrism_3_1(const TriPrism_3_1&);
        virtual ~TriPrism_3_1();

        const std::string name_;
        NumType weight_[8];
        ReferenceGeometry<3>* const refGeoPtr_;
        PointReference<3> gp_[8];
    };

//---------------------------------------------------------------------------
    class TriPrism_5_1
        : public GaussQuadratureRule<3>
    {
    public:
        static TriPrism_5_1& Instance()
            {
                static TriPrism_5_1 theInstance;
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
        TriPrism_5_1();
        TriPrism_5_1(const TriPrism_5_1&);
        virtual ~TriPrism_5_1();

        const std::string name_;
        NumType weight_[21];
        ReferenceGeometry<3>* const refGeoPtr_;
        PointReference<3> gp_[21];
    };

//---------------------------------------------------------------------------
    class TriPrism_7_1
        : public GaussQuadratureRule<3>
    {
    public:
        static TriPrism_7_1& Instance()
            {
                static TriPrism_7_1 theInstance;
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
        TriPrism_7_1();
        TriPrism_7_1(const TriPrism_7_1&);
        virtual ~TriPrism_7_1();

        const std::string name_;
        NumType weight_[64];
        ReferenceGeometry<3>* const refGeoPtr_;
        PointReference<3> gp_[64];
    };

//---------------------------------------------------------------------------
} // close namespace QuadratureRules
#endif

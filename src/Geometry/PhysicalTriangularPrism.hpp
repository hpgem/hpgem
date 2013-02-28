#ifndef PHYSICALTRIANGULARPRISM_HH
#define PHYSICALTRIANGULARPRISM_HH
#include "PhysicalGeometry.hpp"
#include "ReferenceTriangularPrism.hpp"
namespace Geometry
{
    class PhysicalTriangularPrism: public PhysicalGeometry<ThreeD>
    {
    public:

        typedef PhysicalGeometry<ThreeD> PhysicalGeometry3D;
        using PhysicalGeometry3D::VectorOfPointIndexesT;
        using PhysicalGeometry3D::VectorOfPhysicalPointsT;
        using PhysicalGeometry3D::PointIndexT;

    public:

        PhysicalTriangularPrism(
                const VectorOfPointIndexesT&,
                const VectorOfPhysicalPointsT&,
                const ReferenceTriangularPrism* const cube);

        ~PhysicalTriangularPrism() {}

        virtual std::string     getName() const { return "PhysicalTriangularPrism";}

        void                    getGlobalFaceNodeIndices(const PointIndexT, VectorOfPointIndexesT&) const;

        void                    getLocalFaceNodeIndices(const PointIndexT, VectorOfPointIndexesT&) const;
    };
}
#endif

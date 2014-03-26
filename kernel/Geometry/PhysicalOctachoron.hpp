#ifndef PHYSICALOCTACHORON_HH
#define PHYSICALOCTACHORON_HH
#include "PhysicalGeometry.hpp"
namespace Geometry
{
    class PhysicalHexahedron: public PhysicalGeometry
    {
        public:
            typedef PhysicalGeometry PhysicalGeometry4D;
            using PhysicalGeometry3D::VectorOfPointIndexesT;
            using PhysicalGeometry3D::VectorOfPhysicalPointsT;
            using PhysicalGeometry3D::PointIndexT;

        public:

            PhysicalHexahedron(
                    const VectorOfPointIndexesT&,
                    const VectorOfPhysicalPointsT&,
                    const ReferenceHypercube* const);

            ~PhysicalHexahedron() {}

            /// Returns the name of this geometry.
            virtual std::string             getName() const { return "PhysicalOctachron";}

            virtual void getGlobalFaceNodeIndices(const PointIndexT, VectorOfPointIndexesT&) const;

            virtual void getLocalFaceNodeIndices(const PointIndexT, VectorOfPointIndexesT&) const
            
        private:
               
    };
}
#endif

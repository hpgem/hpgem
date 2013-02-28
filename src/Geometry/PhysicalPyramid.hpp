#ifndef PHYSICALPYRAMID_HH
#define PHYSICALPYRAMID_HH

#include "PhysicalGeometry.hpp"
#include "ReferencePyramid.hpp"

namespace Geometry
{
    class PhysicalPyramid: public PhysicalGeometry<ThreeD>
    {
        public:
            typedef PhysicalGeometry<ThreeD> PhysicalGeometry3D;
            using PhysicalGeometry3D::VectorOfPointIndexesT;
            using PhysicalGeometry3D::VectorOfPhysicalPointsT;
            using PhysicalGeometry3D::PointIndexT;

        public:

            PhysicalPyramid(
                    const VectorOfPointIndexesT&,
                    const VectorOfPhysicalPointsT&,
                    const ReferencePyramid* const);

            ~PhysicalPyramid() {}

            /// Returns the name of this geometry.
            virtual std::string             getName() const { return "PhysicalPyramid";}

            virtual void getGlobalFaceNodeIndices(const PointIndexT, VectorOfPointIndexesT&) const;

            virtual void getLocalFaceNodeIndices(const PointIndexT, VectorOfPointIndexesT&) const;

            //PhysicalGeometries getId() const {return id_;}

        private:
            //            const ReferencePyramid* pyramid_;
    };
}
#endif

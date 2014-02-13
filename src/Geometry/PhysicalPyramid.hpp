#ifndef PHYSICALPYRAMID_HH
#define PHYSICALPYRAMID_HH

#include "PhysicalGeometry.hpp"
#include "ReferencePyramid.hpp"

namespace Geometry
{
    class PhysicalPyramid: public PhysicalGeometry
    {
        public:
            typedef PhysicalGeometry PhysicalGeometry3D;
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

            unsigned int getNrOfFaces() const {return 5;}

            //PhysicalGeometries getId() const {return id_;}

        private:
            //            const ReferencePyramid* pyramid_;
    };
}
#endif

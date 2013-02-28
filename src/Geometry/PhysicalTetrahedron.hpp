#ifndef PHYSICALTETRAHEDRON_HH
#define PHYSICALTETRAHEDRON_HH
#include "PhysicalGeometry.hpp"
#include "ReferenceCube.hpp"
namespace Geometry
{
    class PhysicalTetrahedron: public PhysicalGeometry<ThreeD>
    {
        public:

            typedef PhysicalGeometry<ThreeD> PhysicalGeometry3D;
            using PhysicalGeometry3D::VectorOfPointIndexesT;
            using PhysicalGeometry3D::VectorOfPhysicalPointsT;
            using PhysicalGeometry3D::PointIndexT;

        public:

            PhysicalTetrahedron(
                    const VectorOfPointIndexesT&,
                    const VectorOfPhysicalPointsT&,
                    const ReferenceCube* const cube);

            ~PhysicalTetrahedron() {}

            virtual std::string             getName() const { return "PhysicalTetrahedron";}

            void getGlobalFaceNodeIndices(const PointIndexT, VectorOfPointIndexesT&) const;

            void getLocalFaceNodeIndices(const PointIndexT, VectorOfPointIndexesT&) const;
        
            unsigned int getNrOfFaces() const {return refGeometry_->getNrOfCodim1Entities();}
    };
}
#endif

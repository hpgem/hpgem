#ifndef PHYSICALHEXAHEDRON_HH
#define PHYSICALHEXAHEDRON_HH
#include "PhysicalGeometry.hpp"
#include "ReferenceCube.hpp"
namespace Geometry
{
    class PhysicalHexahedron: public PhysicalGeometry<3>
    {
        public:

            typedef PhysicalGeometry<ThreeD> PhysicalGeometry3D;
            using PhysicalGeometry3D::VectorOfPointIndexesT;
            using PhysicalGeometry3D::VectorOfPhysicalPointsT;
            using PhysicalGeometry3D::PointIndexT;

        public:

            PhysicalHexahedron(
                    const VectorOfPointIndexesT&,
                    const VectorOfPhysicalPointsT&,
                    const ReferenceCube* const);

            ~PhysicalHexahedron() {}

            virtual std::string             getName() const { return "PhysicalHexahedron";}

            virtual void getGlobalFaceNodeIndices(const PointIndexT, VectorOfPointIndexesT&) const;

            virtual void getLocalFaceNodeIndices(const PointIndexT, VectorOfPointIndexesT&) const;
        
            unsigned int getNrOfFaces() const {return refGeometry_->getNrOfCodim1Entities();}
    };
}
#endif

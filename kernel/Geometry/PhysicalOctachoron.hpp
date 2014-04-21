#ifndef PHYSICALOCTACHORON_HH
#define PHYSICALOCTACHORON_HH
#include "PhysicalGeometry.hpp"
namespace Geometry
{
	class ReferenceHypercube;

    class PhysicalOctachoron: public PhysicalGeometry
    {
        public:
            typedef PhysicalGeometry PhysicalGeometry4D;
            using PhysicalGeometry4D::VectorOfPointIndexesT;
            using PhysicalGeometry4D::VectorOfPhysicalPointsT;
            using PhysicalGeometry4D::PointIndexT;

        public:

            PhysicalOctachoron(
                    const VectorOfPointIndexesT&,
                    const VectorOfPhysicalPointsT&,
                    const ReferenceHypercube* const);

            ~PhysicalOctachoron() {}

            /// Returns the name of this geometry.
            virtual std::string             getName() const { return "PhysicalOctachron";}

            virtual void getGlobalFaceNodeIndices(const PointIndexT, VectorOfPointIndexesT&) const;

            virtual void getLocalFaceNodeIndices(const PointIndexT, VectorOfPointIndexesT&) const;
            
            virtual unsigned int getNrOfFaces() const {return getRefGeometry()->getNrOfCodim1Entities();}

        private:
               
    };
}
#endif

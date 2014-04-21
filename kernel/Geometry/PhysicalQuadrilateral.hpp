#ifndef PHYSICALQUADRILATERAL_HH
#define PHYSICALQUADRILATERAL_HH

#include "PhysicalGeometry.hpp"
#include "ReferenceSquare.hpp"

namespace Geometry
{
    class PhysicalQuadrilateral: public PhysicalGeometry
    {
        public:

            typedef PhysicalGeometry PhysicalGeometry2D;
            using PhysicalGeometry2D::VectorOfPointIndexesT;
            using PhysicalGeometry2D::VectorOfPhysicalPointsT;
            using PhysicalGeometry2D::PointIndexT;

        public:

            PhysicalQuadrilateral(
                    const VectorOfPointIndexesT&,
                    const VectorOfPhysicalPointsT&,
                    const ReferenceSquare* const);

            ~PhysicalQuadrilateral() {}

            virtual std::string             getName() const { return "PhysicalQuadrilateral";}

            void getGlobalFaceNodeIndices(const PointIndexT, VectorOfPointIndexesT&) const;

            void getLocalFaceNodeIndices(const PointIndexT, VectorOfPointIndexesT&) const;
        
            unsigned int getNrOfFaces() const {return refGeometry_->getNrOfCodim1Entities();}
    };
}
#endif

#ifndef PHYSICALTRIANGLE_HH
#define PHYSICALTRIANGLE_HH
#include "PhysicalGeometry.hpp"
#include "ReferenceTriangle.hpp"
namespace Geometry
{
    class PhysicalTriangle: public PhysicalGeometry
    {
        public:

            typedef PhysicalGeometry PhysicalGeometry2D;
            using PhysicalGeometry2D::VectorOfPointIndexesT;
            using PhysicalGeometry2D::VectorOfPhysicalPointsT;
            using PhysicalGeometry2D::PointIndexT;

        public:

            PhysicalTriangle(
                    const VectorOfPointIndexesT&,
                    const VectorOfPhysicalPointsT&,
                    const ReferenceTriangle* const);

            ~PhysicalTriangle() {}

            virtual std::string             getName() const { return "PhysicalTriangle";}

            virtual void getGlobalFaceNodeIndices(const PointIndexT, VectorOfPointIndexesT&) const;

            virtual void getLocalFaceNodeIndices(const PointIndexT, VectorOfPointIndexesT&) const;
        
            unsigned int getNrOfFaces() const {return refGeometry_->getNrOfCodim1Entities();}
    };
}
#endif

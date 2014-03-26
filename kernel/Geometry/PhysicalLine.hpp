#ifndef PHYSICALLINE_HH
#define PHYSICALLINE_HH
#include "PhysicalGeometry.hpp"
#include "ReferenceLine.hpp"
namespace Geometry
{
    class PhysicalLine: public PhysicalGeometry
    {

    public:

        typedef PhysicalGeometry PhysicalGeometry1D;
        using PhysicalGeometry1D::VectorOfPointIndexesT;
        using PhysicalGeometry1D::VectorOfPhysicalPointsT;
        using PhysicalGeometry1D::PointIndexT;

    public:

        PhysicalLine(
                const VectorOfPointIndexesT&,
                const VectorOfPhysicalPointsT&,
                const ReferenceLine* const);

        ~PhysicalLine() {}

        virtual std::string             getName() const { return "PhysicalLine";}

        void getGlobalFaceNodeIndices(const PointIndexT, VectorOfPointIndexesT&) const;

        void getLocalFaceNodeIndices(const PointIndexT, VectorOfPointIndexesT&) const;
    
        unsigned int getNrOfFaces() const {return refGeometry_->getNrOfCodim1Entities();}
    };
}
#endif

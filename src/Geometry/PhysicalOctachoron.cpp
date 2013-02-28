#include "PhysicalOctachoron.hpp"

namespace Geometry
{
    PhysicalOctachoron::PhysicalOctachoron(
            const VectorOfPointIndexesT& globalNodeIndexes,
            const VectorOfPhysicalPointsT& nodes,
            const ReferenceHypercube* const cube) :
            PhysicalGeometry4D(globalNodeIndexes,nodes, cube)
    {
    }

    void PhysicalHexahedron::getGlobalFaceNodeIndices(const PointIndexT face, VectorOfPointIndexesT& indexes) const
    {
        indexes.resize(6);
        indexes[0] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,0)];
        indexes[1] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,1)];
        indexes[2] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,2)];
        indexes[3] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,3)];
        indexes[4] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,4)];
        indexes[5] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,5)];
        
    }

    void PhysicalHexahedron::getLocalFaceNodeIndices(const PointIndexT face, VectorOfPointIndexesT& indexes) const
    {
        indexes.resize(6);
        indexes[0] = refGeometry_->getLocalNodeIndex(face,0);
        indexes[1] = refGeometry_->getLocalNodeIndex(face,1);
        indexes[2] = refGeometry_->getLocalNodeIndex(face,2);
        indexes[3] = refGeometry_->getLocalNodeIndex(face,3);
        indexes[4] = refGeometry_->getLocalNodeIndex(face,4);
        indexes[5] = refGeometry_->getLocalNodeIndex(face,5);
    }
}

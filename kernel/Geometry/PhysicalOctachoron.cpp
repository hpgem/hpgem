#include "PhysicalOctachoron.hpp"
#include "ReferenceHypercube.hpp"

namespace Geometry
{
    PhysicalOctachoron::PhysicalOctachoron(
            const VectorOfPointIndexesT& globalNodeIndexes,
            const VectorOfPhysicalPointsT& nodes,
            const ReferenceHypercube* const cube) :
            PhysicalGeometry4D(globalNodeIndexes,nodes, cube)
    {
    }

    void PhysicalOctachoron::getGlobalFaceNodeIndices(const PointIndexT face, VectorOfPointIndexesT& indexes) const
    {
        indexes.resize(8);
        indexes[0] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,0)];
        indexes[1] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,1)];
        indexes[2] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,2)];
        indexes[3] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,3)];
        indexes[4] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,4)];
        indexes[5] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,5)];
        indexes[6] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,6)];
        indexes[7] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,7)];
    }

    void PhysicalOctachoron::getLocalFaceNodeIndices(const PointIndexT face, VectorOfPointIndexesT& indexes) const
    {
        indexes.resize(8);
        indexes[0] = refGeometry_->getLocalNodeIndex(face,0);
        indexes[1] = refGeometry_->getLocalNodeIndex(face,1);
        indexes[2] = refGeometry_->getLocalNodeIndex(face,2);
        indexes[3] = refGeometry_->getLocalNodeIndex(face,3);
        indexes[4] = refGeometry_->getLocalNodeIndex(face,4);
        indexes[5] = refGeometry_->getLocalNodeIndex(face,5);
        indexes[6] = refGeometry_->getLocalNodeIndex(face,6);
        indexes[7] = refGeometry_->getLocalNodeIndex(face,7);
    }
}

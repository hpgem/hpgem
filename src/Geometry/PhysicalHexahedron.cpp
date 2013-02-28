/*
 * PhysicalHexahedron.hpp
 *
 *  Created on: Feb 9, 2013
 *      Author: nicorivas
 */
#include "PhysicalHexahedron.hpp"

namespace Geometry
{
    PhysicalHexahedron::PhysicalHexahedron(
            const VectorOfPointIndexesT& globalNodeIndexes,
            const VectorOfPhysicalPointsT& nodes,
            const ReferenceCube* const cube) :
            PhysicalGeometry<3>(globalNodeIndexes,nodes,cube)
    {
    }

    void PhysicalHexahedron::getGlobalFaceNodeIndices(const PointIndexT face, VectorOfPointIndexesT& indexes) const
    {
        indexes.resize(4);
        indexes[0] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,0)];
        indexes[1] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,1)];
        indexes[2] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,2)];
        indexes[3] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,3)];
    }

    void PhysicalHexahedron::getLocalFaceNodeIndices(const PointIndexT face, VectorOfPointIndexesT& indexes) const
    {
        indexes.resize(4);
        indexes[0] = refGeometry_->getLocalNodeIndex(face,0);
        indexes[1] = refGeometry_->getLocalNodeIndex(face,1);
        indexes[2] = refGeometry_->getLocalNodeIndex(face,2);
        indexes[3] = refGeometry_->getLocalNodeIndex(face,3);
    }
}

/*
 * PhysicalTetrahedron.cpp
 *
 *  Created on: Feb 9, 2013
 *      Author: nicorivas
 */
#include "PhysicalTetrahedron.hpp"
#include <vector>

namespace Geometry
{
    PhysicalTetrahedron::PhysicalTetrahedron(
        const VectorOfPointIndexesT& globalNodeIndexes,
        const VectorOfPhysicalPointsT& nodes,
        const ReferenceCube* const cube) :
        PhysicalGeometry<3>(globalNodeIndexes,nodes, cube)
    {
    }

    void PhysicalTetrahedron::getGlobalFaceNodeIndices(const PointIndexT face, VectorOfPointIndexesT& indexes) const
    {
        indexes.resize(3);
        indexes[0] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,0)];
        indexes[1] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,1)];
        indexes[2] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,2)];
    }

    void PhysicalTetrahedron::getLocalFaceNodeIndices(const PointIndexT face, VectorOfPointIndexesT& indexes) const
    {
        indexes.resize(3);
        indexes[0] = refGeometry_->getLocalNodeIndex(face,0);
        indexes[1] = refGeometry_->getLocalNodeIndex(face,1);
        indexes[2] = refGeometry_->getLocalNodeIndex(face,2);
    }
}

/*
 * PhysicalQuadrilateral.cpp
 *
 *  Created on: Feb 9, 2013
 *      Author: nicorivas
 */
#include "PhysicalQuadrilateral.hpp"
#include <vector>

namespace Geometry
{
    PhysicalQuadrilateral::PhysicalQuadrilateral(
        const VectorOfPointIndexesT& globalNodeIndexes,
        const VectorOfPhysicalPointsT& nodes,
        const ReferenceSquare* const square) :
        PhysicalGeometry<2>(globalNodeIndexes,nodes, square)
    {
    }

    void PhysicalQuadrilateral::getGlobalFaceNodeIndices(const PointIndexT face, VectorOfPointIndexesT& indexes) const
    {
        indexes.resize(2);
        indexes[0] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,0)];
        indexes[1] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,1)];
    }

    void PhysicalQuadrilateral::getLocalFaceNodeIndices(const PointIndexT face, VectorOfPointIndexesT& indexes) const
    {
        indexes.resize(2);
        indexes[0] = refGeometry_->getLocalNodeIndex(face,0);
        indexes[1] = refGeometry_->getLocalNodeIndex(face,1);
    }
}

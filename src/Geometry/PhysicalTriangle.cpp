/*
 * PhysicalCube.cpp
 *
 *  Created on: Feb 9, 2013
 *      Author: nicorivas
 */
#include "PhysicalTriangle.hpp"
#include <vector>

namespace Geometry
{
    PhysicalTriangle::PhysicalTriangle(
            const VectorOfPointIndexesT& globalNodeIndexes,
            const VectorOfPhysicalPointsT& nodes,
            const ReferenceTriangle* const triangle) :
        PhysicalGeometry(globalNodeIndexes,nodes, triangle)
    {
    }

    void PhysicalTriangle::getGlobalFaceNodeIndices(const PointIndexT face, VectorOfPointIndexesT& indexes) const
    {
        indexes.resize(2);
        indexes[0] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,0)];
        indexes[1] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,1)];
    }

    void PhysicalTriangle::getLocalFaceNodeIndices(const PointIndexT face, VectorOfPointIndexesT& indexes) const
    {
        indexes.resize(2);
        indexes[0] = refGeometry_->getLocalNodeIndex(face,0);
        indexes[1] = refGeometry_->getLocalNodeIndex(face,1);
    }
}

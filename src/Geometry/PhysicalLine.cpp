/*
 * PhysicalCube.cpp
 *
 *  Created on: Feb 9, 2013
 *      Author: nicorivas
 */
#include "PhysicalLine.hpp"
#include <vector>

namespace Geometry
{
    PhysicalLine::PhysicalLine(
        const VectorOfPointIndexesT& globalNodeIndexes,
        const VectorOfPhysicalPointsT& nodes,
        const ReferenceLine* const line) :
        PhysicalGeometry<OneD>(globalNodeIndexes,nodes, line)
    {
    }

    void PhysicalLine::getGlobalFaceNodeIndices(const PointIndexT face, VectorOfPointIndexesT& indexes) const
    {
        indexes.resize(1);
        indexes[0] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,0)];
    }

    void PhysicalLine::getLocalFaceNodeIndices(const PointIndexT face, VectorOfPointIndexesT& indexes) const
    {
        indexes.resize(1);
        indexes[0] = refGeometry_->getLocalNodeIndex(face,0);
    }
}

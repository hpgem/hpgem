/*
 * PhysicalPyramid.cpp
 *
 *  Created on: Feb 9, 2013
 *      Author: nicorivas
 */
#include "PhysicalPyramid.hpp"
#include <vector>

namespace Geometry
{
    PhysicalPyramid::PhysicalPyramid(
        const VectorOfPointIndexesT& globalNodeIndexes,
        const VectorOfPhysicalPointsT& nodes,
        const ReferencePyramid* const refPyramid) :
        PhysicalGeometry<3>(globalNodeIndexes,nodes, refPyramid)
    {
    }

    void PhysicalPyramid::getGlobalFaceNodeIndices(const PointIndexT face, VectorOfPointIndexesT& indexes) const
    {
        indexes.resize(5);
        for (int i = 0; i < 5; ++i)
        {
            indexes[i] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,i)];
        }
    }

    void PhysicalPyramid::getLocalFaceNodeIndices(const PointIndexT face, VectorOfPointIndexesT& indexes) const
    {
        indexes.resize(5);
        for (int i = 0; i < 5; ++i)
        {
            indexes[i] = refGeometry_->getLocalNodeIndex(face,i);
        }
    }
}

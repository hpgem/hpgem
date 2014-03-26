/*
 * PhysicalPyramid.cpp
 *
 *  Created on: Feb 9, 2013
 *      Author: nicorivas
 */
#include "PhysicalPyramid.hpp"

namespace Geometry
{
    PhysicalPyramid::PhysicalPyramid(
        const VectorOfPointIndexesT& globalNodeIndexes,
        const VectorOfPhysicalPointsT& nodes,
        const ReferencePyramid* const refPyramid) :
        PhysicalGeometry(globalNodeIndexes,nodes, refPyramid)
    {
    }

    void PhysicalPyramid::getGlobalFaceNodeIndices(const PointIndexT face, VectorOfPointIndexesT& indexes) const
    {
        indexes.resize(4);
        if (face==0)
        {
            for (int i = 0; i < 4; ++i)
            {
                indexes[i] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,i)];
            }
        }
        else
        {
            for (int i = 0; i < 3; ++i)
            {
                indexes[i] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,i)];
            }
        }
    }

    void PhysicalPyramid::getLocalFaceNodeIndices(const PointIndexT face, VectorOfPointIndexesT& indexes) const
    {
        indexes.resize(4);
        if (face==0)
        {
            for (int i = 0; i < 4; ++i)
            {
                 indexes[i] = refGeometry_->getLocalNodeIndex(face,i);
            }
        }
        else
        {
            for (int i = 0; i < 3; ++i)
            {
                indexes[i] = refGeometry_->getLocalNodeIndex(face,i);
            }
        }
    }
}

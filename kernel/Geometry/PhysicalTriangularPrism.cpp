/*
 * PhysicalTriangularPrism.hpp
 *
 *  Created on: Feb 9, 2013
 *      Author: nicorivas
 */
#include "PhysicalTriangularPrism.hpp"

namespace Geometry
{
    PhysicalTriangularPrism::PhysicalTriangularPrism(
            const VectorOfPointIndexesT& globalNodeIndexes,
            const VectorOfPhysicalPointsT& nodes,
            const ReferenceTriangularPrism* const triangularPrism) :
            PhysicalGeometry(globalNodeIndexes, nodes, triangularPrism)
    {
    }

    void PhysicalTriangularPrism::getGlobalFaceNodeIndices(const PointIndexT face, VectorOfPointIndexesT& indexes) const
    {
        if (face <= 4)
        {
            if (face <= 1)
            {
                indexes.resize(3);
                indexes[0] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,0)];
                indexes[1] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,1)];
                indexes[2] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,2)];
            }
            else
            {
                indexes.resize(4);
                indexes[0] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,0)];
                indexes[1] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,1)];
                indexes[2] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,2)];
                indexes[3] = globalNodeIndexes_[refGeometry_->getLocalNodeIndex(face,3)];
            }
        }
    }

    void PhysicalTriangularPrism::getLocalFaceNodeIndices(const PointIndexT face, VectorOfPointIndexesT& indexes) const
    {
        if (face <= 4)
        {
            if (face <= 1)
            {
                indexes.resize(3);
                indexes[0] = refGeometry_->getLocalNodeIndex(face,0);
                indexes[1] = refGeometry_->getLocalNodeIndex(face,1);
                indexes[2] = refGeometry_->getLocalNodeIndex(face,2);
            }
            else
            {
                indexes.resize(4);
                indexes[0] = refGeometry_->getLocalNodeIndex(face,0);
                indexes[1] = refGeometry_->getLocalNodeIndex(face,1);
                indexes[2] = refGeometry_->getLocalNodeIndex(face,2);
                indexes[3] = refGeometry_->getLocalNodeIndex(face,3);
            }
        }
    }
}

/*
 * FaceCacheData.cpp
 *
 *  Created on: Feb 3, 2014
 *      Author: brinkf
 */

#include "FaceCacheData.hpp"
#include "Face.hpp"
#include "Geometry/PointReference.hpp"

void Base::FaceCacheData::operator ()(const Base::Face& fa, const Geometry::PointReference& p)
{
    fa.getNormalVector(p, Normal);
    L2Normal = Base::L2Norm(Normal);
}

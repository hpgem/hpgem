/*
 * MappingToRefPointToPoint.cpp
 *
 *  Created on: Jan 21, 2014
 *      Author: brinkf
 */

#include "Geometry/Mappings/MappingToRefPointToPoint.hpp"

namespace Geometry {

MappingToRefPointToPoint::MappingToRefPointToPoint() {
}

MappingToRefPointToPoint::MappingToRefPointToPoint(const MappingToRefPointToPoint&) {

}

MappingToRefPointToPoint::~MappingToRefPointToPoint() {
}

const MappingToRefPointToPoint& MappingToRefPointToPoint::Instance(){
	static const MappingToRefPointToPoint theInstance;
	return theInstance;
}

void MappingToRefPointToPoint::transform(const Geometry::PointReference& p1, Geometry::PointReference& p2) const
{

}

void MappingToRefPointToPoint::calcJacobian(const Geometry::PointReference& p, Geometry::Jacobian& jacobean) const
{

}


} /* namespace Geometry */

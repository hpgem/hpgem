/*
 * MappingToRefPointToPoint.hpp
 *
 *  Created on: Jan 21, 2014
 *      Author: brinkf
 */

#ifndef MAPPINGTOREFPOINTTOPOINT_HPP_
#define MAPPINGTOREFPOINTTOPOINT_HPP_

#include "Geometry/Mappings/MappingReferenceToReference.hpp"

namespace Geometry {

	/*
	 * maps a point to itself
	 *
	 * there is only one possible mapping
	 *
	 * this class is provided for a unified treatment of faces in the 1D case (compared to higher dimensions)
	 *
	 */

	class MappingToRefPointToPoint: public Geometry::MappingReferenceToReference {
	public:
		static const MappingToRefPointToPoint& Instance();
		virtual void transform(const Geometry::PointReference& p1,
									 Geometry::PointReference& p2) const;
		virtual void calcJacobian(const Geometry::PointReference&,
										Geometry::Jacobian&) const;
        virtual int getTargetDimension() const {return 0;}
	private:
		MappingToRefPointToPoint();
		MappingToRefPointToPoint(const MappingToRefPointToPoint&);
		virtual ~MappingToRefPointToPoint();
	};

} /* namespace Geometry */
#endif /* MAPPINGTOREFPOINTTOPOINT_HPP_ */

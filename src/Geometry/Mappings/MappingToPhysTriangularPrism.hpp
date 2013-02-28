/*
 * MappingToPhysTriangularPrism.hpp
 *
 *  Created on: Feb 5, 2013
 *      Author: nicorivas
 */

#ifndef TRIANGULARPRISM_H_
#define TRIANGULARPRISM_H_

#include <vector>
#include "../Jacobian.hpp"
#include "MappingReferenceToPhysical.hpp"

namespace Geometry
{
    /*! \brief Map from the reference (triangular) prism to physical space.
     *
     *  Implementation is based on Maple/triangularprismMapping.mws. For the
     *  purpose of individual methods see the documentation of the base classes,
     *  Ref2PhysSpaceMapping and Mapping. */

    class MappingToPhysTriangularPrism: public MappingReferenceToPhysical<3,3>
    {
        public:
            typedef Geometry::PhysicalGeometry<3> PhysicalGeometryT;
            typedef Geometry::PointReference<3> PointReferenceT;
            typedef Geometry::PointPhysical<3> PointPhysicalT;
            typedef Geometry::Jacobian<3,3> JacobianT;

        public:
            MappingToPhysTriangularPrism(const PhysicalGeometryT*const);

            virtual void transform(const PointReferenceT&, PointPhysicalT&) const;

            virtual void calcJacobian(const PointReferenceT&, JacobianT&) const;

            virtual void reinit(const PhysicalGeometryT*const);

            bool isValidPoint(const PointReferenceT&) const;

        private:
            PointPhysicalT a0, a1, a2, a3, a4, a5;
        
            //PhysicalGeometryT* physicalGeometry_;

    };
};
#endif

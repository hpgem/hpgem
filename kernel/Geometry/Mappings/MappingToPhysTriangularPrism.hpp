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

    class MappingToPhysTriangularPrism: public MappingReferenceToPhysical
    {
        public:
            typedef Geometry::PhysicalGeometry PhysicalGeometryT;
            typedef Geometry::PointReference PointReferenceT;
            typedef Geometry::PointPhysical PointPhysicalT;
            typedef Geometry::Jacobian JacobianT;

        public:
            MappingToPhysTriangularPrism(const PhysicalGeometryT*const);

            virtual void transform(const PointReferenceT&, PointPhysicalT&) const;

            virtual void calcJacobian(const PointReferenceT&, JacobianT&) const;

            virtual void reinit(const PhysicalGeometryT*const);

            bool isValidPoint(const PointReferenceT&) const;
            virtual int getTargetDimension() const {return 3;}

        private:
            PointPhysicalT a0, a1, a2, a3, a4, a5;
            int globalNodeIndices_[6];
            //PhysicalGeometryT* physicalGeometry_;

    };
};
#endif

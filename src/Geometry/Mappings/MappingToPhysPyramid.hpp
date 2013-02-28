/*
 * MappingToPhysPyramid.hpp
 *
 *  Created on: Feb 5, 2013
 *      Author: nicorivas
 */

#ifndef MAPPINGPYRAMID_H_
#define MAPPINGPYRAMID_H_

#include "MappingReferenceToPhysical.hpp"
#include "../GlobalNamespaceGeometry.hpp"
#include "../Jacobian.hpp"
#include <vector>

namespace Geometry
{
    /*!
     * "In geometry, a pyramid is a polyhedron formed by connecting a polygonal base and a point,
     *  called the apex. Each base edge and apex form a triangle." -Wikipedia.
     *
     * This class defines the mappings between reference and physical pyramids.
     */

    class MappingToPhysPyramid: public MappingReferenceToPhysical<3,3>
    {
        private:
            typedef Geometry::PhysicalGeometry<3> PhysicalGeometryT;
            typedef Geometry::PointReference<3> PointReferenceT;
            typedef Geometry::PointPhysical<3> PointPhysicalT;
            typedef Geometry::Jacobian<3,3> JacobianT;

        public:
            MappingToPhysPyramid(const PhysicalGeometryT*const physicalGeometry);

            virtual void transform(const PointReferenceT&, PointPhysicalT&) const;

            virtual void calcJacobian(const PointReferenceT&, JacobianT&) const;

            virtual void reinit(const PhysicalGeometryT*const);

            bool isValidPoint(const PointReferenceT&) const;

        private:
            // ~OC~ Only undefined version is tested
            #undef SAVECOEFFS
            #ifdef SAVECOEFFS
                PointPhysicalT a0, a1, a2, a3, a4;
            #else
                std::vector<int> globalNodeIndices_;
            #endif
    };
};
#endif

/*
 * MeshMoverBase.hpp
 *
 *  Created on: Feb 21, 2013
 *      Author: nicorivas
 */

#ifndef MESHMOVERBASE_HPP_
#define MESHMOVERBASE_HPP_

#include "../Geometry/PointPhysical.hpp"

namespace Base {

    template <unsigned int DIM>
    class MeshMoverBase {

    public:

        typedef Geometry::PointPhysical<DIM> PointPhysicalT;

    public:

        MeshMoverBase() {};

        virtual ~MeshMoverBase() {};

        /// This pure virtual function should be implemented in the Problem that needs moving meshes,
        /// and is called by MeshManipulator on every point.
        virtual void movePoint(PointPhysicalT& point) const {};
    };
}

#endif

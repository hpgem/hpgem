/*
 * MeshMover.hpp
 *
 *  Created on: Feb 21, 2013
 *      Author: nicorivas
 */

#ifndef MESHMOVER_HPP_
#define MESHMOVER_HPP_

#include "Base/MeshMoverBase.hpp"

namespace Base
{
    template <unsigned int DIM>
    class MeshMover : public MeshMoverBase<DIM>
    {

        public:

            typedef Geometry::PointPhysical<DIM> PointPhysicalT;

        public:

            MeshMover() {}

            virtual ~MeshMover() {};

            void movePoint(PointPhysicalT* point) const
            {
                std::cout << "we are moving" << std::endl;
                point->operator[](0) = (0.5-fabs(point->operator[](0)-0.5))/10.0;
                point->operator[](1) = (0.5-fabs(point->operator[](1)-0.5))/10.0;
            }
    };
};

#endif /* MESHMOVER_HPP_ */

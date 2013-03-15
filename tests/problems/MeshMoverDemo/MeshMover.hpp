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

            void movePoint(PointPhysicalT& point) const
            {
                std::cout << "we are moving" << std::endl;
               
                point[0]=point[0];
                point[1]=point[1]*(point[0]+0.5);
            }
    };
};

#endif /* MESHMOVER_HPP_ */

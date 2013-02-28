/*
 * PointPhysical.hpp
 *
 *  Created on: Feb 5, 2013
 *      Author: nicorivas
 */

#ifndef POINTPHYSICAL_HPP_
#define POINTPHYSICAL_HPP_

#include "Point.hpp"

namespace Geometry
{
    template<unsigned int DIM>
    class PointPhysical: public Point<DIM>
    {

    public:

        typedef Point<DIM>                          PointT;
        typedef PointPhysical<DIM>                  PointPhysicalT;
        typedef double                              CoordTypeT;
        typedef typename PointT::VectorOfCoordsT    VectorOfCoordsT;

    public:

        PointPhysical(): PointT() {}

        PointPhysical(const PointT& p): PointT(p) {}

        PointPhysical(const VectorOfCoordsT& coord):Point<DIM>(coord) {}
    
        PointPhysical  operator* (double right) const
                {return PointPhysical(PointT::coordinates_ * right);}

        PointPhysical  operator* (double right)
                {return PointPhysical(PointT::coordinates_ * right);}
        
        PointPhysical  operator+ (const PointPhysical& right)
                {return PointPhysical(PointT::coordinates_  + right.coordinates_);}
        
        PointPhysical  operator+ (const PointPhysical& right) const
                {return PointPhysical(PointT::coordinates_  + right.coordinates_);}
        
        PointPhysical  operator- (const PointPhysical& right)
                {return PointPhysical(PointT::coordinates_  - right.coordinates_);}
        
        PointPhysical  operator- (const PointPhysical& right) const
                {return PointPhysical(PointT::coordinates_  - right.coordinates_);}
        
        PointPhysical& operator= (const PointPhysical& right)
                {PointT::coordinates_ = right.coordinates_; return *this;}
        
//        friend PointT operator*(const double& left, const PointT& right){return PointPhysical(right.coordinates_*left);}
    };
};

#endif /* POINTPHYSICAL_HPP_ */

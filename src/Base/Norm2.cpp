//
//  Norm2.cpp
//  
//
//  Created by Shavarsh Nurijanyan on 6/21/13.
//
//

#include "Norm2.hpp"

namespace Utilities
{
    template <>
    double norm2<1>(const Geometry::PointPhysical<1>& p)
    {
        return std::abs(p[0]);
    }
    
    template <>
    double norm2<2>(const Geometry::PointPhysical<2>& p)
    {
        return std::sqrt(p[0] * p[0] + p[1] * p[1]);
    }
    
    template <>
    double norm2<3>(const Geometry::PointPhysical<3>& p)
    {
        return std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
    }
    
    template <>
    double norm2<4>(const Geometry::PointPhysical<4>& p)
    {
        return std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2] + p[3] * p[3]);
    }
}
    //------------------------------------------------------------------------------
    // Local variables:
    // mode:c++
    // comment-column: 48
    // End:

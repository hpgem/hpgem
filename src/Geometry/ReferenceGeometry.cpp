//
//  ReferenceGeometry_Impl.hpp
//  
//
//  Created by Shavarsh Nurijanyan on 2/4/13.
//
//

#include "ReferenceGeometry.hpp"

#ifndef _ReferenceGeometry_Impl_hpp
#define _ReferenceGeometry_Impl_hpp

namespace Geometry
{
    
    ReferenceGeometry::ReferenceGeometry(const TypeOfReferenceGeometry& geoT):
        geometryType_(geoT)
    {
        
    }
    
    ReferenceGeometry::ReferenceGeometry(unsigned int numberOfNodes, unsigned int DIM, const TypeOfReferenceGeometry& geoT):
        points_(numberOfNodes,DIM),
        geometryType_(geoT)
    {
        
    }
    ReferenceGeometry::ReferenceGeometry(const ReferenceGeometry& other):
        points_(other.points_),
        geometryType_(other.geometryType_)
    {
        
    }
};
#endif

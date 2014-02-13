//
//  RectangularMeshDescriptor.h
//  
//
//  Created by Shavarsh Nurijanyan on 4/8/13.
//
//

#ifndef ____RectangularMeshDescriptor__
#define ____RectangularMeshDescriptor__

#include <iostream>
#include <Geometry/PointPhysical.hpp>
#include <vector>

namespace Base
{
    
    struct RectangularMeshDescriptor
    {
        RectangularMeshDescriptor(unsigned int DIM):
            bottomLeft_(DIM),
            topRight_(DIM),
            numElementsInDIM_(DIM),
            boundaryConditions_(DIM)
        {
        }
        enum Boundary{SOLID_WALL=0, PERIODIC=1};
        
        Geometry::PointPhysical   bottomLeft_;
        Geometry::PointPhysical   topRight_;
        
        std::vector<unsigned int>      numElementsInDIM_;
        std::vector<unsigned int>      boundaryConditions_;//x,y,z// according to this order!!!
    };
    
};

#endif /* defined(____RectangularMeshDescriptor__) */

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
    
    template <unsigned int DIM>
    struct RectangularMeshDescriptor
    {
        RectangularMeshDescriptor():
            bottomLeft_(),
            topLeft_(),
            numElementsInDIM_(DIM),
            boundaryConditions_(DIM)
        {
        }
        enum Boundary{SOLID_WALL=0, PERIODIC=1};
        
        Geometry::PointPhysical<DIM>   bottomLeft_;
        Geometry::PointPhysical<DIM>   topLeft_;
        
        std::vector<unsigned int>      numElementsInDIM_;
        std::vector<unsigned int>      boundaryConditions_;//x,y,z// according to this order!!!
    };
    
};

#endif /* defined(____RectangularMeshDescriptor__) */

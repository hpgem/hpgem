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
            numElementsInDIM_(DIM)
        {
        }
        
        Geometry::PointPhysical<DIM>   bottomLeft_;
        Geometry::PointPhysical<DIM>   topLeft_;
        
        std::vector<unsigned int>      numElementsInDIM_;
    };
    
};

#endif /* defined(____RectangularMeshDescriptor__) */

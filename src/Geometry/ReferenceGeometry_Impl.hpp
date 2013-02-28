//
//  ReferenceGeometry_Impl.hpp
//  
//
//  Created by Shavarsh Nurijanyan on 2/4/13.
//
//

#ifndef _ReferenceGeometry_Impl_hpp
#define _ReferenceGeometry_Impl_hpp

namespace Geometry
{
    template<unsigned int DIM>
    class ReferenceGeometry;
    
    template<unsigned int DIM>
    ReferenceGeometry<DIM>::ReferenceGeometry(const TypeOfReferenceGeometry& geoT):
        geometryType_(geoT)
    {
        
    }
    
    template<unsigned int DIM>
    ReferenceGeometry<DIM>::ReferenceGeometry(unsigned int numberOfNodes, const TypeOfReferenceGeometry& geoT):
        points_(numberOfNodes),
        geometryType_(geoT)
    {
        
    }
    template<unsigned int DIM>
    ReferenceGeometry<DIM>::ReferenceGeometry(const ReferenceGeometry& other):
        points_(other.points_),
        geometryType_(other.geometryType_)
    {
        
    }
};
#endif

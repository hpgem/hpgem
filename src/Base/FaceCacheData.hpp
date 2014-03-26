//------------------------------------------------------------------------------
// File: ElementCacheDataBase.hpp
// MT Julianto, Sun Feb 17 10:32:14 WET 2013
//------------------------------------------------------------------------------
#ifndef FaceCacheData_hpp
#define FaceCacheData_hpp

#include "Base/L2Norm.hpp"
#include "Geometry/PointReference.hpp"
#include "Geometry/PointPhysical.hpp"

namespace Base
{
    // forward declaration
    class Face;

  
    struct FaceCacheData
    {
    	FaceCacheData(unsigned int DIM):Normal(DIM){}
      // cache data
      LinearAlgebra::NumericalVector Normal;
      double L2Normal;

      // calculate the cache data
      void operator()(const Base::Face& fa,
                      const Geometry::PointReference& p);

    };
}

#endif

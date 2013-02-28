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
    template <unsigned int DIM>
    class Face;

  
    template <unsigned int dim>
    struct FaceCacheData
    {
      // cache data
      Geometry::PointPhysical<dim> Normal;
      double L2Normal;

      // calculate the cache data
      void operator()(const Base::Face<dim>& fa,
                      const Geometry::PointReference<dim-1>& p)
      {
          fa.getNormalVector(p, Normal);
          L2Normal = Base::L2Norm<dim>(Normal);
      }

    };
}

#endif

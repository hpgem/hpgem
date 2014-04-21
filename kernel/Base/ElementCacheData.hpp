//------------------------------------------------------------------------------
// File: ElementCacheDataBase.hpp
// MT Julianto, Sun Feb 17 10:32:14 WET 2013
//------------------------------------------------------------------------------
#ifndef ElementCacheData_hpp
#define ElementCacheData_hpp

#include "Geometry/Jacobian.hpp"
#include "Geometry/PointReference.hpp"

namespace Base
{
    // forward declaration
    class Element;

    class ElementCacheData
    {
        public:

        // calculate the cache data
        void operator()(const Base::Element* el, const Geometry::PointReference& p);

        // cache data
        double absDetJac_;
    };
}

#endif

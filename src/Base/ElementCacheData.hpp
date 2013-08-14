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
    template <unsigned int DIM>
    class Element;

    template <unsigned int dim>
    class ElementCacheData
    {
        public:

        // calculate the cache data
        void operator()(const Base::Element<dim>* el, const Geometry::PointReference<dim>& p)
        {
            Geometry::Jacobian<dim,dim> jac;
            el->calcJacobian(p, jac);
            absDetJac_ = std::abs(jac.determinant());
        }

        // cache data
        double absDetJac_;
    };
}

#endif

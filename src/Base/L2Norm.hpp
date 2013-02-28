//------------------------------------------------------------------------------
// File: norm2.hh
// Compute the 2-norm of an object of type Vector<NumType, dim>::Type.
// Lars Pesch, 27/08/2004
//------------------------------------------------------------------------------
#ifndef L2Norm_hpp
#define L2Norm_hpp
//------------------------------------------------------------------------------
// Package includes:
#include "../LinearAlgebra/NumericalVector.hpp"
#include "Geometry/PointPhysical.hpp"
//------------------------------------------------------------------------------
namespace Base
{
    /*! Compute the 2 norm of a vector. */
    template<unsigned int dim>
    double L2Norm(const LinearAlgebra::NumericalVector&);

    template<unsigned int dim>
    double L2Norm(const Geometry::PointPhysical<dim>&);
}
#endif
//------------------------------------------------------------------------------
// Local variables:
// mode:c++
// comment-column: 48
// End:

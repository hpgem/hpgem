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

///\todo merge with duplicate file Norm2.hpp
namespace Base
{
    /*! Compute the 2 norm of a vector. */
    double L2Norm(const LinearAlgebra::NumericalVector&);

    double L2Norm(const Geometry::PointPhysical&);
}
#endif
//------------------------------------------------------------------------------
// Local variables:
// mode:c++
// comment-column: 48
// End:

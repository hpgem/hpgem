//------------------------------------------------------------------------------
// File: L2Norm.cc
//------------------------------------------------------------------------------
// System includes and names imported from them:
#include <cmath>
//------------------------------------------------------------------------------
#include "L2Norm.hpp"
//------------------------------------------------------------------------------
namespace Base
{
    /*template <>
    double L2Norm<1>(const LinearAlgebra::NumericalVector& v)
    {
        return std::abs(v[0]);
    }
    
    template <>
    double L2Norm<2>(const LinearAlgebra::NumericalVector& v)
    {
        return std::sqrt(v[0] * v[0] + v[1] * v[1]);
    }

    template <>
    double L2Norm<3>(const LinearAlgebra::NumericalVector& v)
    {
        return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }*/

    double L2Norm(const LinearAlgebra::NumericalVector& v)
    {
    	double retSquared(0);
    	for(int i=0;i<v.size();++i){
    		retSquared+=v[i]*v[i];
    	}
    	return std::sqrt(retSquared);
    }
    
    /*template <>
    double L2Norm<1>(const Geometry::PointPhysical<1>& v)
    {
        return std::abs(v[0]);
    }
    
    template <>
    double L2Norm<2>(const Geometry::PointPhysical<2>& v)
    {
        return std::sqrt(v[0] * v[0] + v[1] * v[1]);
    }

    template <>
    double L2Norm<3>(const Geometry::PointPhysical<3>& v)
    {
        return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }*/

    ///\deprecated
    double L2Norm(const Geometry::PointPhysical& v)
    {
    	double retSquared(0);
    	for(int i=0;i<v.size();++i){
    		retSquared+=v[i]*v[i];
    	}
    	return std::sqrt(retSquared);
    }
};


//
//  ReferenceGeometry_Impl.hpp
//  
//
//  Created by Shavarsh Nurijanyan on 2/4/13.
//
//

#include "ReferenceGeometry.hpp"

#ifndef _ReferenceGeometry_Impl_hpp
#define _ReferenceGeometry_Impl_hpp

namespace Geometry
{
    
    ReferenceGeometry::ReferenceGeometry(const TypeOfReferenceGeometry& geoT):
        geometryType_(geoT)
    {
        
    }
    
    ReferenceGeometry::ReferenceGeometry(unsigned int numberOfNodes, unsigned int DIM, const TypeOfReferenceGeometry& geoT):
        points_(numberOfNodes,DIM),
        geometryType_(geoT)
    {
        
    }
    ReferenceGeometry::ReferenceGeometry(const ReferenceGeometry& other):
        points_(other.points_),
        geometryType_(other.geometryType_)
    {
        
    }

    double
    ReferenceGeometry::getBasisFunctionValue(const Base::BaseBasisFunction* function,const PointReference& p)
    {//not faster than just computing (until we do high precision computations)
    /*	try{
    		return basisfunctionValues_[function].at(p);
    	}catch(std::out_of_range&){
    		basisfunctionValues_[function][p]=function->eval(p);
    		return basisfunctionValues_[function].at(p);
    	}*/
    	return function->eval(p);
    }

    void
    ReferenceGeometry::getBasisFunctionDerivative(const Base::BaseBasisFunction* function,const PointReference& p, NumericalVector& ret)
    {
    //	try{
    //		ret=basisfunctionDerivatives_[function].at(p);
    //	}catch(std::out_of_range&){
    //		basisfunctionDerivatives_[function][p].resize(ret.size());
    		function->evalDeriv(p,ret);
    //		basisfunctionDerivatives_[function].at(p)=ret;
    //	}
    }
};
#endif

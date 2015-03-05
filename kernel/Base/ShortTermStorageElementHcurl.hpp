//
//  ShortTermStorageElementHcurl.h
//
//
//  Created by Devashish  on 01/06/14.
//
//this file contains the basis functions and vector manipulations

#ifndef _ShortTermStorageElementHcurl_hpp
#define _ShortTermStorageElementHcurl_hpp


#include "Base/ShortTermStorageElementBase.hpp"
#include "Geometry/Jacobian.hpp"
#include "Base/Element.hpp"

namespace Base{
    // H curl conforming space
    class ShortTermStorageElementBase;
    
    class ShortTermStorageElementHcurl:public ShortTermStorageElementBase{
        
    public:
        ShortTermStorageElementHcurl(std::size_t dimension) : ShortTermStorageElementBase(dimension){}
        
        virtual void computeData();
        
        
        //virtual double basisFunction(unsigned int i, const PointReferenceT& p);
        //virtual double basisFunction(unsigned int i, const PointReferenceT& p) const;
         
        virtual void basisFunction(std::size_t i, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret);
        virtual void basisFunction(std::size_t i, const PointReferenceT& P, LinearAlgebra::NumericalVector& ret) const;
        
        
        virtual void basisFunctionCurl(std::size_t i, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret);
        virtual void basisFunctionCurl(std::size_t i, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret) const;
        

    };
    
    
    
    
    
    
}


#endif /*ShortTermStorageElementHcurl.hpp*/




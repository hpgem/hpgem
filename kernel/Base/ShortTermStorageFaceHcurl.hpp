//
//  ShortTermStorageFaceHcurl.h
//  
//
//  Created by Devashish  on 01/06/14.
//
//

#ifndef ____ShortTermStorageFaceHcurl__hpp
#define ____ShortTermStorageFaceHcurl__hpp

#include "Base/ShortTermStorageFaceBase.hpp"
#include "Geometry/Jacobian.hpp"
#include "Base/Face.hpp"

namespace Base{
    
    //H curl conforming transformation for face
    
    
    class ShortTermStorageFaceHcurl:public ShortTermStorageFaceBase{
    
    public:
        ShortTermStorageFaceHcurl(std::size_t dimension):ShortTermStorageFaceBase(dimension){}
        
        virtual void computeData();
        
        
        virtual void   basisFunction(std::size_t i, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret);
		virtual void   basisFunction(std::size_t i, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret) const;
        
		virtual void   basisFunctionNormal(std::size_t i, const LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret);
		virtual void   basisFunctionNormal(std::size_t i, const LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret) const;
        

        
        virtual void basisFunctionCurl(std::size_t i, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret);
        virtual void basisFunctionCurl(std::size_t i, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret) const;
        
        
    };
}

#endif /* ShortTermStorageFaceHcurl_hpp */
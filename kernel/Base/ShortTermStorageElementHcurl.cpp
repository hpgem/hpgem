//
//  ShortTermStorageElementHcurl.cpp
//  
//
//  Created by Devashish  on 01/06/14.
//
//

#include "Base/ShortTermStorageElementHcurl.hpp"
#include "ElementCacheData.hpp"
#include "Element.hpp"
#include "Geometry/PointReference.hpp"
#include "Geometry/ReferenceGeometry.hpp"
#include "Geometry/PointPhysical.hpp"

void Base::ShortTermStorageElementHcurl::computeData(){
    
    ShortTermStorageElementBase::computeData();// calculates the Jacobian in jac_
    int DIM = currentPoint_.size();
    Geometry::Jacobian  jacobianT(DIM, DIM), jacobian(DIM, DIM);
    basisFunctionValues_.resize(element_->getNrOfBasisFunctions());
    basisFunctionCurlValues_.resize(element_->getNrOfBasisFunctions());
    jacobian = jac_;
    jac_.inverse(jac_);
    for (int i = 0;i < DIM ; ++i)
    {
        for(int j = 0; j < DIM; ++j)
        {
            jacobianT(i, j) = jac_(j, i);
        }
    }
    /*
     Inverse(0,0) = jac_(1,1) * jac_(2,2) - jac_(1,2) * jac_(2,1);
     Inverse(1,0) = jac_(2,1) * jac_(0,2) - jac_(2,2) * jac_(0,1);
     Inverse(2,0) = jac_(0,1) * jac_(1,2) - jac_(0,2) * jac_(1,1);
     Inverse(0,1) = jac_(1,2) * jac_(2,0) - jac_(1,0) * jac_(2,2);
     Inverse(1,1) = jac_(2,2) * jac_(0,0) - jac_(2,0) * jac_(0,2);
     Inverse(2,1) = jac_(0,2) * jac_(1,0) - jac_(0,0) * jac_(1,2);
     Inverse(0,2) = jac_(1,0) * jac_(2,1) - jac_(1,1) * jac_(2,0);
     Inverse(1,2) = jac_(2,0) * jac_(0,1) - jac_(2,1) * jac_(0,0);
     Inverse(2,2) = jac_(0,0) * jac_(1,1) - jac_(0,1) * jac_(1,0);
     
     */
    
    for(int i = 0; i < element_->getNrOfBasisFunctions(); ++i)
    {
        basisFunctionValues_[i].resize(3);
        element_->Element::basisFunction(i, currentPoint_, basisFunctionValues_[i]); //basisFunctionValues_[i][0] = element_->basisFunction(i, currentPoint_);
        basisFunctionValues_[i] = jacobianT  * basisFunctionValues_[i];  //basisFunctionValues_[i][0] = jac_ * basisFunctionValues_[i][0];
        
        basisFunctionCurlValues_[i].resize(3);
        element_->Element::basisFunctionCurl(i, currentPoint_, basisFunctionCurlValues_[i]);
        basisFunctionCurlValues_[i] = (jacobian / (std::abs(jacobian.determinant()))) * basisFunctionCurlValues_[i];
    }
    
}



void Base::ShortTermStorageElementHcurl::basisFunction(std::size_t i, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret)
{
    //std::cout<<"Basis Function called in Hcurl"<<std::endl;
    if(!(p==currentPoint_))
    {
        currentPoint_=p;
        computeData();
    }
    ret = basisFunctionValues_[i];
}


void Base::ShortTermStorageElementHcurl::basisFunction(std::size_t i, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret) const
{
    //std::cout<<"Basis Function called in Hcurl const"<<std::endl;
    ret = basisFunctionValues_[i];
    if(!(p == currentPoint_))
    {
        std::cout<<"Warning : The operator being used is slow";
        element_->basisFunction(i, p, ret);
    }
}

void Base::ShortTermStorageElementHcurl::basisFunctionCurl(std::size_t i, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret)
{
    if(!(p==currentPoint_))
    {
        currentPoint_ = p;
        computeData();
    }
    ret = basisFunctionCurlValues_[i]; //define basisFunctionCurlValues_ as basisFunctionDerivatives_ from H1 function
    
}


void Base::ShortTermStorageElementHcurl::basisFunctionCurl(std::size_t i, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret) const
{
    ret = basisFunctionCurlValues_[i];
    if(!(p == currentPoint_))
    {
        std::cout<<"Warning: The operator being used is slow";
        element_->basisFunctionCurl(i, p, ret); // removed the argument this
    }

}


















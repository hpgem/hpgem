/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found below.


 Copyright (c) 2014, University of Twente
 All rights reserved.

 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "Base/ShortTermStorageElementHCurl.h"
#include "ElementCacheData.h"
#include "Element.h"
#include "Geometry/PointReference.h"
#include "Geometry/ReferenceGeometry.h"
#include "Geometry/PointPhysical.h"

void Base::ShortTermStorageElementHcurl::computeData()
{
    
    ShortTermStorageElementBase::computeData(); // calculates the Jacobian in jac_
    std::size_t DIM = currentPoint_.size();
    Geometry::Jacobian jacobianT(DIM, DIM), jacobian(DIM, DIM);
    basisFunctionValues_.resize(element_->getNrOfBasisFunctions());
    basisFunctionCurlValues_.resize(element_->getNrOfBasisFunctions());
    jacobian = jac_;
    jac_ = jac_.inverse();
    for (std::size_t i = 0; i < DIM; ++i)
    {
        for (std::size_t j = 0; j < DIM; ++j)
        {
            jacobianT(i, j) = jac_(j, i);
        }
    }

    for (std::size_t i = 0; i < element_->getNrOfBasisFunctions(); ++i)
    {
        basisFunctionValues_[i].resize(3);
        //Compute the values of basis function i on the current point. Note that we return in the last parameter.
        element_->Element::basisFunction(i, currentPoint_, basisFunctionValues_[i]);
        //Multiplies the values of basis function i with the transposed inverse Jacobian.
        basisFunctionValues_[i] = jacobianT * basisFunctionValues_[i];
                
        basisFunctionCurlValues_[i].resize(3);
        basisFunctionCurlValues_[i] = element_->Element::basisFunctionCurl(i, currentPoint_);
        basisFunctionCurlValues_[i] = (jacobian / (std::abs(jacobian.determinant()))) * basisFunctionCurlValues_[i];
    }
    
}

void Base::ShortTermStorageElementHcurl::basisFunction(std::size_t i, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret)
{
    logger(DEBUG, "basisFunction called in Hcurl");
    if (!(p == currentPoint_))
    {
        currentPoint_ = p;
        computeData();
    }
    ret = basisFunctionValues_[i];
}

void Base::ShortTermStorageElementHcurl::basisFunction(std::size_t i, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret) const
{
    logger(DEBUG, "Basis Function called in Hcurl const");
    ret = basisFunctionValues_[i];
    if (!(p == currentPoint_))
    {
        logger(WARN, "WARNING: you are using a slow operator");
        element_->basisFunction(i, p, ret);
    }
}

void Base::ShortTermStorageElementHcurl::basisFunctionCurl(std::size_t i, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret)
{
    if (!(p == currentPoint_))
    {
        currentPoint_ = p;
        computeData();
    }
    //define basisFunctionCurlValues_ as basisFunctionDerivatives_ from H1 function
    ret = basisFunctionCurlValues_[i];
            
}

void Base::ShortTermStorageElementHcurl::basisFunctionCurl(std::size_t i, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret) const
{
    ret = basisFunctionCurlValues_[i];
    if (!(p == currentPoint_))
    {
        logger(WARN, "WARNING: you are using a slow operator");
        ret = element_->basisFunctionCurl(i, p);
    }
    
}


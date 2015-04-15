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

#include "Base/ShortTermStorageFaceH1.h"
#include "Base/ShortTermStorageElementH1.h"
#include "Geometry/PhysicalGeometry.h"
#include "L2Norm.h"
#include "Geometry/PointPhysical.h"
#include "FaceCacheData.h"
#include "ElementCacheData.h"

void Base::ShortTermStorageFaceH1::computeData()
{
    ShortTermStorageFaceBase::computeData();
    std::size_t n = face_->getNrOfBasisFunctions();
    basisFunctionValues_.resize(n);
    basisFunctionDerivatives_.resize(n);
    basisFunctionsTimesNormal_.resize(n);
    ShortTermStorageElementBase* elementwrapper = new ShortTermStorageElementH1(currentPoint_.size() + 1);
    std::size_t leftFunctions = getPtrElementLeft()->getNrOfBasisFunctions();
    *elementwrapper = *getPtrElementLeft();
    Geometry::PointReference pElement = mapRefFaceToRefElemL(currentPoint_);
    double norm = Base::L2Norm(normal_);
    for (std::size_t i = 0; i < leftFunctions; ++i)
    {
        basisFunctionValues_[i].resize(1);
        basisFunctionValues_[i][0] = elementwrapper->basisFunction(i, pElement);
        basisFunctionsTimesNormal_[i].resize(currentPoint_.size() + 1);
        basisFunctionsTimesNormal_[i] = normal_;
        basisFunctionsTimesNormal_[i] *= basisFunctionValues_[i][0] / norm;
        basisFunctionDerivatives_[i] = elementwrapper->basisFunctionDeriv(i, pElement);
    }
    if (n > leftFunctions)
    {
        *elementwrapper = *getPtrElementRight();
        pElement = mapRefFaceToRefElemR(currentPoint_);
    }
    for (std::size_t i = leftFunctions; i < n; ++i)
    {
        basisFunctionValues_[i].resize(1);
        basisFunctionValues_[i][0] = elementwrapper->basisFunction(i - leftFunctions, pElement);
        basisFunctionsTimesNormal_[i].resize(currentPoint_.size() + 1);
        basisFunctionsTimesNormal_[i] = normal_;
        basisFunctionsTimesNormal_[i] *= -basisFunctionValues_[i][0] / norm;
        basisFunctionDerivatives_[i] = elementwrapper->basisFunctionDeriv(i - leftFunctions, pElement);
    }
    delete elementwrapper;
}

double Base::ShortTermStorageFaceH1::basisFunction(std::size_t i, const Geometry::PointReference& p)
{
    if (!(currentPoint_ == p))
    {
        currentPoint_ = p;
        computeData();
    }
    return basisFunctionValues_[i][0];
}

double Base::ShortTermStorageFaceH1::basisFunction(std::size_t i, const Geometry::PointReference& p) const
{
    if (!(currentPoint_ == p))
    {
        logger(WARN, "Warning: you are using slow data access");
        return face_->basisFunction(i, p);
    }
    return basisFunctionValues_[i][0];
}

void Base::ShortTermStorageFaceH1::basisFunction(std::size_t i, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret)
{
    if (!(currentPoint_ == p))
    {
        currentPoint_ = p;
        computeData();
    }
    ret = basisFunctionValues_[i];
}

void Base::ShortTermStorageFaceH1::basisFunction(std::size_t i, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret) const
{
    ret = basisFunctionValues_[i];
    if (!(currentPoint_ == p))
    {
        logger(WARN, "Warning: you are using slow data access");
        face_->basisFunction(i, p, ret);
    }
}

LinearAlgebra::NumericalVector Base::ShortTermStorageFaceH1::basisFunctionNormal(std::size_t i, const LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p)
{
    if (!(currentPoint_ == p))
    {
        currentPoint_ = p;
        computeData();
    }
    //normal is passed only because in usual situations you do not want to recompute it for every function
    return basisFunctionsTimesNormal_[i];
}

LinearAlgebra::NumericalVector Base::ShortTermStorageFaceH1::basisFunctionNormal(std::size_t i, const LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p) const
{
    if (!(currentPoint_ == p))
    {
        logger(WARN, "Warning: you are using slow data access");
        return face_->basisFunctionNormal(i, normal, p);
    }
    return basisFunctionsTimesNormal_[i];
}

LinearAlgebra::NumericalVector Base::ShortTermStorageFaceH1::basisFunctionDeriv(std::size_t i, const Geometry::PointReference& p)
{
    if (!(currentPoint_ == p))
    {
        currentPoint_ = p;
        computeData();
    }
    return basisFunctionDerivatives_[i];
}

LinearAlgebra::NumericalVector Base::ShortTermStorageFaceH1::basisFunctionDeriv(std::size_t i, const Geometry::PointReference& p) const
{
    if (!(currentPoint_ == p))
    {
        logger(WARN, "Warning: you are using slow data access");
        return face_->basisFunctionDeriv(i, p);
    }
    return basisFunctionDerivatives_[i];
}

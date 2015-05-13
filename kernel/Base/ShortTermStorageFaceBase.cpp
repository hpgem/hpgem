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

#include <limits>

#include "ShortTermStorageFaceBase.h"
#include "Integration/QuadratureRules/GaussQuadratureRule.h"
#include "FaceCacheData.h"
#include "Logger.h"

Base::Face& Base::ShortTermStorageFaceBase::operator=(Base::Face& face)
{
    logger.assert(this != &face, "Trying to assign a Face of the type ShortTermStorageFaceBase to itself.");
    face_ = &face;
    if (currentPoint_->size() == 0)
    {
        computeData();
    }
    else
    {
        currentPoint_ = Geometry::PointReferenceFactory::instance()->makePoint({std::numeric_limits<double>::quiet_NaN()});
    }
    currentPointIndex_ = -1;
    return *this;
}

void Base::ShortTermStorageFaceBase::computeData()
{
    if (useCache_)
    {
        std::vector<FaceCacheData>& cache = face_->getVecCacheData();
        if (recomputeCache_ || (cache.size() != getGaussQuadratureRule()->nrOfPoints()))
        {
            recomputeCacheOff();
            std::size_t n = getGaussQuadratureRule()->nrOfPoints();
            for (std::size_t i = 0; i < n; ++i)
            {
                cache[i](*face_, getGaussQuadratureRule()->getPoint(i));
            }
        }
        currentPointIndex_++;
        normal_ = face_->getNormalVector(*currentPoint_);
    }
    else
    {
        normal_ = face_->getNormalVector(*currentPoint_);
    }
}

LinearAlgebra::NumericalVector Base::ShortTermStorageFaceBase::getNormalVector(const ReferencePointT& pRefFace) const
{
    if (!(currentPoint_ == &pRefFace))
    {
        logger(WARN, "WARNING: you are using slow data access.");
        return face_->getNormalVector(pRefFace);
    }
    return normal_;
}

LinearAlgebra::NumericalVector Base::ShortTermStorageFaceBase::getNormalVector(const ReferencePointT& pRefFace)
{
    if (!(currentPoint_ == &pRefFace))
    {
        currentPoint_ = &pRefFace;
        computeData();
    }
    return normal_;
}

Geometry::PointPhysical Base::ShortTermStorageFaceBase::referenceToPhysical(const Geometry::PointReference& pointReference) const
{
    return face_->referenceToPhysical(pointReference);
}

void Base::ShortTermStorageFaceBase::cacheOn()
{
    useCache_ = true;
}

void Base::ShortTermStorageFaceBase::cacheOff()
{
    useCache_ = false;
}

void Base::ShortTermStorageFaceBase::recomputeCacheOn()
{
    recomputeCache_ = true;
}

void Base::ShortTermStorageFaceBase::recomputeCacheOff()
{
    recomputeCache_ = false;
}

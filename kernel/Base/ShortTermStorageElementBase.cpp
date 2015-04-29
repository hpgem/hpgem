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

#include "ShortTermStorageElementBase.h"

#include "Integration/QuadratureRules/GaussQuadratureRule.h"
#include "ElementCacheData.h"

#include <limits>

Base::Element& Base::ShortTermStorageElementBase::operator=(Base::Element& element)
{
    logger.assert(this != &element, "Trying to assign an Element of the type ShortTermStorageElementBase to itself.");

    element_ = &element;

    currentPoint_ = Geometry::PointReferenceFactory::instance()->makePoint({std::numeric_limits<double>::quiet_NaN()});
    currentPointIndex_ = -1;
    return *this;
}

Geometry::PointPhysical Base::ShortTermStorageElementBase::referenceToPhysical(const PointReferenceT& pointReference) const
{
    return element_->referenceToPhysical(pointReference);
}

Geometry::Jacobian Base::ShortTermStorageElementBase::calcJacobian(const PointReferenceT& pointReference) const
{
    if (!(currentPoint_ == &pointReference))
    {
        logger(WARN, "WARNING: you are using slow data access by using ShortTermStorageElementBase::calcJacobian const, use the non-const version instead.");
        return element_->calcJacobian(pointReference);
    }
    return jac_;
}

Geometry::Jacobian Base::ShortTermStorageElementBase::calcJacobian(const PointReferenceT& pointReference)
{
    if (!(currentPoint_ == &pointReference))
    {
        currentPoint_ = &pointReference;
        computeData();
    }
    return jac_;
}

void Base::ShortTermStorageElementBase::computeData()
{
    if (useCache_)
    {
        std::vector<ElementCacheData>& cache = element_->getVecCacheData();
        if (recomputeCache_ || (cache.size() != getGaussQuadratureRule()->nrOfPoints()))
        {
            recomputeCacheOff();
            std::size_t n = getGaussQuadratureRule()->nrOfPoints();
            for (std::size_t i = 0; i < n; ++i)
            {
                cache[i](element_, getGaussQuadratureRule()->getPoint(i));
            }
        }
        currentPointIndex_++;
        jac_ = element_->calcJacobian(*currentPoint_);
    }
    else
    {
        jac_ = element_->calcJacobian(*currentPoint_);
    }
}

void Base::ShortTermStorageElementBase::cacheOn()
{
    useCache_ = true;
}

void Base::ShortTermStorageElementBase::cacheOff()
{
    useCache_ = false;
}

void Base::ShortTermStorageElementBase::recomputeCacheOn()
{
    recomputeCache_ = true;
}

void Base::ShortTermStorageElementBase::recomputeCacheOff()
{
    recomputeCache_ = false;
}

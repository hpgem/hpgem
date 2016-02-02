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

#include "Geometry/PointReference.h"
#include "Geometry/Jacobian.h"

#ifndef MAPPINGSIMPLEXLINEAR_CPP_
#define MAPPINGSIMPLEXLINEAR_CPP_

/*! The mapping is linear in every reference space dimension, with the offset to
 *  the origin of the dim-simplex. */
template<std::size_t DIM>
Geometry::PointPhysical<DIM> Geometry::MappingToPhysSimplexLinear<DIM>::transform(const PointReference<DIM>& pointReference) const
{
    logger.assert(pointReference.size()==DIM, "Reference point has the wrong dimension");
    Geometry::PointPhysical<DIM> pointPhysical = geometry->getLocalNodeCoordinates(0);
    for (std::size_t i = 1; i <= DIM; ++i)
    {
        Geometry::PointPhysical<DIM> next = geometry->getLocalNodeCoordinates(i);
        next = next - geometry->getLocalNodeCoordinates(0);
        pointPhysical += pointReference[i - 1] * next;
    }
    return pointPhysical;
}

template<std::size_t DIM>
Geometry::PointReference<DIM> Geometry::MappingToPhysSimplexLinear<DIM>::inverseTransform(const PointPhysical<DIM>& pointPhysical) const
{
    LinearAlgebra::SmallVector<DIM> offSet = (pointPhysical - geometry->getLocalNodeCoordinates(0)).getCoordinates();
    calcJacobian({}).solve(offSet);
    return PointReference<DIM>{offSet};
}

/*! The Jacobian results from the transform function by symbolic derivation. */
template<std::size_t DIM>
Geometry::Jacobian<DIM, DIM> Geometry::MappingToPhysSimplexLinear<DIM>::calcJacobian(const PointReference<DIM>& pointReference) const
{
    logger.assert(pointReference.size()==DIM, "Reference point has the wrong dimension");
    Jacobian<DIM, DIM> jacobian;
    const Geometry::PointPhysical<DIM>& first = geometry->getLocalNodeCoordinates(0);
    for (std::size_t i = 0; i < DIM; i++)
    {
        const PointPhysical<DIM>& point = geometry->getLocalNodeCoordinates(i + 1);
        for (std::size_t j = 0; j < DIM; j++)
        {
            jacobian(j, i) = point[j] - first[j];
        }
    }
    return jacobian;
}

/*! For the simplices it actually saves computations to save the coefficients.*/
template<std::size_t DIM>
void Geometry::MappingToPhysSimplexLinear<DIM>::reinit()
{
    //we get a free reinit because transform and jacobian use the nodes directly
}
#endif /* MAPPINGSIMPLECUBENLINEAR_H_ */

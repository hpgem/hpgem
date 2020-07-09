/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2014, University of Twente
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors
 may be used to endorse or promote products derived from this software without
 specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef HPGEM_KERNEL_MAPPINGINTERFACE_H
#define HPGEM_KERNEL_MAPPINGINTERFACE_H

#include <cstdlib>
#include "Logger.h"
#include "Geometry/Jacobian.h"

namespace Geometry {
template <std::size_t DIM>
class PointReference;
template <std::size_t dimFrom, std::size_t dimTo>
class Jacobian;

/*!
 \brief (OC): Base class that defines the functionality of a mapping between two
 coordinate systems.

 Mappings transform coordinates of a point with respect to one system to
 coordinates in another system. This can be to obtain physical coordinates
 for points in a reference geometry, or to map points from a face to one of
 the neighbouring elements' coordinate systems, etc. <BR>

 Using mappings, calculations can be carried out in reference space;
 e.g. quadrature rules are defined on the reference geometry, while some of
 the functions that are to be integrated live in physical space.<BR>

 Apart from mapping points to (physical) space, a Mapping gives information
 about the relative layout of the connected spaces, namely in its Jacobian.
 This information also has to be provided by implementation classes.<BR>

 End of interface description.

 Some remarks:
 <ul>
 <li> At the moment I expect the assignment of a certain Mapping
 implementation to an Element to be done by a function
 makeMapping() in ElementFactory. This will have to be
 user-implemented as the choice of the order, use of curved
 elements near the boundary etc. are problem-dependent.
 <li> In general one might expect a Mapping to be fully
 characterized by the dimension of the space(s) and the order.
 In that respect one may wonder that the Mapping implemenations
 are adapted to the reference geometries. This however allows
 them to be tailored to their actual purpose. Connected to this
 is the choice NOT to implement the transform and calcJacobian
 functions using shape functions. While they could make the
 representation (especially for higher order cases) simpler
 it would probably lead to slower computation: The
 \f$ \sum_{i=1}^N b_i(\xi) \vec{x}_i\f$ would need \f$ N\f$
 calls to basis function objects. (Note that such an
 implementation would require a Mapping to store the
 NodeList-indices of the vertices (and additional points) it uses
 for its mapping, so it would have to be rebuild after a
 compaction step, but would look up the \f$\vec{x}_i\f$
 coordinates for each evaluation, so it would not matter if these
 coordinates changed; with the current implementation
 it is just the other way round: it does not care about the
 compaction step, but it has to be reinitialized after changing
 vertex coordinates. I would expect most problems to work with
 fixed elements, so I deem this the more efficient approach.
 Additionally I expect the performance to be better since the
 necessary vectors are stored in Mapping, hence no need to
 interfere with a list that is somewhere else in memory; another
 minor pro is that the transformation formulae can be written
 (partially) in Horner form, so that some multiplications are
 saved. This effect is small though, as these are multiplications
 of scalars, while a MUCH bigger proportion will go into the
 operations on Point objects. I would expect an
 improvement of runtime performance using expression templates
 (especially if higher than first order is used).
 <li> Also about the previous item: possibly the call to a ValidPoint
 function of the reference geometry could hence be avoided, but I
 am not sure at the current stage (depends on whether some
 mappings can be used for different kinds of geometries).
 <li> This is templated on the change in dimension, where a positive numbers
 means mapping from a small-dimensional space (reference face) to a
 larger-dimensional space (reference element)
 </ul> */
template <int codim>
class MappingInterface {
   public:
    /*! (OC): Jacobian has a gradient in each line, hence as many lines as
     target space (DIM2) and as many columns as original space (DIM1),
     \frac{\partial x_i}{\partial \xi_j}. */
    virtual Jacobian<0, 0 + codim> calcJacobian(
        const PointReference<0> &) const {
        logger(ERROR, "Passed a point of the wrong dimension");
        return Jacobian<0, 0 + codim>();
    }

    virtual Jacobian<1, 1 + codim> calcJacobian(
        const PointReference<1> &) const {
        logger(ERROR, "Passed a point of the wrong dimension");
        return Jacobian<1, 1 + codim>();
    }

    virtual Jacobian<2, 2 + codim> calcJacobian(
        const PointReference<2> &) const {
        logger(ERROR, "Passed a point of the wrong dimension");
        return Jacobian<2, 2 + codim>();
    }

    virtual Jacobian<3, 3 + codim> calcJacobian(
        const PointReference<3> &) const {
        logger(ERROR, "Passed a point of the wrong dimension");
        return Jacobian<3, 3 + codim>();
    }

    virtual Jacobian<4, 4 + codim> calcJacobian(
        const PointReference<4> &) const {
        logger(ERROR, "Passed a point of the wrong dimension");
        return Jacobian<4, 4 + codim>();
    }

    MappingInterface() = default;
    MappingInterface(const MappingInterface &other) = default;  // does nothing

    virtual ~MappingInterface() = default;

    /// concatenated mapping needs to know what kind of intermediate point to
    /// create
    virtual std::size_t getTargetDimension() const = 0;
};
}  // namespace Geometry
#endif // HPGEM_KERNEL_MAPPINGINTERFACE_H

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

#ifndef REFERENCETOPHYSICALM_H_
#define REFERENCETOPHYSICALM_H_

#include "MappingInterface.h"
#include "Logger.h"
#include "Geometry/PhysicalGeometryBase.h"
#include <vector>

namespace Geometry {
template <std::size_t DIM>
class PhysicalGeometry;
template <std::size_t DIM>
class PointPhysical;

/*! ~OC~
 Second layer abstract base class (derived from Mapping) for mappings
 that go from a reference element to physical space (hence different point
 types in the two spaces, cf. member function transform).

 In the current design of mappings, these do not store the point
 coordinates, but this can neither be enforced, nor am I sure that this is
 essential. This should be reconsidered at a later stage. The access to
 physical space information is via a reference to the NodeContainer;
 mappings can index into that one with the global node numbers of the
 elements.

 The reinit-function is meant to alert an object of a change of the layout
 of the Element in physical space. Since the reference geometry of the
 Element does not change, but rather only (some of) the vertex positions,
 the Mapping can be adjusted to the new layout.

 At the moment only reference to physical maps from reference elements are
 supported so there in no need to template this class
*/

class MappingReferenceToPhysical : public MappingInterface<0> {

   public:
    MappingReferenceToPhysical(const PhysicalGeometryBase* target)
        : MappingInterface(), geometry(target) {}

    // Note that the memory of nodes is managed by Mesh, so do not make a deep
    // copy.
    MappingReferenceToPhysical(const MappingReferenceToPhysical& other)
        : MappingInterface(other), geometry(other.geometry) {}

    // Methods.
    //! ~OC~ Transform a point from reference space to physical space.
    virtual PointPhysical<1> transform(const PointReference<1>&) const {
        logger(ERROR, "Passed a point of the wrong dimension");
        return PointPhysical<1>();
    }

    virtual PointPhysical<2> transform(const PointReference<2>&) const {
        logger(ERROR, "Passed a point of the wrong dimension");
        return PointPhysical<2>();
    }

    virtual PointPhysical<3> transform(const PointReference<3>&) const {
        logger(ERROR, "Passed a point of the wrong dimension");
        return PointPhysical<3>();
    }

    virtual PointPhysical<4> transform(const PointReference<4>&) const {
        logger(ERROR, "Passed a point of the wrong dimension");
        return PointPhysical<4>();
    }

    virtual PointReference<1> inverseTransform(const PointPhysical<1>&) const {
        logger(ERROR, "Passed a point of the wrong dimension");
        return {};
    }

    virtual PointReference<2> inverseTransform(const PointPhysical<2>&) const {
        logger(ERROR, "Passed a point of the wrong dimension");
        return {};
    }

    virtual PointReference<3> inverseTransform(const PointPhysical<3>&) const {
        logger(ERROR, "Passed a point of the wrong dimension");
        return {};
    }

    virtual PointReference<4> inverseTransform(const PointPhysical<4>&) const {
        logger(ERROR, "Passed a point of the wrong dimension");
        return {};
    }

    /// Recompute mapping after physical nodes have moved.
    /// Note that this typically has to be done for all elements, so make sure
    /// to use the global iterator to get all the elements when using parallel
    /// computing.
    virtual void reinit() = 0;

   protected:
    const PhysicalGeometryBase* geometry;  /// Pointer to the physical geometry
                                           /// (for reinitialisation)
};

}  // namespace Geometry
#endif /* REFERENCETOPHYSICALM_H_ */

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

#ifndef PHYSICALGEOMETRY_H_
#define PHYSICALGEOMETRY_H_

#include <vector>
#include <iostream>

#include "PointPhysical.h"
#include "PointReference.h"
#include "ReferenceGeometry.h"
#include "PhysicalGeometryBase.h"

namespace Geometry {
class ReferenceGeometry;
template <std::size_t DIM>
class PointPhysical;

/*!\class PhysicalGeometry
 * \brief PhysicalGeometry describes an actual physical shape in real space.
 * \details
 * You shouldn't create a particular Physical<Shape> (although it is possible),
 * but a PhysicalGeometry. The reference geometry passed to the constructor
 * selects what shape this PhysicalGeometry has
 *
 * It contains only the global indexes of its points in globalNodeIndexes_.
 * These global indexes refer to the global node container, of which every
 * PhysicalGeometry has a reference: nodeCoordinates_.
 *
 * It also contains a reference to the corresponding referenceGeometry.
 *
 * ~ Point is the name of a class.                          ~
 * ~ Node is a point that belongs to a geometry.            ~
 * ~ nodeCoordinate is the Point where the Node is located. ~
 */
template <std::size_t DIM>
class PhysicalGeometry : public PhysicalGeometryBase {
   public:
    /// \brief Constructor gets indexes of the nodes, a reference to the node
    /// container, and a pointer to the corresponding reference geometry.
    PhysicalGeometry(const std::vector<std::size_t>& globalNodeIndexes,
                     std::vector<PointPhysical<DIM> >& nodes,
                     const ReferenceGeometry* const refG)
        : PhysicalGeometryBase(globalNodeIndexes, refG),
          nodeCoordinates_(nodes) {
        logger.assert_debug(refG != nullptr,
                            "Invalid reference geometry passed");
    }

    PhysicalGeometry(const PhysicalGeometry& other) = delete;

    /// \brief Returns a pointer to the global container of nodes.
    std::vector<PointPhysical<DIM> >& getNodeCoordinates() {
        return nodeCoordinates_;
    }

    /// \brief Returns a constant pointer of the global container of nodes
    const std::vector<PointPhysical<DIM> >& getNodeCoordinates() const {
        return nodeCoordinates_;
    }

    /// \brief Given a global index, returns a pointer to the corresponding
    /// point.
    PointPhysicalBase* getNodeCoordinatePtr(const std::size_t globalIndex) {
        logger.assert_debug(globalIndex < nodeCoordinates_.size(),
                            "This mesh does not contain a node with index %",
                            globalIndex);
        return &(nodeCoordinates_[globalIndex]);
    }

    /// \brief Given a global index, returns a pointer to the corresponding
    /// point.
    const PointPhysicalBase* getNodeCoordinatePtr(
        const std::size_t globalIndex) const {
        logger.assert_debug(globalIndex < nodeCoordinates_.size(),
                            "This mesh does not contain a node with index %",
                            globalIndex);
        return &(nodeCoordinates_[globalIndex]);
    }

    /// \deprecated Not consistent with naming convention, please use
    /// getNodeCoordinates()
    std::vector<PointPhysical<DIM> >& getNodes() { return nodeCoordinates_; }

    /// \deprecated Not consistent with naming convention, please use
    /// getNodeCoordinates()
    std::vector<PointPhysical<DIM> >& getNodes() const {
        return nodeCoordinates_;
    }

    /// \deprecated Not consistent with naming convention, please use
    /// getNodeCoordinatePtr
    const PointPhysicalBase* getNodePtr(const std::size_t globalIndex) const {
        logger.assert_debug(globalIndex < nodeCoordinates_.size(),
                            "This mesh does not contain a node with index %",
                            globalIndex);
        return &(nodeCoordinates_[globalIndex]);
    }

    ///\deprecated Not consistent with naming convention, please use
    ///getNodeCoordinatePtr
    PointPhysicalBase* getNodePtr(const std::size_t globalIndex) {
        logger.assert_debug(globalIndex < nodeCoordinates_.size(),
                            "This mesh does not contain a node with index %",
                            globalIndex);
        return &(nodeCoordinates_[globalIndex]);
    }

    /// \brief Given a local index, return the physical coordinates of the
    /// corresponding point.
    const PointPhysicalBase& getLocalNodeCoordinates(
        const std::size_t localIndex) const;

    /// \brief Given a global index, return the physical coordinates of the
    /// corresponding point.
    const PointPhysicalBase& getGlobalNodeCoordinates(
        const std::size_t globalIndex) const;

    /// \brief Given a local index, return the physical coordinates of the
    /// corresponding point.
    PointPhysicalBase& getLocalNodeCoordinates(const std::size_t localIndex);

    /// \brief Given a global index, return the physical coordinates of the
    /// corresponding point.
    PointPhysicalBase& getGlobalNodeCoordinates(const std::size_t globalIndex);

   protected:
    /// Reference to the global node container.
    std::vector<PointPhysical<DIM> >& nodeCoordinates_;
};

}  // namespace Geometry

#include "PhysicalGeometry_Impl.h"

#endif /* PHYSICALGEOMETRY_H_ */

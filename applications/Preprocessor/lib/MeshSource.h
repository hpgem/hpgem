/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2020, University of Twente
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
#ifndef HPGEM_MESHSOURCE_H
#define HPGEM_MESHSOURCE_H

#include <vector>
#include "customIterator.h"

namespace Preprocessor {

/**
 * A generic source that can provide the mesh as input for the Preprocessor.
 *
 * This mesh source described the mesh that is the input for the Preprocessor.
 * It needs a listing of the nodes in the mesh (with coordinates) and the
 * elements to which they attach. This source is agnostic with respect to what
 * the actual origin is of this information.
 *
 */
class MeshSource {
   public:
    virtual ~MeshSource() = default;

    struct Node {
        /**
         * Periodically connected coordinates of this node. For example for a
         * mesh for the unit line segment with two elements one would have three
         * coordinates (0, 0.5 and 1). Due to periodicity the coordinates at
         * x=0 and 1 belong to the same node. For such a mesh the MeshSource
         * would thus get two Nodes, one for the middle node: {{0.5}} and one
         * for the periodic outer nodes: {{0}, {1}}.
         */
        std::vector<std::vector<double>> coordinates;

        bool operator==(const Node& o) const {
            return coordinates == o.coordinates;
        }
    };

    struct Element {

        Element() : coordinateIds(), zoneName("Main"){};

        /**
         * Indices of the coordinates of the vertices of this element. The Ids
         * follow that of the coordinates, not of the nodes (that may bin
         * periodically connected coordinates).
         *
         * The shape of the element is determined by the number of coordinates
         * and the dimension of the mesh. The ordering of the nodes should
         * follow the standard hpgem ordering.
         */

        std::vector<std::size_t> coordinateIds;

        /**
         * Zone specifier of this element, for example for differentiating zones
         * with different materials.
         */
        std::string zoneName;

        bool operator==(const Element& o) const {
            return coordinateIds == o.coordinateIds;
        }
    };

    /**
     * Range over the nodes.
     * @return The range providing the nodes in the mesh
     */
    virtual Range<Node> getNodeCoordinates() = 0;

    /**
     * Range over all elements in the mesh.
     * @return The range providing the elements in the mesh
     */
    virtual Range<Element> getElements() = 0;

    /**
     * Get the dimension of the described mesh
     * @return The dimension
     */
    virtual std::size_t getDimension() const = 0;
};

/**
 * A generic source that can provide the mesh as input for the Preprocessor.
 * (replacing MeshSource)
 *
 * This mesh source described the mesh that is the input for the Preprocessor.
 * It needs a listing of the coordinates in the mesh (with associated node) and
 * the elements to which they attach. This source is agnostic with respect to
 * what the actual origin is of this information.
 */
class MeshSource2 {
   public:
    virtual ~MeshSource2() = default;

    struct Coord {
        /// Actual physical coordinate of the node
        std::vector<double> coordinate;
        /**
         * Identifier for the node to which this coordinate belongs. Multiple
         * coordinates could associated by the same node, for example in the
         * case of periodic boundary conditions. As example, consider the unit
         * line segment with two line segments and three points
         *   0 ---- 0.5 --- 1
         * For a periodic line segment the points at 0 and 1 should be
         * connected by giving them the same node id..
         */
        std::size_t nodeId;
    };

    struct Element {
        /**
         * Indices of the coordinates of the vertices of this element. The Ids
         * follow that of the coordinates.
         *
         * The shape of the element is determined by the number of coordinates
         * and the dimension of the mesh. The ordering of the nodes should
         * follow the standard hpgem ordering.
         */
        std::vector<std::size_t> coordinateIds;

        /**
         * Zone specifier of this element, for example for differentiating zones
         * with different materials.
         */
        std::string zoneName;

        /**
         * Dimension of the Element, this is not used in the generation of the
         * hpgem intermediate format, but maybe useful for parsing
         */
        size_t dimension;
        /**
         * Index of the Element, this is not used in the generation of the hpgem
         * intermediate format, but maybe useful for parsing
         */
        size_t id;

        bool operator==(const Element& o) const {
            return coordinateIds == o.coordinateIds && zoneName == o.zoneName;
        }
    };

    /**
     * Reference to a list of coordinates in the mesh.
     * @return The list of coordinates.
     */
    virtual const std::vector<Coord>& getCoordinates() = 0;
    /**
     * The elements in the mesh.
     * @return A list of the elements.
     */
    virtual const std::vector<Element>& getElements() = 0;

    /**
     * Get the dimension of the described mesh
     * @return The dimension
     */
    virtual std::size_t getDimension() const = 0;
};

}  // namespace Preprocessor

#endif  // HPGEM_MESHSOURCE_H

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

#include <map>
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
        Coord() = default;
        Coord(const Coord&) = default;
        Coord(Coord& other) = default;

        Coord(std::size_t nodeId, std::vector<double> coordinate)
            : nodeId(nodeId), coordinate(std::move(coordinate)){};

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
     * List of merges that should be applied to the mesh after initial
     * construction.
     *
     * Each entry in the vector corresponds to a single merge operation. Where
     * in each operation the coordinates identified by the indices on 'left' of
     * the map are merged with the coordinates identified by the indices on the
     * 'right' of the map. Each of these maps is fed to
     * MergePlan#computeMergePlan for the actual merging.
     *
     * Note that unlike identifying coordinates beforehand (using Coord#nodeId),
     * doing a merge via this route will prevent several problems in mesh
     * construction. For example take two triangles in the configuration:
     *
     *  1 -- 0 -- 2
     *    \  |  /
     *     \ | /
     *       3
     * Identifying coordinates 1 & 2 during construction[1] will merge edges:
     *  - 0-1 with 0-2 and
     *  - 1-3 with 2-3
     *  The second one clearly being unintended. By performing the merge {0->0,
     * 1->2} the coordinate 3 is not included. Hence, only the edges 0-1 and 0-2
     * will be collapsed to a single edge, leaving edges 1-3 and 2-3 untouched.
     *
     * [1] This configuration would happen when you have a periodic boundary as
     * result of a rotation. In this case a rotation with node 0 as the centre
     * of a 180 degree rotation.
     */
    virtual const std::vector<std::map<std::size_t, std::size_t>>& getMerges() {
        static std::vector<std::map<std::size_t, std::size_t>> dummy;
        return dummy;
    }

    /**
     * Get the dimension of the described mesh
     * @return The dimension
     */
    virtual std::size_t getDimension() const = 0;
};

}  // namespace Preprocessor

#endif  // HPGEM_MESHSOURCE_H

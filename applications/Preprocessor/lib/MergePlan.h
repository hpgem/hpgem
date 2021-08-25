/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2021, University of Twente
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
#ifndef HPGEM_MERGEPLAN_H
#define HPGEM_MERGEPLAN_H

#include "mesh/Mesh.h"
#include "utils/SparseUnionFind.h"

#include <array>
#include <map>
#include <memory>
#include <set>
#include <vector>

namespace Preprocessor {

/**
 * Description of merging several MeshEntities in a Mesh, for example to create
 * a periodic boundary.
 *
 * @tparam dimension The dimension of the Mesh.
 */
template <std::size_t dimension>
class MergePlan {
   private:
    /**
     * Helper structure for computeMergePlan(). This stores the elements that
     * use the left or right coordinates.
     */
    struct MergeElements {
        MergeElements(const Mesh<dimension>* mesh,
                      const std::map<CoordId, CoordId>& pairing)
            : mesh_(mesh) {
            leftElements = computeAdjacentElements(mesh, pairing, true);
            rightElements = computeAdjacentElements(mesh, pairing, false);
        }

        const Mesh<dimension>* mesh_;

        /**
         * Set of elements for which at least one node uses the left coordinate.
         */
        std::set<EntityGId> leftElements;
        /**
         * Set of elements for which at least one node uses the right
         * coordinate.
         */
        std::set<EntityGId> rightElements;

        /**
         * For each loop over all the MeshEntities that are adjacent to either
         * the left or right elements. Note that this includes both entities
         * that have coordinates in the left/right coordinate pairing, or none
         * at all.
         *
         * @tparam d The dimension of the MeshEntity to iterate over
         * @param left Whether to go over the left=true elements or the
         * right=false ones.
         * @param f The function that is invoked for each MeshEntity. The
         * signature is void(const Element<dimension>&, EntityLId) corresponding
         * to the i-th meshEntity of that dimension the element.
         */
        template <std::size_t d>
        void forEachBoundingEntity(
            bool left,
            std::function<void(const Element<dimension>&, EntityLId)> f) const {

            const std::set<EntityGId>& elements =
                left ? leftElements : rightElements;
            for (const EntityGId& gid : elements) {
                const Element<dimension>& element = mesh_->getElement(gid);
                std::size_t numEntities =
                    element.template getNumberOfIncidentEntities<d>();
                for (std::size_t i = 0; i < numEntities; ++i) {
                    f(element, i);
                }
            }
        }

       private:
        /**
         * Compute the elements that are adjacent to either the left or right
         * side of the coordinate pairing.
         *
         * @param mesh The mesh for the coordinate pairing and to take the
         * elements from.
         * @param pairing The pairing from left to right coordinates
         * @param left Whether to compute for the left (true) or right (false)
         * side of the pairing.
         * @return The set of elements adjacent to the choosen side of the
         * pairing.
         */
        std::set<EntityGId> computeAdjacentElements(
            const Mesh<dimension>* mesh,
            const std::map<CoordId, CoordId>& pairing, bool left);
    };
    /**
     * A pair of entities that will be merged.
     */
    using MergePair = std::pair<EntityGId, EntityGId>;

    /**
     * A list of entities that are going to be merged
     */
    using EntityMergeGroup = std::vector<EntityGId>;

   public:
    MergePlan() : mesh_(nullptr), merges_(){};

    /**
     * Execute this merge plan
     */
    void executeMerge() { executeMerge(tag<0>{}); }

    /**
     * Exposes the planed merges
     * @return The planned merges
     */
    const std::array<std::vector<EntityMergeGroup>, dimension>& getMerges()
        const {
        return merges_;
    }

    /**
     * \brief Create a merge plan based on pairing coordinates.
     *
     * Create a merge plan based on a pairing of coordinates. Each coordinate
     * pair specifies that the coordinate of the 'left' entry should be
     * connected to the coordinate of the 'right' entry (which may be the same).
     * For nodes this means that the associated nodes for the left and right
     * coordinate are merged.
     *
     * For higher dimensional mesh entities the merge is based on the adjacent
     * elements. Each MeshEntity is associated by some nodes (e.g. the endpoints
     * of a line or corners of a triangle). For each adjacent element these
     * associated nodes have coordinates. Two entities are merged if these
     * coordinates obey the coordinate mapping. Thus formally a pair of
     * MeshEntities eleft and eright is merged if:
     *  1. There exist a pair of elements Eleft and Eright adjacent to
     *     eleft/eright
     *  2. The set coordinates of the nodes of eleft/eright on these elements is
     *     CSleft and CSright
     *  3. The coordinate pairing maps CSleft to CSright
     * Note that this is only applied to MeshEntities < dimension, elements
     * themselves are not merged.
     *
     * @param mesh The mesh used for the pairing
     * @param pairing The pairing between left/right nodes
     * @return The resulting MergePlan
     */
    static MergePlan computeMergePlan(
        Mesh<dimension>* mesh, const std::map<CoordId, CoordId>& pairing);

   private:
    MergePlan(Mesh<dimension>* mesh,
              std::array<std::vector<EntityMergeGroup>, dimension> merges)
        : mesh_(mesh), merges_(merges){};

    /**
     * Recursive level to execute the merge at a single level and then move up.
     * @tparam d The level to merge
     */
    template <std::size_t d>
    void executeMerge(tag<d>);
    // Base case (working from 0 up)
    void executeMerge(tag<dimension>){};

    /**
     * The mesh for which this merge is computed
     */
    Mesh<dimension>* mesh_;

    /**
     * For each dimension of MeshEntities the a vector of groups of entities
     * that need to merged.
     */
    std::array<std::vector<EntityMergeGroup>, dimension> merges_;

    // For building the MergePlan

    /**
     * Recursive step in computeMergePlan, targeting a single MeshEntity
     * dimension before recursing to higher dimensions.
     *
     * @tparam d The meshEntity dimension to merge
     * @param mesh The mesh
     * @param pairing The pairing of coordinates
     * @param elements The left/right elements that use the coordinates of the
     * pairing.
     * @param merges The result of computing the merge plan
     */
    template <std::size_t d>
    static void computeMergePlanLevels(
        const Mesh<dimension>* mesh, const std::map<CoordId, CoordId>& pairing,
        const MergeElements& elements,
        std::array<std::vector<EntityMergeGroup>, dimension>& merges, tag<d>);
    // Base case
    static void computeMergePlanLevels(
        const Mesh<dimension>* mesh, const std::map<CoordId, CoordId>& pairing,
        const MergeElements& elements,
        std::array<std::vector<EntityMergeGroup>, dimension>& merges,
        tag<dimension>){};

    /**
     * Helper for computeMergePlanLevels finding all candidate left MeshEntities
     * that could be merged if a suitable matching right partner is found.
     *
     * A MeshEntity is eligible if all its nodes have coordinates in the left
     * side of the coordinate pairing[1]. The coordinates used for the
     * MeshEntity specify also the coordinates that a matching right partner
     * should have. To facilitate finding the right partner the set of
     * candidates is indexed by the CoordIds that a matching right partner
     * should have.
     *
     * [1] Technical note: The coordinates are not associated directly with the
     * MeshEntity but depend on the element. So the correcter statement would
     * be: If the MeshEntity is on the boundary of an Element E, such that the
     * nodes the of the MeshEntity have, on element E, coordinates from the left
     * set. This distinction is important, because it means that the same
     * MeshEntity could occur multiple times in the result. This can happen when
     * a boundary node of the MeshEntity has multiple coordinate ids that are in
     * the left set.
     *
     * @tparam d The dimension of entities to find
     * @param mesh The mesh
     * @param pairing The pairing of left -> right coordinates
     * @param elements The elements from the mesh with suitable entities
     * @return The id's of the MeshEntities on the left that could be merged.
     *   Indexed by the coordinate ids of the nodes as th
     */
    template <std::size_t d>
    static std::map<std::set<CoordId>, EntityGId> findLeftCandidates(
        const std::map<CoordId, CoordId>& pairing,
        const MergeElements& elements);

    /**
     * Create MergeGroups based on pairs of entities that need to be merged.
     *
     * There are several issues that can happen when the using just pairs of
     * entities:
     *  1. Duplicate entries
     *  2. Self merges e.g. (e1, e1)
     *  3. Transitive merges, e.g. (e1, e2) and (e3, e2)
     * Especially the latter requires that one does not use pairs but instead
     * uses lists of entities that need to be merged. This function transforms
     * the pairs of entities to merge into lists of entities to merge. Handling
     * the three cases as follows:
     *  1. Remove duplicates
     *  2. Removed completely
     *  3. Create one list with the elements in some arbitrary order
     *
     * @param mergePairs The pairs to merge
     * @return The lists of entities to merge
     */
    static std::vector<EntityMergeGroup> groupMerges(
        std::vector<MergePair> mergePairs);
};

template <std::size_t dimension>
std::vector<typename MergePlan<dimension>::EntityMergeGroup>
    MergePlan<dimension>::groupMerges(std::vector<MergePair> mergePairs) {
    // Approach:
    //   1. Use a UnionFind to merge the entity ids of the mergePairs giving
    //      the sets of entities that are merged.
    //   2. Output the sets/merge groups with more than 1 entity
    SparseUnionFind groups;
    for (auto const& pair : mergePairs) {
        groups.unionSets(pair.first.id, pair.second.id);
    }

    // Output a vector of merge
    std::vector<EntityMergeGroup> merges;
    // Each merge group in 'groups' is associated with the id of single
    // representative MeshEntity. This map takes that representative id and maps
    // it to the corresponding index of the EntityMergeGroup in 'groups'.
    // Using an index allows resizing merges without invalidating the
    // 'references'.
    std::map<EntityGId, std::size_t> mergeGroup;
    for (const std::size_t rawEntityId : groups) {
        EntityGId entityId = EntityGId(rawEntityId);
        EntityGId representative = EntityGId(groups.findSet(rawEntityId));

        auto iter = mergeGroup.find(representative);
        if (iter != mergeGroup.end()) {
            // Additional entity to add to the merge group
            merges[iter->second].push_back(entityId);
        } else {
            // First entity of the EntityMergeGroup
            std::size_t index = merges.size();
            // Add singleton element as start of the group
            merges.push_back({entityId});
            mergeGroup[representative] = index;
        }
    }
    return merges;
}

template <std::size_t dimension>
MergePlan<dimension> MergePlan<dimension>::computeMergePlan(
    Mesh<dimension>* mesh, const std::map<CoordId, CoordId>& pairing) {
    logger.assert_always(mesh != nullptr, "Null mesh");

    std::array<std::vector<EntityMergeGroup>, dimension> merges;

    if (pairing.empty()) {
        logger(INFO, "Empty pairing");
        return MergePlan<dimension>(mesh, merges);
    }

    // Zero dimension MeshEntities 'nodes' are special in that each coordinate
    // corresponds to a single node. Hence, merging coordinates is equivalent to
    // merging the associated nodes. Which is much less complicated than for
    // higher dimensional MeshEntities.
    {
        std::vector<MergePair> nodePairs;
        nodePairs.reserve(pairing.size());
        const auto& coordinates = mesh->getNodeCoordinates();
        for (auto const& pair : pairing) {
            nodePairs.emplace_back(coordinates[pair.first.id].nodeIndex,
                                   coordinates[pair.second.id].nodeIndex);
        }
        merges[0] = groupMerges(nodePairs);
    }

    // For merging higher dimensional MeshEntities we need to know the elements
    // that are adjacent to at least one of the coordinates for both the left
    // and right coordinates.
    MergeElements elements(mesh, pairing);
    logger(INFO, "Found %/% left/right elements adjacent to a merge coordinate",
           elements.leftElements.size(), elements.rightElements.size());

    // Recurse over the dimensions for the MeshEntities of higher dimension.
    computeMergePlanLevels(mesh, pairing, elements, merges, tag<1>{});

    // Print some debugging information
    logger(INFO, "Merge statistics:");
    for (std::size_t i = 0; i < dimension; ++i) {
        std::map<std::size_t, std::size_t> mergeSizes;
        for (const auto& merge : merges[i]) {
            mergeSizes[merge.size()]++;
        }
        logger(INFO, "\tmerges for entity dimension %", i);
        for (const auto& sizePair : mergeSizes) {
            logger(INFO, "\t\tsize %: %", sizePair.first, sizePair.second);
        }
    }
    return MergePlan<dimension>(mesh, merges);
}

template <std::size_t dimension>
template <std::size_t d>
void MergePlan<dimension>::computeMergePlanLevels(
    const Mesh<dimension>* mesh, const std::map<CoordId, CoordId>& pairing,
    const MergeElements& elements,
    std::array<std::vector<EntityMergeGroup>, dimension>& merges, tag<d>) {
    static_assert(d < dimension, "Can't merge elements or higher dimension");

    // Computing the merge plan for a single level done in two steps
    //  1. Find all eligible left-entities.
    //  2. List all possible right-entities, and check if there is a matching
    //     left-entity.

    // Eligible left-entities, indexed by the mapped coordinates (i.e.
    // coordinates in the right set).
    std::map<std::set<CoordId>, EntityGId> leftCandidates =
        findLeftCandidates<d>(pairing, elements);

    logger(DEBUG, "Found % left candidates of dimension %",
           leftCandidates.size(), d);
    // Phase 2: Find candidate right-entities and try to match them with a
    // left-entity.
    std::vector<MergePair> levelMerges;
    elements.template forEachBoundingEntity<d>(
        false, [&](const Element<dimension>& rightElement, EntityLId rid) {
            // Find the global coordinates for the mesh entity with respect to
            // the rightElement. Note that this may include CoordIds that are
            // not part of the pairing, but this is no problem, as then no match
            // will be found.
            std::vector<std::size_t> boundaryNodes =
                rightElement.getReferenceGeometry()
                    ->template getNodesOfEntity<d>(rid);
            std::set<CoordId> rightCoords;
            std::transform(
                boundaryNodes.begin(), boundaryNodes.end(),
                std::inserter(rightCoords, rightCoords.begin()),
                [&rightElement](const std::size_t& nodeNumber) -> CoordId {
                    return rightElement.getCoordinateIndex(nodeNumber);
                });

            auto candidateIter = leftCandidates.find(rightCoords);
            if (candidateIter != leftCandidates.end()) {
                // A matching left element is found.
                EntityGId leftEntity = candidateIter->second;
                EntityGId rightEntity =
                    rightElement.template getIncidentEntityIndex<d>(rid);
                levelMerges.emplace_back(leftEntity, rightEntity);
            }
        });
    logger(DEBUG, "Matched % right candidates of dimension %",
           levelMerges.size(), d);

    // Convert pairs of merges into merge vectors and store the result.
    merges[d] = groupMerges(levelMerges);

    // Recurse to higher dimensions
    computeMergePlanLevels(mesh, pairing, elements, merges, tag<d + 1>{});
}

template <std::size_t dimension>
template <std::size_t d>
std::map<std::set<CoordId>, EntityGId> MergePlan<dimension>::findLeftCandidates(
    const std::map<CoordId, CoordId>& pairing, const MergeElements& elements) {
    // Eligible left-entities, indexed by the mapped coordinates (i.e.
    // coordinates in the right set).
    std::map<std::set<CoordId>, EntityGId> leftCandidates;

    elements.template forEachBoundingEntity<d>(
        true, [&](const Element<dimension>& leftElement, EntityLId lid) {
            // This is the lid-th MeshEntity on the border of the leftElement.
            // This is a possible candidate, but only if the coordinates for all
            // the nodes are in the left set of the pairing.

            // Listing of which nodes of the Element the entity uses, e.g.
            // the edge uses nodes 0 and 2 of the triangle (nodes 0-2).
            std::vector<std::size_t> boundaryNodes =
                leftElement.getReferenceGeometry()
                    ->template getNodesOfEntity<d>(lid);

            // For each node that is used, find the corresponding coordinate for
            // the element and try and translate it. If that succeeds for all
            // nodes then we have a possible left-candidate-MeshEntity.
            std::set<CoordId> rightCoords;
            for (const std::size_t& boundaryNode : boundaryNodes) {
                CoordId leftCoord =
                    leftElement.getCoordinateIndex(boundaryNode);
                auto iter = pairing.find(leftCoord);
                if (iter == pairing.end()) {
                    // Can't translate the coordinate for this node -> Not a
                    // candidate
                    return;
                } else {
                    rightCoords.insert(iter->second);
                }
            }
            if (boundaryNodes.size() != rightCoords.size()) {
                // Backup check
                logger.assert_always(
                    false,
                    "Multiple nodes of the MeshEntity are merged onto the same "
                    "coordinate. This is not supported by hpgem.");
                return;
            }
            leftCandidates[rightCoords] =
                leftElement.template getIncidentEntityIndex<d>(lid);
        });
    return leftCandidates;
}

template <std::size_t dimension>
template <std::size_t d>
void MergePlan<dimension>::executeMerge(tag<d>) {
    const std::vector<EntityMergeGroup>& levelMerges = merges_[d];

    for (const EntityMergeGroup& mergeGroup : levelMerges) {
        if (mergeGroup.size() <= 1) {
            // Nothing to do
            continue;
        }
        // Use the first entity as the entity that remains
        MeshEntity<d, dimension>& target =
            mesh_->template getEntity<d>(mergeGroup[0]);
        for (std::size_t i = 1; i < mergeGroup.size(); ++i) {
            target.merge(mergeGroup[i]);
        }
    }
    executeMerge(tag<d + 1>{});
}

template <std::size_t dimension>
std::set<EntityGId>
    MergePlan<dimension>::MergeElements::computeAdjacentElements(
        const Mesh<dimension>* mesh, const std::map<CoordId, CoordId>& pairing,
        bool left) {

    // Approach:
    // Go over all coordinates:
    //   - Find the corresponding node
    //   - Of this node take all adjacent Elements
    //     - Check if the element uses the given coordinate index for the node

    std::set<EntityGId> result;
    for (const auto& pair : pairing) {
        CoordId coordId = left ? pair.first : pair.second;
        const MeshEntity<0, dimension>& node =
            mesh->getNode(mesh->getNodeCoordinates()[coordId.id].nodeIndex);
        for (std::size_t i = 0; i < node.getNumberOfElements(); ++i) {
            const Element<dimension>& element = node.getElement(i);
            EntityLId localId = node.getLocalIndex(i);
            CoordId elementCoord = element.getCoordinateIndex(localId);
            if (elementCoord == coordId) {
                result.insert(element.getGlobalIndex());
            }
        }
    }
    return result;
}

}  // namespace Preprocessor

#endif  // HPGEM_MERGEPLAN_H

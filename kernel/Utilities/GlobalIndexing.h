/*
This file forms part of hpGEM. This package has been developed over a number of
years by various people at the University of Twente and a full list of
contributors can be found at http://hpgem.org/about-the-code/team

This code is distributed using BSD 3-Clause License. A copy of which can found
below.


Copyright (c) 2018, Univesity of Twenete
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef HPGEM_KERNEL_GLOBALINDEXING_H
#define HPGEM_KERNEL_GLOBALINDEXING_H

#include "Base/Edge.h"
#include "Base/Element.h"
#include "Base/Face.h"
#include "Base/MeshManipulatorBase.h"
#include "Base/Node.h"

#include "Base/MeshEntity.h"

namespace hpgem {
namespace Utilities {

/// \brief Mapping between basis functions on geometrical objects (Element,
/// Face, etc.) to the global indices.
///
/// This is the mapping between the local indices on geometrical objects used
/// for basis functions and the globally unique indices. The latter are usually
/// used for constructing the global matrices and vectors. Optionally only some
/// of the unknowns can be included in the index.
///
/// This class defines an equivalence between three ways of addressing the same
/// basis function:
///
/// - Using the triplet (O, U, I), which refers to the I-th basis function of
///   unknown U, which is *locally* owned by geometrical object O.
/// - Global index I, a unique number from the range [0, N), where N is the
///   number of basis functions for all unknowns and processors combined.
/// - Local index Il, a unique number from the range [0, Nl), where Nl is the
///   number of basis functions for the geometrical objects owned by this
///   processor. This is only available for basis functions owned by this
///   processor. These are such that they follow the same order as the global
///   index I.
///
/// As the layout of the basis functions on the same geometrical object is
/// sequential, we only need to store the offset for the global index of the
/// 0-th basis function of each unknown. This has several benefits:
///
///  - Only one call per geometrical object and per unknown is needed to the
///    GlobalIndexing instead of one per local basis function. Thus increasing
///    the efficiency of the code using the indexing.
///  - It is valid to ask for the global index of an unknown for a object
///    without any basis function. That index will not be used, but it might
///    save on performance. (The result however is an arbitrary value)
///
/// The current implementation supports three ways of assigning the global
/// indices to the basis functions. The general layout for both looks like:
///
/// - basis functions for element e1 (sequentially)
/// - for each face f associated with e1
///   - basis functions for f (sequentially)
/// - for each edge e associated with e1
///   - basis functions for e (sequentially)
/// - for each node n associated with e1
///   - basis functions for n (sequentially)
/// [next element e2, etc.]
///
/// Where we associated the face with the left element and the edges and nodes
/// with the first element in their element listing.
///
/// The difference between the three layouts is in where the basis functions for
/// different unknowns and processors are placed. This is best illustrated by an
/// example, for when there are 2 processors (0,1) and 2 unknowns (U, P).
///
/// SEQUENTIAL UP0 (interleaved), UP1 (interleaved)
/// BLOCKED_PROCESSOR U0, P0, U1, P1
/// BLOCKED_GLOBAL U0, U1, P0, P1
///
/// Thus with the SEQUENTIAL layout all the basis functions for the same
/// geometrical object are consecutively numbered, i.e. the basis functions of
/// unknown 1 for an element directly follow the basis functions of unknown 0 of
/// the same element. Both blocked layouts do the opposite and sequentially
/// number all the basis functions of each unknown, thus first all basis
/// functions for unknown 0, then followed by those of unknown 1, etc. The
/// BLOCKED_GLOBAL applies this method to all basis functions on all the
/// processors, thereby creating indices for block matrices and vectors. The
/// BLOCKED_PROCESSOR only applies this process to the basis functions on a
/// single processor. The basis functions for different processors are then laid
/// out sequentially, thus all basis functions of processor 0 before those of
/// processor 1.
///
class GlobalIndexing {
   private:
    struct Offsets;
    struct ConstOffsetVisitor;
    struct OffsetVisitor;

   public:
    enum Layout { SEQUENTIAL, BLOCKED_PROCESSOR, BLOCKED_GLOBAL };

    GlobalIndexing();
    explicit GlobalIndexing(Base::MeshManipulatorBase* mesh,
                            Layout layout = BLOCKED_PROCESSOR,
                            const std::vector<std::size_t>* unknowns = nullptr);

    /// \brief Recreate the mapping for a new mesh.
    ///
    /// \param mesh The new  mesh
    /// \param blocked Whether to construct a blocked mapping.
    /// \param unknowns The subset of unknowns to include (nullptr = all)
    void reset(Base::MeshManipulatorBase* mesh, Layout layout,
               const std::vector<std::size_t>* unknowns = nullptr);

    /// \brief Lookup the global index of the 0-th basis function local to an
    /// MeshEntity.
    ///
    /// \param entity The entity that owns the basis functiosn
    /// \param unknown The unknown for which this basis function is.
    /// \return The global index of the basis function.
    int getGlobalIndex(const Base::MeshEntity* entity,
                       std::size_t unknown) const {
        logger.assert_debug(unknown < totalNumberOfUnknowns_,
                            "No such unknown %", unknown);
        const Offsets& offset = offsets_[unknown];
        ConstOffsetVisitor visitor(offset);
        entity->accept(visitor);
        return visitor.entityOffset_;
    }

    /// \brief Lookup the local index of a basis function local to a MeshEntity
    /// that is owned by this processor.
    ///
    /// \param entity The MeshEntity for which to look it up, this processor
    /// should own the basis functions for the element.
    /// \param unknown The unknown for which to look up the index.
    /// \return The local index of the 0-th local basis function on the
    /// given element for the given unknown.
    int getProcessorLocalIndex(const Base::MeshEntity* entity,
                               std::size_t unknown) const {
        logger.assert_debug(unknown < totalNumberOfUnknowns_,
                            "No such unknown %", unknown);
        const Offsets& offset = offsets_[unknown];
        ConstOffsetVisitor visitor(offset);
        entity->accept(visitor);
        return visitor.entityOffset_ - offset.blockStart_ + offset.localOffset_;
    }

    /// Convert a global index of a basis function into a local index
    ///
    /// \param globalIndex The global index to convert
    /// \return The local index or -1 if it is not locally owned.
    int globalToProcessorLocalIndex(int globalIndex) const {
        for (std::size_t unknown = 0; unknown < totalNumberOfUnknowns_;
             ++unknown) {
            const Offsets& offset = offsets_[unknown];
            if (offset.owns(globalIndex)) {
                return globalIndex - offset.blockStart_ + offset.localOffset_;
            }
        }
        return -1;
    }

    std::size_t getNumberOfLocalBasisFunctions() const {
        return localNumberOfBasisFunctions_;
    }

    /// The total number of unknowns in the mesh. For partial indices this
    /// includes those unknowns that are not included.
    /// \return The total number of unknowns.
    std::size_t getTotalNumberOfUnknowns() const {
        return totalNumberOfUnknowns_;
    }

    std::size_t getNumberOfIncludedUnknowns() const {
        return includedUnknowns_.size();
    }

    // Note, non const references as we do not own the Mesh, we just hold a
    // reference to it.
    /// \brief The mesh that for which this GlobalIndex is build
    ///
    /// \return The mesh, or null if none.
    Base::MeshManipulatorBase* getMesh() const { return mesh_; }

    /// \brief Retrieve the global indices associated with the basis
    /// functions that are associated with element.
    ///
    /// This maps the element local indexing of basis functions to the
    /// global indices of the basis functions. The indices are in the
    /// standard order for the element:
    /// - basis functions local to the element for unknown 0
    /// - for each face f
    ///   - basis functions local to f for unknown 0
    /// - for each edge e
    ///   - basis functions local to e for unknown 0
    /// - for each node n (dim > 1)
    ///   - basis functions local to n for unknown 0
    /// [repeated for unknown 1,2,...]
    ///
    /// \param element The element to construct the mapping for.
    /// \param indices The mapping.
    void getGlobalIndices(const Base::Element* element,
                          std::vector<int>& indices) const {
        std::size_t size = getGlobalIndices(element, 0, indices);
        indices.resize(size);
    }

    /// \brief Retrieve the global indices associated with the basis
    /// functions that are associated with element.
    ///
    /// \see getGlobalIndices(const Base::Element, std::vector<int>)
    std::vector<int> getGlobalIndices(const Base::Element* element) const {
        std::vector<int> result;
        getGlobalIndices(element, result);
        return result;
    }

    /// \brief Retrieve the global indices associated with the basis
    /// functions that are associated with the elements sharing a face.
    ///
    /// This mapping is the mapping of the left element flowed by the
    /// mapping of the right element (if it exists). Note that this is the
    /// same order as used in the FaceMatrix.
    /// \param face The face to construct the mapping for.
    /// \param indices The mapping.
    void getGlobalIndices(const Base::Face* face,
                          std::vector<int>& indices) const;

    /// Verify that the index is complete (i.e. all global indices are
    /// available) Only works when assert_debug is available
    void verifyCompleteIndex() const;

    /// Check whether a certain unknown is included in the index.
    /// \param unknown The unknown to check
    /// \return Whether it is included
    bool isIncludedInIndex(std::size_t unknown) const {
        logger.assert_debug(unknown < totalNumberOfUnknowns_,
                            "No such unknown %", unknown);
        return offsets_[unknown].includedInIndex_;
    }

    /// A vector of the unknowns that are included in the index.
    const std::vector<std::size_t>& getIncludedUnknowns() const {
        return includedUnknowns_;
    }

   private:
    /// Construct the indexing with a blocked layout
    /// \param global Whether to use BlockedGlobal layout (true) or
    /// BlockedProcessor (false)
    void constructBlocked(bool global);
    void constructUnblocked();

    /// \brief Communicate transfer the global indices for locally owned
    /// basis functions with the push/pull partners.
    ///
    /// This transfers the knowledge of the global indices for the locally
    /// owned basis functions to the neighbour processors, so that they also
    /// know the global indices for the basis functions on their shadow
    /// elements (element owned by the neighbouring processor).
    void communicatePushPullElements();

#ifdef HPGEM_USE_MPI
    // Helper methods for communicatePushPullElements

    /// \brief Write message for the first round
    ///
    /// \param elements The elements from the push list
    /// \param message The vector to write the message in
    /// \param targetProcessor The processor to which the message will be
    ///        send
    /// \param secondRoundTags The tags for the elements to send in the
    ///        second round.
    void createInitialMessage(const std::vector<Base::Element*>& elements,
                              std::vector<std::size_t>& message,
                              std::size_t targetProcessor,
                              std::set<std::size_t>& secondRoundTags) const;
    /// \brief Write the message for the second round
    /// \param tags The tags to send
    /// \param message The vector to write the message to.
    void createSecondMessage(const std::set<std::size_t>& tags,
                             std::vector<std::size_t>& message) const;
    void elementMessage(std::size_t elementId,
                        std::vector<std::size_t>& message) const;
    void faceMessage(std::size_t faceId,
                     std::vector<std::size_t>& message) const;
    void edgeMessage(std::size_t edgeId,
                     std::vector<std::size_t>& message) const;
    void nodeMessage(std::size_t nodeId,
                     std::vector<std::size_t>& message) const;

    /// Process a message
    ///
    /// \param message Storage with the message
    /// \param count The length of the message (in element s of message)
    void processMessage(const std::vector<std::size_t>& message,
                        std::size_t count);

    /// Use MPI_Probe and MPI_Recv to receive a variable size message
    /// \param receiveMessage The storage for receiving the message, may
    ///        grow to accommodate the message. Will not be shrunk to
    ///        prevent unnecessary (de)allocations
    /// \return The length of the message.
    std::size_t probeAndReceive(std::vector<std::size_t>& receiveMessage) const;
#endif

    /// \brief Same as getGlobalIndices(const Base::Element *, std::vector<int>)
    /// but offsetting the place in the indices vector used for the output.
    /// \param element The element to construct the mapping for
    /// \param offset The index where to place the first global id in indices.
    /// \param indices The mapping.
    /// \return The number of global indices in the `indices` vector that are
    /// used.
    std::size_t getGlobalIndices(const Base::Element* element,
                                 std::size_t offset,
                                 std::vector<int>& indices) const;

    /// \brief Helper structure to store the information for a single unknown.
    ///
    /// This helper structure stores two pieces of information for included
    /// unknowns:
    ///
    ///  - The offsets for the first basis function owned by geometrical objects
    ///    (as indexed by their ID). Thus given the id I of a element we can
    ///    lookup the global index of the first local basis function on that
    ///    element in elementOffsets_[I] and analogously for the other
    ///    geometrical objects.
    ///  - Information about the block of indices used for this unknown on this
    ///    processor. This can be used to translate between local and global
    ///    indices and to determine whether a certain basis function is owned by
    ///    this processor.
    ///
    /// The block information consists of three numbers (G, L, N). G, the
    /// blockStart, is the global index of the first basis function in the
    /// block, L, the localOffset, the processor local index of the first basis
    /// function in the block and N, the numberOfBasisFunctionsInBlock, is the
    /// number of basis functions that are part of the block.
    ///
    /// For an example how these numbers we look at an example with 2 processors
    /// (0,1) and 2 variables (U, P) with 10 and 5, respectively, degrees of
    /// freedom on each processor.
    ////
    /// - SEQUENTIAL, the indices of U0 and P0 are interleaved and thus part of
    ///   the same block. Hence we have
    ///   (G0, L0, N0) = (0,  0, 15) (for both U, P)
    ///   (G1, L1, N1) = (15, 0, 15) (for both U, P)
    /// - BLOCKED_PROCESSOR, the order of indices is U0, P0, U1, P1, thus the U
    ///   and P variables can be in separate blocks
    ///   (G0, L0, N0) = ( 0,  0, 10) for U0,
    ///   (G0, L0, N0) = (10, 10,  5) for P0,
    ///   (G1, L1, N1) = (15,  0, 10) for U1,
    ///   (G1, L1, N1) = (25, 10,  5) for P1,
    /// - BLOCKED_GLOBAL the order of the indices is U0, U1, P0, P1, thus the U
    ///   and P variables must be in different blocks. The ordering is thus
    ///   (G0, L0, N0) = ( 0,  0, 10) (U0)
    ///   (G1, L1, N1) = (10,  0, 10) (U1)
    ///   (G0, L0, N0) = (20, 10,  5) (P0)
    ///   (G0, L0, N0) = (25, 10,  5) (P1)
    ///
    /// Note that in all of these three cases the blockStart G increases by the
    /// number of basis functions of the previous block (N). Similarly when
    /// looking at the localOffset L for a single processor it will also
    /// increase by the number of basis functions in the previous block N of the
    /// previous block of the same processor.
    struct Offsets {
        /// The global index G, the start of the block with indices.
        int blockStart_ = 0;
        /// The local offset L
        int localOffset_ = 0;
        /// Number of basis functions N in this block
        int numberOfBasisFunctionsInBlock_ = 0;
        /// Whether or not this unknown is included in the GlobalIndexing
        bool includedInIndex_ = false;
        /// Mapping between the element id and the global id of the first basis
        /// function, note these contain both offsets for owned geometrical
        /// objects and for objects on the boundary.
        std::map<std::size_t, int> elementOffsets_;
        std::map<std::size_t, int> faceOffsets_;
        std::map<std::size_t, int> edgeOffsets_;
        std::map<std::size_t, int> nodeOffsets_;
        /// Set the block parameters and correct the offsets
        /// \param globalOffset The offset for the block start
        /// \param localOffset The local offset
        /// \param numberOfBasisFunctions The number of basis functions in the
        /// block
        void setOffset(int globalOffset, int localOffset,
                       int numberOfBasisFunctions);
        /// Check if the global index falls within the block described by this
        /// Offset. Note that there may be multiple blocks in which the same
        /// global index falls when using SEQUENTIAL mode. If the Offset is not
        /// included in the index then the result is false.
        bool owns(int globalIndex) const {
            return includedInIndex_ && blockStart_ <= globalIndex &&
                   globalIndex < blockStart_ + numberOfBasisFunctionsInBlock_;
        }
    };

    /**
     * Visitor to access the *Offset_ members (e.g. elementOffsets_) with the
     * possibility of editing the offset. If there is no offset entry
     * available for the MeshEntity then it will be created (and initialized to
     * the default 0).
     */
    struct OffsetVisitor : public Base::ConstMeshEntityVisitor {
        OffsetVisitor(Offsets& offsets)
            : visitOffsets_(offsets), entityOffset_(nullptr){};

        void visit(const Base::Element& element) final {
            entityOffset_ = &visitOffsets_.elementOffsets_[element.getID()];
        }
        void visit(const Base::Face& face) final {
            entityOffset_ = &visitOffsets_.faceOffsets_[face.getID()];
        }
        void visit(const Base::Edge& edge) final {
            entityOffset_ = &visitOffsets_.edgeOffsets_[edge.getID()];
        }
        void visit(const Base::Node& node) final {
            entityOffset_ = &visitOffsets_.nodeOffsets_[node.getID()];
        }

        Offsets& visitOffsets_;
        /// Pointer to the entry, possibly newly initialized
        int* entityOffset_;
    };

    /**
     * Visitor to access the *Offsets_ members (e.g. elementOffsets_) of the
     * Offsets class. It is undefined behaviour to accept a MeshEntity for which
     * there is no offset entry available.
     */
    struct ConstOffsetVisitor : public Base::ConstMeshEntityVisitor {
        ConstOffsetVisitor(const Offsets& offsets) : visitOffsets_(offsets){};
        void visit(const Base::Element& element) final {
            const auto basisStart =
                visitOffsets_.elementOffsets_.find(element.getID());
            logger.assert_debug(
                basisStart != visitOffsets_.elementOffsets_.end(),
                "No offset for Element %", element.getID());
            entityOffset_ = basisStart->second;
        }

        void visit(const Base::Face& face) final {
            const auto basisStart =
                visitOffsets_.faceOffsets_.find(face.getID());
            logger.assert_debug(basisStart != visitOffsets_.faceOffsets_.end(),
                                "No offset for Face %", face.getID());
            entityOffset_ = basisStart->second;
        }

        void visit(const Base::Edge& edge) final {
            const auto basisStart =
                visitOffsets_.edgeOffsets_.find(edge.getID());
            logger.assert_debug(basisStart != visitOffsets_.edgeOffsets_.end(),
                                "No offset for Edge %", edge.getID());
            entityOffset_ = basisStart->second;
        }

        void visit(const Base::Node& node) final {
            const auto basisStart =
                visitOffsets_.nodeOffsets_.find(node.getID());
            logger.assert_debug(basisStart != visitOffsets_.nodeOffsets_.end(),
                                "No offset for Node %", node.getID());
            entityOffset_ = basisStart->second;
        }

        const Offsets& visitOffsets_;
        int entityOffset_;
    };

    /// Offsets for each of the unknowns
    std::vector<Offsets> offsets_;
    /// The total number of unknowns in the mesh
    std::size_t totalNumberOfUnknowns_;
    /// Total number of local basis functions over all the unknowns.
    std::size_t localNumberOfBasisFunctions_;
    /// Total number of basis functions over all unknowns and processors;
    std::size_t globalNumberOfBasisFunctions_;
    /// Ordered list of the unknowns in use
    std::vector<std::size_t> includedUnknowns_;
    /// The mesh, non const because of the push/pull elements in the mesh.
    Base::MeshManipulatorBase* mesh_;
};
}  // namespace Utilities

}  // namespace hpgem

#endif  // HPGEM_KERNEL_GLOBALINDEXING_H

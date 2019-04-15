/*
This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
http://hpgem.org/about-the-code/team

This code is distributed using BSD 3-Clause License. A copy of which can found below.


Copyright (c) 2018, Univesity of Twenete
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef HPGEM_GLOBALINDEXING_H
#define HPGEM_GLOBALINDEXING_H

#include "../Base/Edge.h"
#include "../Base/Element.h"
#include "../Base/Face.h"
#include "../Base/MeshManipulatorBase.h"
#include "../Base/Node.h"

namespace Utilities
{

    /// \brief Mapping between basis functions on geometrical objects (Element,
    /// Face, etc.) to the global indices.
    ///
    /// This is the mapping between the local indices on geometrical objects used
    /// for basis functions and the globally unique indices. The latter are usually
    /// used for constructing the global matrices and vectors.
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
    /// BLOCKED_PROCESSOR only applies this proces to the basis functions on a single
    /// processor. The basis functions for different processors are then laid out
    /// sequentially, thus all basis functions of processor 0 before those of
    /// processor 1.
    ///
    class GlobalIndexing
    {
    public:
        enum Layout
        {
            SEQUENTIAL,
            BLOCKED_PROCESSOR,
            BLOCKED_GLOBAL
        };


        GlobalIndexing();
        GlobalIndexing(Base::MeshManipulatorBase* mesh, Layout layout);

        /// \brief Recreate the mapping for a new mesh.
        ///
        /// \param mesh The new  mesh
        /// \param blocked Whether to construct a blocked mapping.
        void reset(Base::MeshManipulatorBase* mesh, Layout layout);

        /// \brief Lookup the global index of the 0-th basis function local to an element.
        /// \param element The element of the basis function
        /// \param unknown The unknown for which this basis function is.
        /// \return The global index of the basis function.
        int getGlobalIndex(const Base::Element* element, std::size_t unknown) const
        {
            logger.assert_debug(unknown < numberOfUnknowns_, "No such unknown %", unknown);
            // TODO: Remove the basis function argument from all functions, as we lay them out sequentially.
            const auto basisStart = offsets[unknown].elementOffsets_.find(element->getID());
            logger.assert_debug(basisStart != offsets[unknown].elementOffsets_.end(),
                                "No indices known for element %", element->getID());
            return basisStart->second;
        }
        int getGlobalIndex(const Base::Face* face, std::size_t unknown) const
        {
            logger.assert_debug(unknown < numberOfUnknowns_, "No such unknown %", unknown);
            const auto basisStart = offsets[unknown].faceOffsets_.find(face->getID());
            logger.assert_debug(basisStart != offsets[unknown].faceOffsets_.end(),
                                "No indices known for face %", face->getID());
            return basisStart->second;
        }
        int getGlobalIndex(const Base::Edge* edge, std::size_t unknown) const
        {
            logger.assert_debug(unknown < numberOfUnknowns_, "No such unknown %", unknown);
            const auto basisStart = offsets[unknown].edgeOffsets_.find(edge->getID());
            logger.assert_debug(basisStart != offsets[unknown].edgeOffsets_.end(),
                                "No indices known for edge %", edge->getID());
            return basisStart->second;
        }
        int getGlobalIndex(const Base::Node* node, std::size_t unknown) const
        {
            logger.assert_debug(unknown < numberOfUnknowns_, "No such unknown %", unknown);
            const auto basisStart = offsets[unknown].nodeOffsets_.find(node->getID());
            logger.assert_debug(basisStart != offsets[unknown].nodeOffsets_.end(),
                                "No indices known for node %", node->getID());
            return basisStart->second;
        }

        /// \brief Lookup the local index of a basis function local to an element
        /// that is owned by this processor.
        ///
        /// \param element The element for which to look it up, this processor
        /// should own the basis functions for the element.
        /// \param unknown The unknown for which to look up the index.
        /// \return The local index of the 0-th local basis function on the
        /// given element for the given unknown.
        int getProcessorLocalIndex(const Base::Element *element, std::size_t unknown) const
        {
            logger.assert_debug(unknown < numberOfUnknowns_, "No such unknown %", unknown);
            const Offsets& offset = offsets[unknown];
            const auto basisStart = offset.elementOffsets_.find(element->getID());
            logger.assert_debug(basisStart != offset.elementOffsets_.end(),
                                "No indices known for element %", element->getID());
            return basisStart->second - offset.blockStart_ + offset.localOffset_;
        }

        int getProcessorLocalIndex(const Base::Face *face, std::size_t unknown) const
        {
            logger.assert_debug(unknown < numberOfUnknowns_, "No such unknown %", unknown);
            const Offsets& offset = offsets[unknown];
            const auto basisStart = offset.faceOffsets_.find(face->getID());
            logger.assert_debug(basisStart != offset.faceOffsets_.end(),
                                "No indices known for face %", face->getID());
            return basisStart->second - offset.blockStart_ + offset.localOffset_;
        }
        int getProcessorLocalIndex(const Base::Edge *edge, std::size_t unknown) const
        {
            logger.assert_debug(unknown < numberOfUnknowns_, "No such unknown %", unknown);
            const Offsets& offset = offsets[unknown];
            const auto basisStart = offset.edgeOffsets_.find(edge->getID());
            logger.assert_debug(basisStart != offset.edgeOffsets_.end(),
                                "No indices known for edge %", edge->getID());
            return basisStart->second - offset.blockStart_ + offset.localOffset_;
        }
        int getProcessorLocalIndex(const Base::Node *node, std::size_t unknown) const
        {
            logger.assert_debug(unknown < numberOfUnknowns_, "No such unknown %", unknown);
            const Offsets& offset = offsets[unknown];
            const auto basisStart = offset.nodeOffsets_.find(node->getID());
            logger.assert_debug(basisStart != offset.nodeOffsets_.end(),
                                "No indices known for node %", node->getID());
            return basisStart->second - offset.blockStart_ + offset.localOffset_;
        }

        /// Convert a global index of a basis function into a local index
        ///
        /// \param globalIndex The global index to convert
        /// \return The local index or -1 if it is not locally owned.
        int globalToProcessorLocalIndex(int globalIndex) const
        {
            for (std::size_t unknown = 0; unknown < numberOfUnknowns_; ++unknown)
            {
                const Offsets& offset = offsets[unknown];
                if (offset.owns(globalIndex))
                {
                    return globalIndex - offset.blockStart_ + offset.localOffset_;
                }
            }
            return -1;
        }

        bool isLocallyOwned(const Base::Element* element, std::size_t unknown) const
        {
            logger.assert_debug(unknown < numberOfUnknowns_, "No such unknown %", unknown);
            const Offsets& offset = offsets[unknown];
            const auto basisStart = offset.elementOffsets_.find(element->getID());
            logger.assert_debug(basisStart != offset.elementOffsets_.end(),
                                "No indices available for element %", element->getID());
            int globalId = basisStart->second;
            return offset.owns(globalId);
        }

        bool isLocallyOwned(const Base::Face* face, std::size_t unknown) const
        {
            //TODO: The implementation of these functions mimics what GlobalMatrix
            // did. However, this seems a crazy waste of energy as we should be able
            // to determine if a face/edge/node is locally owned or not without all
            // the extra information about the unknown.
            logger.assert_debug(unknown < numberOfUnknowns_, "No such unknown %", unknown);
            const Offsets& offset = offsets[unknown];
            const auto basisStart = offset.faceOffsets_.find(face->getID());
            if (basisStart == offset.faceOffsets_.end())
                return false;
            int globalId = basisStart->second;
            return offset.owns(globalId);
        }

        bool isLocallyOwned(const Base::Edge* edge, std::size_t unknown) const
        {
            logger.assert_debug(unknown < numberOfUnknowns_, "No such unknown %", unknown);
            const Offsets& offset = offsets[unknown];
            const auto basisStart = offset.edgeOffsets_.find(edge->getID());
            if(basisStart == offset.edgeOffsets_.end())
                return false;
            int globalId = basisStart->second;
            return offset.owns(globalId);
        }

        bool isLocallyOwned(const Base::Node* node, std::size_t unknown) const
        {
            logger.assert_debug(unknown < numberOfUnknowns_, "No such unknown %", unknown);
            const Offsets& offset = offsets[unknown];
            const auto basisStart = offset.nodeOffsets_.find(node->getID());
            if(basisStart == offset.nodeOffsets_.end())
                return false;
            int globalId = basisStart->second;
            return offset.owns(globalId);
        }

        std::size_t getNumberOfLocalBasisFunctions() const
        {
            return localBasisFunctions_;
        }

        std::size_t getNumberOfUnknowns() const
        {
            return numberOfUnknowns_;
        }

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
        void getGlobalIndices(const Base::Element *element, std::vector<int>& indices)
        {
            getGlobalIndices(element, 0, indices);
        }

        /// \brief Retrieve the global indices associated with the basis
        /// functions that are associated with element.
        ///
        /// \see getGlobalIndices(const Base::Element, std::vector<int>)
        std::vector<int> getGlobalIndices(const Base::Element *element)
        {
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
        void getGlobalIndices(const Base::Face *face, std::vector<int>& indices);

    private:

        /// Construct the indexing with a blocked layout
        /// \param mesh The mesh to base the index on
        /// \param global Whether to use BlockedGlobal layout (true) or
        /// BlockedProcessor (false)
        void constructBlocked(Base::MeshManipulatorBase& mesh, bool global);
        void constructUnblocked(Base::MeshManipulatorBase& mesh);

        /// \brief Communicate transfer the global indices for locally owned
        /// basis functions with the push/pull partners.
        ///
        /// This transfers the knowledge of the global indices for the locally
        /// owned basis functions to the neighbour processors, so that they also
        /// know the global indices for the basis functions on their shadow
        /// elements (element owned by the neighbouring processor).
        /// \param mesh The mesh, needed for the push/pull partners.
        void communicatePushPullElements(Base::MeshManipulatorBase& mesh);

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
        void createInitialMessage(
                const std::vector<Base::Element*>& elements,
                std::vector<std::size_t>& message,
                std::size_t targetProcessor,
                std::set<std::size_t>& secondRoundTags)  const;
        /// \brief Write the message for the second round
        /// \param tags The tags to send
        /// \param message The vector to write the message to.
        void createSecondMessage(const std::set<std::size_t>& tags, std::vector<std::size_t>& message)  const;
        void elementMessage(std::size_t elementId, std::vector<std::size_t>& message) const;
        void faceMessage(std::size_t faceId, std::vector<std::size_t>& message) const;
        void edgeMessage(std::size_t edgeId, std::vector<std::size_t>& message) const;
        void nodeMessage(std::size_t nodeId, std::vector<std::size_t>& message) const;

        /// Process a message
        ///
        /// \param message Storage with the message
        /// \param count The length of the message (in element s of message)
        void processMessage(const std::vector<std::size_t>& message, std::size_t count);

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
        void getGlobalIndices(const Base::Element *element, std::size_t offset, std::vector<int> &indices);

        /// \brief Helper structure to store the information for a single unknown.
        ///
        /// This helper structure stores two pieces of information:
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
        struct Offsets
        {
            /// The global index G, the start of the block with indices.
            int blockStart_;
            /// The local offset L
            int localOffset_;
            /// Number of basis functions N in this block
            int numberOfBasisFunctionsInBlock_;
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
            /// \param numberOfBasisFunctions The number of basis functions in the block
            void setOffset(int globalOffset, int localOffset, int numberOfBasisFunctions);
            bool owns(int globalIndex) const
            {
                return blockStart_ <= globalIndex && globalIndex < blockStart_ + numberOfBasisFunctionsInBlock_;
            }
        };
        /// Offsets for each of the unknowns
        std::vector<Offsets> offsets;
        std::size_t numberOfUnknowns_;
        /// Total number of local basis functions over all the unknowns.
        std::size_t localBasisFunctions_;
        /// Dimension of the mesh.
        std::size_t meshDimension;
    };
}

#endif //HPGEM_GLOBALINDEXING_H

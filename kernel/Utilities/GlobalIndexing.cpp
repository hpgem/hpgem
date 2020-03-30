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

#include "GlobalIndexing.h"

#include "Base/MpiContainer.h"

namespace Utilities
{

    GlobalIndexing::GlobalIndexing()
            : mesh_ (nullptr)
            , localBasisFunctions_ (0)
            , offsets_ (0)
            , numberOfUnknowns_ (0)
    {}

    GlobalIndexing::GlobalIndexing(Base::MeshManipulatorBase *mesh, Layout layout, const std::vector<std::size_t>* unknowns)
        : GlobalIndexing()
    {
        // Actual initialization
        reset(mesh, layout, unknowns);
    }

    std::size_t GlobalIndexing::getGlobalIndices(const Base::Element *element, std::size_t offset, std::vector<int> &indices) const
    {
        logger.assert_debug(element != nullptr, "Null pointer as element");
        const std::size_t numberOfUnknowns = element->getNumberOfUnknowns();
        // Make sure we have enough space
        std::size_t totalBasisFunctions = 0;
        for (std::size_t i = 0; i < numberOfUnknowns; ++i)
        {
            if (offsets_[i].includedInIndex_)
            {
                totalBasisFunctions += element->getNumberOfBasisFunctions(i);
            }
        }

        if(offset + totalBasisFunctions > indices.size())
        {
            indices.resize(offset + totalBasisFunctions);
        }

        std::size_t localBasisIndex = offset;
        for (std::size_t unknown = 0; unknown < numberOfUnknowns; ++unknown)
        {
            if (!offsets_[unknown].includedInIndex_)
            {
                continue;
            }
            std::size_t numberOfElementBasisFunctions = element->getLocalNumberOfBasisFunctions(unknown);
            int elementBasis0 = getGlobalIndex(element, unknown);
            for (std::size_t i = 0; i < numberOfElementBasisFunctions; ++i)
            {
                indices[localBasisIndex++] = elementBasis0 + i;
            }

            std::size_t numberOfFaces = element->getPhysicalGeometry()->getNumberOfFaces();
            for (std::size_t face = 0; face < numberOfFaces; face++)
            {
                std::size_t numberOfFaceBasisFunctions = element->getFace(face)->getLocalNumberOfBasisFunctions(unknown);
                int faceBasis0 = getGlobalIndex(element->getFace(face), unknown);
                for (std::size_t i = 0; i < numberOfFaceBasisFunctions; ++i)
                {
                    indices[localBasisIndex++] = faceBasis0 + i;
                }
            }

            std::size_t numberOfEdges = element->getNumberOfEdges();
            for (std::size_t edge = 0; edge < numberOfEdges; edge++)
            {
                std::size_t numberOfEdgeBasisFunctions = element->getEdge(edge)->getLocalNumberOfBasisFunctions(unknown);
                int edgeBasis0 = getGlobalIndex(element->getEdge(edge), unknown);
                for (std::size_t i = 0; i < numberOfEdgeBasisFunctions; ++i)
                {
                    indices[localBasisIndex++] = edgeBasis0 + i;
                }
            }

            if (mesh_->dimension() > 1)
            {
                std::size_t numberOfNodes = element->getNumberOfNodes();
                for (std::size_t node = 0; node < numberOfNodes; node++)
                {
                    std::size_t numberOfNodeBasisFunctions = element->getNode(node)->getLocalNumberOfBasisFunctions(unknown);
                    int nodeBasis0 = getGlobalIndex(element->getNode(node), unknown);
                    for (std::size_t i = 0; i < numberOfNodeBasisFunctions; ++i)
                    {
                        indices[localBasisIndex++] = nodeBasis0 + i;
                    }
                }
            }
        }

        std::size_t expectedNumberOfBasisFunctions = 0;
        for (std::size_t unknown : includedUknowns_)
            expectedNumberOfBasisFunctions += element->getNumberOfBasisFunctions(unknown);

        logger.assert_debug(localBasisIndex == offset + expectedNumberOfBasisFunctions,
                "Not all basis functions have been assigned an index.");
        return localBasisIndex;
    }

    void GlobalIndexing::getGlobalIndices(const Base::Face *face, std::vector<int> &indices) const
    {
        std::size_t size;
        size = getGlobalIndices(face->getPtrElementLeft(), 0, indices);
        if (face->isInternal())
        {
            size = getGlobalIndices(face->getPtrElementRight(), size, indices);
        }
        // Note this will resize when using mixed geometry faces and when
        // switching from boundary to internal faces. This is an unfortunate
        // consequence from having an easy function signature for use.
        indices.resize(size);
    }

    void GlobalIndexing::reset(Base::MeshManipulatorBase *mesh, Layout layout, const std::vector<std::size_t>* unknowns)
    {
        mesh_ = mesh;
        if (mesh_ != nullptr)
        {
            // We do not support empty meshes as we do not know the number of unknowns,
            // nor is it sensible to create an empty indexing_ for it.
            // But we do support meshes that have empty submeshes
#ifndef HPGEM_USE_MPI
            logger.assert_debug(mesh_->elementColBegin() != mesh_->elementColEnd(), "Empty mesh not supported.");
            // Note as seen below do we assume that the number of unknowns is the same for each element, face, etc.
            numberOfUnknowns_ = (*mesh_->elementColBegin())->getNumberOfUnknowns();
#else
            if(mesh_->elementColBegin() == mesh_->elementColEnd()) {
                numberOfUnknowns_ = 0;
            } else {
                numberOfUnknowns_ = (*mesh_->elementColBegin())->getNumberOfUnknowns();
            }
            MPI_Allreduce(MPI_IN_PLACE, &numberOfUnknowns_, 1, Base::Detail::toMPIType(numberOfUnknowns_), MPI_MAX, Base::MPIContainer::Instance().getComm());
#endif

            offsets_.clear();
            offsets_.resize(numberOfUnknowns_);
            // Setup which unknowns are used
            for (std::size_t i = 0; i < numberOfUnknowns_; ++i)
            {
                // Default to true if no subset is requested.
                offsets_[i].includedInIndex_ = (unknowns == nullptr);
            }
            // Set the subset
            if (unknowns != nullptr)
            {
                logger.assert_always(unknowns->size() > 0, "Building index with no unknowns is not (yet) supported");
                for (std::size_t unknown : *unknowns) {
                    logger.assert_always(unknown < numberOfUnknowns_, "Can not include unknown % as there are only %",
                                         unknown, numberOfUnknowns_);
                    offsets_[unknown].includedInIndex_ = true;
                }
                includedUknowns_ = *unknowns;
                std::sort(includedUknowns_.begin(), includedUknowns_.end());
            }
            else
            {
                includedUknowns_.resize(numberOfUnknowns_);
                for (std::size_t i = 0; i < numberOfUnknowns_; ++i)
                {
                    includedUknowns_[i] = i;
                }
            }
            switch (layout)
            {
                case SEQUENTIAL:
                    constructUnblocked();
                    break;
                case BLOCKED_GLOBAL:
                    constructBlocked(true);
                    break;
                case BLOCKED_PROCESSOR:
                    constructBlocked(false);
                    break;
                default:
                    logger.assert_debug(false, "Unknown index layout %", layout);
            }
        }
        else
        {
            localBasisFunctions_ = 0;
            numberOfUnknowns_ = 0;
            offsets_.clear();
        }
    }

    void GlobalIndexing::constructUnblocked()
    {
        // Construct local ordering.
        std::size_t index = 0;
        for (Base::Element *element : mesh_->getElementsList())
        {
            for (std::size_t unknown : includedUknowns_)
            {
                offsets_[unknown].elementOffsets_[element->getID()] = index;
                index += element->getLocalNumberOfBasisFunctions(unknown);
            }

            for (Base::Face *face : element->getFacesList())
            {
                if (face->getPtrElementLeft() == element)
                {
                    for (std::size_t unknown : includedUknowns_)
                    {
                        offsets_[unknown].faceOffsets_[face->getID()] = index;
                        index += face->getLocalNumberOfBasisFunctions(unknown);
                    }
                }
            }

            for (Base::Edge *edge : element->getEdgesList())
            {
                if (edge->getElement(0) == element)
                {
                    for (std::size_t unknown : includedUknowns_)
                    {
                        offsets_[unknown].edgeOffsets_[edge->getID()] = index;
                        index += edge->getLocalNumberOfBasisFunctions(unknown);
                    }
                }
            }
            // Faces and nodes are the same in 1D
            if (mesh_->dimension() > 1)
            {
                for (Base::Node *node : element->getNodesList())
                {
                    if (node->getElement(0) == element)
                    {
                        for (std::size_t unknown : includedUknowns_)
                        {
                            offsets_[unknown].nodeOffsets_[node->getID()] = index;
                            index += node->getLocalNumberOfBasisFunctions(unknown);
                        }
                    }
                }
            }
        }
        localBasisFunctions_ = index;

        // Offset of the first basis function with MPI, for no MPI 0.
        std::size_t mpiOffset = 0;
#ifdef HPGEM_USE_MPI
        // Compute the blocks in global index space used by each processor
        auto &mpiInstance = Base::MPIContainer::Instance();

        std::size_t n = mpiInstance.getNumberOfProcessors();
        std::size_t rank = mpiInstance.getProcessorID();

        // Compute a the offsets needed to make the local ids unique
        std::vector<std::size_t> globalOffset(n + 1);
        globalOffset[0] = 0;
        globalOffset[rank + 1] = localBasisFunctions_;
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                      globalOffset.data() + 1, 1, Base::Detail::toMPIType(n), mpiInstance.getComm());
        std::partial_sum(globalOffset.begin(), globalOffset.end(), globalOffset.begin());

        mpiOffset = globalOffset[rank];
#endif
        for (std::size_t unknown = 0; unknown < numberOfUnknowns_; ++unknown)
        {
            Offsets& offset = offsets_[unknown];
            if (!offset.includedInIndex_)
                continue;
            offset.setOffset(mpiOffset, 0, localBasisFunctions_);
        }
        communicatePushPullElements();
    }

    void GlobalIndexing::constructBlocked(bool global)
    {
        std::vector<size_t> numberOfBasisFunctions(numberOfUnknowns_);
        localBasisFunctions_ = 0;

        // Number the local basis functions for each unknown by 0...Nu-1, where Nu
        // is the number of local basis functions for this unknown.
        //
        for (std::size_t unknown : includedUknowns_)
        {
            Offsets& offset = offsets_[unknown];
            std::size_t index = 0;
            for (Base::Element *element : mesh_->getElementsList())
            {
                offset.elementOffsets_[element->getID()] = index;
                index += element->getLocalNumberOfBasisFunctions(unknown);

                for (Base::Face *face : element->getFacesList())
                {
                    if (face->getPtrElementLeft() == element)
                    {
                        offset.faceOffsets_[face->getID()] = index;
                        index += face->getLocalNumberOfBasisFunctions(unknown);
                    }
                }

                for (Base::Edge *edge : element->getEdgesList())
                {
                    if (edge->getElement(0) == element)
                    {
                        offset.edgeOffsets_[edge->getID()] = index;
                        index += edge->getLocalNumberOfBasisFunctions(unknown);
                    }
                }
                // Faces and nodes are the same in 1D
                if (mesh_->dimension() > 1)
                {
                    for (Base::Node *node : element->getNodesList())
                    {
                        if (node->getElement(0) == element)
                        {
                            offset.nodeOffsets_[node->getID()] = index;
                            index += node->getLocalNumberOfBasisFunctions(unknown);
                        }
                    }
                }
            }
            numberOfBasisFunctions[unknown] = index;
            localBasisFunctions_ += index;
        }


#ifdef HPGEM_USE_MPI
        auto &mpiInstance = Base::MPIContainer::Instance();

        std::size_t n = mpiInstance.getNumberOfProcessors();
        std::size_t rank = mpiInstance.getProcessorID();
        // Offset due to the basis functions for previous unknowns.
        std::size_t baseOffset = 0;
        std::size_t localBaseOffset = 0;
        if (global)
        {
            for (std::size_t unknown : includedUknowns_)
            {
                // Compute a the offsets needed to make the local ids unique
                std::vector<std::size_t> globalOffset(n + 1);
                globalOffset[0] = baseOffset;
                globalOffset[rank + 1] = numberOfBasisFunctions[unknown];
                MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                              globalOffset.data() + 1, 1, Base::Detail::toMPIType(n), mpiInstance.getComm());
                std::partial_sum(globalOffset.begin(), globalOffset.end(), globalOffset.begin());

                std::size_t mpiOffset = globalOffset[rank];

                offsets_[unknown].setOffset(mpiOffset, localBaseOffset, numberOfBasisFunctions[unknown]);
                // globalOffset[n] contains the sum of all the basis functions for
                // this unknown plus those for previous unknowns.
                baseOffset = globalOffset[n];
                localBaseOffset += numberOfBasisFunctions[unknown];
            }
        }
        else
        {
            std::vector<std::size_t> globalOffset(n + 1);
            globalOffset[0] = 0;
            globalOffset[rank+1] = localBasisFunctions_;
            MPI_Allgather(
                    MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                    globalOffset.data() + 1, 1, Base::Detail::toMPIType(n),
                    mpiInstance.getComm());
            std::partial_sum(globalOffset.begin(), globalOffset.end(), globalOffset.begin());
            std::size_t localOffset = 0;
            for (std::size_t unknown : includedUknowns_)
            {
                offsets_[unknown].setOffset(globalOffset[rank] + localOffset, localOffset,
                                            numberOfBasisFunctions[unknown]);
                localOffset += numberOfBasisFunctions[unknown];
            }
        }
#else
        // Compute the offsets for all the unknowns
        std::size_t offset = 0;
        for (std::size_t unknown : includedUknowns_)
        {
            offsets_[unknown].setOffset(offset, offset, numberOfBasisFunctions[unknown]);
            offset += numberOfBasisFunctions[unknown];
        }
#endif
        communicatePushPullElements();
    }

    void GlobalIndexing::verifyCompleteIndex() const
    {
        // Check basic consistency for unused unknowns
        for (std::size_t unknown = 0; unknown < numberOfUnknowns_; ++unknown)
        {
            if (offsets_[unknown].includedInIndex_)
                continue;
            logger.assert_debug(offsets_[unknown].blockStart_ == 0, "Non zero block start for excluded unknown");
            logger.assert_debug(offsets_[unknown].numberOfBasisFunctionsInBlock_ == 0, "Non zero number of basis functions for excluded unknown");
            logger.assert_debug(offsets_[unknown].elementOffsets_.empty(), "Element offsets for excluded unknown.");
            logger.assert_debug(offsets_[unknown].faceOffsets_.empty(), "Face offsets for excluded unknown.");
            logger.assert_debug(offsets_[unknown].edgeOffsets_.empty(), "Edge offsets for excluded unknown.");
            logger.assert_debug(offsets_[unknown].nodeOffsets_.empty(), "Node offsets for excluded unknown.");
        }

        // Verify the used unknowns
        for( const Base::Element *element : mesh_->getElementsList(Base::IteratorType::GLOBAL))
        {
            for(std::size_t unknown : includedUknowns_)
            {
                // Internally checks if it is present
                getGlobalIndex(element, unknown);
            }
        }

        for( const Base::Face *face : mesh_->getFacesList(Base::IteratorType::GLOBAL))
        {
            for(std::size_t unknown : includedUknowns_)
            {
                // Internally checks if it is present
                getGlobalIndex(face, unknown);
            }
        }

        for( const Base::Edge *edge : mesh_->getEdgesList(Base::IteratorType::GLOBAL))
        {
            for(std::size_t unknown : includedUknowns_)
            {
                // Internally checks if it is present
                getGlobalIndex(edge, unknown);
            }
        }

        if (mesh_->dimension() > 1)
        {
            for (const Base::Node *node : mesh_->getNodesList(Base::IteratorType::GLOBAL))
            {
                for (std::size_t unknown : includedUknowns_)
                {
                    // Internally checks if it is present
                    getGlobalIndex(node, unknown);
                }
            }
        }
    }


#ifdef HPGEM_USE_MPI

    void GlobalIndexing::createInitialMessage(
            const std::vector<Base::Element*>& elements,
            std::vector<std::size_t> &message,
            std::size_t targetProcessor,
            std::set<std::size_t>& secondRoundTags) const
    {
        // To prevent sending information for faces, edges and nodes more than
        // once, we store their tags. Alternatively we could do something
        // similar by asking face->getOwningElement() == element, except that
        // this fails when the owning element is not a ghost element of the
        // targetProcessor. Instead of compensating for the case, we take the
        // easier approach of storing the message tags of all faces, edges and
        // nodes that have already been send.
        std::set<std::size_t> sendTags;

        // For format see processMessage method
        for(Base::Element *element : elements)
        {
            elementMessage(element->getID(), message);

            for (auto face : element->getFacesList())
            {
                std::size_t tag = 4*face->getID() + 1;
                if (face->isOwnedByCurrentProcessor())
                {
                    // Prevent sending it for both sides of the face by only
                    // sending it if the tag was actually inserted
                    if(sendTags.insert(tag).second)
                    {
                        faceMessage(face->getID(), message);
                    }
                }
                // As the pushElements are all owned by the current processor
                // and the face is not owned by the current processor, we know
                // that the left element must exist and belong to a different
                // processor.
                else if(face->getPtrElementLeft()->getOwner() != targetProcessor)
                {
                    secondRoundTags.insert(tag);
                }
            }

            // Same for edge, nodes
            for (auto edge : element->getEdgesList())
            {
                std::size_t tag = 4*edge->getID() + 2;
                if (edge->isOwnedByCurrentProcessor())
                {
                    if(sendTags.insert(tag).second)
                    {
                        edgeMessage(edge->getID(), message);
                    }
                }
                else
                {
                    bool targetIsNeighbour = false;
                    for(Base::Element *edgeElem : edge->getElements())
                    {
                        if(edgeElem->getOwner() == targetProcessor)
                        {
                            targetIsNeighbour = true;
                            break;
                        }
                    }
                    if(!targetIsNeighbour)
                    {
                        secondRoundTags.insert(tag);
                    }
                }
            }

            if (mesh_->dimension() > 1)
            {
                for (auto node : element->getNodesList())
                {
                    std::size_t tag = 4*node->getID()+3;
                    if (node->isOwnedByCurrentProcessor())
                    {
                        if(sendTags.insert(tag).second)
                        {
                            nodeMessage(node->getID(), message);
                        }
                    }
                    else
                    {
                        bool targetIsNeighbour = false;
                        for(Base::Element* nodeElem : node->getElements())
                        {
                            if(nodeElem->getOwner() == targetProcessor)
                            {
                                targetIsNeighbour = true;
                                break;
                            }
                        }
                        if(!targetIsNeighbour)
                        {
                            secondRoundTags.insert(tag);
                        }
                    }

                }
            }
        }
    }

    void GlobalIndexing::createSecondMessage(const std::set<std::size_t> &tags, std::vector<std::size_t> &message) const
    {
        // For format see the processMessage method
        message.clear();
        for(std::size_t tag : tags)
        {
            switch (tag % 4)
            {
                case 1:
                    faceMessage(tag/4, message);
                    break;
                case 2:
                    edgeMessage(tag/4, message);
                    break;
                case 3:
                    nodeMessage(tag/4, message);
                    break;
                default:
                    logger.assert_always(false, "Invalid tag");
            }
        }
    }

    void GlobalIndexing::elementMessage(std::size_t elementId, std::vector<std::size_t> &message) const
    {
        message.emplace_back(4*elementId);
        for(std::size_t unknown : includedUknowns_)
        {
            message.emplace_back(offsets_.at(unknown).elementOffsets_.at(elementId));
        }
    }

    void GlobalIndexing::faceMessage(std::size_t faceId, std::vector<std::size_t> &message) const
    {
        message.emplace_back(4*faceId+1);
        for(std::size_t unknown : includedUknowns_)
        {
            message.emplace_back(offsets_.at(unknown).faceOffsets_.at(faceId));
        }
    }

    void GlobalIndexing::edgeMessage(std::size_t edgeId, std::vector<std::size_t> &message) const
    {
        message.emplace_back(4*edgeId+2);
        for(std::size_t unknown : includedUknowns_)
        {
            if (!offsets_[unknown].includedInIndex_)
                continue;
            message.emplace_back(offsets_.at(unknown).edgeOffsets_.at(edgeId));
        }
    }

    void GlobalIndexing::nodeMessage(std::size_t nodeId, std::vector<std::size_t> &message) const
    {
        message.emplace_back(4*nodeId+3);
        for(std::size_t unknown : includedUknowns_)
        {
            if (!offsets_[unknown].includedInIndex_)
                continue;
            message.emplace_back(offsets_.at(unknown).nodeOffsets_.at(nodeId));
        }
    }

    void GlobalIndexing::processMessage(const std::vector<std::size_t> &message, std::size_t count)
    {
        // The body of the message in the first and second round share the same
        // format. They consist of multiple parts, where each part consists of
        //  - A tag, identifying the geometric object (element, face, etc.)
        //    which the next information is about.
        //  - numberOfUnknowns_ global indices for the geometric object in
        //    order of the unknown.
        // The tags are four times the id of the object plus 0 for elements 1
        // for faces, 2 for edges and 3 for nodes. This ensures that all tags
        // uniquely identify a geometric object without having to communicate
        // before hand.

        logger.assert_debug(count % (1+numberOfUnknowns_) == 0,
                            "Incorrect message of size %", count);
        for(std::size_t offset = 0; offset < count; offset+=(1+numberOfUnknowns_))
        {
            std::size_t tag = message[offset];
            std::size_t id = tag / 4;
            switch (tag % 4)
            {
                case 0:
                    for(std::size_t unknown : includedUknowns_)
                        offsets_[unknown].elementOffsets_[id] = message[offset + 1 + unknown];
                    break;
                case 1:
                    for(std::size_t unknown : includedUknowns_)
                        offsets_[unknown].faceOffsets_[id] = message[offset + 1 + unknown];
                    break;
                case 2:
                    for(std::size_t unknown : includedUknowns_)
                        offsets_[unknown].edgeOffsets_[id] = message[offset + 1 + unknown];
                    break;
                case 3:
                    for(std::size_t unknown : includedUknowns_)
                        offsets_[unknown].nodeOffsets_[id] = message[offset + 1 + unknown];
                    break;
                default:
                    logger.assert_always(false, "Error invalid tag %", tag);
            }
        }
    }

    std::size_t GlobalIndexing::probeAndReceive(std::vector<std::size_t> &receiveMessage) const
    {
        auto &mpiInstance = Base::MPIContainer::Instance();
        // Sample type
        MPI_Datatype mpiType;
        {
            std::size_t sample = 1;
            mpiType = Base::Detail::toMPIType(sample);
        }

        // Probe the message size
        MPI_Status status;
        MPI_Probe(MPI_ANY_SOURCE, 0, mpiInstance.getComm(), &status);
        // Make place for the messages
        int count;
        MPI_Get_count(&status, mpiType, & count);
        logger.assert_debug(count >= 0, "expected to receive a message of more than % integers", count);

        // To prevent unnecessary (de)allocation we only grow the message
        // vector, its size is therefore only an upper bound on the message
        // size.
        if(static_cast<std::size_t>(count) > receiveMessage.size())
            receiveMessage.resize(count);
        // Actual receive
        MPI_Recv(receiveMessage.data(), count, mpiType,
                 status.MPI_SOURCE, 0, mpiInstance.getComm(), MPI_STATUS_IGNORE);
        return count;
    }

#endif


    void GlobalIndexing::communicatePushPullElements()
    {
#ifdef HPGEM_USE_MPI
        // The global indices on ghost elements, faces, etc. are only known on
        // their owning processor. This method communicates these to all the
        // neighbouring processors. This is achieved in two rounds.
        //
        // In the first round, we communicate the global indices of all elements
        // in the pushElements list, along with all the global indices of the
        // faces, edges and nodes that border such elements and are available
        // (i.e. are owned by the current processor). We then receive the same
        // information from other processor for our own ghost elements.
        //
        // After this first round we all the information about the ghost
        // elements, but the information of the enclosing faces, edges and nodes
        // can be missing. Consider as example the following line
        //
        //    [0]---0---[0]---1---[1]---2---[2]
        //
        // Where for each of the three line segments (elements) and each node
        // (in brackets) the owning processor is marked. Processor two has only
        // part of the line:
        //
        //              [0]---1---[1]---2---[2]
        //
        // Where the element owned by processor 1 is a ghost element. After the
        // first transfer round the globalIndices for the element and node
        // belonging to processor 1 are communicated to processor 2, but the
        // index for the node owned by processor 0 is not, as it does not belong
        // to processor 1, nor is there a ghost element from processor 0
        // communicating the information. Hence we need a second round.
        //
        // In the second round we communicate the global indices of faces, edges
        // and nodes that:
        //  - border an element in the pushElements list;
        //  - are not owned by the current processor;
        //  - are not adjacent to an element that is owned by the target
        //    processor.
        // As we own all the elements in the pushElements list, we do have this
        // information after round one. Furthermore, we reduce the amount of
        // messages as much as possible with the second two conditions, for which
        // the target processor has already received the information in round one.
        // After sending the information we again receive and process the
        // information send from the other processors.


        // Setup //
        ///////////
        auto &mpiInstance = Base::MPIContainer::Instance();
        // Number of messages in the first round
        std::size_t firstRoundSends = mesh_->getPushElements().size();
        std::vector<MPI_Request> sendRequests (firstRoundSends);
        // As we do asynchronous sending we need to keep the message available
        // till after the send finishes. We therefore store a vector of the raw
        // messages here.
        std::vector<std::vector<std::size_t>> sendMessages (firstRoundSends);
        // For each target processor, the set of tags that need to be send in
        // the second round.
        std::map<std::size_t, std::set<std::size_t>> secondRound;

        MPI_Datatype mpiType;
        {
            std::size_t sample = 1;
            mpiType = Base::Detail::toMPIType(sample);
        }

        // First round - sending //
        ///////////////////////////
        std::size_t index = 0;
        for(const auto& entry : mesh_->getPushElements())
        {
            int targetProcessor = entry.first;

            auto& message = sendMessages[index];
            auto& secondRoundTags = secondRound[targetProcessor];

            createInitialMessage(entry.second, message, targetProcessor, secondRoundTags);
            bool needsSecondRound = !secondRoundTags.empty();
            // Append to the message whether we need a second round.
            if (needsSecondRound)
            {
                message.emplace_back(1);
            }
            else
            {
                message.emplace_back(0);
                // Prevent second round message
                secondRound.erase(targetProcessor);
            }
            // Actual sending
            MPI_Request& request = sendRequests[index];
            MPI_Isend(message.data(), message.size(), mpiType,
                    targetProcessor, 0, mpiInstance.getComm(), &request);
            ++index;
        }

        // First round - receive //
        ///////////////////////////
        std::size_t numReceives = mesh_->getPullElements().size();
        // Storage for the receiving message
        std::vector<std::size_t> receiveMessage;
        // Counter the number of messages to receive in the second round
        std::size_t secondRoundReceives = 0;
        for(std::size_t receiveId = 0; receiveId < numReceives; ++receiveId)
        {
            int count = probeAndReceive(receiveMessage);
            // Process and remove the tag for whether or not a second round will be send
            if(receiveMessage[count - 1])
                ++secondRoundReceives;
            --count;
            // Process the message
            processMessage(receiveMessage, count);
        }

        // First round - closing //
        ///////////////////////////

        // Await finish of all requests to prevent early cleanup of
        // the sendMessages vector.
        MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE);
        sendRequests.clear();
        // We need to ensure that all processes have finished receiving.
        // Otherwise a race condition might occur where a message from the
        // second round is received by another processor expecting first round
        // messages.
        MPI_Barrier(mpiInstance.getComm());


        // Second round - sending //
        ////////////////////////////
        index = 0;
        for(auto& entry : secondRound)
        {
            // Reuse the message storage space from first round, preventing
            // extra allocation.
            std::vector<std::size_t>& sendMessage = sendMessages[index++];
            sendMessage.clear();
            createSecondMessage(entry.second, sendMessage);
            sendRequests.emplace_back();
            MPI_Request& request = sendRequests[sendRequests.size() - 1];
            MPI_Isend(sendMessage.data(), sendMessage.size(), mpiType,
                      entry.first, 0, mpiInstance.getComm(), &request);
        }

        // Second round - receiving //
        //////////////////////////////
        for(std::size_t receiveId = 0; receiveId < secondRoundReceives; ++receiveId)
        {

            int count = probeAndReceive(receiveMessage);
            processMessage(receiveMessage, count);


        }

        // Second round - closing //
        ////////////////////////////

        // Await the completion sending all messages
        MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE);
        sendRequests.clear();
        // Prevent interference by ensuring that all communication has finished.
        MPI_Barrier(mpiInstance.getComm());

        // Verify that everything is complete
        verifyCompleteIndex();
#endif
    }


    void GlobalIndexing::Offsets::setOffset(int globalOffset, int localOffset, int numberOfBasisFunctions)
    {
        logger.assert_debug(localOffset >= 0, "Negative local offset");
        logger.assert_debug(globalOffset >= 0, "Negative global offset");
        this->localOffset_ = localOffset;
        this->blockStart_ = globalOffset;
        this->numberOfBasisFunctionsInBlock_ = numberOfBasisFunctions;
        if (globalOffset != 0)
        {
            {
                auto iterEnd = elementOffsets_.end();
                for (auto iter = elementOffsets_.begin(); iter != iterEnd; ++iter)
                {
                    iter->second += globalOffset;
                }
            }
            {
                auto iterEnd = faceOffsets_.end();
                for (auto iter = faceOffsets_.begin(); iter != iterEnd; ++iter)
                {
                    iter->second += globalOffset;
                }
            }
            {
                auto iterEnd = edgeOffsets_.end();
                for (auto iter = edgeOffsets_.begin(); iter != iterEnd; ++iter)
                {
                    iter->second += globalOffset;
                }
            }
            {
                auto iterEnd = nodeOffsets_.end();
                for (auto iter = nodeOffsets_.begin(); iter != iterEnd; ++iter)
                {
                    iter->second += globalOffset;
                }
            }
        }

    }
}

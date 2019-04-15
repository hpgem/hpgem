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
            : numberOfUnknowns_(0), offsets(0), localBasisFunctions_(0), meshDimension(0)
    {}

    GlobalIndexing::GlobalIndexing(Base::MeshManipulatorBase *mesh, Layout layout)
    {
        reset(mesh, layout);
    }

    void GlobalIndexing::getGlobalIndices(const Base::Element *element, std::size_t offset, std::vector<int> &indices)
    {
        logger.assert_debug(element != nullptr, "Null pointer as element");
        indices.resize(offset + element->getTotalNumberOfBasisFunctions());

        std::size_t localBasisIndex = offset;
        for (std::size_t unknown = 0; unknown < element->getNumberOfUnknowns(); ++unknown)
        {
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

            if (meshDimension > 1)
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

        logger.assert_debug(localBasisIndex == indices.size(), "Not all basis functions have been assigned an index.");
    }

    void GlobalIndexing::getGlobalIndices(const Base::Face *face, std::vector<int> &indices)
    {
        getGlobalIndices(face->getPtrElementLeft(), indices);
        if (face->isInternal())
        {
            getGlobalIndices(face->getPtrElementRight(), indices.size(), indices);
        }
    }

    void GlobalIndexing::reset(Base::MeshManipulatorBase *mesh, Layout layout)
    {
        if (mesh != nullptr)
        {
            meshDimension = mesh->dimension();
            // We do not support empty meshes as we do not know the number of unknowns,
            // nor is it sensible to create an empty indexing_ for it.
            logger.assert_debug(mesh->elementColBegin() != mesh->elementColEnd(), "Empty mesh not supported.");
            // Note as seen below do we assume that the number of unknowns is the same for each element, face, etc.
            numberOfUnknowns_ = (*mesh->elementColBegin())->getNumberOfUnknowns();
            offsets.clear();
            offsets.resize(numberOfUnknowns_);
            switch (layout)
            {
                case SEQUENTIAL:
                    constructUnblocked(*mesh);
                    break;
                case BLOCKED_GLOBAL:
                    constructBlocked(*mesh, true);
                    break;
                case BLOCKED_PROCESSOR:
                    constructBlocked(*mesh, false);
                    break;
                default:
                    logger.assert_debug(false, "Unknown index layout %", layout);
            }
        }
        else
        {
            meshDimension = 0;
            numberOfUnknowns_ = 0;
            offsets.clear();
        }
    }

    void GlobalIndexing::constructUnblocked(Base::MeshManipulatorBase &mesh)
    {
        // Construct local ordering.
        std::size_t index = 0;
        for (Base::Element *element : mesh.getElementsList())
        {
            for (std::size_t unknown = 0; unknown < numberOfUnknowns_; ++unknown)
            {
                offsets[unknown].elementOffsets_[element->getID()] = index;
                index += element->getLocalNumberOfBasisFunctions(unknown);
            }

            for (Base::Face *face : element->getFacesList())
            {
                if (face->getPtrElementLeft() == element)
                {
                    for (std::size_t unknown = 0; unknown < numberOfUnknowns_; ++unknown)
                    {
                        offsets[unknown].faceOffsets_[face->getID()] = index;
                        index += face->getLocalNumberOfBasisFunctions(unknown);
                    }
                }
            }

            for (Base::Edge *edge : element->getEdgesList())
            {
                if (edge->getElement(0) == element)
                {
                    for (std::size_t unknown = 0; unknown < numberOfUnknowns_; ++unknown)
                    {
                        offsets[unknown].edgeOffsets_[edge->getID()] = index;
                        index += edge->getLocalNumberOfBasisFunctions(unknown);
                    }
                }
            }
            // Faces and nodes are the same in 1D
            if (meshDimension > 1)
            {
                for (Base::Node *node : element->getNodesList())
                {
                    if (node->getElement(0) == element)
                    {
                        for (std::size_t unknown = 0; unknown < numberOfUnknowns_; ++unknown)
                        {
                            offsets[unknown].nodeOffsets_[node->getID()] = index;
                            index += node->getLocalNumberOfBasisFunctions(unknown);
                        }
                    }
                }
            }
        }
        localBasisFunctions_ = index;

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

        std::size_t mpiOffset = globalOffset[rank];

        for (std::size_t unknown = 0; unknown < numberOfUnknowns_; ++unknown)
        {
            offsets[unknown].setOffset(mpiOffset, 0, localBasisFunctions_);
        }
#else
        for (std::size_t unknown = 0; unknown < numberOfUnknowns_; ++unknown)
        {
            offsets[unknown].setOffset(0, 0, localBasisFunctions_);
        }
#endif
        communicatePushPullElements(mesh);
    }

    void GlobalIndexing::constructBlocked(Base::MeshManipulatorBase &mesh, bool global)
    {
        std::vector<size_t> numberOfBasisFunctions(numberOfUnknowns_);
        localBasisFunctions_ = 0;

        // Number the local basis functions for each unknown by 0...Nu-1, where Nu
        // is the number of local basis functions for this unknown.
        //
        for (std::size_t unknown = 0; unknown < numberOfUnknowns_; ++unknown)
        {
            std::size_t index = 0;
            for (Base::Element *element : mesh.getElementsList())
            {
                offsets[unknown].elementOffsets_[element->getID()] = index;
                index += element->getLocalNumberOfBasisFunctions(unknown);

                for (Base::Face *face : element->getFacesList())
                {
                    if (face->getPtrElementLeft() == element)
                    {
                        offsets[unknown].faceOffsets_[face->getID()] = index;
                        index += face->getLocalNumberOfBasisFunctions(unknown);
                    }
                }

                for (Base::Edge *edge : element->getEdgesList())
                {
                    if (edge->getElement(0) == element)
                    {
                        offsets[unknown].edgeOffsets_[edge->getID()] = index;
                        index += edge->getLocalNumberOfBasisFunctions(unknown);
                    }
                }
                // Faces and nodes are the same in 1D
                if (meshDimension > 1)
                {
                    for (Base::Node *node : element->getNodesList())
                    {
                        if (node->getElement(0) == element)
                        {
                            offsets[unknown].nodeOffsets_[node->getID()] = index;
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
            for (std::size_t unknown = 0; unknown < numberOfUnknowns_; ++unknown)
            {
                // Compute a the offsets needed to make the local ids unique
                std::vector<std::size_t> globalOffset(n + 1);
                globalOffset[0] = baseOffset;
                globalOffset[rank + 1] = numberOfBasisFunctions[unknown];
                MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                              globalOffset.data() + 1, 1, Base::Detail::toMPIType(n), mpiInstance.getComm());
                std::partial_sum(globalOffset.begin(), globalOffset.end(), globalOffset.begin());

                std::size_t mpiOffset = globalOffset[rank];

                offsets[unknown].setOffset(mpiOffset, localBaseOffset, numberOfBasisFunctions[unknown]);
                // globalOffset[n+1] contains the sum of all the basis functions for
                // this unknown plus those for previous unknowns.
                baseOffset = globalOffset[n + 1];
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
            for (std::size_t unknown = 0; unknown < numberOfUnknowns_; ++unknown)
            {
                offsets[unknown].setOffset(globalOffset[rank] + localOffset, localOffset,
                        numberOfBasisFunctions[unknown]);
                localOffset += numberOfBasisFunctions[unknown];
            }
        }
#else
        // Compute the offsets for all the unknowns
        std::size_t offset = 0;
        for (std::size_t unknown = 0; unknown < numberOfUnknowns_; ++unknown)
        {
            offsets[unknown].setOffset(offset, offset, numberOfBasisFunctions[unknown]);
            offset += numberOfBasisFunctions[unknown];
        }
#endif
        communicatePushPullElements(mesh);
    }

#ifdef HPGEM_USE_MPI

    void GlobalIndexing::createInitialMessage(
            const std::vector<Base::Element*>& elements,
            std::vector<std::size_t> &message,
            std::size_t targetProcessor,
            std::vector<std::size_t>& secondRoundTags)
    {
        for(Base::Element *element : elements)
        {
            elementMessage(element->getID(), message);

            for (auto face : element->getFacesList())
            {
                if (face->isOwnedByCurrentProcessor())
                {
                    faceMessage(face->getID(), message);
                }
                else if(face->getPtrElementLeft()->getOwner() != targetProcessor)
                {
                    secondRoundTags.emplace_back(4*face->getID() + 1);
                }
            }

            for (auto edge : element->getEdgesList())
            {
                if (edge->isOwnedByCurrentProcessor())
                {
                    edgeMessage(edge->getID(), message);
                }
                else
                {
                    bool targetIsNeighbour = false;
                    for(Base::Element *edgeElem : edge->getElements())
                    {
                        targetIsNeighbour |= edgeElem->getOwner() == targetProcessor;
                    }
                    if(!targetIsNeighbour)
                    {
                        secondRoundTags.emplace_back(4*edge->getID() + 2);
                    }
                }
            }

            if (meshDimension > 1)
            {
                for (auto node : element->getNodesList())
                {
//                    if(node->getID() == 43)
//                    {
//                        auto &mpiInstance = Base::MPIContainer::Instance();
//                        std::cout << "Node 43 " << mpiInstance.getProcessorID() << "->" << targetProcessor;
//                    }
                    if (node->isOwnedByCurrentProcessor())
                    {
                        if(node->getID() == 43)
                        {
//                            std::cout << " Directly" << std::endl;
//                            std::cout << "Neighbouring owners";
//                            for(auto *nodeElem : node->getElements())
//                            {
//                                std::cout << " (" << nodeElem->getID()<< ":"<< nodeElem->getOwner() << ")";
//                            }
//                            std::cout << std::endl;
                        }
                        nodeMessage(node->getID(), message);
                    }
                    else
                    {
                        bool targetIsNeighbour = false;
                        for(Base::Element* nodeElem : node->getElements())
                        {
                            bool neigh = nodeElem->getOwner() == targetProcessor;
//                            if(neigh && node->getID() == 43)
//                            {
//                                std::cout << " Skipped due to element " << nodeElem->getID() << " owned by " << nodeElem->getOwner() << std::endl;
//                            }
                            targetIsNeighbour |= nodeElem->getOwner() == targetProcessor;
                            if(neigh)
                                break;
                        }
                        if(!targetIsNeighbour)
                        {
//                            if(node->getID() == 43)
//                            {
//                                std::cout << " Second" << std::endl;
//                            }
                            secondRoundTags.emplace_back(4*node->getID() + 3);
                        }
                    }

                }
            }
        }
    }

    void GlobalIndexing::elementMessage(std::size_t elementId, std::vector<std::size_t> &message)
    {
        message.emplace_back(4*elementId);
        for(std::size_t unknown = 0; unknown < numberOfUnknowns_; ++unknown)
        {
            message.emplace_back(offsets[unknown].elementOffsets_[elementId]);
        }
    }

    void GlobalIndexing::faceMessage(std::size_t faceId, std::vector<std::size_t> &message)
    {
        message.emplace_back(4*faceId+1);
        for(std::size_t unknown = 0; unknown < numberOfUnknowns_; ++unknown)
        {
            message.emplace_back(offsets[unknown].faceOffsets_[faceId]);
        }
    }

    void GlobalIndexing::edgeMessage(std::size_t edgeId, std::vector<std::size_t> &message)
    {
        message.emplace_back(4*edgeId+2);
        for(std::size_t unknown = 0; unknown < numberOfUnknowns_; ++unknown)
        {
            message.emplace_back(offsets[unknown].edgeOffsets_[edgeId]);
        }
    }

    void GlobalIndexing::nodeMessage(std::size_t nodeId, std::vector<std::size_t> &message)
    {
        message.emplace_back(4*nodeId+3);
        for(std::size_t unknown = 0; unknown < numberOfUnknowns_; ++unknown)
        {
            message.emplace_back(offsets[unknown].nodeOffsets_[nodeId]);
        }
    }

    void GlobalIndexing::processMessage(std::vector<std::size_t> &message, std::size_t count)
    {
        logger.assert_debug(count % (1+numberOfUnknowns_) == 0,
                            "Incorrect message of size %", count);
        for(std::size_t offset = 0; offset < count; offset+=(1+numberOfUnknowns_))
        {
            std::size_t tag = message[offset];
            std::size_t id = tag / 4;
            switch (tag % 4)
            {
                case 0:
                    for(std::size_t unknown = 0; unknown < numberOfUnknowns_; ++unknown)
                        offsets[unknown].elementOffsets_[id] = message[offset + 1 + unknown];
                    break;
                case 1:
                    for(std::size_t unknown = 0; unknown < numberOfUnknowns_; ++unknown)
                        offsets[unknown].faceOffsets_[id] = message[offset + 1 + unknown];
                    break;
                case 2:
                    for(std::size_t unknown = 0; unknown < numberOfUnknowns_; ++unknown)
                        offsets[unknown].edgeOffsets_[id] = message[offset + 1 + unknown];
                    break;
                case 3:
                    for(std::size_t unknown = 0; unknown < numberOfUnknowns_; ++unknown)
                        offsets[unknown].nodeOffsets_[id] = message[offset + 1 + unknown];
                    break;
                default:
                    logger.assert_always(false, "Error invalid tag %", tag);
            }
        }
    }

#endif


    void GlobalIndexing::communicatePushPullElements(Base::MeshManipulatorBase &mesh)
    {
#ifdef HPGEM_USE_MPI
        auto &mpiInstance = Base::MPIContainer::Instance();
        std::size_t sends = mesh.getPushElements().size();
        std::vector<MPI_Request> requests (sends);
        // Storage space for sending all the messages
        std::vector<std::vector<std::size_t>> sendMessages (sends);
        std::map<std::size_t, std::vector<std::size_t>> secondRound;

        MPI_Datatype mpiType;
        {
            std::size_t sample = 1;
            mpiType = Base::Detail::toMPIType(sample);
        }

        std::size_t index = 0;
        for(const auto& entry : mesh.getPushElements())
        {
            int targetProcessor = entry.first;
            // Note reference to keep the message alive till the end of the
            // method
            std::vector<std::size_t>& message = sendMessages[index];
            std::vector<std::size_t>& seconds = secondRound[targetProcessor];

            createInitialMessage(entry.second, message, targetProcessor, seconds);
            // Add extra note whether a second round is needed or not.
            bool needsSecondRound = !seconds.empty();
            message.emplace_back(needsSecondRound ? 1 : 0);
            if (!needsSecondRound)
                secondRound.erase(targetProcessor);
            // Send message
            MPI_Request& request = requests[index];
            std::cout << mpiInstance.getProcessorID() << " Sending to "
                << targetProcessor << " count: " << message.size() << ", seconds: "
                << seconds.size() << std::endl;
            MPI_Isend(message.data(), message.size(), mpiType,
                    targetProcessor, 0, mpiInstance.getComm(), &request);
            ++index;
        }

        std::size_t receives = mesh.getPullElements().size();
        std::vector<std::size_t> receiveMessage;
        std::size_t secondRounds = 0;
        for(std::size_t receiveId = 0; receiveId < receives; ++receiveId)
        {

            MPI_Status status;
            MPI_Probe(MPI_ANY_SOURCE, 0, mpiInstance.getComm(), &status);
            // Make place for the messages
            int count;
            MPI_Get_count(&status, mpiType, & count);
            // To prevent unnecessary (de)allocation we only grow the message
            // vector, its size is therefore only an upper bound on the message
            // size.
            if(count > receiveMessage.size())
                receiveMessage.resize(count);
            // Actual receive
            MPI_Recv(receiveMessage.data(), count, mpiType,
                    status.MPI_SOURCE, 0, mpiInstance.getComm(), MPI_STATUS_IGNORE);
            // Process and remove the second round tag
            if(receiveMessage[count - 1])
                ++secondRounds;
            --count;
            // Process the message
            processMessage(receiveMessage, count);
        }
        // Await finish of all requests to prevent early cleanup of
        // the sendMessages vector.
        MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
        requests.clear();
        std::cout << mpiInstance.getProcessorID() << " Finished first round" << std::endl;
        // Prevent interference between first and second round
        MPI_Barrier(mpiInstance.getComm());
        //
        // Second rounds.
        index = 0;
        for(auto& entry : secondRound)
        {
            // repurpose the messages vector
            std::vector<std::size_t>& sendMessage = sendMessages[index++];
            sendMessage.clear();
            for(std::size_t tag : entry.second)
            {
                switch (tag % 4)
                {
                    case 1:
                        faceMessage(tag/4, sendMessage);
                        break;
                    case 2:
                        edgeMessage(tag/4, sendMessage);
                        break;
                    case 3:
                        nodeMessage(tag/4, sendMessage);
                        break;
                    default:
                        logger.assert_always(false, "Invalid tag");
                }
            }
            requests.emplace_back();
            MPI_Request& request = requests[requests.size() - 1];
            std::cout << mpiInstance.getProcessorID() << " Sending second message to "
                      << entry.first << " count: " << sendMessage.size() << std::endl;
            MPI_Isend(sendMessage.data(), sendMessage.size(), mpiType,
                      entry.first, 0, mpiInstance.getComm(), &request);
        }
        // Receive, second rounds

        std::cout << mpiInstance.getProcessorID() << " Finished sending second rounds, receiving "
                << secondRounds << std::endl;
        for(std::size_t receiveId = 0; receiveId < secondRounds; ++receiveId)
        {

            MPI_Status status;
            MPI_Probe(MPI_ANY_SOURCE, 0, mpiInstance.getComm(), &status);
            // Make place for the messages
            int count = 0;
            MPI_Get_count(&status, mpiType, & count);
            // To prevent unnecessary (de)allocation we only grow the message
            // vector, its size is therefore only an upper bound on the message
            // size.
            if(count > receiveMessage.size())
                receiveMessage.resize(count);

            // Actual receive
            MPI_Recv(receiveMessage.data(), count, mpiType,
                     status.MPI_SOURCE, 0, mpiInstance.getComm(), MPI_STATUS_IGNORE);
            std::cout << mpiInstance.getProcessorID() << " Received second message from "
                      << status.MPI_SOURCE << " size " << count << std::endl;
            // Process the message
            processMessage(receiveMessage, count);
        }
        // Await finish of all requests to prevent early cleanup of
        // the sendMessages vector.
        MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
        requests.clear();
        // Prevent interfering messages
        MPI_Barrier(mpiInstance.getComm());

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

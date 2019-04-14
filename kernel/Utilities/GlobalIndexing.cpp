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

    void GlobalIndexing::communicatePushPullElements(Base::MeshManipulatorBase &mesh)
    {
#ifdef HPGEM_USE_MPI
        auto &mpiInstance = Base::MPIContainer::Instance();
        // Push/pull boundary elements
        for (std::size_t unknown = 0; unknown < numberOfUnknowns_; ++unknown)
        {

            for (auto entry : mesh.getPullElements())
            {
                int sourceProcessor = entry.first;
                for (Base::Element *element : entry.second)
                {
                    //the id's of elements, faces edges and nodes are interleaved so all information can be communicated simulateously without tag collisions

                    mpiInstance.receive(offsets[unknown].elementOffsets_[element->getID()], sourceProcessor,
                                        4 * element->getID());
                    for (auto face : element->getFacesList())
                    {
                        if (offsets[unknown].faceOffsets_.count(face->getID()) == 0)
                        {
                            mpiInstance.receive(offsets[unknown].faceOffsets_[face->getID()], sourceProcessor,
                                                4 * face->getID() + 1);
                        }
                    }
                    for (auto edge : element->getEdgesList())
                    {
                        if (offsets[unknown].edgeOffsets_.count(edge->getID()) == 0)
                        {
                            mpiInstance.receive(offsets[unknown].edgeOffsets_[edge->getID()], sourceProcessor,
                                                4 * edge->getID() + 2);
                        }
                    }
                    if (meshDimension > 1)
                    {
                        for (auto node : element->getNodesList())
                        {
                            if (offsets[unknown].nodeOffsets_.count(node->getID()) == 0)
                            {
                                mpiInstance.receive(offsets[unknown].nodeOffsets_[node->getID()], sourceProcessor,
                                                    4 * node->getID() + 3);
                            }
                        }
                    }
                }
            }
            for (auto entry : mesh.getPushElements())
            {
                int targetProcessor = entry.first;
                for (Base::Element *element : entry.second)
                {
                    //the id's of elements, faces edges and nodes are interleaved so all information can be communicated simultaneously without tag collisions
                    mpiInstance.send(offsets[unknown].elementOffsets_[element->getID()], targetProcessor,
                                     4 * element->getID());
                    for (auto face : element->getFacesList())
                    {
                        if (offsets[unknown].faceOffsets_.count(face->getID()) == 1)
                        {
                            mpiInstance.send(offsets[unknown].faceOffsets_[face->getID()], targetProcessor,
                                             4 * face->getID() + 1);
                        }
                    }
                    for (auto edge : element->getEdgesList())
                    {
                        if (offsets[unknown].edgeOffsets_.count(edge->getID()) == 1)
                        {
                            mpiInstance.send(offsets[unknown].edgeOffsets_[edge->getID()], targetProcessor,
                                             4 * edge->getID() + 2);
                        }
                    }
                    if (meshDimension > 1)
                    {
                        for (auto node : element->getNodesList())
                        {
                            if (offsets[unknown].nodeOffsets_.count(node->getID()) == 1)
                            {
                                mpiInstance.send(offsets[unknown].nodeOffsets_[node->getID()], targetProcessor,
                                                 4 * node->getID() + 3);
                            }
                        }
                    }
                }
            }
            mpiInstance.sync();
        }
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

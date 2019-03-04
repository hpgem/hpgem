/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2017, University of Twente
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "output.h"
#include "LinearAlgebra/SmallVector.h"
#include <fstream>
#include <set>
#include <map>

namespace Detail {
    template<std::size_t d>
    struct tag{};

    template<std::size_t dimension>
    void printOtherEntityCounts(std::ofstream& output, const Preprocessor::Mesh<dimension>& mesh, tag<0>){}

    template<std::size_t d, std::size_t dimension>
    void printOtherEntityCounts(std::ofstream& output, const Preprocessor::Mesh<dimension>& mesh, tag<d>) {
        output << " " << mesh.template getNumberOfEntities<d>();
        printOtherEntityCounts(output, mesh, tag<d - 1>{});
    }

    template<std::size_t dimension, typename indexType>
    void printOtherEntities(std::ofstream& output, const Preprocessor::Mesh<dimension>& mesh, Preprocessor::MeshData<indexType, dimension, dimension>& partitions, tag<0>){}

    template<std::size_t d, std::size_t dimension, typename indexType>
    void printOtherEntities(std::ofstream& output, const Preprocessor::Mesh<dimension>& mesh, Preprocessor::MeshData<indexType, dimension, dimension>& partitions, tag<d>) {
        for(auto entity : mesh.template getEntities<d>()) {
            std::set<std::size_t> localPartitions;
            output << entity.getNumberOfElements() << " ";
            for(std::size_t i = 0; i < entity.getNumberOfElements(); ++i) {
                auto element = entity.getElement(i);
                output << element.getGlobalIndex() << " " << entity.getLocalIndex(i) << " ";
                for(auto node : element.getNodesList()) {
                    for(auto neighbour : node.getElementsList()){
                        localPartitions.insert(partitions[neighbour]);
                    }
                }
            }
            output << localPartitions.size() << " ";
            for(auto partition : localPartitions) {
                output << partition << " ";
            }
            output << std::endl;
        }
        printOtherEntities(output, mesh, partitions, tag<d - 1>{});
    }
}

template<typename indexType, std::size_t dimension>
void Preprocessor::outputMesh(Mesh<dimension>& mesh, MeshData<indexType, dimension, dimension> partitions, std::size_t numberOfPartitions) {
    if(mesh.getNumberOfNodes() == 0) {
        logger(WARN, "outputting empty mesh");
    }
    std::ofstream output(outputFileName.getValue());
    output << std::hexfloat;
    output << "mesh 1" << std::endl;
    output << mesh.getNumberOfNodes() << " " << mesh.getNumberOfElements() << " " << dimension;
    printOtherEntityCounts(output, mesh, ::Detail::tag<dimension - 1>{});
    output << std::endl;
    std::size_t reservedSpace = std::max(std::log10(mesh.getNumberOfNodes()), 0.) + 2;
    std::string whiteSpace(reservedSpace, ' ');
    output << numberOfPartitions << " ";
    auto partitionInformation = output.tellp();
    //reserve space
    for(std::size_t i = 0; i < numberOfPartitions; ++i) {
        output << whiteSpace;
    }
    output << std::endl;
    std::vector<std::size_t> partitionData(numberOfPartitions, 0);
    MeshData<std::vector<std::size_t>, dimension, 0> coordinateIndices(&mesh);
    for(auto node : mesh.getNodes()) {
        std::set<std::size_t> nodePartitions;
        std::set<LinearAlgebra::SmallVector<dimension>> nodeCoordinates;
        for(std::size_t i = 0; i < node.getNumberOfElements(); ++i) {
            auto element = node.getElement(i);
            nodePartitions.insert(partitions[element]);
            //make sure the shadow elements have all their nodes reside in all required partitions
            auto list = element.getNodesList();
            for(auto otherNode : list) {
                auto otherList = otherNode.getElementsList();
                for(auto otherElement : otherList) {
                    if(partitions[otherElement] != partitions[element]) {
                        nodePartitions.insert(partitions[otherElement]);
                    }
                }
            }
            if(!nodeCoordinates.insert(element.getCoordinate(node.getLocalIndex(i))).second) {
                logger(DEBUG, "node % is linked to coordinate % multiple times", node.getGlobalIndex(), node.getElement(i).getCoordinate(node.getLocalIndex(i)));
            }
        }
        logger(DEBUG, "%", nodePartitions.size());
        coordinateIndices[node].resize(node.getNumberOfElements());
        for(std::size_t i = 0; i < node.getNumberOfElements(); ++i) {
            coordinateIndices[node][i] = std::distance(nodeCoordinates.begin(),
                                                       nodeCoordinates.find(node.getElement(i).getCoordinate(node.getLocalIndex(i))));
        }
        output << nodePartitions.size() << " ";
        for(auto partition : nodePartitions) {
            partitionData[partition]++;
            output << partition << " ";
        }
        output << std::endl << nodeCoordinates.size() << " ";
        for(auto coordinate : nodeCoordinates) {
            for(std::size_t i = 0; i < dimension; ++i) {
                output << coordinate[i] << " ";
            }
        }
        output << std::endl;
    }
    for(auto element : mesh.getElements()) {
        output << element.getNumberOfNodes() << " ";
        for(auto node : element.getNodesList()) {
            output << node.getGlobalIndex() << " ";
            output << coordinateIndices[node][node.getElementIndex(element)] << " ";

        }
        output << partitions[element] << " ";
        std::set<std::size_t> shadowPartitions;
        for(auto node : element.getNodesList()) {
            for(auto otherElement : node.getElementsList()) {
                if(partitions[otherElement] != partitions[element]) {
                    shadowPartitions.insert(partitions[otherElement]);
                }
            }
        }
        output << shadowPartitions.size() << " ";
        for(auto index : shadowPartitions) {
            output << index << " ";
        }
        output << std::endl;
    }
    printOtherEntities(output, mesh, partitions, ::Detail::tag<dimension - 1>{});
    output.seekp(partitionInformation);
    for(std::size_t i = 0; i < numberOfPartitions; ++i) {
        output << partitionData[i] << " ";
    }
    output.close();
}

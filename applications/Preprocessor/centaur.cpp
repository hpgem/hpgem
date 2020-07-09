/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2017, University of Twente
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

#include <algorithm>
#include <iostream>

#include "centaur.h"
#include "Logger.h"

using namespace hpgem;


namespace Preprocessor {

UnstructuredInputStream<std::istringstream> CentaurReader::readLine() {
    // todo figure out ifstream.get(streambuf)
    logger.assert_always(!!centaurFile,
                         "There is expected to be an extra line, but the file "
                         "is not ready for reading");
    std::uint32_t lineSize, referenceSize;
    centaurFile >> lineSize;
    std::string buffer(lineSize, '\0');
    // workaround: string::data() returns const char* before c++17
    centaurFile.read(&buffer.front(), lineSize);
    centaurFile >> referenceSize;
    if (HPGEM_LOGLEVEL >= Log::DEBUG) {
        // swap the comment if you would like to see doubles instead
        // std::cout << buffer;
        char toHex[] = "0123456789abcdef";
        for (unsigned char input : buffer) {
            std::cout << toHex[input / 16] << toHex[input % 16] << " ";
        }
        // doubleBuffer buffers doubles, it has nothing to do with double
        // buffering note that this might break if your hardware is picky about
        // the alignment of its data
        /*const double* doubleBuffer = reinterpret_cast<const
        double*>(buffer.data()); for(std::size_t i = 0; i < buffer.size() /
        sizeof(double); ++i) { std::cout << doubleBuffer[i] << " ";
        }*/
        std::cout << std::endl;
    }
    logger.assert_always(lineSize == referenceSize,
                         "read error in centaur file");
    return std::istringstream(buffer, std::ios_base::binary);
}

CentaurReader::CentaurReader(std::string filename) {
    centaurFile.open(filename, std::ios::binary);
    logger.assert_always(centaurFile.is_open(),
                         "Cannot open Centaur meshfile.");
    logger.assert_always(centaurFile.good(),
                         "Something is not so good about this mesh");
    auto currentLine = readLine();

    float version;
    currentLine >> version >> centaurFileType;
    logger(INFO, "This mesh is in Centaur version % format", version);

    nodeStart = centaurFile.tellg();
    numberOfNodes = skipGroup();  // nodes

    elementStart = centaurFile.tellg();
    numberOfElements = 0;
    if (centaurFileType > 1) numberOfElements += skipGroup();  // hexahedra
    numberOfElements += skipGroup();  // triangles/prisms
    if (centaurFileType > 1) numberOfElements += skipGroup();  // pyramids
    numberOfElements += skipGroup();  // quadrilaterals/tetrahedra
    logger(DEBUG, "done with skipping elements");
    if (centaurFileType < 0) skipGroup();  // redundant boundary nodes
    if (centaurFileType > 3 || centaurFileType < 0)
        skipGroup(2);
    else
        skipGroup(1);                      // boundary faces
    if (centaurFileType < 0) skipGroup();  // segments
    skipGroup(1, false);                   // face/segment to group connections
    skipGroup(2, false);                   // group to b.c. type connections
    logger(VERBOSE, "done with skipping information");
    readNodeConnections();
}

std::function<void(int&)> difference_generator(
    const std::vector<std::size_t>& input) {
    std::size_t i = 1;
    return [=](int& result) mutable {
        result = static_cast<long long>(input[i + 1]) -
                 static_cast<long long>(input[i]);
        i++;
    };
}

Range<std::vector<std::vector<double>>> CentaurReader::getNodeCoordinates() {
    centaurFile.seekg(nodeStart);
    std::size_t dimension = 3;
    if (centaurFileType < 0) dimension = 2;
    std::vector<double> coordinate(dimension);
    std::uint32_t nodesOnLine;
    auto currentLine = readLine();
    currentLine >> nodesOnLine;  // actually available nodes, but for old
                                 // centaur files that is the same
    if (centaurFileType > 4) currentLine >> nodesOnLine;
    logger(VERBOSE, "Processing % nodes per line", nodesOnLine);
    std::uint32_t remainderThisLine = nodesOnLine;
    std::size_t index = 1;
    currentLine = readLine();
    auto increment = [=, currentLine = std::move(currentLine)](
                         std::vector<std::vector<double>>& next) mutable {
        if (remainderThisLine == 0) {
            remainderThisLine = nodesOnLine;
            currentLine = readLine();
        }
        auto position = boundaryConnections.find(index);
        // skip over boundary nodes that are already processed
        while (position != boundaryConnections.end() &&
               position->second[0] != index) {
            position = boundaryConnections.find(++index);
            for (auto& value : coordinate) {
                currentLine >> value;
            }
            if (--remainderThisLine == 0) {
                remainderThisLine = nodesOnLine;
                currentLine = readLine();
            }
        }
        next.resize(1);
        for (auto& value : coordinate) {
            currentLine >> value;
        }
        next[0] = coordinate;
        remainderThisLine--;
        index++;
        if (position != boundaryConnections.end() &&
            position->second.size() > 1) {
            std::istringstream::pos_type currentFilePosition =
                centaurFile.tellg();
            auto oldRemainder = remainderThisLine;
            std::istringstream::pos_type currentLinePosition =
                currentLine.tellg();
            auto oldLine = std::move(currentLine);
            currentLine = std::istringstream(oldLine.str());
            currentLine.seekg(currentLinePosition);
            for (auto offset :
                 Range<int>{static_cast<int>(position->second[1]) -
                                static_cast<int>(position->second[0]),
                            difference_generator(position->second),
                            position->second.size() - 1}) {
                logger.assert_always(offset > 0, "Need to skip % positions",
                                     offset);
                while (static_cast<uint32_t>(offset) > remainderThisLine) {
                    currentLine = readLine();
                    offset -= remainderThisLine;
                    remainderThisLine = nodesOnLine;
                }
                while (offset-- > 0) {
                    --remainderThisLine;
                    for (auto& value : coordinate) {
                        currentLine >> value;
                    }
                }
                next.push_back(coordinate);
            }
            currentLine = std::move(oldLine);
            centaurFile.seekg(currentFilePosition);
            remainderThisLine = oldRemainder;
        }
        if (HPGEM_LOGLEVEL >= Log::VERBOSE) {
            for (auto outputCoordinate : next) {
                std::cout << "(" << outputCoordinate[0];
                for (std::size_t i = 1; i < outputCoordinate.size(); ++i) {
                    std::cout << ", " << outputCoordinate[i];
                }
                std::cout << ")" << std::endl;
            }
        }
        logger(VERBOSE, "");
    };
    std::vector<std::vector<double>> coordinates;
    increment(coordinates);
    return Range<std::vector<std::vector<double>>>{
        coordinates, std::move(increment), numberOfNodes};
}

Range<std::vector<std::size_t>> CentaurReader::getElements() {
    centaurFile.seekg(elementStart);
    std::size_t currentGroupRemainder = 0;
    std::size_t entitiesOnLine = 0;
    std::size_t remainderThisLine = 0;
    std::size_t groupsProcessed = 0;
    UnstructuredInputStream<std::istringstream> currentLine;
    auto increment = [=, currentLine = std::move(currentLine)](
                         std::vector<std::size_t>& next) mutable {
        while (currentGroupRemainder == 0) {
            if (centaurFileType < 2)
                groupsProcessed += 2;
            else
                groupsProcessed++;
            currentLine = readLine();
            currentLine >> currentGroupRemainder;
            entitiesOnLine = currentGroupRemainder;
            if (centaurFileType > 4) currentLine >> entitiesOnLine;
            remainderThisLine = entitiesOnLine;
            if (groupsProcessed == 1)
                next.resize(8);
            else if (groupsProcessed == 3)
                next.resize(5);
            else if (groupsProcessed == 4)
                next.resize(4);
            else if (centaurFileType < 0)
                next.resize(3);
            else
                next.resize(6);
            currentLine = readLine();
        }
        if (remainderThisLine == 0) {
            currentLine = readLine();
            remainderThisLine = entitiesOnLine;
        }
        for (auto& index : next) {
            uint32_t input;
            currentLine >> input;
            index = toHpgemNumbering[input];
        }
        remainderThisLine--;
        currentGroupRemainder--;
    };
    std::vector<std::size_t> first;
    increment(first);
    return {first, std::move(increment), numberOfElements};
}

std::uint32_t CentaurReader::skipGroup(std::size_t linesPerEntity,
                                       bool multiline) {
    std::uint32_t numberOfEntities;
    auto currentLine = readLine();
    currentLine >> numberOfEntities;
    logger(VERBOSE, "skipping % entities", numberOfEntities);
    std::uint32_t numberOfEntitiesPerLine = numberOfEntities;
    if (centaurFileType > 4 && multiline)
        currentLine >> numberOfEntitiesPerLine;
    logger(VERBOSE, "there are % entities per line", numberOfEntitiesPerLine);
    if (numberOfEntities == 0) readLine();
    for (std::size_t i = 0; i < numberOfEntities;
         i += numberOfEntitiesPerLine) {
        for (std::size_t j = 0; j < linesPerEntity; ++j) readLine();
    }
    return numberOfEntities;
}

void CentaurReader::readNodeConnections() {
    std::uint32_t numberOfPeriodicTransformations = 1;
    if (centaurFileType > 3) {
        auto currentLine = readLine();
        currentLine >> numberOfPeriodicTransformations;
    }
    for (std::size_t i = 0; i < numberOfPeriodicTransformations; ++i) {
        if (centaurFileType > 3) {
            readLine();
            readLine();  // rotation/scaling information
        }
        std::uint32_t numberOfPeriodicNodes;
        auto currentLine = readLine();
        currentLine >> numberOfPeriodicNodes;
        currentLine = readLine();
        for (std::size_t j = 0; j < numberOfPeriodicNodes; ++j) {
            std::uint32_t matchingNodes[2];
            currentLine >> matchingNodes[0] >> matchingNodes[1];
            std::vector<std::size_t>& first =
                boundaryConnections[matchingNodes[0]];
            std::vector<std::size_t>& second =
                boundaryConnections[matchingNodes[1]];
            std::vector<std::size_t> temp;
            if (matchingNodes[0] > matchingNodes[1])
                std::swap(matchingNodes[0], matchingNodes[1]);
            std::set_union(first.begin(), first.end(), matchingNodes,
                           matchingNodes + 2, std::back_inserter(temp));
            first.clear();
            std::set_union(temp.begin(), temp.end(), second.begin(),
                           second.end(), std::back_inserter(first));
            second = first;
            if (HPGEM_LOGLEVEL >= Log::DEBUG) {
                for (auto k : first) {
                    std::cout << k << " ";
                }
                std::cout << "(" << matchingNodes[0] << ", " << matchingNodes[1]
                          << ")" << std::endl;
            }
        }
    }
    std::size_t hpgemIndex = 0;
    constexpr std::size_t iMax = std::numeric_limits<std::size_t>::max();
    toHpgemNumbering.resize(numberOfNodes + 1, iMax);
    numberOfNodes = 0;
    // the indexing in hpgem has the nodes regrouped such that periodically
    // connected coordinates are together indices in centaur start at 1
    for (std::size_t i = 1; i < toHpgemNumbering.size(); ++i) {
        if (toHpgemNumbering[i] == iMax) {
            numberOfNodes++;
            auto position = boundaryConnections.find(i);
            if (position == boundaryConnections.end()) {
                logger(DEBUG, "% -> %", i, hpgemIndex);
                toHpgemNumbering[i] = hpgemIndex++;
            } else {
                for (auto nodeIndex : position->second) {
                    logger(DEBUG, "% -> %", nodeIndex, hpgemIndex);
                    toHpgemNumbering[nodeIndex] = hpgemIndex++;
                }
            }
        }
    }
}
}  // namespace Preprocessor

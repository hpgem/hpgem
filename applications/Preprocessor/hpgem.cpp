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

#include <array>
#include "hpgem.h"
#include "Logger.h"

namespace Preprocessor {

    //this is clearly not a bool so std::vector<Bool> will not behave as badly as std::vector<bool>
    struct Bool {
        bool data;
        operator bool&(){
            return data;
        }
        operator const bool&() const{
            return data;
        }
    };

    class StructuredReader : public PrivateReader {
    public:
        StructuredReader(std::ifstream&& hpgemFile) : hpgemFile(std::move(hpgemFile)) {
            std::size_t versionNumber;
            this->hpgemFile >> std::boolalpha;
            this->hpgemFile >> versionNumber;
            logger.assert_always(versionNumber == 1, "mesh file is too new for this version of hpGEM");
            std::size_t dimension;
            this->hpgemFile >> dimension;
            numberOfElements.resize(dimension);
            numberOfNodes.resize(dimension);
            lowerLeft.resize(dimension);
            size.resize(dimension);
            isPeriodic.resize(dimension);
            for(std::size_t i = 0; i < dimension; ++i) {
                this->hpgemFile >> numberOfElements[i];
            }
            for(std::size_t i = 0; i < dimension; ++i) {
                this->hpgemFile >> lowerLeft[i];
            }
            for(std::size_t i = 0; i < dimension; ++i) {
                this->hpgemFile >> size[i];
                size[i] -= lowerLeft[i];
            }
            std::string triangles;
            this->hpgemFile >> triangles;
            if(triangles == "triangular") {
                hasTriangles = true;
            } else {
                logger.assert_always(triangles == "rectangular", "mesh must be either triangular or rectangular, but it is %", triangles);
                hasTriangles = false;
            }
            for(std::size_t i = 0; i < dimension; ++i) {
                this->hpgemFile >> isPeriodic[i];
            }
            for(std::size_t i = 0; i < dimension; ++i) {
                numberOfNodes[i] = numberOfElements[i] + (isPeriodic[i]?0:1);
            }
        }

        Range<std::vector<std::vector<double>>> getNodeCoordinates() override {
            std::vector<std::size_t> loopIndices(numberOfNodes.size());
            std::size_t nodeIndex = -1;
            std::size_t cumulativeIndex = 0;
            std::size_t totalNumberOfNodes = 1;
            for(std::size_t i = 0; i <  loopIndices.size(); ++i) {
                totalNumberOfNodes *= numberOfNodes[i];
            }
            coordinateIDs.resize(totalNumberOfNodes);
            auto increment = [=](std::vector<std::vector<double>>& next) mutable {
                next.resize(1);
                next[0].resize(loopIndices.size());
                coordinateIDs[++nodeIndex].push_back(cumulativeIndex++);
                for(std::size_t i = 0; i < loopIndices.size(); ++i) {
                    next[0][i] = lowerLeft[i] + loopIndices[i] * size[i] / numberOfElements[i];
                }
                for(std::size_t i = 0; i < loopIndices.size(); ++i) {
                    if(isPeriodic[i] && loopIndices[i] == 0) {
                        std::size_t currentSize = next.size();
                        for(std::size_t j = 0; j < currentSize; ++j) {
                            std::vector<double> newCoordinate = next[j];
                            newCoordinate[i] += size[i];
                            next.push_back(newCoordinate);
                            coordinateIDs[nodeIndex].push_back(cumulativeIndex++);
                        }
                    }
                }
                //emulates nested loops over all dimensions, for an unknown amount of dimensions
                for(std::size_t i = 0; i < loopIndices.size(); ++i) {
                    if(++loopIndices[i] < numberOfNodes[i]) return;
                    loopIndices[i] = 0;
                }
            };
            std::vector<std::vector<double>> first;
            increment(first);
            return {first, std::move(increment), totalNumberOfNodes};
        }

        Range<std::vector<std::size_t>> getElements() override {
            std::vector<std::size_t> loopIndices(numberOfNodes.size());
            std::vector<std::size_t> multipliedNodeCounts = {1};
            std::vector<std::vector<std::size_t>> cornersOfTriangles;
            std::size_t triangleIndex = 0;
            std::vector<std::size_t> rotate = {2, 0, 3, 1, 6, 4, 7, 5};
            if(loopIndices.size() == 2) {
                cornersOfTriangles = {{0, 1, 2}, {1, 3, 2}};
            } else if(loopIndices.size() == 3) {
                cornersOfTriangles = {{1, 2, 4, 7}, {0, 1, 2, 4}, {1, 3, 2, 7}, {1, 4, 5, 7}, {2, 6, 4, 7}};
            } else {
                logger.assert_always(!hasTriangles, "I don't know how to generate triangles for this dimension");
            }

            for(std::size_t i = 1; i < loopIndices.size(); ++i) {
                multipliedNodeCounts.push_back(multipliedNodeCounts.back() * numberOfNodes[i - 1]);
            }
            auto increment = [=](std::vector<std::size_t>& next) mutable {
                next = {0};
                for(std::size_t i = 0; i < loopIndices.size(); ++i) {
                    std::size_t currentSize = next.size();
                    for(std::size_t j = 0; j < currentSize; ++j) {
                        next[j] += loopIndices[i] * multipliedNodeCounts[i];
                        next.push_back(next[j] + multipliedNodeCounts[i]);
                    }
                }
                for(std::size_t j = 0; j < next.size(); ++j) {
                    std::size_t subIndex = 0;
                    std::size_t numberOfPreviousPeriodicDimensions = 0;
                    for(std::size_t i = 0; i < loopIndices.size(); ++i) {
                        if(((j & (1 << i)) != 0) && (loopIndices[i] == numberOfNodes[i] - 1)) {
                            next[j] -= (numberOfNodes[i] * multipliedNodeCounts[i]);
                            subIndex += (1 << numberOfPreviousPeriodicDimensions++);
                        }
                        if(isPeriodic[i] && (loopIndices[i] + (j & (1 << i)) == 0)) numberOfPreviousPeriodicDimensions++;
                    }
                    next[j] = coordinateIDs[next[j]][subIndex];
                }
                if(hasTriangles) {
                    auto temp = next;
                    for(std::size_t i = 0; i < cornersOfTriangles[triangleIndex].size(); ++i) {
                        temp[i] = next[cornersOfTriangles[triangleIndex][i]];
                    }
                    next = std::move(temp);
                    next.resize(cornersOfTriangles[triangleIndex].size());
                    if(++triangleIndex == cornersOfTriangles.size()) triangleIndex = 0;
                }
                if(triangleIndex == 0) {
                    //emulates nested loops over all dimensions, for an unknown amount of dimensions
                    for(std::size_t i = 0; i < loopIndices.size(); ++i) {
                        for(std::size_t j = 0; j < cornersOfTriangles.size(); ++j) {
                            for(std::size_t k = 0; k < cornersOfTriangles[j].size(); ++k) {
                                cornersOfTriangles[j][k] = rotate[cornersOfTriangles[j][k]];
                            }
                        }
                        if(++loopIndices[i] < numberOfElements[i]) return;
                        for(std::size_t j = 0; j < cornersOfTriangles.size(); ++j) {
                            for(std::size_t k = 0; k < cornersOfTriangles[j].size(); ++k) {
                                if(loopIndices[i] % 2 == 1)
                                cornersOfTriangles[j][k] = rotate[cornersOfTriangles[j][k]];
                            }
                        }
                        loopIndices[i] = 0;
                    }
                }
            };
            std::vector<std::size_t> first;
            increment(first);
            std::size_t totalNumberOfElements = 1;
            for(std::size_t i = 0; i < loopIndices.size(); ++i) {
                totalNumberOfElements *= numberOfElements[i];
            }
            if(hasTriangles) totalNumberOfElements *= cornersOfTriangles.size();
            return {first, std::move(increment), totalNumberOfElements};
        }

        std::size_t getDimension() override {
            return isPeriodic.size();
        }

        std::size_t getTargetProcessorCount() override {
            return 1;
        }

        Range<std::size_t> getProcessorBindings() override {
            std::size_t numberOfElements = 1;
            for(auto count : this->numberOfElements) {
                numberOfElements *= count;
            }
            return {0UL, [](std::size_t& next){next = 0UL;}, numberOfElements};
        }

    private:
        std::ifstream hpgemFile;
        std::vector<std::size_t> numberOfElements;
        std::vector<std::size_t> numberOfNodes;
        std::vector<double> lowerLeft;
        std::vector<double> size;
        bool hasTriangles;
        std::vector<Bool> isPeriodic;
        std::vector<std::vector<std::size_t>> coordinateIDs;
    };

    class PreprocessorReader : public PrivateReader {
    public:
        PreprocessorReader(std::ifstream&& hpgemFile) : hpgemFile(std::move(hpgemFile)) {
            std::size_t versionNumber;
            this->hpgemFile >> versionNumber;
            logger.assert_always(versionNumber == 1, "mesh file is too new for this version of hpGEM");
            this->hpgemFile >> numberOfNodes >> numberOfElements >> dimension;
            coordinateIDs.resize(numberOfNodes);
            //we are not interested in information about faces, edges or entities of dimension 2
            this->hpgemFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            this->hpgemFile >> targetProcessorCount;
            //we are not interested in how the nodes are split across the partitions
            this->hpgemFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            nodesStart = this->hpgemFile.tellg();
            for(std::size_t i = 0; i < numberOfNodes; ++i) {
                this->hpgemFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                this->hpgemFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
            elementsStart = this->hpgemFile.tellg();
        }

        Range<std::vector<std::vector<double>>> getNodeCoordinates() override {
            hpgemFile.seekg(nodesStart);
            std::size_t currentNode = 0;
            std::size_t processedNodes = 0;
            auto increment = [=](std::vector<std::vector<double>>& next) mutable {
                hpgemFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                std::size_t numberOfCoordinates;
                hpgemFile >> numberOfCoordinates;
                next.resize(numberOfCoordinates, std::vector<double>(dimension));
                for(auto& coordinate : next) {
                    coordinateIDs[currentNode].push_back(processedNodes++);
                    for(std::size_t i = 0; i < dimension; ++i) {
                        coordinate[i] = readDouble(hpgemFile);
                    }
                }
                currentNode++;
                hpgemFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            };
            std::vector<std::vector<double>> start;
            increment(start);
            return {start, std::move(increment), numberOfNodes};
        }

        Range<std::vector<std::size_t>> getElements() override {
            hpgemFile.seekg(elementsStart);
            auto increment = [this](std::vector<std::size_t>& next) mutable {
                std::size_t numberOfNodes;
                std::size_t nodeID, localCoordinateID;
                hpgemFile >> numberOfNodes;
                next.resize(numberOfNodes);
                for(auto& coordinateID : next) {
                    hpgemFile >> nodeID >> localCoordinateID;
                    coordinateID = coordinateIDs[nodeID][localCoordinateID];
                }
                hpgemFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            };
            std::vector<std::size_t> start;
            increment(start);
            return {start, std::move(increment), numberOfElements};
        }

        std::size_t getDimension() override {
            return dimension;
        }

        std::size_t getTargetProcessorCount() override {
            return targetProcessorCount;
        }

        Range<std::size_t> getProcessorBindings() override {
            hpgemFile.seekg(elementsStart);
            auto increment = [this](std::size_t& next) mutable {
                std::size_t numberOfNodes, ignore;
                hpgemFile >> numberOfNodes;
                for(std::size_t i = 0; i < numberOfNodes; ++i) {
                    hpgemFile >> ignore >> ignore;
                }
                hpgemFile >> next;
                hpgemFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            };
            std::size_t start;
            increment(start);
            return {start, std::move(increment), numberOfElements};
        }

    private:

        double readDouble(std::istream& input) const
        {
            const int BUFFER_LENGTH = 50;
            char buffer[BUFFER_LENGTH];

            while (isspace(input.peek()))
                input.get();
            // Read the input.
            int c = input.get();
            int bufferIndex = 0;
            // bufferIndex+1 to keep place to add an \0
            while (c != -1 && !isspace(c) && bufferIndex+1 < BUFFER_LENGTH)
            {
                buffer[bufferIndex] = (char)c;
                bufferIndex++;
                c = input.get();
            }
            // The last read character was not part of the number, so put
            // it back. Of course at EOF (-1) we should not put it back.
            if (c != -1)
            {
                input.unget();
            }
            // Add the null-terminator to get a valid c-string
            buffer[bufferIndex] = '\0';
            double result = 0;
            if(bufferIndex > 2 && buffer[0] == '0' && (buffer[1] == 'x' || buffer[1] == 'X'))
            {
                // Resort to sscanf for hex float parsing, as it is (one of) the only
                // functions that can do this.
                sscanf(buffer, "%la", &result);
            }
            else
            {
                result = strtod(buffer, nullptr);
            }
            return result;
        }

        std::ifstream hpgemFile;
        std::ifstream::pos_type nodesStart;
        std::ifstream::pos_type elementsStart;
        std::size_t numberOfNodes;
        std::size_t numberOfElements;
        std::size_t dimension;
        std::size_t targetProcessorCount;
        std::vector<std::vector<std::size_t>> coordinateIDs;
    };

    HpgemReader::HpgemReader(std::string filename) {
        std::ifstream hpgemFile{filename};
        std::string fileType;
        logger.assert_always(hpgemFile.is_open(), "failed to open %", filename);
        logger.assert_always(hpgemFile.good(), "something is not so good about this file");
        hpgemFile >> fileType;
        if(fileType == "mesh") {
            impl = std::make_shared<PreprocessorReader>(std::move(hpgemFile));
        } else if(fileType == "structured") {
            impl = std::make_shared<StructuredReader>(std::move(hpgemFile));
        }
        logger.assert_always(impl != nullptr, "File type % not recognized", fileType);
    }
}

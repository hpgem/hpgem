/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2014, University of Twente
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

#include "MeshManipulator.h"

#include "Geometry/PhysicalGeometry.h"
#include "Geometry/ReferenceTriangle.h"
#include "FE/BaseBasisFunction.h"
#include "FE/BasisFunctionSet.h"
#include "FE/OrientedBasisFunctionSet.h"
#include "Edge.h"
#include "ConfigurationData.h"
#include "Element.h"
#include "Face.h"
#include "MeshMoverBase.h"
#include "Geometry/PointPhysical.h"
#include "Geometry/ReferenceSquare.h"
#include "Geometry/GlobalNamespaceGeometry.h"
#include "Geometry/PointReference.h"
#include "Geometry/Mappings/MappingReferenceToPhysical.h"
#include "ElementFactory.h"
#include "FaceFactory.h"
#include "L2Norm.h"
#include "Geometry/Jacobian.h"
#include "Geometry/ReferenceGeometry.h"
#include "FE/BasisFunctions1DH1ConformingLine.h"
#include "FE/BasisFunctions2DH1ConformingSquare.h"
#include "FE/BasisFunctions2DH1ConformingTriangle.h"
#include "FE/BasisFunctions2DNedelec.h"
#include "FE/BasisFunctions3DH1ConformingCube.h"
#include "FE/BasisFunctions3DH1ConformingPrism.h"
#include "FE/BasisFunctions3DH1ConformingPyramid.h"
#include "FE/BasisFunctions3DH1ConformingTetrahedron.h"
#include "FE/BasisFunctions3DNedelec.h"
#include "FE/BasisFunctions3DAinsworthCoyle.h"
#include "FE/BasisFunctionsMonomials.h"
#include "Logger.h"

#include <algorithm>
#include <cctype>
#include <iostream>
#include <unordered_set>
#include <array>
#include <vector>
#include <numeric>
#include <type_traits>
#include <typeinfo>

//(crude) fix for pre-c++14 limitations of std::hash: just use the std::hash of
// the underlying type
template <typename T>
class EnumHash {
    // cause a compile error if someone uses this for non-enum types
    using onlyForEnums =
        typename std::enable_if<std::is_enum<T>::value, T>::type;

   public:
    std::size_t operator()(const T &t) const {
        return std::hash<typename std::underlying_type<T>::type>()(
            static_cast<typename std::underlying_type<T>::type>(t));
    }
};
namespace hpgem {
namespace Base {
template <std::size_t DIM>
void MeshManipulator<DIM>::useMonomialBasisFunctions(std::size_t order) {
    FE::BasisFunctionSet *bFset1 = new FE::BasisFunctionSet(order);
    switch (DIM) {
        case 1:
            FE::assembleMonomialBasisFunctions1D(*bFset1, order);
            break;
        case 2:
            FE::assembleMonomialBasisFunctions2D(*bFset1, order);
            break;
        case 3:

            FE::assembleMonomialBasisFunctions3D(*bFset1, order);
            break;
        case 4:
            FE::assembleMonomialBasisFunctions4D(*bFset1, order);
            break;
        default:
            logger(ERROR, "No basisfunctions exist in this dimension");
    }
    collBasisFSet_.resize(1);
    collBasisFSet_[0] = std::shared_ptr<const FE::BasisFunctionSet>(bFset1);
    // Register them all for usage
    for (Element *element : getElementsList(IteratorType::GLOBAL)) {
        for (std::size_t unknown = 0; unknown < configData_->numberOfUnknowns_;
             ++unknown) {
            element->setDefaultBasisFunctionSet(0, unknown);
        }
    }
}

template <std::size_t DIM>
void MeshManipulator<DIM>::useDefaultDGBasisFunctions(std::size_t order) {
    collBasisFSet_.clear();
    std::unordered_map<Geometry::ReferenceGeometryType, std::size_t,
                       EnumHash<Geometry::ReferenceGeometryType>>
        shapeToIndex;
    for (Element *element : getElementsList(IteratorType::GLOBAL)) {
        try {
            element->setDefaultBasisFunctionSet(shapeToIndex.at(
                element->getReferenceGeometry()->getGeometryType()));
        } catch (std::out_of_range &) {
            switch (element->getReferenceGeometry()->getGeometryType()) {
                case Geometry::ReferenceGeometryType::LINE:
                    shapeToIndex[Geometry::ReferenceGeometryType::LINE] =
                        collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        FE::createDGBasisFunctionSet1DH1Line(order));
                    break;
                case Geometry::ReferenceGeometryType::SQUARE:
                    shapeToIndex[Geometry::ReferenceGeometryType::SQUARE] =
                        collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        FE::createDGBasisFunctionSet2DH1Square(order));
                    break;
                case Geometry::ReferenceGeometryType::TRIANGLE:
                    shapeToIndex[Geometry::ReferenceGeometryType::TRIANGLE] =
                        collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        FE::createDGBasisFunctionSet2DH1Triangle(order));
                    break;
                case Geometry::ReferenceGeometryType::CUBE:
                    shapeToIndex[Geometry::ReferenceGeometryType::CUBE] =
                        collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        FE::createDGBasisFunctionSet3DH1Cube(order));
                    break;
                case Geometry::ReferenceGeometryType::TETRAHEDRON:
                    shapeToIndex[Geometry::ReferenceGeometryType::TETRAHEDRON] =
                        collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        FE::createDGBasisFunctionSet3DH1Tetrahedron(order));
                    break;
                case Geometry::ReferenceGeometryType::TRIANGULARPRISM:
                    shapeToIndex
                        [Geometry::ReferenceGeometryType::TRIANGULARPRISM] =
                            collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        FE::createDGBasisFunctionSet3DH1ConformingPrism(order));
                    break;
                case Geometry::ReferenceGeometryType::PYRAMID:
                    shapeToIndex[Geometry::ReferenceGeometryType::PYRAMID] =
                        collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        FE::createDGBasisFunctionSet3DH1ConformingPyramid(
                            order));
                    break;
                case Geometry::ReferenceGeometryType::HYPERCUBE:
                    logger(ERROR,
                           "No well-conditioned basis functions have been "
                           "implemented for %s",
                           element->getReferenceGeometry()->getName());
                    break;
                case Geometry::ReferenceGeometryType::POINT:
                    logger(ERROR,
                           "A point is not a valid geometry for an Element!");
                    break;
                default:
                    logger(
                        ERROR,
                        "A new geometry has been implemented, please add it to "
                        "the cases in "
                        "MeshManipulator::useDefaultDGBasisFunctions and "
                        "MeshManipulator::useDefaultConformingBasisFunctions");
            }
            element->setDefaultBasisFunctionSet(shapeToIndex.at(
                element->getReferenceGeometry()->getGeometryType()));
        }
    }
}

template <std::size_t DIM>
void MeshManipulator<DIM>::useDefaultDGBasisFunctions(std::size_t order,
                                                      std::size_t unknown) {
    // collBasisFSet_.clear();
    logger.assert_debug(
        unknown > 0,
        "useDefaultDGBasisFunctions(std::size_t unknown) will not clear "
        "collBasisFSet_ of the default basis functions. Use "
        "useDefaultBasisFunctions() instead for the zeroth unknown.");
    std::unordered_map<Geometry::ReferenceGeometryType, std::size_t,
                       EnumHash<Geometry::ReferenceGeometryType>>
        shapeToIndex;
    for (Element *element : getElementsList(IteratorType::GLOBAL)) {
        try {
            element->setDefaultBasisFunctionSet(
                shapeToIndex.at(
                    element->getReferenceGeometry()->getGeometryType()),
                unknown);
        } catch (std::out_of_range &) {
            switch (element->getReferenceGeometry()->getGeometryType()) {
                case Geometry::ReferenceGeometryType::LINE:
                    shapeToIndex[Geometry::ReferenceGeometryType::LINE] =
                        collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        FE::createDGBasisFunctionSet1DH1Line(order));
                    break;
                case Geometry::ReferenceGeometryType::SQUARE:
                    shapeToIndex[Geometry::ReferenceGeometryType::SQUARE] =
                        collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        FE::createDGBasisFunctionSet2DH1Square(order));
                    break;
                case Geometry::ReferenceGeometryType::TRIANGLE:
                    shapeToIndex[Geometry::ReferenceGeometryType::TRIANGLE] =
                        collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        FE::createDGBasisFunctionSet2DH1Triangle(order));
                    break;
                case Geometry::ReferenceGeometryType::CUBE:
                    shapeToIndex[Geometry::ReferenceGeometryType::CUBE] =
                        collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        FE::createDGBasisFunctionSet3DH1Cube(order));
                    break;
                case Geometry::ReferenceGeometryType::TETRAHEDRON:
                    shapeToIndex[Geometry::ReferenceGeometryType::TETRAHEDRON] =
                        collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        FE::createDGBasisFunctionSet3DH1Tetrahedron(order));
                    break;
                case Geometry::ReferenceGeometryType::TRIANGULARPRISM:
                    shapeToIndex
                        [Geometry::ReferenceGeometryType::TRIANGULARPRISM] =
                            collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        FE::createDGBasisFunctionSet3DH1ConformingPrism(order));
                    break;
                case Geometry::ReferenceGeometryType::PYRAMID:
                    shapeToIndex[Geometry::ReferenceGeometryType::PYRAMID] =
                        collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        FE::createDGBasisFunctionSet3DH1ConformingPyramid(
                            order));
                    break;
                case Geometry::ReferenceGeometryType::HYPERCUBE:
                    logger(ERROR,
                           "No well-conditioned basis functions have been "
                           "implemented for %s",
                           element->getReferenceGeometry()->getName());
                    break;
                case Geometry::ReferenceGeometryType::POINT:
                    logger(ERROR,
                           "A point is not a valid geometry for an Element!");
                    break;
                default:
                    logger(
                        ERROR,
                        "A new geometry has been implemented, please add it to "
                        "the cases in "
                        "MeshManipulator::useDefaultDGBasisFunctions and "
                        "MeshManipulator::useDefaultConformingBasisFunctions");
            }
            element->setDefaultBasisFunctionSet(
                shapeToIndex.at(
                    element->getReferenceGeometry()->getGeometryType()),
                unknown);
        }
    }
}

template <std::size_t DIM>
void MeshManipulator<DIM>::useNedelecDGBasisFunctions(std::size_t order) {
    collBasisFSet_.clear();
    std::unordered_map<Geometry::ReferenceGeometryType, std::size_t,
                       EnumHash<Geometry::ReferenceGeometryType>>
        shapeToIndex;
    for (Element *element : getElementsList(IteratorType::GLOBAL)) {
        try {
            element->setDefaultBasisFunctionSet(shapeToIndex.at(
                element->getReferenceGeometry()->getGeometryType()));
        } catch (std::out_of_range &) {
            switch (element->getReferenceGeometry()->getGeometryType()) {
                case Geometry::ReferenceGeometryType::TRIANGLE:
                    shapeToIndex[Geometry::ReferenceGeometryType::TRIANGLE] =
                        collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        FE::createDGBasisFunctionSet2DNedelec(order));
                    break;
                case Geometry::ReferenceGeometryType::TETRAHEDRON:
                    shapeToIndex[Geometry::ReferenceGeometryType::TETRAHEDRON] =
                        collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        FE::createDGBasisFunctionSet3DNedelec(order));
                    break;
                default:
                    logger(ERROR,
                           "No Nedelec basis functions have been implemented "
                           "for %s",
                           element->getReferenceGeometry()->getName());
            }
            element->setDefaultBasisFunctionSet(shapeToIndex.at(
                element->getReferenceGeometry()->getGeometryType()));
        }
    }
}

template <std::size_t DIM>
void MeshManipulator<DIM>::useAinsworthCoyleDGBasisFunctions(
    std::size_t order) {
    collBasisFSet_.clear();
    std::unordered_map<Geometry::ReferenceGeometryType, std::size_t,
                       EnumHash<Geometry::ReferenceGeometryType>>
        shapeToIndex;
    for (Element *element : getElementsList(IteratorType::GLOBAL)) {
        try {
            element->setDefaultBasisFunctionSet(shapeToIndex.at(
                element->getReferenceGeometry()->getGeometryType()));
        } catch (std::out_of_range &) {
            switch (element->getReferenceGeometry()->getGeometryType()) {
                case Geometry::ReferenceGeometryType::TETRAHEDRON:
                    shapeToIndex[Geometry::ReferenceGeometryType::TETRAHEDRON] =
                        collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        FE::createDGBasisFunctionSet3DAinsworthCoyle(order));
                    break;
                default:
                    logger(ERROR,
                           "No Nedelec basis functions have been implemented "
                           "for %s",
                           element->getReferenceGeometry()->getName());
            }
            element->setDefaultBasisFunctionSet(shapeToIndex.at(
                element->getReferenceGeometry()->getGeometryType()));
        }
    }
}

template <std::size_t DIM>
void MeshManipulator<DIM>::useDefaultConformingBasisFunctions(
    std::size_t order) {
    logger.assert_debug(order > 0,
                        "Basis function may not have an empty union of "
                        "supporting elements. Use a DG basis function on a "
                        "single element non-periodic mesh instead");
    collBasisFSet_.clear();
    std::unordered_map<Geometry::ReferenceGeometryType, std::size_t,
                       EnumHash<Geometry::ReferenceGeometryType>>
        shapeToElementIndex;
    std::unordered_map<Geometry::ReferenceGeometryType, std::size_t,
                       EnumHash<Geometry::ReferenceGeometryType>>
        numberOfFaceSets;
    std::unordered_map<Geometry::ReferenceGeometryType, std::size_t,
                       EnumHash<Geometry::ReferenceGeometryType>>
        numberOfEdgeSets;
    std::unordered_map<Geometry::ReferenceGeometryType, std::size_t,
                       EnumHash<Geometry::ReferenceGeometryType>>
        numberOfNodeSets;
    logger.suppressWarnings([&]() {
        for (Element *element : getElementsList(IteratorType::GLOBAL)) {
            try {
                element->setDefaultBasisFunctionSet(shapeToElementIndex.at(
                    element->getReferenceGeometry()->getGeometryType()));
            }
            // there is more relevant code after the huge catch block
            catch (std::out_of_range &) {
                auto type = element->getReferenceGeometry()->getGeometryType();
                std::vector<const FE::BasisFunctionSet *> nodeSet;
                std::vector<const FE::OrientedBasisFunctionSet *> faceSet;
                std::vector<const FE::OrientedBasisFunctionSet *> edgeSet;
                switch (type) {
                    case Geometry::ReferenceGeometryType::LINE:
                        shapeToElementIndex[type] = collBasisFSet_.size();
                        collBasisFSet_.emplace_back(
                            FE::createInteriorBasisFunctionSet1DH1Line(order));
                        faceSet =
                            FE::createVertexBasisFunctionSet1DH1Line(order);
                        for (const FE::BasisFunctionSet *set : faceSet) {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfFaceSets[type] = collBasisFSet_.size() -
                                                 shapeToElementIndex[type] - 1;
                        numberOfEdgeSets[type] = 0;
                        numberOfNodeSets[type] = 0;
                        break;
                    case Geometry::ReferenceGeometryType::SQUARE:
                        shapeToElementIndex[type] = collBasisFSet_.size();
                        collBasisFSet_.emplace_back(
                            FE::createInteriorBasisFunctionSet2DH1Square(
                                order));
                        faceSet =
                            FE::createFaceBasisFunctionSet2DH1Square(order);
                        for (const FE::BasisFunctionSet *set : faceSet) {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfFaceSets[type] = collBasisFSet_.size() -
                                                 shapeToElementIndex[type] - 1;
                        numberOfEdgeSets[type] = 0;
                        nodeSet =
                            FE::createVertexBasisFunctionSet2DH1Square(order);
                        for (const FE::BasisFunctionSet *set : nodeSet) {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfNodeSets[type] =
                            collBasisFSet_.size() - shapeToElementIndex[type] -
                            numberOfFaceSets[type] - numberOfEdgeSets[type] - 1;
                        break;
                    case Geometry::ReferenceGeometryType::TRIANGLE:
                        shapeToElementIndex[type] = collBasisFSet_.size();
                        collBasisFSet_.emplace_back(
                            FE::createInteriorBasisFunctionSet2DH1Triangle(
                                order));
                        faceSet =
                            FE::createFaceBasisFunctionSet2DH1Triangle(order);
                        for (const FE::BasisFunctionSet *set : faceSet) {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfFaceSets[type] = collBasisFSet_.size() -
                                                 shapeToElementIndex[type] - 1;
                        numberOfEdgeSets[type] = 0;
                        nodeSet =
                            FE::createVertexBasisFunctionSet2DH1Triangle(order);
                        for (const FE::BasisFunctionSet *set : nodeSet) {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfNodeSets[type] =
                            collBasisFSet_.size() - shapeToElementIndex[type] -
                            numberOfFaceSets[type] - numberOfEdgeSets[type] - 1;
                        break;
                    case Geometry::ReferenceGeometryType::CUBE:
                        shapeToElementIndex[type] = collBasisFSet_.size();
                        collBasisFSet_.emplace_back(
                            FE::createInteriorBasisFunctionSet3DH1Cube(order));
                        faceSet = FE::createFaceBasisFunctionSet3DH1Cube(order);
                        for (const FE::BasisFunctionSet *set : faceSet) {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfFaceSets[type] = collBasisFSet_.size() -
                                                 shapeToElementIndex[type] - 1;
                        edgeSet = FE::createEdgeBasisFunctionSet3DH1Cube(order);
                        for (const FE::BasisFunctionSet *set : edgeSet) {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfEdgeSets[type] = collBasisFSet_.size() -
                                                 shapeToElementIndex[type] -
                                                 numberOfFaceSets[type] - 1;
                        nodeSet =
                            FE::createVertexBasisFunctionSet3DH1Cube(order);
                        for (const FE::BasisFunctionSet *set : nodeSet) {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfNodeSets[type] =
                            collBasisFSet_.size() - shapeToElementIndex[type] -
                            numberOfFaceSets[type] - numberOfEdgeSets[type] - 1;
                        break;
                    case Geometry::ReferenceGeometryType::TETRAHEDRON:
                        shapeToElementIndex[type] = collBasisFSet_.size();
                        collBasisFSet_.emplace_back(
                            FE::createInteriorBasisFunctionSet3DH1Tetrahedron(
                                order));
                        faceSet = FE::createFaceBasisFunctionSet3DH1Tetrahedron(
                            order);
                        for (const FE::BasisFunctionSet *set : faceSet) {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfFaceSets[type] = collBasisFSet_.size() -
                                                 shapeToElementIndex[type] - 1;
                        edgeSet = FE::createEdgeBasisFunctionSet3DH1Tetrahedron(
                            order);
                        for (const FE::BasisFunctionSet *set : edgeSet) {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfEdgeSets[type] = collBasisFSet_.size() -
                                                 shapeToElementIndex[type] -
                                                 numberOfFaceSets[type] - 1;
                        nodeSet =
                            FE::createVertexBasisFunctionSet3DH1Tetrahedron(
                                order);
                        for (const FE::BasisFunctionSet *set : nodeSet) {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfNodeSets[type] =
                            collBasisFSet_.size() - shapeToElementIndex[type] -
                            numberOfFaceSets[type] - numberOfEdgeSets[type] - 1;
                        break;
                    case Geometry::ReferenceGeometryType::TRIANGULARPRISM:
                        shapeToElementIndex[type] = collBasisFSet_.size();
                        collBasisFSet_.emplace_back(
                            FE::createInteriorBasisFunctionSet3DH1ConformingPrism(
                                order));
                        faceSet =
                            FE::createFaceBasisFunctionSet3DH1ConformingPrism(
                                order);
                        for (const FE::BasisFunctionSet *set : faceSet) {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfFaceSets[type] = collBasisFSet_.size() -
                                                 shapeToElementIndex[type] - 1;
                        edgeSet =
                            FE::createEdgeBasisFunctionSet3DH1ConformingPrism(
                                order);
                        for (const FE::BasisFunctionSet *set : edgeSet) {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfEdgeSets[type] = collBasisFSet_.size() -
                                                 shapeToElementIndex[type] -
                                                 numberOfFaceSets[type] - 1;
                        nodeSet =
                            FE::createVertexBasisFunctionSet3DH1ConformingPrism(
                                order);
                        for (const FE::BasisFunctionSet *set : nodeSet) {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfNodeSets[type] =
                            collBasisFSet_.size() - shapeToElementIndex[type] -
                            numberOfFaceSets[type] - numberOfEdgeSets[type] - 1;
                        break;
                    case Geometry::ReferenceGeometryType::PYRAMID:
                        shapeToElementIndex[type] = collBasisFSet_.size();
                        collBasisFSet_.emplace_back(
                            FE::createInteriorBasisFunctionSet3DH1ConformingPyramid(
                                order));
                        nodeSet = FE::
                            createVertexBasisFunctionSet3DH1ConformingPyramid(
                                order);
                        for (const FE::BasisFunctionSet *set : nodeSet) {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfFaceSets[type] = 0;
                        numberOfEdgeSets[type] = 0;
                        numberOfNodeSets[type] = collBasisFSet_.size() -
                                                 shapeToElementIndex[type] - 1;
                        break;
                    case Geometry::ReferenceGeometryType::HYPERCUBE:
                        logger(ERROR,
                               "No well-conditioned basis functions have been "
                               "implemented for %s",
                               element->getReferenceGeometry()->getName());
                        break;
                    case Geometry::ReferenceGeometryType::POINT:
                        logger(
                            ERROR,
                            "A point is not a valid geometry for an Element!");
                        break;
                    default:
                        logger(
                            ERROR,
                            "A new geometry has been implemented, please add "
                            "it to the cases in "
                            "MeshManipulator::useDefaultDGBasisFunctions and "
                            "MeshManipulator::"
                            "useDefaultConformingBasisFunctions");
                }
                element->setDefaultBasisFunctionSet(
                    shapeToElementIndex.at(type));
            }
        }
        for (Face *face : getFacesList(IteratorType::GLOBAL)) {
            std::size_t faceNumber = face->localFaceNumberLeft();
            auto type = face->getPtrElementLeft()
                            ->getReferenceGeometry()
                            ->getGeometryType();
            for (std::size_t i = shapeToElementIndex[type] + 1;
                 i < shapeToElementIndex[type] + numberOfFaceSets[type] + 1;
                 ++i) {
                logger.assert_debug(
                    typeid(*collBasisFSet_[i]) ==
                        typeid(const FE::OrientedBasisFunctionSet),
                    "This is not supposed to happen");
                if (static_cast<const FE::OrientedBasisFunctionSet *>(
                        collBasisFSet_[i].get())
                        ->checkOrientation(0, faceNumber)) {
                    face->getPtrElementLeft()->setFaceBasisFunctionSet(
                        i, faceNumber);
                    // the number of basis functions depends on the shape of the
                    // face, not on the shape of the element
                    face->setLocalNumberOfBasisFunctions(
                        collBasisFSet_[i]->size());
                }
            }
            if (face->isInternal()) {
                faceNumber = face->localFaceNumberRight();
                type = face->getPtrElementRight()
                           ->getReferenceGeometry()
                           ->getGeometryType();
                std::size_t orientation = face->getFaceToFaceMapIndex();
                for (std::size_t i = shapeToElementIndex[type] + 1;
                     i < shapeToElementIndex[type] + numberOfFaceSets[type] + 1;
                     ++i) {
                    logger.assert_debug(
                        typeid(*collBasisFSet_[i]) ==
                            typeid(const FE::OrientedBasisFunctionSet),
                        "This is not supposed to happen");
                    if (static_cast<const FE::OrientedBasisFunctionSet *>(
                            collBasisFSet_[i].get())
                            ->checkOrientation(orientation, faceNumber)) {
                        face->getPtrElementRight()->setFaceBasisFunctionSet(
                            i, faceNumber);
                    }
                }
            }
        }
        for (Edge *edge : getEdgesList(IteratorType::GLOBAL)) {
            for (std::size_t i = 0; i < edge->getNumberOfElements(); ++i) {
                Element *element = edge->getElement(i);
                std::size_t edgeNumber = edge->getEdgeNumber(i);
                std::size_t orientation = edge->getOrientation(i);
                auto type = element->getReferenceGeometry()->getGeometryType();
                for (std::size_t j =
                         shapeToElementIndex[type] + numberOfFaceSets[type] + 1;
                     j < shapeToElementIndex[type] + numberOfFaceSets[type] +
                             numberOfEdgeSets[type] + 1;
                     ++j) {
                    logger.assert_debug(
                        typeid(*collBasisFSet_[j]) ==
                            typeid(const FE::OrientedBasisFunctionSet),
                        "This is not supposed to happen");
                    if (static_cast<const FE::OrientedBasisFunctionSet *>(
                            collBasisFSet_[j].get())
                            ->checkOrientation(orientation, edgeNumber)) {
                        element->setEdgeBasisFunctionSet(j, edgeNumber);
                        edge->setLocalNumberOfBasisFunctions(
                            collBasisFSet_[j]->size());
                    }
                }
            }
        }
        for (Node *node : getNodesList(IteratorType::GLOBAL)) {
            if (DIM > 1) {
                for (std::size_t i = 0; i < node->getNumberOfElements(); ++i) {
                    Element *element = node->getElement(i);
                    std::size_t nodeNumber = node->getNodeNumber(i);
                    auto type =
                        element->getReferenceGeometry()->getGeometryType();
                    element->setVertexBasisFunctionSet(
                        shapeToElementIndex[type] + numberOfFaceSets[type] +
                            numberOfEdgeSets[type] + 1 + nodeNumber,
                        nodeNumber);
                    node->setLocalNumberOfBasisFunctions(
                        collBasisFSet_[shapeToElementIndex[type] +
                                       numberOfFaceSets[type] +
                                       numberOfEdgeSets[type] + 1 + nodeNumber]
                            ->size());
                }
            }
        }
    });
}

template <std::size_t DIM>
void MeshManipulator<DIM>::useDefaultConformingBasisFunctions(
    std::size_t order, std::size_t unknown) {
    logger.assert_debug(order > 0,
                        "Basis function may not have an empty union of "
                        "supporting elements. Use a DG basis function on a "
                        "single element non-periodic mesh instead");
    // collBasisFSet_.clear();
    std::unordered_map<Geometry::ReferenceGeometryType, std::size_t,
                       EnumHash<Geometry::ReferenceGeometryType>>
        shapeToElementIndex;
    std::unordered_map<Geometry::ReferenceGeometryType, std::size_t,
                       EnumHash<Geometry::ReferenceGeometryType>>
        numberOfFaceSets;
    std::unordered_map<Geometry::ReferenceGeometryType, std::size_t,
                       EnumHash<Geometry::ReferenceGeometryType>>
        numberOfEdgeSets;
    std::unordered_map<Geometry::ReferenceGeometryType, std::size_t,
                       EnumHash<Geometry::ReferenceGeometryType>>
        numberOfNodeSets;
    for (Element *element : getElementsList(IteratorType::GLOBAL)) {
        try {
            element->setDefaultBasisFunctionSet(
                shapeToElementIndex.at(
                    element->getReferenceGeometry()->getGeometryType()),
                unknown);
        }
        // there is more relevant code after the huge catch block
        catch (std::out_of_range &) {
            auto type = element->getReferenceGeometry()->getGeometryType();
            std::vector<const FE::BasisFunctionSet *> nodeSet;
            std::vector<const FE::OrientedBasisFunctionSet *> faceSet;
            std::vector<const FE::OrientedBasisFunctionSet *> edgeSet;
            switch (type) {
                case Geometry::ReferenceGeometryType::LINE:
                    shapeToElementIndex[type] = collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        FE::createInteriorBasisFunctionSet1DH1Line(order));
                    faceSet = FE::createVertexBasisFunctionSet1DH1Line(order);
                    for (const FE::BasisFunctionSet *set : faceSet) {
                        collBasisFSet_.emplace_back(set);
                    }
                    numberOfFaceSets[type] =
                        collBasisFSet_.size() - shapeToElementIndex[type] - 1;
                    numberOfEdgeSets[type] = 0;
                    numberOfNodeSets[type] = 0;
                    break;
                case Geometry::ReferenceGeometryType::SQUARE:
                    shapeToElementIndex[type] = collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        FE::createInteriorBasisFunctionSet2DH1Square(order));
                    faceSet = FE::createFaceBasisFunctionSet2DH1Square(order);
                    for (const FE::BasisFunctionSet *set : faceSet) {
                        collBasisFSet_.emplace_back(set);
                    }
                    numberOfFaceSets[type] =
                        collBasisFSet_.size() - shapeToElementIndex[type] - 1;
                    numberOfEdgeSets[type] = 0;
                    nodeSet = FE::createVertexBasisFunctionSet2DH1Square(order);
                    for (const FE::BasisFunctionSet *set : nodeSet) {
                        collBasisFSet_.emplace_back(set);
                    }
                    numberOfNodeSets[type] =
                        collBasisFSet_.size() - shapeToElementIndex[type] -
                        numberOfFaceSets[type] - numberOfEdgeSets[type] - 1;
                    break;
                case Geometry::ReferenceGeometryType::TRIANGLE:
                    shapeToElementIndex[type] = collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        FE::createInteriorBasisFunctionSet2DH1Triangle(order));
                    faceSet = FE::createFaceBasisFunctionSet2DH1Triangle(order);
                    for (const FE::BasisFunctionSet *set : faceSet) {
                        collBasisFSet_.emplace_back(set);
                    }
                    numberOfFaceSets[type] =
                        collBasisFSet_.size() - shapeToElementIndex[type] - 1;
                    numberOfEdgeSets[type] = 0;
                    nodeSet =
                        FE::createVertexBasisFunctionSet2DH1Triangle(order);
                    for (const FE::BasisFunctionSet *set : nodeSet) {
                        collBasisFSet_.emplace_back(set);
                    }
                    numberOfNodeSets[type] =
                        collBasisFSet_.size() - shapeToElementIndex[type] -
                        numberOfFaceSets[type] - numberOfEdgeSets[type] - 1;
                    break;
                case Geometry::ReferenceGeometryType::CUBE:
                    shapeToElementIndex[type] = collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        FE::createInteriorBasisFunctionSet3DH1Cube(order));
                    faceSet = FE::createFaceBasisFunctionSet3DH1Cube(order);
                    for (const FE::BasisFunctionSet *set : faceSet) {
                        collBasisFSet_.emplace_back(set);
                    }
                    numberOfFaceSets[type] =
                        collBasisFSet_.size() - shapeToElementIndex[type] - 1;
                    edgeSet = FE::createEdgeBasisFunctionSet3DH1Cube(order);
                    for (const FE::BasisFunctionSet *set : edgeSet) {
                        collBasisFSet_.emplace_back(set);
                    }
                    numberOfEdgeSets[type] = collBasisFSet_.size() -
                                             shapeToElementIndex[type] -
                                             numberOfFaceSets[type] - 1;
                    nodeSet = FE::createVertexBasisFunctionSet3DH1Cube(order);
                    for (const FE::BasisFunctionSet *set : nodeSet) {
                        collBasisFSet_.emplace_back(set);
                    }
                    numberOfNodeSets[type] =
                        collBasisFSet_.size() - shapeToElementIndex[type] -
                        numberOfFaceSets[type] - numberOfEdgeSets[type] - 1;
                    break;
                case Geometry::ReferenceGeometryType::TETRAHEDRON:
                    shapeToElementIndex[type] = collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        FE::createInteriorBasisFunctionSet3DH1Tetrahedron(
                            order));
                    faceSet =
                        FE::createFaceBasisFunctionSet3DH1Tetrahedron(order);
                    for (const FE::BasisFunctionSet *set : faceSet) {
                        collBasisFSet_.emplace_back(set);
                    }
                    numberOfFaceSets[type] =
                        collBasisFSet_.size() - shapeToElementIndex[type] - 1;
                    edgeSet =
                        FE::createEdgeBasisFunctionSet3DH1Tetrahedron(order);
                    for (const FE::BasisFunctionSet *set : edgeSet) {
                        collBasisFSet_.emplace_back(set);
                    }
                    numberOfEdgeSets[type] = collBasisFSet_.size() -
                                             shapeToElementIndex[type] -
                                             numberOfFaceSets[type] - 1;
                    nodeSet =
                        FE::createVertexBasisFunctionSet3DH1Tetrahedron(order);
                    for (const FE::BasisFunctionSet *set : nodeSet) {
                        collBasisFSet_.emplace_back(set);
                    }
                    numberOfNodeSets[type] =
                        collBasisFSet_.size() - shapeToElementIndex[type] -
                        numberOfFaceSets[type] - numberOfEdgeSets[type] - 1;
                    break;
                case Geometry::ReferenceGeometryType::TRIANGULARPRISM:
                    shapeToElementIndex[type] = collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        FE::createInteriorBasisFunctionSet3DH1ConformingPrism(
                            order));
                    faceSet = FE::createFaceBasisFunctionSet3DH1ConformingPrism(
                        order);
                    for (const FE::BasisFunctionSet *set : faceSet) {
                        collBasisFSet_.emplace_back(set);
                    }
                    numberOfFaceSets[type] =
                        collBasisFSet_.size() - shapeToElementIndex[type] - 1;
                    edgeSet = FE::createEdgeBasisFunctionSet3DH1ConformingPrism(
                        order);
                    for (const FE::BasisFunctionSet *set : edgeSet) {
                        collBasisFSet_.emplace_back(set);
                    }
                    numberOfEdgeSets[type] = collBasisFSet_.size() -
                                             shapeToElementIndex[type] -
                                             numberOfFaceSets[type] - 1;
                    nodeSet =
                        FE::createVertexBasisFunctionSet3DH1ConformingPrism(
                            order);
                    for (const FE::BasisFunctionSet *set : nodeSet) {
                        collBasisFSet_.emplace_back(set);
                    }
                    numberOfNodeSets[type] =
                        collBasisFSet_.size() - shapeToElementIndex[type] -
                        numberOfFaceSets[type] - numberOfEdgeSets[type] - 1;
                    break;
                case Geometry::ReferenceGeometryType::PYRAMID:
                    shapeToElementIndex[type] = collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        FE::createInteriorBasisFunctionSet3DH1ConformingPyramid(
                            order));
                    nodeSet =
                        FE::createVertexBasisFunctionSet3DH1ConformingPyramid(
                            order);
                    for (const FE::BasisFunctionSet *set : nodeSet) {
                        collBasisFSet_.emplace_back(set);
                    }
                    numberOfFaceSets[type] = 0;
                    numberOfEdgeSets[type] = 0;
                    numberOfNodeSets[type] =
                        collBasisFSet_.size() - shapeToElementIndex[type] - 1;
                    break;
                case Geometry::ReferenceGeometryType::HYPERCUBE:
                    logger(ERROR,
                           "No well-conditioned basis functions have been "
                           "implemented for %s",
                           element->getReferenceGeometry()->getName());
                    break;
                case Geometry::ReferenceGeometryType::POINT:
                    logger(ERROR,
                           "A point is not a valid geometry for an Element!");
                    break;
                default:
                    logger(
                        ERROR,
                        "A new geometry has been implemented, please add it to "
                        "the cases in "
                        "MeshManipulator::useDefaultDGBasisFunctions and "
                        "MeshManipulator::useDefaultConformingBasisFunctions");
            }
            element->setDefaultBasisFunctionSet(shapeToElementIndex.at(type),
                                                unknown);
            const_cast<ConfigurationData *>(configData_)
                ->numberOfBasisFunctions_ +=
                collBasisFSet_[shapeToElementIndex.at(type)]->size();
        }
    }
    for (Face *face : getFacesList(IteratorType::GLOBAL)) {
        std::size_t faceNumber = face->localFaceNumberLeft();
        auto type = face->getPtrElementLeft()
                        ->getReferenceGeometry()
                        ->getGeometryType();
        for (std::size_t i = shapeToElementIndex[type] + 1;
             i < shapeToElementIndex[type] + numberOfFaceSets[type] + 1; ++i) {
            logger.assert_debug(typeid(*collBasisFSet_[i]) ==
                                    typeid(const FE::OrientedBasisFunctionSet),
                                "This is not supposed to happen");
            if (static_cast<const FE::OrientedBasisFunctionSet *>(
                    collBasisFSet_[i].get())
                    ->checkOrientation(0, faceNumber)) {
                face->getPtrElementLeft()->setFaceBasisFunctionSet(
                    i, faceNumber, unknown);
                const_cast<ConfigurationData *>(configData_)
                    ->numberOfBasisFunctions_ += collBasisFSet_[i]->size();
                // the number of basis functions depends on the shape of the
                // face, not on the shape of the element
                face->setLocalNumberOfBasisFunctions(collBasisFSet_[i]->size(),
                                                     unknown);
            }
        }
        if (face->isInternal()) {
            faceNumber = face->localFaceNumberRight();
            type = face->getPtrElementRight()
                       ->getReferenceGeometry()
                       ->getGeometryType();
            std::size_t orientation = face->getFaceToFaceMapIndex();
            for (std::size_t i = shapeToElementIndex[type] + 1;
                 i < shapeToElementIndex[type] + numberOfFaceSets[type] + 1;
                 ++i) {
                logger.assert_debug(
                    typeid(*collBasisFSet_[i]) ==
                        typeid(const FE::OrientedBasisFunctionSet),
                    "This is not supposed to happen");
                if (static_cast<const FE::OrientedBasisFunctionSet *>(
                        collBasisFSet_[i].get())
                        ->checkOrientation(orientation, faceNumber)) {
                    face->getPtrElementRight()->setFaceBasisFunctionSet(
                        i, faceNumber, unknown);
                    const_cast<ConfigurationData *>(configData_)
                        ->numberOfBasisFunctions_ += collBasisFSet_[i]->size();
                }
            }
        }
    }
    for (Edge *edge : getEdgesList(IteratorType::GLOBAL)) {
        for (std::size_t i = 0; i < edge->getNumberOfElements(); ++i) {
            Element *element = edge->getElement(i);
            std::size_t edgeNumber = edge->getEdgeNumber(i);
            std::size_t orientation = edge->getOrientation(i);
            auto type = element->getReferenceGeometry()->getGeometryType();
            for (std::size_t j =
                     shapeToElementIndex[type] + numberOfFaceSets[type] + 1;
                 j < shapeToElementIndex[type] + numberOfFaceSets[type] +
                         numberOfEdgeSets[type] + 1;
                 ++j) {
                logger.assert_debug(
                    typeid(*collBasisFSet_[j]) ==
                        typeid(const FE::OrientedBasisFunctionSet),
                    "This is not supposed to happen");
                if (static_cast<const FE::OrientedBasisFunctionSet *>(
                        collBasisFSet_[j].get())
                        ->checkOrientation(orientation, edgeNumber)) {
                    element->setEdgeBasisFunctionSet(j, edgeNumber, unknown);
                    const_cast<ConfigurationData *>(configData_)
                        ->numberOfBasisFunctions_ += collBasisFSet_[j]->size();
                    edge->setLocalNumberOfBasisFunctions(
                        collBasisFSet_[j]->size(), unknown);
                }
            }
        }
    }
    for (Node *node : getNodesList(IteratorType::GLOBAL)) {
        if (DIM > 1) {
            for (std::size_t i = 0; i < node->getNumberOfElements(); ++i) {
                Element *element = node->getElement(i);
                std::size_t nodeNumber = node->getNodeNumber(i);
                auto type = element->getReferenceGeometry()->getGeometryType();
                element->setVertexBasisFunctionSet(
                    shapeToElementIndex[type] + numberOfFaceSets[type] +
                        numberOfEdgeSets[type] + 1 + nodeNumber,
                    nodeNumber, unknown);
                const_cast<ConfigurationData *>(configData_)
                    ->numberOfBasisFunctions_ +=
                    collBasisFSet_[shapeToElementIndex[type] +
                                   numberOfFaceSets[type] +
                                   numberOfEdgeSets[type] + 1 + nodeNumber]
                        ->size();
                node->setLocalNumberOfBasisFunctions(
                    collBasisFSet_[shapeToElementIndex[type] +
                                   numberOfFaceSets[type] +
                                   numberOfEdgeSets[type] + 1 + nodeNumber]
                        ->size(),
                    unknown);
            }
        }
    }
}

template <std::size_t DIM>
MeshManipulator<DIM>::MeshManipulator(const ConfigurationData *config,
                                      std::size_t numberOfElementMatrices,
                                      std::size_t numberOfElementVectors,
                                      std::size_t numberOfFaceMatrices,
                                      std::size_t numberOfFaceVectors)
    : MeshManipulatorBase(config, DIM, numberOfElementMatrices,
                          numberOfElementVectors, numberOfFaceMatrices,
                          numberOfFaceVectors),
      meshMover_(nullptr) {}

template <std::size_t DIM>
MeshManipulator<DIM>::MeshManipulator(const MeshManipulator &other)
    : MeshManipulatorBase(other),
      theMesh_(other.theMesh_),
      meshMover_(other.meshMover_),
      collBasisFSet_(other.collBasisFSet_) {}

template <std::size_t DIM>
MeshManipulator<DIM>::~MeshManipulator() {
    delete meshMover_;
}

template <std::size_t DIM>
void MeshManipulator<DIM>::setDefaultBasisFunctionSet(
    FE::BasisFunctionSet *bFSet) {
    logger.assert_debug(bFSet != nullptr, "Invalid basis function set passed");
    if (collBasisFSet_.size() == 0) collBasisFSet_.resize(1);
    collBasisFSet_[0] = std::shared_ptr<const FE::BasisFunctionSet>(bFSet);
    const_cast<ConfigurationData *>(configData_)->numberOfBasisFunctions_ =
        bFSet->size();
    for (Base::Face *face : getFacesList(IteratorType::GLOBAL)) {
        face->setLocalNumberOfBasisFunctions(0);
    }
    for (Base::Edge *edge : getEdgesList(IteratorType::GLOBAL)) {
        edge->setLocalNumberOfBasisFunctions(0);
    }
    for (Base::Node *node : getNodesList(IteratorType::GLOBAL)) {
        node->setLocalNumberOfBasisFunctions(0);
    }
    for (ElementIterator it = elementColBegin(IteratorType::GLOBAL);
         it != elementColEnd(IteratorType::GLOBAL); ++it) {
        (*it)->setDefaultBasisFunctionSet(0);
    }
}

template <std::size_t DIM>
void MeshManipulator<DIM>::addVertexBasisFunctionSet(
    const std::vector<const FE::BasisFunctionSet *> &bFsets) {
    std::size_t firstNewEntry = collBasisFSet_.size();
    for (const FE::BasisFunctionSet *set : bFsets) {
        logger.assert_debug(set != nullptr,
                            "Invalid basis function set detected");
        collBasisFSet_.emplace_back(set);
    }
    for (Node *node : getNodesList()) {
        for (std::size_t i = 0; i < node->getNumberOfElements(); ++i) {
            node->getElement(i)->setVertexBasisFunctionSet(
                firstNewEntry + node->getNodeNumber(i), node->getNodeNumber(i));
        }
        node->setLocalNumberOfBasisFunctions(bFsets[0]->size());
    }
    const_cast<ConfigurationData *>(configData_)->numberOfBasisFunctions_ +=
        (*elementColBegin())->getNumberOfNodes() * bFsets[0]->size();
}

template <std::size_t DIM>
void MeshManipulator<DIM>::addFaceBasisFunctionSet(
    const std::vector<const FE::OrientedBasisFunctionSet *> &bFsets) {
    std::size_t firstNewEntry = collBasisFSet_.size();
    for (const FE::BasisFunctionSet *set : bFsets) {
        logger.assert_debug(set != nullptr,
                            "Invalid basis function set detected");
        collBasisFSet_.emplace_back(set);
    }
    for (Face *face : getFacesList()) {
        std::size_t faceNumber = face->localFaceNumberLeft();
        for (std::size_t i = 0; i < bFsets.size(); ++i) {
            if (bFsets[i]->checkOrientation(0, faceNumber)) {
                face->getPtrElementLeft()->setFaceBasisFunctionSet(
                    firstNewEntry + i, faceNumber);
            }
        }
        if (face->isInternal()) {
            faceNumber = face->localFaceNumberRight();
            std::size_t orientation = face->getFaceToFaceMapIndex();
            for (std::size_t i = 0; i < bFsets.size(); ++i) {
                if (bFsets[i]->checkOrientation(orientation, faceNumber)) {
                    face->getPtrElementRight()->setFaceBasisFunctionSet(
                        firstNewEntry + i, faceNumber);
                }
            }
        }
        face->setLocalNumberOfBasisFunctions(bFsets[0]->size());
    }
    const_cast<ConfigurationData *>(configData_)->numberOfBasisFunctions_ +=
        (*elementColBegin())->getPhysicalGeometry()->getNumberOfFaces() *
        bFsets[0]->size();
}

template <std::size_t DIM>
void MeshManipulator<DIM>::addEdgeBasisFunctionSet(
    const std::vector<const FE::OrientedBasisFunctionSet *> &bFsets) {
    std::size_t firstNewEntry = collBasisFSet_.size();
    for (const FE::BasisFunctionSet *set : bFsets) {
        logger.assert_debug(set != nullptr,
                            "Invalid basis function set detected");
        collBasisFSet_.emplace_back(set);
    }
    logger(DEBUG, "In MeshManipulator::addEdgeBasisFunctionSet: ");
    for (Edge *edge : getEdgesList()) {
        for (std::size_t i = 0; i < edge->getNumberOfElements(); ++i) {
            for (std::size_t j = 0; j < bFsets.size(); ++j) {
                if (bFsets[j]->checkOrientation(edge->getOrientation(i),
                                                edge->getEdgeNumber(i))) {
                    edge->getElement(i)->setEdgeBasisFunctionSet(
                        firstNewEntry + j, edge->getEdgeNumber(i));
                    logger(DEBUG, "% % %", edge->getOrientation(i),
                           edge->getEdgeNumber(i), bFsets[j]->size());
                }
            }
        }
        edge->setLocalNumberOfBasisFunctions(bFsets[0]->size());
    }
    const_cast<ConfigurationData *>(configData_)->numberOfBasisFunctions_ +=
        (*elementColBegin())->getPhysicalGeometry()->getNumberOfFaces() *
        bFsets[0]->size();
}

template <std::size_t DIM>
Base::Element *MeshManipulator<DIM>::addElement(
    const std::vector<std::size_t> &globalNodeIndexes, std::size_t zoneId,
    std::size_t owner, bool owning) {
    logger.assert_debug(
        [&]() -> bool {
            for (std::size_t i = 0; i < globalNodeIndexes.size(); ++i)
                for (std::size_t j = 0; j < i; ++j)
                    if (globalNodeIndexes[i] == globalNodeIndexes[j])
                        return false;
            return true;
        }(),
        "Trying to pass the same node twice");
    logger.assert_debug(zoneId < getZones().size(), "Invalid zone ID");
    auto result = theMesh_.addElement(globalNodeIndexes, zoneId, owner, owning);
    return result;
}

template <std::size_t DIM>
void MeshManipulator<DIM>::move() {
    for (Geometry::PointPhysical<DIM> &p : theMesh_.getNodeCoordinates()) {
        if (meshMover_ != nullptr) {
            meshMover_->movePoint(p);
        }
    }
}

template <std::size_t DIM>
void MeshManipulator<DIM>::setMeshMover(const MeshMoverBase<DIM> *meshMover) {
    // can be set to nullptr if you don't want to move the mesh anymore
    meshMover_ = meshMover;
}

template <std::size_t DIM>
bool MeshManipulator<DIM>::addFace(Element *leftElementPtr,
                                   std::size_t leftElementLocalFaceNo,
                                   Element *rightElementPtr,
                                   std::size_t rightElementLocalFaceNo,
                                   const Geometry::FaceType &faceType) {
    logger.assert_debug(leftElementPtr != nullptr, "Invalid element passed");
    // rightElementPtr may be nullptr for boundary faces
    auto result =
        theMesh_.addFace(leftElementPtr, leftElementLocalFaceNo,
                         rightElementPtr, rightElementLocalFaceNo, faceType);
    if (result)
        logger.suppressWarnings([this]() {
            theMesh_.getSubmesh().add(
                *(--theMesh_.getFacesList(Base::IteratorType::GLOBAL).end()));
        });
    return result;
}

template <std::size_t DIM>
Edge *MeshManipulator<DIM>::addEdge() {
    auto edge = theMesh_.addEdge();
    theMesh_.getSubmesh().add(edge);
    return edge;
}

template <std::size_t DIM>
void MeshManipulator<DIM>::addNode() {
    theMesh_.addNode();
    logger.suppressWarnings([this]() {
        theMesh_.getSubmesh().add(
            *(--theMesh_.getNodesList(Base::IteratorType::GLOBAL).end()));
    });
}

template <std::size_t DIM>
std::tuple<const Base::Element *, Geometry::PointReference<DIM>>
    MeshManipulator<DIM>::physicalToReference(
        Geometry::PointPhysical<DIM> pointPhysical) const {
    try {
        // straightforward post-order traversal is slow since it doesn't benefit
        // from information on the coarse level
        return physicalToReference_detail(pointPhysical,
                                          getElementsList().getRootEntries());
    } catch (const char *message) {
        ///\todo MPI applications might be unable to use this routine
        /// effectively, should we wrap the message in a proper std::exception
        /// and retrow instead?
        logger(ERROR, message, pointPhysical);
    }
}

template <std::size_t DIM>
template <typename Iterable>
std::tuple<const Base::Element *, Geometry::PointReference<DIM>>
    MeshManipulator<DIM>::physicalToReference_detail(
        Geometry::PointPhysical<DIM> pointPhysical,
        Iterable elementContainer) const {
    for (auto singleEntry : elementContainer) {
        const Base::Element *element = singleEntry->getData();
        Geometry::PointReference<DIM> pointReference =
            element->physicalToReference(pointPhysical);
        if (element->getReferenceGeometry()->isInternalPoint(pointReference)) {
            if (singleEntry->hasChild()) {
                return physicalToReference_detail(pointPhysical,
                                                  singleEntry->getChildren());
            }
            return std::make_tuple(element, pointReference);
        }
    }
    throw "The point % lies outsize the subdomain managed by this thread";
}

template <std::size_t DIM>
void MeshManipulator<DIM>::addMeasurePoint(
    Geometry::PointPhysical<DIM> pointPhysical) {
    try {
        // bypass the API since it spawns an error, and the pointPhysical is
        // probably on another subdomain anyway
        measurePoints_.push_back(physicalToReference_detail(
            pointPhysical, getElementsList().getRootEntries()));
    } catch (const char *) {
        // silently ignore failure to find the point
        ///\todo should we error anyway if there is no MPI running?
    }
}

template <std::size_t DIM>
const std::vector<
    std::tuple<const Base::Element *, Geometry::PointReference<DIM>>>
    &MeshManipulator<DIM>::getMeasurePoints() const {
    return measurePoints_;
}

template <std::size_t DIM>
std::vector<std::tuple<const Base::Element *, Geometry::PointReference<DIM>>>
    &MeshManipulator<DIM>::getMeasurePoints() {
    return measurePoints_;
}
}  // namespace Base

template <std::size_t DIM>
std::ostream &operator<<(std::ostream &os,
                         const Base::MeshManipulator<DIM> &mesh) {
    for (Geometry::PointPhysical<DIM> p : mesh.getNodeCoordinates()) {
        os << "Node "
           << " " << p << std::endl;
    }

    for (Base::Element *element : mesh.getElementsList()) {
        os << "Element " << element->getID() << " " << element << std::endl;
    }
    return os;
}

///\todo refactor out smaller functions until each of them fits on a screen
namespace Base {
template <std::size_t DIM>
void MeshManipulator<DIM>::refine(
    const Geometry::RefinementMapping *refinementMapping,
    std::function<bool(const Element *)> shouldRefine) {
    if (!shouldRefine) {
        shouldRefine = [refinementMapping](const Element *element) {
            return element->getRefinementMap() == nullptr &&
                   refinementMapping->getBigElementReferenceGeometry() ==
                       element->getReferenceGeometry();
        };
    }
    ElementFactory::instance().setCollectionOfBasisFunctionSets(
        &collBasisFSet_);
    ElementFactory::instance().setNumberOfMatrices(numberOfElementMatrices_);
    ElementFactory::instance().setNumberOfVectors(numberOfElementVectors_);
    ElementFactory::instance().setNumberOfTimeLevels(
        configData_->numberOfTimeLevels_);
    ElementFactory::instance().setNumberOfUnknowns(
        configData_->numberOfUnknowns_);
    FaceFactory::instance().setNumberOfFaceMatrices(numberOfFaceMatrices_);
    FaceFactory::instance().setNumberOfFaceVectors(numberOfFaceVectors_);

    // do the refinement
    for (Base::Element *element : getElementsList()) {
        if (shouldRefine(element)) {
            logger.assert_debug(
                element->getRefinementMap() == nullptr,
                "Element was already refined, refine the children instead");
            logger.assert_debug(
                element->getReferenceGeometry() ==
                    refinementMapping->getBigElementReferenceGeometry(),
                "The mapping is intended for %s, but the element is a %",
                refinementMapping->getBigElementReferenceGeometry()->getName(),
                element->getReferenceGeometry()->getName());

            auto newReferenceNodes = refinementMapping->getNewNodeLocations(
                Geometry::PointReference<DIM>{});
            std::vector<std::size_t> globalPointIndices(
                refinementMapping->getNumberOfNewNodes() +
                element->getNumberOfNodes());
            std::vector<Node *> newNodes(
                refinementMapping->getNumberOfNewNodes() +
                    element->getNumberOfNodes(),
                nullptr);

            // insert exisiting data
            for (std::size_t i = 0; i < element->getNumberOfNodes(); ++i) {
                globalPointIndices[i] =
                    element->getPhysicalGeometry()->getNodeIndex(i);
                newNodes[i] = element->getNode(i);
            }

            for (std::size_t i = 0; i < element->getNumberOfFaces(); ++i) {
                const Face *face = element->getFace(i);
                // if the face was already refined in the other element and it
                // is also refined in this element
                if (face->getPositionInTree()->hasChild() &&
                    refinementMapping->getCodim1RefinementMaps()[i]
                            ->getNumberOfNewNodes() > 0) {
                    auto elementIndices =
                        refinementMapping->getCodim1LocalNodeIndices(i);
                    auto childIterator =
                        ++face->getPositionInTree()->getIterator(
                            Base::TreeTraversalMethod::PREORDER);
                    childIterator.setAllLevelTraversal();
                    for (std::size_t j = 0;
                         j < face->getPositionInTree()->getNumberOfChildren();
                         ++j) {
                        // optimization: doing ++childIterator after the last
                        // child of this face is done will incur a linear cost
                        // since the iterator has to check every leaf to see if
                        // there is one on a lower level
                        if (j > 0) ++childIterator;
                        Base::Face *subFace = *childIterator;
                        for (std::size_t k = 0;
                             k < subFace->getReferenceGeometry()
                                     ->getNumberOfNodes();
                             ++k) {
                            // make the existing subnodes available for the new
                            // element
                            std::size_t elementIndex = elementIndices
                                [refinementMapping->getCodim1RefinementMaps()[i]
                                     ->getSubElementLocalNodeIndices(j)[k]];
                            std::size_t indexOfFaceNode =
                                subFace->getPtrElementLeft()
                                    ->getReferenceGeometry()
                                    ->getLocalNodeIndexFromFaceAndIndexOnFace(
                                        subFace->localFaceNumberLeft(), k);
                            globalPointIndices[elementIndex] =
                                subFace->getPtrElementLeft()
                                    ->getPhysicalGeometry()
                                    ->getNodeIndex(indexOfFaceNode);
                            newNodes[elementIndex] =
                                subFace->getPtrElementLeft()->getNode(
                                    indexOfFaceNode);
                        }
                    }
                }
            }

            std::size_t offset = element->getNumberOfNodes();
            // generate new data
            for (std::size_t i = 0;
                 i < refinementMapping->getNumberOfNewNodes(); ++i) {
                if (newNodes[i + offset] == nullptr) {
                    addNode();
                    newNodes[i + offset] =
                        theMesh_.getNodesList(IteratorType::GLOBAL).back();
                    Geometry::PointPhysical<DIM> nodeCoordinate =
                        element->referenceToPhysical(newReferenceNodes[i]);
                    globalPointIndices[i + offset] =
                        theMesh_.getNumberOfNodeCoordinates();
                    theMesh_.addNodeCoordinate(nodeCoordinate);
                }
            }

            // create the sub-elements
            std::vector<Base::Element *> subElements;
            for (std::size_t i = 0;
                 i < refinementMapping->getNumberOfSubElements(); ++i) {
                std::vector<std::size_t> subElementNodeIndices;
                subElementNodeIndices.reserve(
                    refinementMapping->getSubElementReferenceGeometry(i)
                        ->getNumberOfNodes());
                for (std::size_t j :
                     refinementMapping->getSubElementLocalNodeIndices(i)) {
                    subElementNodeIndices.push_back(globalPointIndices[j]);
                }
                auto newElement = ElementFactory::instance().makeElement(
                    subElementNodeIndices, theMesh_.getNodeCoordinates(),
                    element->getZone(), element->getOwner(),
                    element->isOwnedByCurrentProcessor());
                subElements.push_back(newElement);
                for (std::size_t j = 0;
                     j < refinementMapping->getSubElementReferenceGeometry(i)
                             ->getNumberOfNodes();
                     ++j) {
                    newNodes[refinementMapping->getSubElementLocalNodeIndices(
                                 i)[j]]
                        ->addElement(newElement, j);
                }
            }
            theMesh_.addSubElements(element, subElements);

            // connect sub-faces
            for (std::size_t i = 0; i < element->getNumberOfFaces(); ++i) {
                const Face *face = element->getFace(i);
                if (face->getPositionInTree()->hasChild()) {
                    auto childIterator =
                        ++face->getPositionInTree()->getIterator(
                            Base::TreeTraversalMethod::PREORDER);
                    childIterator.setAllLevelTraversal();
                    for (std::size_t j = 0;
                         j < refinementMapping->getCodim1RefinementMaps()[i]
                                 ->getNumberOfSubElements();
                         ++j) {
                        // optimization: doing ++childIterator after the last
                        // child of this face is done will incur a linear cost
                        // since the iterator has to check every leaf to see if
                        // there is one on a lower level
                        if (j > 0) ++childIterator;
                        std::size_t elementIndex, faceIndex;
                        std::tie(elementIndex, faceIndex) =
                            refinementMapping->getSubElementAndLocalFaceIndex(
                                i, j);
                        (*childIterator)
                            ->addElement(subElements[elementIndex], faceIndex);
                    }
                } else {
                    std::vector<Base::Face *> subFaces;
                    for (std::size_t j = 0;
                         j < refinementMapping->getCodim1RefinementMaps()[i]
                                 ->getNumberOfSubElements();
                         ++j) {
                        std::size_t elementIndex, faceIndex;
                        std::tie(elementIndex, faceIndex) =
                            refinementMapping->getSubElementAndLocalFaceIndex(
                                i, j);
                        if (face->isInternal()) {
                            subFaces.push_back(FaceFactory::instance().makeFace(
                                subElements[elementIndex], faceIndex,
                                Geometry::FaceType::REFINEMENT_BOUNDARY));
                        } else {
                            subFaces.push_back(FaceFactory::instance().makeFace(
                                subElements[elementIndex], faceIndex,
                                face->getFaceType()));
                        }
                    }
                    theMesh_.addSubFaces(face, subFaces);
                }
            }
        }
    }
    faceFactory();
    edgeFactory();

    // move the measurepoints down the tree
    for (auto &pair : measurePoints_) {
        const Base::Element *element = std::get<0>(pair);
        if (element->getPositionInTree()->hasChild()) {
            pair = physicalToReference_detail(
                element->referenceToPhysical(std::get<1>(pair)),
                element->getPositionInTree()->getChildren());
        }
    }
}
namespace Detail {
inline void readSegmentHeader(std::istream &stream,
                              std::string expectedHeader) {
    std::string line;
    stream >> line;
    logger.assert_always(line == expectedHeader,
                         "Different segment header expected '%' but got '%'",
                         expectedHeader, line);
}
}  // namespace Detail

template <std::size_t DIM>
void MeshManipulator<DIM>::readMesh(const std::string &filename) {
    // set to correct value in case some other meshManipulator changed things
    ElementFactory::instance().setCollectionOfBasisFunctionSets(
        &collBasisFSet_);
    ElementFactory::instance().setNumberOfMatrices(numberOfElementMatrices_);
    ElementFactory::instance().setNumberOfVectors(numberOfElementVectors_);
    ElementFactory::instance().setNumberOfTimeLevels(
        configData_->numberOfTimeLevels_);
    ElementFactory::instance().setNumberOfUnknowns(
        configData_->numberOfUnknowns_);
    FaceFactory::instance().setNumberOfFaceMatrices(numberOfFaceMatrices_);
    FaceFactory::instance().setNumberOfFaceVectors(numberOfFaceVectors_);
    getElementsList().setSingleLevelTraversal(0);
    logger.suppressWarnings([this]() {
        getElementsList(IteratorType::GLOBAL).setSingleLevelTraversal(0);
    });

    using namespace std::string_literals;
    std::ifstream input;
    input.open(filename.c_str());
    logger.assert_always(input.is_open(), "Cannot open input file: %",
                         filename);
    logger.assert_always(input.good(),
                         "Something is not so good about this mesh");
    std::string rawInput;
    std::getline(input, rawInput);
    logger.assert_always(
        rawInput.substr(0, 4) == "mesh",
        "incorrect file type, please use the preprocessor first");
    // Extract the version
    std::size_t version;
    {
        std::istringstream temp(rawInput.substr(4));
        temp >> version;
        logger.assert_always(version >= 1 && version <= 2,
                             "Version % can't be read by this reader.",
                             version);
    }
    // Header //
    ////////////
    std::size_t numberOfNodes, numberOfElements, meshDimension;
    input >> numberOfNodes >> numberOfElements >> meshDimension;
    logger.assert_always(meshDimension == DIM,
                         "The mesh in this input file has the wrong dimension "
                         "(read %, but expected %)",
                         meshDimension, DIM);
    std::size_t numberOfFaces = 0;
    std::size_t numberOfEdges = 0;
    if (DIM > 1) input >> numberOfFaces;
    if (DIM > 2) input >> numberOfEdges;
    std::size_t numberOfPartitions, localNumberOfNodes;
    input >> numberOfPartitions;
    logger.assert_debug(MPIContainer::Instance().getProcessorID() >= 0,
                        "We got assigned processor ID % by MPI (expected >=0)",
                        MPIContainer::Instance().getProcessorID());
    std::size_t processorID = MPIContainer::Instance().getProcessorID();
    for (std::size_t i = 0; i < processorID + 1; ++i) {
        input >> localNumberOfNodes;
    }
    input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    logger.assert_debug(MPIContainer::Instance().getNumberOfProcessors() >= 0,
                        "MPI thinks we are running only on % processors",
                        MPIContainer::Instance().getNumberOfProcessors());
    logger.assert_always(
        numberOfPartitions ==
            static_cast<std::size_t>(
                MPIContainer::Instance().getNumberOfProcessors()),
        "This mesh is targeting % parallel threads, but you are running on % "
        "threads, please rerun the preprocessor first",
        numberOfPartitions, MPIContainer::Instance().getNumberOfProcessors());

    // ZONES //
    ///////////

    // For each zone, the corresponding zoneId in the mesh.
    std::vector<std::size_t> meshZoneIndices;
    if (version == 1) {
        // Zones are only present in version > 1, provide a default zone
        Zone &zone = addZone("Main");
        meshZoneIndices.push_back(zone.getZoneId());
    } else {
        Detail::readSegmentHeader(input, "zones");
        // Number of meshZoneIndices
        std::size_t numZones;
        input >> numZones;
        logger.assert_always(
            numZones > 0 || numberOfElements == 0,
            "Got a file with zero meshZoneIndices but with elements");
        input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        // Zone names
        for (std::size_t i = 0; i < numZones; ++i) {
            std::string zoneName;
            getline(input, zoneName);
            Zone &zone = addZone(zoneName);
            meshZoneIndices.push_back(zone.getZoneId());
        }
    }

    // NODES //
    ///////////

    if (version > 1) {
        Detail::readSegmentHeader(input, "nodes");
    }

    std::size_t processedNodes = 0;
    std::map<std::size_t, std::size_t> localNodeIndex;
    std::map<std::size_t, std::size_t> startOfCoordinates;

    // we need some leeway later on to synchronize face indices with node
    // indices later on
    if (DIM == 1) GlobalUniqueIndex::instance().getNodeIndex();

    for (std::size_t i = 0; i < numberOfNodes; ++i) {
        std::size_t nodePartitions;
        input >> nodePartitions;
        bool present = false;
        for (std::size_t j = 0; j < nodePartitions && !present; ++j) {
            std::size_t partition;
            input >> partition;
            // getProcessorID is already nonnegative according to preceding
            // check
            if (partition == processorID) {
                input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                startOfCoordinates[i] = getNumberOfNodeCoordinates();
                localNodeIndex[i] = processedNodes++;
                std::size_t numberOfCoordinates;
                input >> numberOfCoordinates;
                for (std::size_t k = 0; k < numberOfCoordinates; ++k) {
                    double nextValue;
                    LinearAlgebra::SmallVector<DIM> nextCoordinate;
                    for (std::size_t l = 0; l < DIM; ++l) {
                        nextValue = readDouble(input);
                        nextCoordinate[l] = nextValue;
                    }
                    getMesh().addNodeCoordinate(nextCoordinate);
                }
                present = true;
                addNode();
            }
        }
        if (!present) {
            input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            // reserve the index for use on another processor
            GlobalUniqueIndex::instance().getNodeIndex();
        }
    }
    logger.assert_always(
        processedNodes == localNumberOfNodes,
        "The mesh file lies about how many nodes are in this partition");

    // Elements //
    //////////////

    if (version > 1) {
        Detail::readSegmentHeader(input, "elements");
    }

    std::map<std::size_t, Element *> actualElement;

    for (std::size_t i = 0; i < numberOfElements; ++i) {
        // Read all the information about the element
        std::size_t nodesPerElement;
        input >> nodesPerElement;
        // Read node + coordinate indices
        std::vector<std::size_t> coordinateIndices;
        std::vector<std::size_t> nodeNumbers;
        for (std::size_t j = 0; j < nodesPerElement; ++j) {
            std::size_t globalIndex, coordinateOffset;
            input >> globalIndex >> coordinateOffset;
            coordinateIndices.push_back(startOfCoordinates[globalIndex] +
                                        coordinateOffset);
            nodeNumbers.push_back(localNodeIndex[globalIndex]);
        }
        // Owning partition (= MPI process id)
        std::size_t partition;
        input >> partition;
        // Shadow partitions (= MPI process ids where it is in the shadow layer)
        std::size_t numberOfShadowPartitions;
        input >> numberOfShadowPartitions;
        std::vector<std::size_t> shadowParitions(numberOfShadowPartitions);
        // Flag to keep track of whether it is in the shadow layer of this
        // processor
        bool inShadow = false;

        for (std::size_t j = 0; j < numberOfShadowPartitions; ++j) {
            input >> shadowParitions[j];
            inShadow |= processorID == shadowParitions[j];
        }

        std::size_t zoneIndex;
        if (version == 1) {
            // Version 1 did not have any zones, this uses the default
            zoneIndex = 0;
        } else {
            input >> zoneIndex;
        }
        logger.assert_always(zoneIndex < meshZoneIndices.size(),
                             "Zone index % out of bounds", zoneIndex);

        // Create the element if needed
        bool owning = partition == processorID;
        if (owning || inShadow) {
            Base::Element *element =
                addElement(coordinateIndices, meshZoneIndices[zoneIndex],
                           partition, owning);
            actualElement[i] = element;
            if (owning) {
                getMesh().getSubmesh().add(element);
                for (std::size_t &shadowPartition : shadowParitions) {
                    getMesh().getSubmesh().addPush(element, shadowPartition);
                }
            } else {
                getMesh().getSubmesh().addPull(element, partition);
            }
            for (std::size_t j = 0; j < nodeNumbers.size(); ++j) {
                getNodesList()[nodeNumbers[j]]->addElement(element, j);
            }
        } else {
            // reserve the index for use in another processor
            GlobalUniqueIndex::instance().getElementIndex();
        }
    }

    logger.suppressWarnings([&]() {
        getFacesList(IteratorType::GLOBAL).setPreOrderTraversal();
        getEdgesList(IteratorType::GLOBAL).setPreOrderTraversal();
    });

    // FACES //
    ///////////

    if (version > 1 && DIM > 1) {
        Detail::readSegmentHeader(input, "codim1");
    }

    for (std::size_t i = 0; i < numberOfFaces; ++i) {
        std::size_t localNumberOfFaces;
        input >> localNumberOfFaces;
        std::array<std::size_t, 2> globalElementIndices;
        std::array<std::size_t, 2> localFaceNumbers;
        logger.assert_always(0 < localNumberOfFaces && localNumberOfFaces < 3,
                             "The mesh file thinks there are % faces, but "
                             "there should be 1 or 2",
                             localNumberOfFaces);
        for (std::size_t j = 0; j < localNumberOfFaces; ++j) {
            std::size_t nextElementIndex, nextFaceNumber;
            input >> nextElementIndex >> nextFaceNumber;
            globalElementIndices[j] = nextElementIndex;
            localFaceNumbers[j] = nextFaceNumber;
        }
        bool faceIsInPartition = false;
        std::size_t numberOfPartitions;
        input >> numberOfPartitions;
        for (std::size_t j = 0; j < numberOfPartitions; ++j) {
            std::size_t partition;
            input >> partition;
            if (partition == processorID) {
                faceIsInPartition = true;
                if (localNumberOfFaces == 1) {
                    logger.assert_always(
                        actualElement[globalElementIndices[0]] != nullptr,
                        "local face is bounded by nonlocal element");
                    addFace(actualElement[globalElementIndices[0]],
                            localFaceNumbers[0], nullptr, 0,
                            Geometry::FaceType::WALL_BC);
                } else {
                    if (actualElement[globalElementIndices[0]] == nullptr) {
                        logger.assert_always(
                            actualElement[globalElementIndices[1]] != nullptr,
                            "local face is bounded by nonlocal element");
                        addFace(actualElement[globalElementIndices[1]],
                                localFaceNumbers[1], nullptr, 0,
                                Geometry::FaceType::PARTIAL_FACE);
                    } else if (actualElement[globalElementIndices[1]] ==
                               nullptr) {
                        logger.assert_always(
                            actualElement[globalElementIndices[0]] != nullptr,
                            "local face is bounded by nonlocal element");
                        addFace(actualElement[globalElementIndices[0]],
                                localFaceNumbers[0], nullptr, 0,
                                Geometry::FaceType::PARTIAL_FACE);
                    } else {
                        addFace(actualElement[globalElementIndices[0]],
                                localFaceNumbers[0],
                                actualElement[globalElementIndices[1]],
                                localFaceNumbers[1]);
                    }
                }
            }
        }
        if (!faceIsInPartition) {
            GlobalUniqueIndex::instance().getFaceIndex();
        }
    }

    // EDGES //
    ///////////

    if (version > 1 && DIM > 2) {
        Detail::readSegmentHeader(input, "codim2");
    }

    for (std::size_t i = 0; i < numberOfEdges; ++i) {
        std::size_t localNumberOfEdges;
        input >> localNumberOfEdges;
        std::vector<std::size_t> globalElementIndices(localNumberOfEdges);
        std::vector<std::size_t> localEdgeNumbers(localNumberOfEdges);
        for (std::size_t j = 0; j < localNumberOfEdges; ++j) {
            std::size_t nextElementIndex, nextEdgeNumber;
            input >> nextElementIndex >> nextEdgeNumber;
            globalElementIndices[j] = nextElementIndex;
            localEdgeNumbers[j] = nextEdgeNumber;
        }
        bool edgeIsInPartition = false;
        std::size_t numberOfPartitions;
        input >> numberOfPartitions;
        for (std::size_t j = 0; j < numberOfPartitions; ++j) {
            std::size_t partition;
            input >> partition;
            if (partition == processorID) {
                edgeIsInPartition = true;
                auto edge = addEdge();
                for (std::size_t k = 0; k < globalElementIndices.size(); ++k) {
                    if (actualElement[globalElementIndices[k]]) {
                        edge->addElement(actualElement[globalElementIndices[k]],
                                         localEdgeNumbers[k]);
                    }
                }
            }
        }
        if (!edgeIsInPartition) {
            GlobalUniqueIndex::instance().getEdgeIndex();
        }
    }

    std::size_t lastIndex = GlobalUniqueIndex::instance().getFaceIndex();

    if (DIM == 1) {
        logger.suppressWarnings([&]() {
            for (auto node : getNodesList(Base::IteratorType::GLOBAL)) {
                while (lastIndex++ < node->getID() - 1)
                    GlobalUniqueIndex::instance().getFaceIndex();
                if (node->getNumberOfElements() == 1) {
                    addFace(node->getElement(0), node->getNodeNumber(0),
                            nullptr, 0, Geometry::FaceType::WALL_BC);
                } else {
                    addFace(node->getElement(0), node->getNodeNumber(0),
                            node->getElement(1), node->getNodeNumber(1));
                }
            }
        });
    }

    for (auto pair : getPullElements()) {
        for (Base::Element *element : pair.second) {
            for (Base::Face *face : element->getFacesList()) {
                if (face->isInternal()) {
                    if (face->getFaceType() == Geometry::FaceType::INTERNAL) {
                        face->setFaceType(
                            Geometry::FaceType::SUBDOMAIN_BOUNDARY);
                    } else {
                        face->setFaceType(
                            Geometry::FaceType::PERIODIC_SUBDOMAIN_BC);
                    }
                }
            }
        }
    }
}

/// \bug does not do the bc flags yet
template <std::size_t DIM>
void MeshManipulator<DIM>::faceFactory() {
    getFacesList(IteratorType::GLOBAL).setPreOrderTraversal();
    std::vector<std::size_t> nodeIndices;
    std::vector<Element *> candidates;

    getElementsList(IteratorType::GLOBAL).setPreOrderTraversal();
    for (Element *element : theMesh_.getElementsList(IteratorType::GLOBAL)) {
        for (std::size_t i = 0; i < element->getNumberOfFaces(); ++i) {
            std::vector<const Node *> localNodes;
            // if this face is not there yet
            if (element->getFace(i) == nullptr) {
                localNodes.clear();
                candidates.clear();
                nodeIndices = element->getReferenceGeometry()
                                  ->getCodim1EntityLocalIndices(i);

                candidates = element->getNode(nodeIndices[0])->getElements();
                localNodes.push_back(element->getNode(nodeIndices[0]));
                std::sort(candidates.begin(), candidates.end(),
                          [](Element *left, Element *right) {
                              return left->getID() < right->getID();
                          });

                logger(DEBUG, "Candidates: ");
                for (Element *coutElement : candidates) {
                    logger(DEBUG, "Element %: %", coutElement->getID(),
                           *coutElement);
                }
                for (std::size_t j = 1; j < nodeIndices.size(); ++j) {
                    localNodes.push_back(element->getNode(nodeIndices[j]));
                    std::vector<Element *> temp, nextIndices;
                    nextIndices =
                        element->getNode(nodeIndices[j])->getElements();
                    std::sort(nextIndices.begin(), nextIndices.end(),
                              [](Element *left, Element *right) {
                                  return left->getID() < right->getID();
                              });
                    std::set_intersection(
                        candidates.begin(), candidates.end(),
                        nextIndices.begin(), nextIndices.end(),
                        std::back_inserter(temp),
                        [](Element *left, Element *right) {
                            return left->getID() < right->getID();
                        });
                    candidates = std::move(temp);
                    logger(DEBUG, "Candidates: ");
                    for (Element *coutElement : candidates) {
                        logger(DEBUG, "Element %: %", coutElement->getID(),
                               *coutElement);
                    }
                }

                // the current element does not bound the face or more than two
                // elements bound the face
                logger.assert_always(
                    candidates.size() == 1 || candidates.size() == 2,
                    "Detected % bounding elements for face %, which is "
                    "impossible",
                    candidates.size(),
                    theMesh_.getFacesList(IteratorType::GLOBAL).size() + 1);
                // boundary face
                if (candidates.size() == 1) {
                    logger.assert_debug(candidates[0] == element,
                                        "dropped the original element");
                    addFace(element, i, nullptr, 0,
                            Geometry::FaceType::WALL_BC);
                }
                if (candidates.size() == 2) {
                    Element *other;
                    if (candidates[0] == element) {
                        other = candidates[1];
                    } else {
                        other = candidates[0];
                    }
                    bool matchFound = false;
                    std::vector<std::size_t> otherNodeIndices;
                    for (std::size_t j = 0; j < other->getNumberOfFaces();
                         ++j) {
                        otherNodeIndices = other->getReferenceGeometry()
                                               ->getCodim1EntityLocalIndices(j);
                        bool match = true;
                        for (std::size_t k : otherNodeIndices) {
                            if (std::find(localNodes.begin(), localNodes.end(),
                                          other->getNode(k)) ==
                                localNodes.end()) {
                                match = false;
                            }
                        }
                        if (match) {
                            logger.assert_debug(
                                !matchFound,
                                "Found two opposing faces for face % "
                                " of element % in opposing element %.",
                                i, element->getID(), other->getID());
                            addFace(element, i, other, j);
                            matchFound = true;
                        }
                    }
                    logger.assert_debug(
                        matchFound,
                        "Could not find matching face for face % "
                        "of element % in opposing element %.",
                        i, element->getID(), other->getID());
                }
            }
        }
    }

    logger(VERBOSE, "Total number of Faces: %",
           getFacesList(IteratorType::GLOBAL).size());
}

// the algorithm for the edge factory is based on that of the face factory
// with some minor adaptation to account for the fact that there may be
// more than two elements per edge
///\bug does not do 4D yet
template <std::size_t DIM>
void MeshManipulator<DIM>::edgeFactory() {
    getEdgesList(IteratorType::GLOBAL).setPreOrderTraversal();
    //'edges' in DIM 2 are actually nodes
    if (DIM != 2) {
        std::vector<std::size_t> nodeList, otherNodeList;

        const Node *nodes[2];

        getElementsList(IteratorType::GLOBAL).setPreOrderTraversal();
        for (Element *element :
             theMesh_.getElementsList(IteratorType::GLOBAL)) {
            for (std::size_t i = 0; i < element->getNumberOfEdges(); ++i) {
                if (element->getEdge(i) == nullptr) {
                    nodeList = element->getReferenceGeometry()
                                   ->getCodim2EntityLocalIndices(i);
                    std::vector<Element *> candidates(0);
                    auto &leftElements =
                        element->getNode(nodeList[0])->getElements();
                    auto &rightElements =
                        element->getNode(nodeList[1])->getElements();
                    std::set_intersection(
                        leftElements.begin(), leftElements.end(),
                        rightElements.begin(), rightElements.end(),
                        std::back_inserter(candidates),
                        [](Element *a, Element *b) {
                            return a->getID() < b->getID();
                        });
                    logger.assert_debug(
                        candidates.size() > 0,
                        "current element is not adjacent to its own edges");
                    addEdge();
                    nodes[0] = element->getNode(nodeList[0]);
                    nodes[1] = element->getNode(nodeList[1]);
                    Edge *newEdge =
                        *(--theMesh_.getEdgesList(IteratorType::GLOBAL).end());
                    newEdge->addElement(element, i);
                    for (std::size_t j = 1; j < candidates.size(); ++j) {
                        Element *other = candidates[j];
                        for (std::size_t k = 0; k < other->getNumberOfEdges();
                             ++k) {
                            otherNodeList =
                                other->getReferenceGeometry()
                                    ->getCodim2EntityLocalIndices(k);
                            if ((other->getNode(otherNodeList[0]) == nodes[0] ||
                                 other->getNode(otherNodeList[0]) ==
                                     nodes[1]) &&
                                (other->getNode(otherNodeList[1]) == nodes[0] ||
                                 other->getNode(otherNodeList[1]) ==
                                     nodes[1])) {
                                newEdge->addElement(other, k);
                            }
                        }
                    }
                }
            }
        }
    }
}

template <std::size_t DIM>
Mesh<DIM> &MeshManipulator<DIM>::getMesh() {
    return theMesh_;
}

template <std::size_t DIM>
const Mesh<DIM> &MeshManipulator<DIM>::getMesh() const {
    return theMesh_;
}

template <std::size_t DIM>
double MeshManipulator<DIM>::readDouble(std::istream &input) const {
    const int BUFFER_LENGTH = 50;
    char buffer[BUFFER_LENGTH];

    while (isspace(input.peek())) input.get();
    // Read the input.
    int c = input.get();
    int bufferIndex = 0;
    // bufferIndex+1 to keep place to add an \0
    while (c != -1 && !isspace(c) && bufferIndex + 1 < BUFFER_LENGTH) {
        buffer[bufferIndex] = (char)c;
        bufferIndex++;
        c = input.get();
    }
    // The last read character was not part of the number, so put
    // it back. Of course at EOF (-1) we should not put it back.
    if (c != -1) {
        input.unget();
    }
    // Add the null-terminator to get a valid c-string
    buffer[bufferIndex] = '\0';
    double result = 0;
    if (bufferIndex > 2 && buffer[0] == '0' &&
        (buffer[1] == 'x' || buffer[1] == 'X')) {
        // Resort to sscanf for hex float parsing, as it is (one of) the only
        // functions that can do this.
        sscanf(buffer, "%la", &result);
    } else {
        result = strtod(buffer, nullptr);
    }
    return result;
}
}  // namespace Base
}  // namespace hpgem
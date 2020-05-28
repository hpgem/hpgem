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

#ifdef HPGEM_USE_QHULL
#include "libqhullcpp/QhullFacet.h"
#include "libqhullcpp/QhullQh.h"
#include "libqhullcpp/RboxPoints.h"
#include "libqhullcpp/QhullError.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullFacetSet.h"
#include "libqhullcpp/QhullLinkedList.h"
#include "libqhullcpp/QhullVertex.h"
#include "libqhullcpp/QhullVertexSet.h"
#include "libqhullcpp/Qhull.h"
#endif

#include "MeshManipulator.h"

#include "Geometry/PhysicalGeometry.h"
#include "Geometry/ReferenceTriangle.h"
#include "Edge.h"
#include "Base/BasisFunctionSet.h"
#include "ConfigurationData.h"
#include "Element.h"
#include "Face.h"
#include "MeshMoverBase.h"
#include "OrientedBasisFunctionSet.h"
#include "Geometry/PointPhysical.h"
#include "Geometry/ReferenceSquare.h"
#include "Geometry/GlobalNamespaceGeometry.h"
#include "Geometry/PointReference.h"
#include "BaseBasisFunction.h"
#include "Geometry/Mappings/MappingReferenceToPhysical.h"
#include "ElementFactory.h"
#include "FaceFactory.h"
#include "L2Norm.h"
#include "Geometry/Jacobian.h"
#include "Geometry/ReferenceGeometry.h"
#include "Utilities/BasisFunctions1DH1ConformingLine.h"
#include "Utilities/BasisFunctions2DH1ConformingSquare.h"
#include "Utilities/BasisFunctions2DH1ConformingTriangle.h"
#include "Utilities/BasisFunctions2DNedelec.h"
#include "Utilities/BasisFunctions3DH1ConformingCube.h"
#include "Utilities/BasisFunctions3DH1ConformingPrism.h"
#include "Utilities/BasisFunctions3DH1ConformingPyramid.h"
#include "Utilities/BasisFunctions3DH1ConformingTetrahedron.h"
#include "Utilities/BasisFunctions3DNedelec.h"
#include "Utilities/BasisFunctions3DAinsworthCoyle.h"
#include "Utilities/BasisFunctionsMonomials.h"
#include "Logger.h"

#include <algorithm>
#include <ctype.h>
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

namespace Base {
template <std::size_t DIM>
void MeshManipulator<DIM>::useMonomialBasisFunctions(std::size_t order) {
    Base::BasisFunctionSet *bFset1 = new Base::BasisFunctionSet(order);
    switch (DIM) {
        case 1:
            Utilities::assembleMonomialBasisFunctions1D(*bFset1, order);
            break;
        case 2:
            Utilities::assembleMonomialBasisFunctions2D(*bFset1, order);
            break;
        case 3:

            Utilities::assembleMonomialBasisFunctions3D(*bFset1, order);
            break;
        case 4:
            Utilities::assembleMonomialBasisFunctions4D(*bFset1, order);
            break;
        default:
            logger(ERROR, "No basisfunctions exist in this dimension");
    }
    collBasisFSet_.resize(1);
    collBasisFSet_[0] = std::shared_ptr<const BasisFunctionSet>(bFset1);
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
                        Utilities::createDGBasisFunctionSet1DH1Line(order));
                    break;
                case Geometry::ReferenceGeometryType::SQUARE:
                    shapeToIndex[Geometry::ReferenceGeometryType::SQUARE] =
                        collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        Utilities::createDGBasisFunctionSet2DH1Square(order));
                    break;
                case Geometry::ReferenceGeometryType::TRIANGLE:
                    shapeToIndex[Geometry::ReferenceGeometryType::TRIANGLE] =
                        collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        Utilities::createDGBasisFunctionSet2DH1Triangle(order));
                    break;
                case Geometry::ReferenceGeometryType::CUBE:
                    shapeToIndex[Geometry::ReferenceGeometryType::CUBE] =
                        collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        Utilities::createDGBasisFunctionSet3DH1Cube(order));
                    break;
                case Geometry::ReferenceGeometryType::TETRAHEDRON:
                    shapeToIndex[Geometry::ReferenceGeometryType::TETRAHEDRON] =
                        collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        Utilities::createDGBasisFunctionSet3DH1Tetrahedron(
                            order));
                    break;
                case Geometry::ReferenceGeometryType::TRIANGULARPRISM:
                    shapeToIndex
                        [Geometry::ReferenceGeometryType::TRIANGULARPRISM] =
                            collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        Utilities::createDGBasisFunctionSet3DH1ConformingPrism(
                            order));
                    break;
                case Geometry::ReferenceGeometryType::PYRAMID:
                    shapeToIndex[Geometry::ReferenceGeometryType::PYRAMID] =
                        collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        Utilities::
                            createDGBasisFunctionSet3DH1ConformingPyramid(
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
                        Utilities::createDGBasisFunctionSet1DH1Line(order));
                    break;
                case Geometry::ReferenceGeometryType::SQUARE:
                    shapeToIndex[Geometry::ReferenceGeometryType::SQUARE] =
                        collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        Utilities::createDGBasisFunctionSet2DH1Square(order));
                    break;
                case Geometry::ReferenceGeometryType::TRIANGLE:
                    shapeToIndex[Geometry::ReferenceGeometryType::TRIANGLE] =
                        collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        Utilities::createDGBasisFunctionSet2DH1Triangle(order));
                    break;
                case Geometry::ReferenceGeometryType::CUBE:
                    shapeToIndex[Geometry::ReferenceGeometryType::CUBE] =
                        collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        Utilities::createDGBasisFunctionSet3DH1Cube(order));
                    break;
                case Geometry::ReferenceGeometryType::TETRAHEDRON:
                    shapeToIndex[Geometry::ReferenceGeometryType::TETRAHEDRON] =
                        collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        Utilities::createDGBasisFunctionSet3DH1Tetrahedron(
                            order));
                    break;
                case Geometry::ReferenceGeometryType::TRIANGULARPRISM:
                    shapeToIndex
                        [Geometry::ReferenceGeometryType::TRIANGULARPRISM] =
                            collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        Utilities::createDGBasisFunctionSet3DH1ConformingPrism(
                            order));
                    break;
                case Geometry::ReferenceGeometryType::PYRAMID:
                    shapeToIndex[Geometry::ReferenceGeometryType::PYRAMID] =
                        collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        Utilities::
                            createDGBasisFunctionSet3DH1ConformingPyramid(
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
                        Utilities::createDGBasisFunctionSet2DNedelec(order));
                    break;
                case Geometry::ReferenceGeometryType::TETRAHEDRON:
                    shapeToIndex[Geometry::ReferenceGeometryType::TETRAHEDRON] =
                        collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        Utilities::createDGBasisFunctionSet3DNedelec(order));
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
                        Utilities::createDGBasisFunctionSet3DAinsworthCoyle(
                            order));
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
                std::vector<const Base::BasisFunctionSet *> nodeSet;
                std::vector<const Base::OrientedBasisFunctionSet *> faceSet;
                std::vector<const Base::OrientedBasisFunctionSet *> edgeSet;
                switch (type) {
                    case Geometry::ReferenceGeometryType::LINE:
                        shapeToElementIndex[type] = collBasisFSet_.size();
                        collBasisFSet_.emplace_back(
                            Utilities::createInteriorBasisFunctionSet1DH1Line(
                                order));
                        faceSet =
                            Utilities::createVertexBasisFunctionSet1DH1Line(
                                order);
                        for (const Base::BasisFunctionSet *set : faceSet) {
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
                            Utilities::createInteriorBasisFunctionSet2DH1Square(
                                order));
                        faceSet =
                            Utilities::createFaceBasisFunctionSet2DH1Square(
                                order);
                        for (const Base::BasisFunctionSet *set : faceSet) {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfFaceSets[type] = collBasisFSet_.size() -
                                                 shapeToElementIndex[type] - 1;
                        numberOfEdgeSets[type] = 0;
                        nodeSet =
                            Utilities::createVertexBasisFunctionSet2DH1Square(
                                order);
                        for (const Base::BasisFunctionSet *set : nodeSet) {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfNodeSets[type] =
                            collBasisFSet_.size() - shapeToElementIndex[type] -
                            numberOfFaceSets[type] - numberOfEdgeSets[type] - 1;
                        break;
                    case Geometry::ReferenceGeometryType::TRIANGLE:
                        shapeToElementIndex[type] = collBasisFSet_.size();
                        collBasisFSet_.emplace_back(
                            Utilities::
                                createInteriorBasisFunctionSet2DH1Triangle(
                                    order));
                        faceSet =
                            Utilities::createFaceBasisFunctionSet2DH1Triangle(
                                order);
                        for (const Base::BasisFunctionSet *set : faceSet) {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfFaceSets[type] = collBasisFSet_.size() -
                                                 shapeToElementIndex[type] - 1;
                        numberOfEdgeSets[type] = 0;
                        nodeSet =
                            Utilities::createVertexBasisFunctionSet2DH1Triangle(
                                order);
                        for (const Base::BasisFunctionSet *set : nodeSet) {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfNodeSets[type] =
                            collBasisFSet_.size() - shapeToElementIndex[type] -
                            numberOfFaceSets[type] - numberOfEdgeSets[type] - 1;
                        break;
                    case Geometry::ReferenceGeometryType::CUBE:
                        shapeToElementIndex[type] = collBasisFSet_.size();
                        collBasisFSet_.emplace_back(
                            Utilities::createInteriorBasisFunctionSet3DH1Cube(
                                order));
                        faceSet = Utilities::createFaceBasisFunctionSet3DH1Cube(
                            order);
                        for (const Base::BasisFunctionSet *set : faceSet) {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfFaceSets[type] = collBasisFSet_.size() -
                                                 shapeToElementIndex[type] - 1;
                        edgeSet = Utilities::createEdgeBasisFunctionSet3DH1Cube(
                            order);
                        for (const Base::BasisFunctionSet *set : edgeSet) {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfEdgeSets[type] = collBasisFSet_.size() -
                                                 shapeToElementIndex[type] -
                                                 numberOfFaceSets[type] - 1;
                        nodeSet =
                            Utilities::createVertexBasisFunctionSet3DH1Cube(
                                order);
                        for (const Base::BasisFunctionSet *set : nodeSet) {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfNodeSets[type] =
                            collBasisFSet_.size() - shapeToElementIndex[type] -
                            numberOfFaceSets[type] - numberOfEdgeSets[type] - 1;
                        break;
                    case Geometry::ReferenceGeometryType::TETRAHEDRON:
                        shapeToElementIndex[type] = collBasisFSet_.size();
                        collBasisFSet_.emplace_back(
                            Utilities::
                                createInteriorBasisFunctionSet3DH1Tetrahedron(
                                    order));
                        faceSet = Utilities::
                            createFaceBasisFunctionSet3DH1Tetrahedron(order);
                        for (const Base::BasisFunctionSet *set : faceSet) {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfFaceSets[type] = collBasisFSet_.size() -
                                                 shapeToElementIndex[type] - 1;
                        edgeSet = Utilities::
                            createEdgeBasisFunctionSet3DH1Tetrahedron(order);
                        for (const Base::BasisFunctionSet *set : edgeSet) {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfEdgeSets[type] = collBasisFSet_.size() -
                                                 shapeToElementIndex[type] -
                                                 numberOfFaceSets[type] - 1;
                        nodeSet = Utilities::
                            createVertexBasisFunctionSet3DH1Tetrahedron(order);
                        for (const Base::BasisFunctionSet *set : nodeSet) {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfNodeSets[type] =
                            collBasisFSet_.size() - shapeToElementIndex[type] -
                            numberOfFaceSets[type] - numberOfEdgeSets[type] - 1;
                        break;
                    case Geometry::ReferenceGeometryType::TRIANGULARPRISM:
                        shapeToElementIndex[type] = collBasisFSet_.size();
                        collBasisFSet_.emplace_back(
                            Utilities::
                                createInteriorBasisFunctionSet3DH1ConformingPrism(
                                    order));
                        faceSet = Utilities::
                            createFaceBasisFunctionSet3DH1ConformingPrism(
                                order);
                        for (const Base::BasisFunctionSet *set : faceSet) {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfFaceSets[type] = collBasisFSet_.size() -
                                                 shapeToElementIndex[type] - 1;
                        edgeSet = Utilities::
                            createEdgeBasisFunctionSet3DH1ConformingPrism(
                                order);
                        for (const Base::BasisFunctionSet *set : edgeSet) {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfEdgeSets[type] = collBasisFSet_.size() -
                                                 shapeToElementIndex[type] -
                                                 numberOfFaceSets[type] - 1;
                        nodeSet = Utilities::
                            createVertexBasisFunctionSet3DH1ConformingPrism(
                                order);
                        for (const Base::BasisFunctionSet *set : nodeSet) {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfNodeSets[type] =
                            collBasisFSet_.size() - shapeToElementIndex[type] -
                            numberOfFaceSets[type] - numberOfEdgeSets[type] - 1;
                        break;
                    case Geometry::ReferenceGeometryType::PYRAMID:
                        shapeToElementIndex[type] = collBasisFSet_.size();
                        collBasisFSet_.emplace_back(
                            Utilities::
                                createInteriorBasisFunctionSet3DH1ConformingPyramid(
                                    order));
                        nodeSet = Utilities::
                            createVertexBasisFunctionSet3DH1ConformingPyramid(
                                order);
                        for (const Base::BasisFunctionSet *set : nodeSet) {
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
                logger.assert_debug(typeid(*collBasisFSet_[i]) ==
                                        typeid(const OrientedBasisFunctionSet),
                                    "This is not supposed to happen");
                if (static_cast<const OrientedBasisFunctionSet *>(
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
                            typeid(const OrientedBasisFunctionSet),
                        "This is not supposed to happen");
                    if (static_cast<const OrientedBasisFunctionSet *>(
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
                            typeid(const OrientedBasisFunctionSet),
                        "This is not supposed to happen");
                    if (static_cast<const OrientedBasisFunctionSet *>(
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
            std::vector<const Base::BasisFunctionSet *> nodeSet;
            std::vector<const Base::OrientedBasisFunctionSet *> faceSet;
            std::vector<const Base::OrientedBasisFunctionSet *> edgeSet;
            switch (type) {
                case Geometry::ReferenceGeometryType::LINE:
                    shapeToElementIndex[type] = collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        Utilities::createInteriorBasisFunctionSet1DH1Line(
                            order));
                    faceSet =
                        Utilities::createVertexBasisFunctionSet1DH1Line(order);
                    for (const Base::BasisFunctionSet *set : faceSet) {
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
                        Utilities::createInteriorBasisFunctionSet2DH1Square(
                            order));
                    faceSet =
                        Utilities::createFaceBasisFunctionSet2DH1Square(order);
                    for (const Base::BasisFunctionSet *set : faceSet) {
                        collBasisFSet_.emplace_back(set);
                    }
                    numberOfFaceSets[type] =
                        collBasisFSet_.size() - shapeToElementIndex[type] - 1;
                    numberOfEdgeSets[type] = 0;
                    nodeSet = Utilities::createVertexBasisFunctionSet2DH1Square(
                        order);
                    for (const Base::BasisFunctionSet *set : nodeSet) {
                        collBasisFSet_.emplace_back(set);
                    }
                    numberOfNodeSets[type] =
                        collBasisFSet_.size() - shapeToElementIndex[type] -
                        numberOfFaceSets[type] - numberOfEdgeSets[type] - 1;
                    break;
                case Geometry::ReferenceGeometryType::TRIANGLE:
                    shapeToElementIndex[type] = collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        Utilities::createInteriorBasisFunctionSet2DH1Triangle(
                            order));
                    faceSet = Utilities::createFaceBasisFunctionSet2DH1Triangle(
                        order);
                    for (const Base::BasisFunctionSet *set : faceSet) {
                        collBasisFSet_.emplace_back(set);
                    }
                    numberOfFaceSets[type] =
                        collBasisFSet_.size() - shapeToElementIndex[type] - 1;
                    numberOfEdgeSets[type] = 0;
                    nodeSet =
                        Utilities::createVertexBasisFunctionSet2DH1Triangle(
                            order);
                    for (const Base::BasisFunctionSet *set : nodeSet) {
                        collBasisFSet_.emplace_back(set);
                    }
                    numberOfNodeSets[type] =
                        collBasisFSet_.size() - shapeToElementIndex[type] -
                        numberOfFaceSets[type] - numberOfEdgeSets[type] - 1;
                    break;
                case Geometry::ReferenceGeometryType::CUBE:
                    shapeToElementIndex[type] = collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        Utilities::createInteriorBasisFunctionSet3DH1Cube(
                            order));
                    faceSet =
                        Utilities::createFaceBasisFunctionSet3DH1Cube(order);
                    for (const Base::BasisFunctionSet *set : faceSet) {
                        collBasisFSet_.emplace_back(set);
                    }
                    numberOfFaceSets[type] =
                        collBasisFSet_.size() - shapeToElementIndex[type] - 1;
                    edgeSet =
                        Utilities::createEdgeBasisFunctionSet3DH1Cube(order);
                    for (const Base::BasisFunctionSet *set : edgeSet) {
                        collBasisFSet_.emplace_back(set);
                    }
                    numberOfEdgeSets[type] = collBasisFSet_.size() -
                                             shapeToElementIndex[type] -
                                             numberOfFaceSets[type] - 1;
                    nodeSet =
                        Utilities::createVertexBasisFunctionSet3DH1Cube(order);
                    for (const Base::BasisFunctionSet *set : nodeSet) {
                        collBasisFSet_.emplace_back(set);
                    }
                    numberOfNodeSets[type] =
                        collBasisFSet_.size() - shapeToElementIndex[type] -
                        numberOfFaceSets[type] - numberOfEdgeSets[type] - 1;
                    break;
                case Geometry::ReferenceGeometryType::TETRAHEDRON:
                    shapeToElementIndex[type] = collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        Utilities::
                            createInteriorBasisFunctionSet3DH1Tetrahedron(
                                order));
                    faceSet =
                        Utilities::createFaceBasisFunctionSet3DH1Tetrahedron(
                            order);
                    for (const Base::BasisFunctionSet *set : faceSet) {
                        collBasisFSet_.emplace_back(set);
                    }
                    numberOfFaceSets[type] =
                        collBasisFSet_.size() - shapeToElementIndex[type] - 1;
                    edgeSet =
                        Utilities::createEdgeBasisFunctionSet3DH1Tetrahedron(
                            order);
                    for (const Base::BasisFunctionSet *set : edgeSet) {
                        collBasisFSet_.emplace_back(set);
                    }
                    numberOfEdgeSets[type] = collBasisFSet_.size() -
                                             shapeToElementIndex[type] -
                                             numberOfFaceSets[type] - 1;
                    nodeSet =
                        Utilities::createVertexBasisFunctionSet3DH1Tetrahedron(
                            order);
                    for (const Base::BasisFunctionSet *set : nodeSet) {
                        collBasisFSet_.emplace_back(set);
                    }
                    numberOfNodeSets[type] =
                        collBasisFSet_.size() - shapeToElementIndex[type] -
                        numberOfFaceSets[type] - numberOfEdgeSets[type] - 1;
                    break;
                case Geometry::ReferenceGeometryType::TRIANGULARPRISM:
                    shapeToElementIndex[type] = collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        Utilities::
                            createInteriorBasisFunctionSet3DH1ConformingPrism(
                                order));
                    faceSet = Utilities::
                        createFaceBasisFunctionSet3DH1ConformingPrism(order);
                    for (const Base::BasisFunctionSet *set : faceSet) {
                        collBasisFSet_.emplace_back(set);
                    }
                    numberOfFaceSets[type] =
                        collBasisFSet_.size() - shapeToElementIndex[type] - 1;
                    edgeSet = Utilities::
                        createEdgeBasisFunctionSet3DH1ConformingPrism(order);
                    for (const Base::BasisFunctionSet *set : edgeSet) {
                        collBasisFSet_.emplace_back(set);
                    }
                    numberOfEdgeSets[type] = collBasisFSet_.size() -
                                             shapeToElementIndex[type] -
                                             numberOfFaceSets[type] - 1;
                    nodeSet = Utilities::
                        createVertexBasisFunctionSet3DH1ConformingPrism(order);
                    for (const Base::BasisFunctionSet *set : nodeSet) {
                        collBasisFSet_.emplace_back(set);
                    }
                    numberOfNodeSets[type] =
                        collBasisFSet_.size() - shapeToElementIndex[type] -
                        numberOfFaceSets[type] - numberOfEdgeSets[type] - 1;
                    break;
                case Geometry::ReferenceGeometryType::PYRAMID:
                    shapeToElementIndex[type] = collBasisFSet_.size();
                    collBasisFSet_.emplace_back(
                        Utilities::
                            createInteriorBasisFunctionSet3DH1ConformingPyramid(
                                order));
                    nodeSet = Utilities::
                        createVertexBasisFunctionSet3DH1ConformingPyramid(
                            order);
                    for (const Base::BasisFunctionSet *set : nodeSet) {
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
                                    typeid(const OrientedBasisFunctionSet),
                                "This is not supposed to happen");
            if (static_cast<const OrientedBasisFunctionSet *>(
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
                logger.assert_debug(typeid(*collBasisFSet_[i]) ==
                                        typeid(const OrientedBasisFunctionSet),
                                    "This is not supposed to happen");
                if (static_cast<const OrientedBasisFunctionSet *>(
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
                logger.assert_debug(typeid(*collBasisFSet_[j]) ==
                                        typeid(const OrientedBasisFunctionSet),
                                    "This is not supposed to happen");
                if (static_cast<const OrientedBasisFunctionSet *>(
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
void MeshManipulator<DIM>::setDefaultBasisFunctionSet(BasisFunctionSet *bFSet) {
    logger.assert_debug(bFSet != nullptr, "Invalid basis function set passed");
    if (collBasisFSet_.size() == 0) collBasisFSet_.resize(1);
    collBasisFSet_[0] = std::shared_ptr<const BasisFunctionSet>(bFSet);
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
    const std::vector<const BasisFunctionSet *> &bFsets) {
    std::size_t firstNewEntry = collBasisFSet_.size();
    for (const BasisFunctionSet *set : bFsets) {
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
    const std::vector<const OrientedBasisFunctionSet *> &bFsets) {
    std::size_t firstNewEntry = collBasisFSet_.size();
    for (const BasisFunctionSet *set : bFsets) {
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
    const std::vector<const OrientedBasisFunctionSet *> &bFsets) {
    std::size_t firstNewEntry = collBasisFSet_.size();
    for (const BasisFunctionSet *set : bFsets) {
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
    const std::vector<std::size_t> &globalNodeIndexes, std::size_t owner,
    bool owning) {
    logger.assert_debug(
        [&]() -> bool {
            for (std::size_t i = 0; i < globalNodeIndexes.size(); ++i)
                for (std::size_t j = 0; j < i; ++j)
                    if (globalNodeIndexes[i] == globalNodeIndexes[j])
                        return false;
            return true;
        }(),
        "Trying to pass the same node twice");
    auto result = theMesh_.addElement(globalNodeIndexes, owner, owning);
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
                    subElementNodeIndices, theMesh_.getNodeCoordinates());
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
    logger.assert_always(
        rawInput == "mesh 1"s,
        "This file is too new to be read by the current version of hpGEM");
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

    std::map<std::size_t, Element *> actualElement;

    for (std::size_t i = 0; i < numberOfElements; ++i) {
        std::size_t nodesPerElement;
        input >> nodesPerElement;
        std::vector<std::size_t> coordinateIndices;
        std::vector<std::size_t> nodeNumbers;
        for (std::size_t j = 0; j < nodesPerElement; ++j) {
            std::size_t globalIndex, coordinateOffset;
            input >> globalIndex >> coordinateOffset;
            coordinateIndices.push_back(startOfCoordinates[globalIndex] +
                                        coordinateOffset);
            nodeNumbers.push_back(localNodeIndex[globalIndex]);
        }
        std::size_t partition;
        input >> partition;
        Base::Element *element = nullptr;
        if (partition == processorID) {
            element = addElement(coordinateIndices, partition, true);
            actualElement[i] = element;
            getMesh().getSubmesh().add(element);
        }
        std::size_t numberOfShadowPartitions;
        input >> numberOfShadowPartitions;
        for (std::size_t j = 0; j < numberOfShadowPartitions; ++j) {
            std::size_t shadowPartition;
            input >> shadowPartition;
            if (shadowPartition == processorID) {
                element = addElement(coordinateIndices, partition, false);
                actualElement[i] = element;
                getMesh().getSubmesh().addPull(element, partition);
            }
            if (partition == processorID) {
                getMesh().getSubmesh().addPush(element, shadowPartition);
            }
        }
        if (element != nullptr) {
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

#ifdef HPGEM_USE_QHULL
template <std::size_t DIM>
void MeshManipulator<DIM>::createUnstructuredMesh(
    Geometry::PointPhysical<DIM> BottomLeft,
    Geometry::PointPhysical<DIM> TopRight, std::size_t TotalNoNodes,
    std::function<double(Geometry::PointPhysical<DIM>)> domainDescription,
    std::vector<Geometry::PointPhysical<DIM>> fixedPoints,
    std::function<double(Geometry::PointPhysical<DIM>)> relativeEdgeLength,
    double growFactor) {
    // impossible to create a mesh with more fixed nodes that total nodes
    // note that when equality is met, this will only do a delaunay
    // triangulation
    logger.assert_debug(
        fixedPoints.size() <= TotalNoNodes,
        "Cannot create a mesh with more fixed nodes than total nodes");

    // set to correct value in case some other meshmanipulator changed things
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

    // periodic unstructured mesh generation not yet implemented
    logger.assert_debug(
        !(periodicX_ || periodicY_ || periodicZ_),
        "Unstructured mesh generator does not support periodic boundaries");

    // guess the required distance between two nodes
    double dist = std::pow(double(TotalNoNodes), -1. / double(dimension()));

    for (std::size_t i = 0; i < dimension(); ++i) {
        dist *= std::pow(TopRight[i] - BottomLeft[i], 1. / double(dimension()));
    }

    std::vector<Geometry::PointPhysical<DIM>> hpGEMCoordinates = fixedPoints;

    // seed approximately N points inside the bounding box (total amount will be
    // tweaked later)
    Geometry::PointPhysical<DIM> nextPoint = BottomLeft;
    for (std::size_t i = 0; i < fixedPoints.size(); ++i) {
        logger.assert_debug(domainDescription(fixedPoints[i]) < 1e-10,
                            "One of the fixed points is outside of the domain");
        theMesh_.addNode();
        theMesh_.addNodeCoordinate(fixedPoints[i]);
    }
    // cant do nested for loops for generic dimension
    while (nextPoint[DIM - 1] < TopRight[DIM - 1] - 1e-10) {
        std::size_t incrementalDimension = 0;
        // if the point is already far enough to the right, reset and continue
        // with the next dimension
        for (; nextPoint[incrementalDimension] >
               TopRight[incrementalDimension] + 1e-10;
             ++incrementalDimension) {
            nextPoint[incrementalDimension] = BottomLeft[incrementalDimension];
        }
        nextPoint[incrementalDimension] += dist;
        if (domainDescription(nextPoint) < 0) {
            hpGEMCoordinates.push_back(nextPoint);
            theMesh_.addNode();
            theMesh_.addNodeCoordinate(nextPoint);
        }
    }
    std::size_t nFixedPoints = fixedPoints.size();
    // there are not enough points to do a triangulation
    logger.assert_debug(
        DIM < nFixedPoints,
        "Could not construct enough points for the initial triangulation");
    // there is inherent rounding down in the gridding and some nodes are
    // outside the domain (so they are discarded)
    logger.assert_debug(hpGEMCoordinates.size() <= TotalNoNodes,
                        "Constructed too many nodes");

    while (hpGEMCoordinates.size() < TotalNoNodes) {
        // start of QHull magic to create a triangulation
        orgQhull::RboxPoints qHullCoordinates;
        qHullCoordinates.setDimension(DIM);
        qHullCoordinates.reserveCoordinates(DIM * hpGEMCoordinates.size());

        for (Geometry::PointPhysical<DIM> point : hpGEMCoordinates) {
            qHullCoordinates.append(DIM, point.data());
        }

        // create the triangulation, pass "d" for delaunay
        //"QJ" because there are likely to be groups of more that (d+1)
        // cocircular nodes in a regular grid, so joggle them up a bit
        orgQhull::Qhull triangulation;
        triangulation.runQhull(qHullCoordinates, "d Qbb Qx Qc Qt");

        for (orgQhull::QhullFacet triangle : triangulation.facetList()) {
            if (triangle.isGood() && !triangle.isUpperDelaunay()) {
                Geometry::PointPhysical<DIM> center;
                std::vector<std::size_t> pointIndices;
                for (auto vertex : triangle.vertices()) {
                    logger.assert_debug(
                        vertex.point().id() >= 0,
                        "QHull breaks our assumptions on indexes");
                    center += hpGEMCoordinates[vertex.point().id()];
                    pointIndices.push_back(vertex.point().id());
                }
                center = center / pointIndices.size();
                if (domainDescription(center) < 0) {
                    auto newElement = addElement(pointIndices);
                    for (std::size_t i = 0; i < pointIndices.size(); ++i) {
                        theMesh_
                            .getNodesList(IteratorType::GLOBAL)[pointIndices[i]]
                            ->addElement(newElement, i);
                    }
                }
            }
        }
        // end of QHull magic

        // extract connectivity information
        faceFactory();
        edgeFactory();

        // compute current and expected (relative) edge length
        std::vector<double> expectedLength;
        std::multimap<double, std::size_t> knownLengths;
        std::vector<double> currentLength;
        // for proper scaling
        double totalcurrentLength = 0;
        double totalexpectedLength = 0;
        bool needsExpansion = false;

        // compute the expected relative length at the coordinates
        for (std::size_t i = 0; i < hpGEMCoordinates.size(); ++i) {
            double newLength = relativeEdgeLength(hpGEMCoordinates[i]);
            expectedLength.push_back(newLength);
            if (std::isnan(newLength) || std::isinf(newLength)) {
                needsExpansion |= true;
            } else {
                // cannot deliberately construct tangled meshes
                logger.assert_debug(
                    newLength > 0,
                    "Found an edge that is supposed to have a negative length");
                knownLengths.insert({newLength, i});
            }
        }

        // if the desired relative edge length is not known everywhere, slowly
        // make them larger because apparently the user is not interested in
        // controlling edge lengths for this part but sudden enlargement leads
        // to a bad mesh quality
        if (needsExpansion) {
            // iterate over all nodes, sorted by edge lengths
            for (std::pair<double, std::size_t> entry : knownLengths) {
                Node *current =
                    theMesh_.getNodesList(IteratorType::GLOBAL)[entry.second];
                for (Element *element : current->getElements()) {
                    for (std::size_t i = 0; i < element->getNumberOfNodes();
                         ++i) {
                        if (std::isnan(
                                expectedLength[element->getNode(i)->getID()]) ||
                            std::isinf(
                                expectedLength[element->getNode(i)->getID()])) {
                            expectedLength[element->getNode(i)->getID()] =
                                growFactor * entry.first;

                            // inserting does not invalidate the iterators;
                            // new node has a larger edge length, so it is
                            // guaranteed to be visited later on
                            knownLengths.insert(knownLengths.end(),
                                                {growFactor * entry.first,
                                                 element->getNode(i)->getID()});
                        }
                    }
                }
            }
        }
        // iterate over all edges to compute total length and scaling factor
        // the volume scales with (total edge length)^dimension
        // the total volume filled by the edges should be constant
        // so scale appropriately
        if (DIM == 1) {
            // the algorithm is mostly dimension independent, but the data type
            // it operates on is not
            for (Element *element :
                 theMesh_.getElementsList(IteratorType::GLOBAL)) {
                Geometry::PointPhysical<DIM> firstNode =
                    element->getPhysicalGeometry()->getLocalNodeCoordinates(0);
                Geometry::PointPhysical<DIM> secondNode =
                    element->getPhysicalGeometry()->getLocalNodeCoordinates(1);
                currentLength.push_back(L2Norm(firstNode - secondNode));
                totalcurrentLength += currentLength.back();
                totalexpectedLength +=
                    expectedLength[element->getNode(0)->getID()] / 2;
                totalexpectedLength +=
                    expectedLength[element->getNode(1)->getID()] / 2;
            }
        } else if (DIM == 2) {
            for (Face *face : theMesh_.getFacesList(IteratorType::GLOBAL)) {
                Geometry::PointPhysical<DIM> firstNode, secondNode;
                std::vector<std::size_t> nodeIndices =
                    face->getPtrElementLeft()
                        ->getReferenceGeometry()
                        ->getCodim1EntityLocalIndices(
                            face->localFaceNumberLeft());
                firstNode = face->getPtrElementLeft()
                                ->getPhysicalGeometry()
                                ->getLocalNodeCoordinates(nodeIndices[0]);
                secondNode = face->getPtrElementLeft()
                                 ->getPhysicalGeometry()
                                 ->getLocalNodeCoordinates(nodeIndices[1]);
                currentLength.push_back(L2Norm(firstNode - secondNode));
                totalcurrentLength +=
                    currentLength.back() * currentLength.back();
                totalexpectedLength +=
                    std::pow(expectedLength[face->getPtrElementLeft()
                                                ->getNode(nodeIndices[0])
                                                ->getID()] +
                                 expectedLength[face->getPtrElementLeft()
                                                    ->getNode(nodeIndices[1])
                                                    ->getID()],
                             2.) /
                    4.;
            }
        } else {
            for (Edge *edge : theMesh_.getEdgesList(IteratorType::GLOBAL)) {
                Geometry::PointPhysical<DIM> firstNode, secondNode;
                std::vector<std::size_t> nodeIndices =
                    edge->getElement(0)
                        ->getReferenceGeometry()
                        ->getCodim2EntityLocalIndices(edge->getEdgeNumber(0));
                firstNode = edge->getElement(0)
                                ->getPhysicalGeometry()
                                ->getLocalNodeCoordinates(nodeIndices[0]);
                secondNode = edge->getElement(0)
                                 ->getPhysicalGeometry()
                                 ->getLocalNodeCoordinates(nodeIndices[1]);
                currentLength.push_back(L2Norm(firstNode - secondNode));
                totalcurrentLength += currentLength.back() *
                                      currentLength.back() *
                                      currentLength.back();
                totalexpectedLength +=
                    std::pow(expectedLength[edge->getElement(0)
                                                ->getNode(nodeIndices[0])
                                                ->getID()] +
                                 expectedLength[edge->getElement(0)
                                                    ->getNode(nodeIndices[1])
                                                    ->getID()],
                             3.) /
                    8.;
            }
        }

        // all regions of the domain where elements are allowed to be as large
        // as possible must be connected to regions where relativeEdgeLength
        // provides a limitation
        logger.assert_debug(
            !std::isnan(totalexpectedLength) &&
                !std::isinf(totalexpectedLength),
            "could not infer edge sizes for the entirety of the domain");

        // sort the centers of the edges such that the centers of the large
        // edges are indexed first note that in this case the inverse measure is
        // computed, because that will result in a more natural force
        // computation later on
        std::multimap<double, Geometry::PointPhysical<DIM>> centerPoints;
        if (DIM == 1) {
            // the algorithm is mostly dimension independent, but the data type
            // it operates on is not
            for (Element *element :
                 theMesh_.getElementsList(IteratorType::GLOBAL)) {
                Geometry::PointPhysical<DIM> firstNode, secondNode;
                firstNode =
                    element->getPhysicalGeometry()->getLocalNodeCoordinates(0);
                secondNode =
                    element->getPhysicalGeometry()->getLocalNodeCoordinates(1);
                // length is scaled in case somebody hasty decides to add
                // smoothing at this point all edges should be squeezed a little
                // if the algorithm is to work correctly so pretend the volume
                // is 1.5 times as large remember to scale back from a volume
                // measure to a length measure
                double length = (expectedLength[element->getNode(0)->getID()] +
                                 expectedLength[element->getNode(1)->getID()]) /
                                currentLength[centerPoints.size()] * 2 *
                                totalcurrentLength / totalexpectedLength;
                centerPoints.insert({length, (firstNode + secondNode) / 2.});
            }
        } else if (DIM == 2) {
            for (Face *face : theMesh_.getFacesList(IteratorType::GLOBAL)) {
                Geometry::PointPhysical<DIM> firstNode, secondNode;
                std::vector<std::size_t> nodeIndices =
                    face->getPtrElementLeft()
                        ->getReferenceGeometry()
                        ->getCodim1EntityLocalIndices(
                            face->localFaceNumberLeft());
                firstNode = face->getPtrElementLeft()
                                ->getPhysicalGeometry()
                                ->getLocalNodeCoordinates(nodeIndices[0]);
                secondNode = face->getPtrElementLeft()
                                 ->getPhysicalGeometry()
                                 ->getLocalNodeCoordinates(nodeIndices[1]);
                // length is scaled in case somebody hasty decides to add
                // smoothing at this point all edges should be squeezed a little
                // if the algorithm is to work correctly so pretend the volume
                // is 1.5 times as large remember to scale back from a volume
                // measure to a length measure
                double length =
                    (expectedLength[face->getPtrElementLeft()
                                        ->getNode(nodeIndices[0])
                                        ->getID()] +
                     expectedLength[face->getPtrElementLeft()
                                        ->getNode(nodeIndices[1])
                                        ->getID()]) /
                    currentLength[centerPoints.size()] *
                    std::pow(2 * totalcurrentLength / totalexpectedLength,
                             1. / 2.);
                centerPoints.insert({length, (firstNode + secondNode) / 2});
            }
        } else {
            for (Edge *edge : theMesh_.getEdgesList(IteratorType::GLOBAL)) {
                Geometry::PointPhysical<DIM> firstNode, secondNode;
                std::vector<std::size_t> nodeIndices =
                    edge->getElement(0)
                        ->getReferenceGeometry()
                        ->getCodim2EntityLocalIndices(edge->getEdgeNumber(0));
                firstNode = edge->getElement(0)
                                ->getPhysicalGeometry()
                                ->getLocalNodeCoordinates(nodeIndices[0]);
                secondNode = edge->getElement(0)
                                 ->getPhysicalGeometry()
                                 ->getLocalNodeCoordinates(nodeIndices[1]);
                // length is scaled in case somebody hasty decides to add
                // smoothing at this point all edges should be squeezed a little
                // if the algorithm is to work correctly so pretend the volume
                // is 1.5 times as large remember to scale back from a volume
                // measure to a length measure
                double length =
                    (expectedLength[edge->getElement(0)
                                        ->getNode(nodeIndices[0])
                                        ->getID()] +
                     expectedLength[edge->getElement(0)
                                        ->getNode(nodeIndices[1])
                                        ->getID()]) /
                    currentLength[centerPoints.size()] *
                    std::pow(2 * totalcurrentLength / totalexpectedLength,
                             1. / 3.);
                centerPoints.insert({length, (firstNode + secondNode) / 2});
            }
        }
        // insert nodes in the longest edges only, in an attempt to make them
        // all equally long cannot use range based loop because the centerpoint
        // is not marked
        auto it = centerPoints.begin();
        std::size_t nNewNodes = std::min(
            centerPoints.size() / 2, TotalNoNodes - hpGEMCoordinates.size());
        for (std::size_t i = 0; i < nNewNodes && it != centerPoints.end();
             ++i, ++it) {
            if (domainDescription(it->second) < 0) {
                hpGEMCoordinates.push_back(it->second);
            } else {
                --i;
            }
        }

        // cleanest solution, but not the fastest
        theMesh_.clear();
        for (Geometry::PointPhysical<DIM> point : hpGEMCoordinates) {
            theMesh_.addNode();
            theMesh_.addNodeCoordinate(point);
        }
    }
    // start of QHull magic to create a triangulation
    orgQhull::RboxPoints qHullCoordinates;
    qHullCoordinates.setDimension(DIM);
    qHullCoordinates.reserveCoordinates(DIM * hpGEMCoordinates.size());

    for (Geometry::PointPhysical<DIM> point : hpGEMCoordinates) {
        qHullCoordinates.append(DIM, point.data());
    }

    // create the triangulation, pass "d" for delaunay
    // no "QJ" because periodic boundary nodes must keep the current exact
    // distaince (including direction)
    orgQhull::Qhull triangulation(qHullCoordinates, "d QbB Qx Qc Qt");

    for (orgQhull::QhullFacet triangle : triangulation.facetList()) {
        if (triangle.isGood() && !triangle.isUpperDelaunay()) {
            Geometry::PointPhysical<DIM> center;
            std::vector<std::size_t> pointIndices;
            for (auto vertexIt1 = triangle.vertices().begin();
                 vertexIt1 != triangle.vertices().end(); ++vertexIt1) {
                center += hpGEMCoordinates[(*vertexIt1).point().id()];
                pointIndices.push_back((*vertexIt1).point().id());
            }
            center = center / pointIndices.size();
            if (domainDescription(center) < 0) {
                auto newElement = addElement(pointIndices);
                for (std::size_t i = 0; i < pointIndices.size(); ++i) {
                    theMesh_
                        .getNodesList(IteratorType::GLOBAL)[pointIndices[i]]
                        ->addElement(newElement, i);
                }
            }
        }
    }
    // end of QHull magic

    edgeFactory();
    faceFactory();

    std::vector<std::size_t> fixedPointIdxs;
    for (std::size_t i = 0; i < nFixedPoints; ++i) {
        fixedPointIdxs.push_back(i);
    }

    updateMesh(domainDescription, fixedPointIdxs, relativeEdgeLength,
               growFactor);
}

///\bug Assumes a DG basis function set is used. (Workaround: set the basis
/// function set again after calling this routine if you are using something
/// conforming)
template <std::size_t DIM>
void MeshManipulator<DIM>::updateMesh(
    std::function<double(Geometry::PointPhysical<DIM>)> domainDescription,
    std::vector<std::size_t> fixedPointIdxs,
    std::function<double(Geometry::PointPhysical<DIM>)> relativeEdgeLength,
    double growFactor,
    std::function<bool(Geometry::PointPhysical<DIM>)> isOnPeriodic,
    std::function<Geometry::PointPhysical<DIM>(Geometry::PointPhysical<DIM>)>
        duplicatePeriodic,
    std::function<bool(Geometry::PointPhysical<DIM>)> isOnOtherPeriodic,
    std::function<Geometry::PointPhysical<DIM>(Geometry::PointPhysical<DIM>)>
        safeSpot,
    std::vector<std::size_t> dontConnect) {
    std::sort(fixedPointIdxs.begin(), fixedPointIdxs.end());
    bool needsExpansion = false;
    double totalCurrentLength = 0;
    double oldQuality = 0;
    double worstQuality = 0.5;

    std::set<std::pair<std::size_t, std::size_t>> periodicPairing{};

    // set to correct value in case some other meshmanipulator changed things
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

    for (Node *node : theMesh_.getNodesList(IteratorType::GLOBAL)) {
        Geometry::PointPhysical<DIM> point;
        Geometry::PointPhysical<DIM> compare;
        point =
            node->getElement(0)->getPhysicalGeometry()->getLocalNodeCoordinates(
                node->getNodeNumber(0));
        std::set<std::size_t> equivalentIndices{};
        equivalentIndices.insert(
            node->getElement(0)->getPhysicalGeometry()->getNodeIndex(
                node->getNodeNumber(0)));
        for (std::size_t i = 1; i < node->getNumberOfElements(); ++i) {
            compare = node->getElement(i)
                          ->getPhysicalGeometry()
                          ->getLocalNodeCoordinates(node->getNodeNumber(i));
            if (compare != point) {
                equivalentIndices.insert(
                    node->getElement(i)->getPhysicalGeometry()->getNodeIndex(
                        node->getNodeNumber(i)));
            }
        }
        auto equivalentIterator = equivalentIndices.begin();
        ++equivalentIterator;
        for (; equivalentIterator != equivalentIndices.end();
             ++equivalentIterator) {
            periodicPairing.insert(
                {*equivalentIndices.begin(), *equivalentIterator});
        }
    }

    // compute the lengths of the edges and how far the nodes have moved, to see
    // if the nodes have moved so far that a retriangulation is in order
    double maxShift = 0;
    // except don't bother if a retriangulation is in order anyway
    if (oldNodeLocations_.size() == theMesh_.getNodeCoordinates().size()) {
        std::vector<double> unscaledShift{};
        unscaledShift.reserve(theMesh_.getNumberOfNodes(IteratorType::GLOBAL));
        // compute current and expected (relative) edge length
        std::vector<double> expectedLength{};
        expectedLength.reserve(theMesh_.getNumberOfNodes(IteratorType::GLOBAL));
        std::multimap<double, std::size_t> knownLengths;
        std::vector<double> currentLength{};
        // for proper scaling
        double totalexpectedLength = 0;
        for (Node *node : theMesh_.getNodesList(IteratorType::GLOBAL)) {
            Geometry::PointPhysical<DIM> point =
                node->getElement(0)
                    ->getPhysicalGeometry()
                    ->getLocalNodeCoordinates(node->getNodeNumber(0));
            unscaledShift.push_back(
                L2Norm(oldNodeLocations_[expectedLength.size()] - point));
            expectedLength.push_back(relativeEdgeLength(point));
            if (isnan(expectedLength.back()) || isinf(expectedLength.back())) {
                needsExpansion |= true;
            } else {
                knownLengths.insert(
                    {expectedLength.back(), expectedLength.size() - 1});
            }
        }
        if (needsExpansion) {
            // iterate over all nodes, sorted by edge lengths
            for (std::pair<double, std::size_t> entry : knownLengths) {
                Node *current =
                    theMesh_.getNodesList(IteratorType::GLOBAL)[entry.second];
                for (Element *element : current->getElements()) {
                    for (std::size_t i = 0; i < element->getNumberOfNodes();
                         ++i) {
                        if (std::isnan(
                                expectedLength[element->getNode(i)->getID()]) ||
                            std::isinf(
                                expectedLength[element->getNode(i)->getID()])) {
                            expectedLength[element->getNode(i)->getID()] =
                                growFactor * entry.first;

                            // inserting does not invalidate the iterators;
                            // new node has a larger edge length, so it is
                            // guaranteed to be visited later on

                            knownLengths.insert(knownLengths.end(),
                                                {growFactor * entry.first,
                                                 element->getNode(i)->getID()});
                        }
                    }
                }
            }
        }
        // iterate over all edges to compute total length and scaling factor
        // the volume scales with (total edge length)^dimension
        // the total volume filled by the edges should be constant
        // so scale appropriately
        if (DIM == 1) {
            currentLength.reserve(
                theMesh_.getNumberOfElements(IteratorType::GLOBAL));
            // the algorithm is mostly dimension independent, but the data type
            // it operates on is not
            for (Element *element :
                 theMesh_.getElementsList(IteratorType::GLOBAL)) {
                Geometry::PointPhysical<DIM> firstNode, secondNode;
                firstNode =
                    element->getPhysicalGeometry()->getLocalNodeCoordinates(0);
                secondNode =
                    element->getPhysicalGeometry()->getLocalNodeCoordinates(1);
                currentLength.push_back(L2Norm(firstNode - secondNode));
                totalCurrentLength += currentLength.back();
                totalexpectedLength +=
                    expectedLength[element->getNode(0)->getID()] / 2;
                totalexpectedLength +=
                    expectedLength[element->getNode(1)->getID()] / 2;
                maxShift = std::max(
                    maxShift,
                    L2Norm(firstNode -
                           oldNodeLocations_[element->getPhysicalGeometry()
                                                 ->getNodeIndex(0)]) /
                        currentLength.back());
                maxShift = std::max(
                    maxShift,
                    L2Norm(secondNode -
                           oldNodeLocations_[element->getPhysicalGeometry()
                                                 ->getNodeIndex(1)]) /
                        currentLength.back());
            }
        } else if (DIM == 2) {
            currentLength.reserve(
                theMesh_.getNumberOfFaces(IteratorType::GLOBAL));
            for (Face *face : theMesh_.getFacesList(IteratorType::GLOBAL)) {
                Geometry::PointPhysical<DIM> firstNode, secondNode;
                std::vector<std::size_t> nodeIndices =
                    face->getPtrElementLeft()
                        ->getReferenceGeometry()
                        ->getCodim1EntityLocalIndices(
                            face->localFaceNumberLeft());
                firstNode = face->getPtrElementLeft()
                                ->getPhysicalGeometry()
                                ->getLocalNodeCoordinates(nodeIndices[0]);
                secondNode = face->getPtrElementLeft()
                                 ->getPhysicalGeometry()
                                 ->getLocalNodeCoordinates(nodeIndices[1]);
                currentLength.push_back(L2Norm(firstNode - secondNode));
                totalCurrentLength +=
                    currentLength.back() * currentLength.back();
                totalexpectedLength +=
                    std::pow(expectedLength[face->getPtrElementLeft()
                                                ->getNode(nodeIndices[0])
                                                ->getID()] +
                                 expectedLength[face->getPtrElementLeft()
                                                    ->getNode(nodeIndices[1])
                                                    ->getID()],
                             2.) /
                    4.;
                maxShift = std::max(
                    maxShift,
                    L2Norm(
                        firstNode -
                        oldNodeLocations_[face->getPtrElementLeft()
                                              ->getPhysicalGeometry()
                                              ->getNodeIndex(nodeIndices[0])]) /
                        currentLength.back());
                maxShift = std::max(
                    maxShift,
                    L2Norm(
                        secondNode -
                        oldNodeLocations_[face->getPtrElementLeft()
                                              ->getPhysicalGeometry()
                                              ->getNodeIndex(nodeIndices[1])]) /
                        currentLength.back());
            }
            worstQuality = 1;
            for (Element *element :
                 theMesh_.getElementsList(IteratorType::GLOBAL)) {
                std::array<double, 3> edgeLengths;
                for (std::size_t i = 0; i < 3; ++i) {
                    edgeLengths[i] =
                        currentLength[element->getFace(i)->getID()];
                }
                worstQuality = std::min(
                    worstQuality,
                    (edgeLengths[0] + edgeLengths[1] - edgeLengths[2]) *
                        (edgeLengths[1] + edgeLengths[2] - edgeLengths[0]) *
                        (edgeLengths[2] + edgeLengths[0] - edgeLengths[1]) /
                        edgeLengths[0] / edgeLengths[1] / edgeLengths[2]);
            }
        } else {
            currentLength.reserve(
                theMesh_.getNumberOfEdges(IteratorType::GLOBAL));
            for (Edge *edge : theMesh_.getEdgesList(IteratorType::GLOBAL)) {
                Geometry::PointPhysical<DIM> firstNode, secondNode;
                std::vector<std::size_t> nodeIndices =
                    edge->getElement(0)
                        ->getReferenceGeometry()
                        ->getCodim2EntityLocalIndices(edge->getEdgeNumber(0));
                firstNode = edge->getElement(0)
                                ->getPhysicalGeometry()
                                ->getLocalNodeCoordinates(nodeIndices[0]);
                secondNode = edge->getElement(0)
                                 ->getPhysicalGeometry()
                                 ->getLocalNodeCoordinates(nodeIndices[1]);
                currentLength.push_back(L2Norm(firstNode - secondNode));
                totalCurrentLength += currentLength.back() *
                                      currentLength.back() *
                                      currentLength.back();
                totalexpectedLength +=
                    std::pow(expectedLength[edge->getElement(0)
                                                ->getNode(nodeIndices[0])
                                                ->getID()] +
                                 expectedLength[edge->getElement(0)
                                                    ->getNode(nodeIndices[1])
                                                    ->getID()],
                             3.) /
                    8.;
                maxShift = std::max(
                    maxShift,
                    L2Norm(
                        firstNode -
                        oldNodeLocations_[edge->getElement(0)
                                              ->getPhysicalGeometry()
                                              ->getNodeIndex(nodeIndices[0])]) /
                        currentLength.back());
                maxShift = std::max(
                    maxShift,
                    L2Norm(
                        secondNode -
                        oldNodeLocations_[edge->getElement(0)
                                              ->getPhysicalGeometry()
                                              ->getNodeIndex(nodeIndices[1])]) /
                        currentLength.back());
            }
        }

        // all regions of the domain where elements are allowed to be as large
        // as possible must be connected to regions where relativeEdgeLength
        // provides a limitation
        logger.assert_debug(
            !std::isnan(totalexpectedLength) &&
                !std::isinf(totalexpectedLength),
            "Could not infer edge sizes for the entirety of the domain");
    }
    std::size_t counter = 0;
    double maxMovement = std::numeric_limits<double>::infinity();
    std::map<size_t, double> currentLength{};
    std::map<std::size_t, Geometry::PointPhysical<DIM>> movement;
    // stop after n iterations, or (when the nodes have stopped moving and the
    // mesh is not becoming worse), or when the mesh is great, or when the mesh
    // is decent, but worsening
    while ((counter < 10000 &&
            (maxMovement > 1e-3 || oldQuality - worstQuality > 1e-3) &&
            worstQuality < 0.8 &&
            (worstQuality < 2. / 3. || oldQuality - worstQuality < 0)) ||
           counter < 5) {
        counter++;
        if ((maxShift > 0.1 && (oldQuality - worstQuality) < 5e-3 * maxShift) ||
            (oldNodeLocations_.size() !=
             theMesh_.getNumberOfNodeCoordinates()) ||
            worstQuality < 1e-6) {
            maxShift = 0;

            orgQhull::RboxPoints qHullCoordinates{};
            qHullCoordinates.setDimension(DIM);
            oldNodeLocations_.clear();
            oldNodeLocations_.reserve(theMesh_.getNumberOfNodeCoordinates());
            for (Geometry::PointPhysical<DIM> &point :
                 theMesh_.getNodeCoordinates()) {
                oldNodeLocations_.push_back(point);
            }
            theMesh_.clear();
            auto pairingIterator = periodicPairing.begin();
            for (Geometry::PointPhysical<DIM> point : oldNodeLocations_) {
                theMesh_.addNodeCoordinate(point);
                if (pairingIterator == periodicPairing.end()) {
                    theMesh_.addNode();
                } else {
                    // skip one insertion for each master/slave pair
                    pairingIterator++;
                }
            }

            std::vector<std::size_t> vertexIndex{};
            vertexIndex.resize(theMesh_.getNumberOfNodeCoordinates(),
                               std::numeric_limits<std::size_t>::max());
            pairingIterator = periodicPairing.begin();
            std::size_t currentNodeNumber = 0;
            for (std::size_t i = 0;
                 i < theMesh_.getNumberOfNodeCoordinates();) {
                vertexIndex[i] = currentNodeNumber;
                // see if there are any new boundary nodes
                if (isOnPeriodic(theMesh_.getNodeCoordinates()[i]) &&
                    !(pairingIterator != periodicPairing.end() &&
                      pairingIterator->first == i)) {
                    std::size_t j = pairingIterator->first;
                    periodicPairing.insert(
                        {i, theMesh_.getNumberOfNodeCoordinates()});
                    logger(DEBUG, "periodic pair: % % ", i,
                           theMesh_.getNumberOfNodeCoordinates());
                    pairingIterator = std::find_if(
                        periodicPairing.begin(), periodicPairing.end(),
                        [=](const std::pair<std::size_t, std::size_t> &p)
                            -> bool { return p.first == std::min(i, j); });
                    Geometry::PointPhysical<DIM> newNodeCoordinate =
                        duplicatePeriodic(theMesh_.getNodeCoordinates()[i]);
                    logger(DEBUG, "new periodic pair coordinates: % %",
                           theMesh_.getNodeCoordinates()[i], newNodeCoordinate);
                    theMesh_.addNodeCoordinate(newNodeCoordinate);
                    oldNodeLocations_.push_back(newNodeCoordinate);
                    vertexIndex.resize(theMesh_.getNumberOfNodeCoordinates(),
                                       std::numeric_limits<std::size_t>::max());
                }
                // see if there are any non-boundary nodes that slipped into the
                // boundary
                if (isOnOtherPeriodic(theMesh_.getNodeCoordinates()[i]) &&
                    !(pairingIterator != periodicPairing.end() &&
                      pairingIterator->first == i)) {
                    theMesh_.getNodeCoordinates()[i] =
                        safeSpot(theMesh_.getNodeCoordinates()[i]);
                }
                // assign boundary nodes
                while (pairingIterator != periodicPairing.end() &&
                       pairingIterator->first == i) {
                    logger(DEBUG, "periodic pair: % % ", pairingIterator->first,
                           pairingIterator->second);
                    logger.assert_debug(
                        Base::L2Norm(
                            duplicatePeriodic(
                                theMesh_.getNodeCoordinates()[pairingIterator
                                                                  ->first]) -
                            theMesh_.getNodeCoordinates()[pairingIterator
                                                              ->second]) < 1e-9,
                        "periodic pair is not moving simulateously");
                    vertexIndex[pairingIterator->second] = currentNodeNumber;
                    ++pairingIterator;
                }
                currentNodeNumber++;
                // skip over already set boundary nodes
                while (i < theMesh_.getNumberOfNodeCoordinates() &&
                       vertexIndex[i] <
                           std::numeric_limits<std::size_t>::max()) {
                    ++i;
                }
            }
            logger(DEBUG, "periodic pairs end");

            // all periodic boundary pairs are used
            logger.assert_debug(pairingIterator == periodicPairing.end(),
                                "Somehow missed some periodic pair");
            // the actual amount of vertices and the assigned amount of vertices
            // match
            logger.assert_debug(currentNodeNumber == theMesh_.getNumberOfNodes(
                                                         IteratorType::GLOBAL),
                                "Missed some node indexes");

            qHullCoordinates.reserveCoordinates(
                DIM * theMesh_.getNumberOfNodeCoordinates());
            for (Geometry::PointPhysical<DIM> &point :
                 theMesh_.getNodeCoordinates()) {
                qHullCoordinates.append(DIM, point.data());
            }

            orgQhull::Qhull triangulation(qHullCoordinates,
                                          "d PF1e-10 QbB Qx Qc Qt");

            for (orgQhull::QhullFacet triangle : triangulation.facetList()) {
                if (triangle.isGood() && !triangle.isUpperDelaunay()) {
                    logger(DEBUG, "adding %", triangle);
                    Geometry::PointPhysical<DIM> center;
                    std::vector<std::size_t> pointIndices{};
                    bool shouldConnect = false;
                    for (auto vertexIt1 = triangle.vertices().begin();
                         vertexIt1 != triangle.vertices().end(); ++vertexIt1) {
                        logger.assert_debug(
                            (*vertexIt1).point().id() >= 0,
                            "QHull breaks our assumptions on indexes");
                        logger(DEBUG, "% % %", shouldConnect,
                               vertexIndex[(*vertexIt1).point().id()],
                               *vertexIt1);
                        center += oldNodeLocations_[(*vertexIt1).point().id()];
                        pointIndices.push_back((*vertexIt1).point().id());
                        shouldConnect |=
                            (std::find(
                                 dontConnect.begin(), dontConnect.end(),
                                 vertexIndex[(*vertexIt1).point().id()]) ==
                             dontConnect.end());
                    }
                    center = center / pointIndices.size();
                    if (domainDescription(center) < -1e-10 && shouldConnect) {
                        auto newElement = addElement(pointIndices);
                        for (std::size_t i = 0; i < pointIndices.size(); ++i) {
                            theMesh_
                                .getNodesList(IteratorType::GLOBAL)
                                    [vertexIndex[pointIndices[i]]]
                                ->addElement(newElement, i);
                        }
                    } else {
                        logger(VERBOSE, "external element % ignored", triangle);
                    }
                }
                if (!triangle.isGood() && !triangle.isUpperDelaunay()) {
                    logger(VERBOSE, "small element % ignored", triangle);
                }
            }
            for (Node *node : theMesh_.getNodesList(IteratorType::GLOBAL)) {
                // all of the nodes should be in the interior of the domain or
                // near the boundary of the domain
                if (node->getNumberOfElements() == 0) {
                    for (std::size_t i = 0; i < vertexIndex.size(); ++i) {
                        if (vertexIndex[i] == node->getID()) {
                            logger(DEBUG, "% % %", i,
                                   theMesh_.getNodeCoordinates()[i],
                                   domainDescription(
                                       theMesh_.getNodeCoordinates()[i]));
                        }
                    }
                }
                logger.assert_debug(
                    node->getNumberOfElements() > 0,
                    "There is an node without any elements connected to it");
            }
            edgeFactory();
            faceFactory();
        }
        oldQuality = worstQuality;

        std::map<std::size_t, double> expectedLength{};
        std::multimap<double, std::size_t> knownLengths{};
        // for proper scaling
        totalCurrentLength = 0;
        double totalexpectedLength = 0.;
        for (Node *node : theMesh_.getNodesList(IteratorType::GLOBAL)) {
            Geometry::PointPhysical<DIM> point =
                node->getElement(0)
                    ->getPhysicalGeometry()
                    ->getLocalNodeCoordinates(node->getNodeNumber(0));
            expectedLength[node->getID()] = relativeEdgeLength(point);
            if (!isnan(expectedLength[node->getID()]) &&
                !isinf(expectedLength[node->getID()])) {
                knownLengths.insert(
                    {expectedLength[node->getID()], expectedLength.size() - 1});
            } else {
                needsExpansion |= true;
            }
        }
        if (needsExpansion) {
            // iterate over all nodes, sorted by edge lengths
            for (std::pair<double, std::size_t> entry : knownLengths) {
                Node *current =
                    theMesh_.getNodesList(IteratorType::GLOBAL)[entry.second];
                for (Element *element : current->getElements()) {
                    for (std::size_t i = 0; i < element->getNumberOfNodes();
                         ++i) {
                        if (std::isnan(
                                expectedLength[element->getNode(i)->getID()]) ||
                            std::isinf(
                                expectedLength[element->getNode(i)->getID()])) {
                            expectedLength[element->getNode(i)->getID()] =
                                growFactor * entry.first;

                            // inserting does not invalidate the iterators;
                            // new node has a larger edge length, so it is
                            // guaranteed to be visited later on

                            knownLengths.insert(
                                knownLengths.end(),
                                {growFactor * entry.first,
                                 std::find(
                                     getNodesList(IteratorType::GLOBAL).begin(),
                                     getNodesList(IteratorType::GLOBAL).end(),
                                     element->getNode(i)) -
                                     getNodesList(IteratorType::GLOBAL)
                                         .begin()});
                        }
                    }
                }
            }
        }

        // iterate over all edges to compute total length and scaling factor
        // the volume scales with (total edge length)^dimension
        // the total volume filled by the edges should be constant
        // so scale appropriately
        totalCurrentLength = 0;
        totalexpectedLength = 0;
        if (DIM == 1) {
            // the algorithm is mostly dimension independent, but the data type
            // it operates on is not
            for (Element *element :
                 theMesh_.getElementsList(IteratorType::GLOBAL)) {
                Geometry::PointPhysical<DIM> firstNode, secondNode;
                firstNode =
                    element->getPhysicalGeometry()->getLocalNodeCoordinates(0);
                secondNode =
                    element->getPhysicalGeometry()->getLocalNodeCoordinates(1);
                currentLength[element->getID()] =
                    L2Norm(firstNode - secondNode);
                totalCurrentLength += currentLength[element->getID()];
                totalexpectedLength +=
                    expectedLength[element->getNode(0)->getID()] / 2.;
                totalexpectedLength +=
                    expectedLength[element->getNode(1)->getID()] / 2.;
            }
        } else if (DIM == 2) {
            for (Face *face : theMesh_.getFacesList(IteratorType::GLOBAL)) {
                Geometry::PointPhysical<DIM> firstNode, secondNode;
                std::vector<std::size_t> nodeIndices =
                    face->getPtrElementLeft()
                        ->getReferenceGeometry()
                        ->getCodim1EntityLocalIndices(
                            face->localFaceNumberLeft());
                firstNode = face->getPtrElementLeft()
                                ->getPhysicalGeometry()
                                ->getLocalNodeCoordinates(nodeIndices[0]);
                secondNode = face->getPtrElementLeft()
                                 ->getPhysicalGeometry()
                                 ->getLocalNodeCoordinates(nodeIndices[1]);
                currentLength[face->getID()] = L2Norm(firstNode - secondNode);
                totalCurrentLength +=
                    currentLength[face->getID()] * currentLength[face->getID()];
                totalexpectedLength +=
                    std::pow(expectedLength[face->getPtrElementLeft()
                                                ->getNode(nodeIndices[0])
                                                ->getID()] +
                                 expectedLength[face->getPtrElementLeft()
                                                    ->getNode(nodeIndices[1])
                                                    ->getID()],
                             2.) /
                    4.;
            }
        } else {
            for (Edge *edge : theMesh_.getEdgesList(IteratorType::GLOBAL)) {
                Geometry::PointPhysical<DIM> firstNode, secondNode;
                std::vector<std::size_t> nodeIndices =
                    edge->getElement(0)
                        ->getReferenceGeometry()
                        ->getCodim2EntityLocalIndices(edge->getEdgeNumber(0));
                firstNode = edge->getElement(0)
                                ->getPhysicalGeometry()
                                ->getLocalNodeCoordinates(nodeIndices[0]);
                secondNode = edge->getElement(0)
                                 ->getPhysicalGeometry()
                                 ->getLocalNodeCoordinates(nodeIndices[1]);
                currentLength[edge->getID()] = L2Norm(firstNode - secondNode);
                totalCurrentLength += currentLength[edge->getID()] *
                                      currentLength[edge->getID()] *
                                      currentLength[edge->getID()];
                totalexpectedLength +=
                    std::pow(expectedLength[edge->getElement(0)
                                                ->getNode(nodeIndices[0])
                                                ->getID()] +
                                 expectedLength[edge->getElement(0)
                                                    ->getNode(nodeIndices[1])
                                                    ->getID()],
                             3.) /
                    8.;
            }
        }

        movement.clear();

        if (DIM == 1) {
            for (Element *element :
                 theMesh_.getElementsList(IteratorType::GLOBAL)) {
                Geometry::PointPhysical<DIM> firstNode, secondNode;
                firstNode =
                    element->getPhysicalGeometry()->getLocalNodeCoordinates(0);
                secondNode =
                    element->getPhysicalGeometry()->getLocalNodeCoordinates(1);
                // it is impossible to detect if a node inside the domain should
                // be clipped to the edge instead make sure that the nodes that
                // DO belong dont get pulled into the interior roundoff error
                // should make sure that nodes move away from the boundary if
                // there are too many all edges should be squeezed a little if
                // the algorithm is to work correctly so pretend the volume
                // is 1.4 times as large remember to scale back from a volume
                // measure to a length measure the non-linearity makes
                // everything slightly more robust
                double length = (expectedLength[element->getNode(0)->getID()] +
                                 expectedLength[element->getNode(1)->getID()]) /
                                currentLength[element->getID()] * 1.4 *
                                totalCurrentLength / totalexpectedLength / 2.;
                movement[element->getNode(0)->getID()] +=
                    std::max(length - 1., 0.) * (firstNode - secondNode) *
                    (length + 1.) * 0.5;
                movement[element->getNode(1)->getID()] +=
                    std::max(length - 1., 0.) * (secondNode - firstNode) *
                    (length + 1.) * 0.5;
            }
        } else if (DIM == 2) {
            for (Face *face : theMesh_.getFacesList(IteratorType::GLOBAL)) {
                Geometry::PointPhysical<DIM> firstNode, secondNode;
                std::vector<std::size_t> nodeIndices =
                    face->getPtrElementLeft()
                        ->getReferenceGeometry()
                        ->getCodim1EntityLocalIndices(
                            face->localFaceNumberLeft());
                firstNode = face->getPtrElementLeft()
                                ->getPhysicalGeometry()
                                ->getLocalNodeCoordinates(nodeIndices[0]);
                secondNode = face->getPtrElementLeft()
                                 ->getPhysicalGeometry()
                                 ->getLocalNodeCoordinates(nodeIndices[1]);
                double length =
                    (expectedLength[face->getPtrElementLeft()
                                        ->getNode(nodeIndices[0])
                                        ->getID()] +
                     expectedLength[face->getPtrElementLeft()
                                        ->getNode(nodeIndices[1])
                                        ->getID()]) /
                    currentLength[face->getID()] *
                    std::pow(1.4 * totalCurrentLength / totalexpectedLength,
                             1. / 2.) /
                    2.;
                movement[face->getPtrElementLeft()
                             ->getNode(nodeIndices[0])
                             ->getID()] += std::max(length - 1., 0.) *
                                           (firstNode - secondNode) *
                                           (length + 1.) * 0.5;
                movement[face->getPtrElementLeft()
                             ->getNode(nodeIndices[1])
                             ->getID()] += std::max(length - 1., 0.) *
                                           (secondNode - firstNode) *
                                           (length + 1.) * 0.5;
            }
        } else if (DIM == 3) {
            for (Edge *edge : theMesh_.getEdgesList(IteratorType::GLOBAL)) {
                Geometry::PointPhysical<DIM> firstNode, secondNode;
                std::vector<std::size_t> nodeIndices =
                    edge->getElement(0)
                        ->getReferenceGeometry()
                        ->getCodim2EntityLocalIndices(edge->getEdgeNumber(0));
                firstNode = edge->getElement(0)
                                ->getPhysicalGeometry()
                                ->getLocalNodeCoordinates(nodeIndices[0]);
                secondNode = edge->getElement(0)
                                 ->getPhysicalGeometry()
                                 ->getLocalNodeCoordinates(nodeIndices[1]);
                double length =
                    (expectedLength[edge->getElement(0)
                                        ->getNode(nodeIndices[0])
                                        ->getID()] +
                     expectedLength[edge->getElement(0)
                                        ->getNode(nodeIndices[1])
                                        ->getID()]) /
                    currentLength[edge->getID()] *
                    std::pow(1.4 * totalCurrentLength / totalexpectedLength,
                             1. / 3.) /
                    2.;
                movement
                    [edge->getElement(0)->getNode(nodeIndices[0])->getID()] +=
                    std::max(length - 1., 0.) * (firstNode - secondNode) *
                    (length + 1.) * 0.5;
                movement
                    [edge->getElement(0)->getNode(nodeIndices[1])->getID()] +=
                    std::max(length - 1., 0.) * (secondNode - firstNode) *
                    (length + 1.) * 0.5;
            }
        }

        // forward Euler discretisation of an optimally damped mass-spring
        // system, with time step 0.02 this time step could be 0.1, but there is
        // a stability issue where springs aligned along the periodic boundary
        // are applied twice
        maxMovement = 0;
        auto moveIterator = movement.begin();
        auto fixIterator = fixedPointIdxs.begin();
        for (std::size_t i = 0;
             i < theMesh_.getNumberOfNodes(IteratorType::GLOBAL);
             ++moveIterator, ++i) {
            if (fixIterator != fixedPointIdxs.end() && i == *fixIterator) {
                ++fixIterator;
                std::get<1>(*moveIterator) *= 0;
            } else {
                Node *node = theMesh_.getNodesList(IteratorType::GLOBAL)[i];
                Geometry::PointPhysical<DIM> &point =
                    theMesh_
                        .getNodeCoordinates()[node->getElement(0)
                                                  ->getPhysicalGeometry()
                                                  ->getNodeIndex(
                                                      node->getNodeNumber(0))];
                point += 0.1 * std::get<1>(*moveIterator);
                logger.assert_debug(!(std::isnan(point[0])), "%", i);
                bool isPeriodic = false;
                std::map<std::size_t, bool> hasMoved{};
                hasMoved[node->getElement(0)
                             ->getPhysicalGeometry()
                             ->getNodeIndex(node->getNodeNumber(0))] = true;
                for (std::size_t j = 1; j < node->getNumberOfElements(); ++j) {
                    if (!hasMoved[node->getElement(j)
                                      ->getPhysicalGeometry()
                                      ->getNodeIndex(node->getNodeNumber(j))]) {
                        Geometry::PointPhysical<DIM> &other =
                            theMesh_.getNodeCoordinates()
                                [node->getElement(j)
                                     ->getPhysicalGeometry()
                                     ->getNodeIndex(node->getNodeNumber(j))];
                        other += 0.1 * std::get<1>(*moveIterator);
                        hasMoved[node->getElement(j)
                                     ->getPhysicalGeometry()
                                     ->getNodeIndex(node->getNodeNumber(j))] =
                            true;
                        isPeriodic = true;
                    }
                }
                if (domainDescription(point) > 0 && !isPeriodic) {
                    // the point is outside of the domain, move it back inside
                    double currentValue = domainDescription(point);
                    LinearAlgebra::SmallVector<DIM> gradient;
                    LinearAlgebra::SmallVector<DIM> offset;
                    // one-sided numerical derivative
                    for (std::size_t j = 0; j < DIM; ++j) {
                        offset[j] = 1e-7;
                        gradient[j] =
                            (currentValue - domainDescription(point + offset)) *
                            1e7;
                        offset[j] = 0;
                    }
                    point += currentValue * gradient / L2Norm(gradient);
                    std::get<1>(*moveIterator) +=
                        10 * currentValue * gradient / L2Norm(gradient);
                    currentValue = domainDescription(point);
                    // second step for robustness and accuracy if needed
                    if (currentValue > 0) {
                        for (std::size_t j = 0; j < DIM; ++j) {
                            offset[j] = 1e-7;
                            gradient[j] = (currentValue -
                                           domainDescription(point + offset)) *
                                          1e7;
                            offset[j] = 0;
                        }
                        point += currentValue * gradient / L2Norm(gradient);
                        std::get<1>(*moveIterator) +=
                            10 * currentValue * gradient / L2Norm(gradient);
                        // if two steps are not enough, more are also not likely
                        // to help
                        currentValue = domainDescription(point);
                        if (currentValue > 1e-10) {
                            logger(WARN,
                                   "NOTE: Failed to move point % (%) back into "
                                   "the domain."
                                   "\n Distance from boundary is %. Algorithm "
                                   "may crash.\n Consider fixing "
                                   "points at corners to remedy this issue.",
                                   i, point, currentValue);
                        }
                    }
                }
                if (isPeriodic) {
                    // do a total of four newton iteration before giving up
                    Geometry::PointPhysical<DIM> testPoint;
                    for (std::size_t j = 0; j < 4; ++j) {
                        // make sure the node stays on the periodic boundary, to
                        // prevent faces with 3 or more elements connected to
                        // them
                        for (std::size_t k = 0; k < node->getNumberOfElements();
                             ++k) {
                            testPoint = node->getElement(k)
                                            ->getPhysicalGeometry()
                                            ->getLocalNodeCoordinates(
                                                node->getNodeNumber(k));
                            double currentValue = domainDescription(testPoint);
                            if (currentValue > 0) {
                                LinearAlgebra::SmallVector<DIM> gradient;
                                LinearAlgebra::SmallVector<DIM> offset;
                                for (std::size_t l = 0; l < DIM; ++l) {
                                    offset[l] = 1e-7;
                                    gradient[l] = (currentValue -
                                                   domainDescription(testPoint +
                                                                     offset)) *
                                                  1e7;
                                    offset[l] = 0;
                                }
                                hasMoved.clear();
                                for (std::size_t l = 0;
                                     l < node->getNumberOfElements(); ++l) {
                                    if (!hasMoved[node->getElement(l)
                                                      ->getPhysicalGeometry()
                                                      ->getNodeIndex(
                                                          node->getNodeNumber(
                                                              l))]) {
                                        Geometry::PointPhysical<DIM> &other =
                                            theMesh_.getNodeCoordinates()
                                                [node->getElement(l)
                                                     ->getPhysicalGeometry()
                                                     ->getNodeIndex(
                                                         node->getNodeNumber(
                                                             l))];
                                        other += currentValue * gradient /
                                                 L2Norm(gradient);
                                        hasMoved[node->getElement(l)
                                                     ->getPhysicalGeometry()
                                                     ->getNodeIndex(
                                                         node->getNodeNumber(
                                                             l))] = true;
                                    }
                                }
                                std::get<1>(*moveIterator) +=
                                    10 * currentValue * gradient /
                                    L2Norm(gradient);
                            }
                        }
                    }
                    for (std::size_t j = 0; j < node->getNumberOfElements();
                         ++j) {
                        testPoint = node->getElement(j)
                                        ->getPhysicalGeometry()
                                        ->getLocalNodeCoordinates(
                                            node->getNodeNumber(j));
                        if (domainDescription(testPoint) > 1e-10) {
                            logger(WARN,
                                   "NOTE: Failed to move periodic testPoint % "
                                   "(%) back to the periodic boundary.\n "
                                   "Distance from boundary is %. Algorithm may "
                                   "crash.\n "
                                   "Consider fixing points at corners to "
                                   "remedy this issue.",
                                   i, testPoint, domainDescription(testPoint));
                        }
                    };
                }
            }
        }

        worstQuality = 1;
        if (DIM == 1) {
            // quality measure is not an issue in 1D just create a mesh with the
            // proper lengths
            worstQuality = 0.5;
            // the algorithm is mostly dimension independent, but the data type
            // it operates on is not
            for (Element *element :
                 theMesh_.getElementsList(IteratorType::GLOBAL)) {
                Geometry::PointPhysical<DIM> firstNode, secondNode;
                firstNode =
                    element->getPhysicalGeometry()->getLocalNodeCoordinates(0);
                secondNode =
                    element->getPhysicalGeometry()->getLocalNodeCoordinates(1);
                maxMovement =
                    std::max(maxMovement,
                             L2Norm(movement[element->getNode(0)->getID()]) /
                                 10 / currentLength[element->getID()]);
                maxMovement =
                    std::max(maxMovement,
                             L2Norm(movement[element->getNode(1)->getID()]) /
                                 10 / currentLength[element->getID()]);
                maxShift = std::max(
                    maxShift,
                    L2Norm(firstNode -
                           oldNodeLocations_[element->getPhysicalGeometry()
                                                 ->getNodeIndex(0)]) /
                        currentLength[element->getID()]);
                maxShift = std::max(
                    maxShift,
                    L2Norm(secondNode -
                           oldNodeLocations_[element->getPhysicalGeometry()
                                                 ->getNodeIndex(1)]) /
                        currentLength[element->getID()]);
            }
        } else if (DIM == 2) {
            // ratio between incircle and circumcircle (scaled so equilateral is
            // quality 1 and reference is quality ~.8)
            for (Element *element :
                 theMesh_.getElementsList(IteratorType::GLOBAL)) {
                std::array<double, 3> edgeLengths{};
                for (std::size_t i = 0; i < 3; ++i) {
                    edgeLengths[i] =
                        currentLength[element->getFace(i)->getID()];
                }
                worstQuality = std::min(
                    worstQuality,
                    (edgeLengths[0] + edgeLengths[1] - edgeLengths[2]) *
                        (edgeLengths[1] + edgeLengths[2] - edgeLengths[0]) *
                        (edgeLengths[2] + edgeLengths[0] - edgeLengths[1]) /
                        edgeLengths[0] / edgeLengths[1] / edgeLengths[2]);
            }
            for (Face *face : theMesh_.getFacesList(IteratorType::GLOBAL)) {
                Geometry::PointPhysical<DIM> firstNode, secondNode;
                std::vector<std::size_t> nodeIndices =
                    face->getPtrElementLeft()
                        ->getReferenceGeometry()
                        ->getCodim1EntityLocalIndices(
                            face->localFaceNumberLeft());
                firstNode = face->getPtrElementLeft()
                                ->getPhysicalGeometry()
                                ->getLocalNodeCoordinates(nodeIndices[0]);
                secondNode = face->getPtrElementLeft()
                                 ->getPhysicalGeometry()
                                 ->getLocalNodeCoordinates(nodeIndices[1]);
                maxMovement = std::max(
                    maxMovement, L2Norm(movement[face->getPtrElementLeft()
                                                     ->getNode(nodeIndices[0])
                                                     ->getID()]) /
                                     10 / currentLength[face->getID()]);
                maxMovement = std::max(
                    maxMovement, L2Norm(movement[face->getPtrElementLeft()
                                                     ->getNode(nodeIndices[1])
                                                     ->getID()]) /
                                     10 / currentLength[face->getID()]);
                maxShift = std::max(
                    maxShift,
                    L2Norm(
                        firstNode -
                        oldNodeLocations_[face->getPtrElementLeft()
                                              ->getPhysicalGeometry()
                                              ->getNodeIndex(nodeIndices[0])]) /
                        currentLength[face->getID()]);
                maxShift = std::max(
                    maxShift,
                    L2Norm(
                        secondNode -
                        oldNodeLocations_[face->getPtrElementLeft()
                                              ->getPhysicalGeometry()
                                              ->getNodeIndex(nodeIndices[1])]) /
                        currentLength[face->getID()]);
            }
        } else {
            // ratio between volume and cubed average edge length (scaled so
            // equilateral is quality 1 and reference is quality ~.8)
            for (Element *element :
                 theMesh_.getElementsList(IteratorType::GLOBAL)) {
                std::array<double, 6> edgeLengths{};
                for (std::size_t i = 0; i < 6; ++i) {
                    edgeLengths[i] =
                        currentLength[element->getEdge(i)->getID()];
                }
                double average =
                    std::accumulate(edgeLengths.begin(), edgeLengths.end(), 0) /
                    6;
                const Geometry::PointReference<DIM> &center =
                    element->getReferenceGeometry()->getCenter();
                Geometry::Jacobian<DIM, DIM> jac =
                    element->calcJacobian(center);
                worstQuality = std::min(
                    worstQuality, jac.determinant() / average * std::sqrt(2));
            }
            for (Edge *edge : theMesh_.getEdgesList(IteratorType::GLOBAL)) {
                Geometry::PointPhysical<DIM> firstNode, secondNode;
                std::vector<std::size_t> nodeIndices =
                    edge->getElement(0)
                        ->getReferenceGeometry()
                        ->getCodim2EntityLocalIndices(edge->getEdgeNumber(0));
                firstNode = edge->getElement(0)
                                ->getPhysicalGeometry()
                                ->getLocalNodeCoordinates(nodeIndices[0]);
                secondNode = edge->getElement(0)
                                 ->getPhysicalGeometry()
                                 ->getLocalNodeCoordinates(nodeIndices[1]);
                maxMovement = std::max(
                    maxMovement, L2Norm(movement[edge->getElement(0)
                                                     ->getNode(nodeIndices[0])
                                                     ->getID()]) /
                                     10 / currentLength[edge->getID()]);
                maxMovement = std::max(
                    maxMovement, L2Norm(movement[edge->getElement(0)
                                                     ->getNode(nodeIndices[1])
                                                     ->getID()]) /
                                     10 / currentLength[edge->getID()]);
                maxShift = std::max(
                    maxShift,
                    L2Norm(
                        firstNode -
                        oldNodeLocations_[edge->getElement(0)
                                              ->getPhysicalGeometry()
                                              ->getNodeIndex(nodeIndices[0])]) /
                        currentLength[edge->getID()]);
                maxShift = std::max(
                    maxShift,
                    L2Norm(
                        secondNode -
                        oldNodeLocations_[edge->getElement(0)
                                              ->getPhysicalGeometry()
                                              ->getNodeIndex(nodeIndices[1])]) /
                        currentLength[edge->getID()]);
            }
        }

        // no teleporting nodes in the final iteration
        ///\todo temporarily toggled off for debug reasons
        if (counter % 50 == 1 && false) {
            // the actual sorting is more expensive than computing the lengths
            // and this does not happen very often
            std::multimap<double,
                          std::pair<Geometry::PointPhysical<DIM>, std::size_t>>
                centerPoints{};
            if (DIM == 1) {
                // the algorithm is mostly dimension independent, but the data
                // type it operates on is not
                for (Element *element :
                     theMesh_.getElementsList(IteratorType::GLOBAL)) {
                    Geometry::PointPhysical<DIM> firstNode, secondNode;
                    firstNode =
                        element->getPhysicalGeometry()->getLocalNodeCoordinates(
                            0);
                    secondNode =
                        element->getPhysicalGeometry()->getLocalNodeCoordinates(
                            1);
                    // length is scaled in case somebody hasty decides to add
                    // smoothing at this point all edges should be squeezed a
                    // little if the algorithm is to work correctly so pretend
                    // the volume is 1.5 times as large remember to scale back
                    // from a volume measure to a length measure
                    double length =
                        (expectedLength[element->getNode(0)->getID()] +
                         expectedLength[element->getNode(1)->getID()]) /
                        currentLength[centerPoints.size()] * 2 *
                        totalCurrentLength / totalexpectedLength;
                    centerPoints.insert({length,
                                         {(firstNode + secondNode) / 2,
                                          element->getNode(0)->getID()}});
                }
            } else if (DIM == 2) {
                for (Face *face : theMesh_.getFacesList(IteratorType::GLOBAL)) {
                    Geometry::PointPhysical<DIM> firstNode, secondNode;
                    std::vector<std::size_t> nodeIndices =
                        face->getPtrElementLeft()
                            ->getReferenceGeometry()
                            ->getCodim1EntityLocalIndices(
                                face->localFaceNumberLeft());
                    firstNode = face->getPtrElementLeft()
                                    ->getPhysicalGeometry()
                                    ->getLocalNodeCoordinates(nodeIndices[0]);
                    secondNode = face->getPtrElementLeft()
                                     ->getPhysicalGeometry()
                                     ->getLocalNodeCoordinates(nodeIndices[1]);
                    // length is scaled in case somebody hasty decides to add
                    // smoothing at this point all edges should be squeezed a
                    // little if the algorithm is to work correctly so pretend
                    // the volume is 1.5 times as large remember to scale back
                    // from a volume measure to a length measure
                    double length =
                        (expectedLength[face->getPtrElementLeft()
                                            ->getNode(nodeIndices[0])
                                            ->getID()] +
                         expectedLength[face->getPtrElementLeft()
                                            ->getNode(nodeIndices[1])
                                            ->getID()]) /
                        currentLength[centerPoints.size()] *
                        std::pow(2 * totalCurrentLength / totalexpectedLength,
                                 1. / 2.);
                    centerPoints.insert({length,
                                         {(firstNode + secondNode) / 2,
                                          face->getPtrElementLeft()
                                              ->getNode(nodeIndices[0])
                                              ->getID()}});
                }
            } else {
                for (Edge *edge : theMesh_.getEdgesList(IteratorType::GLOBAL)) {
                    Geometry::PointPhysical<DIM> firstNode, secondNode;
                    std::vector<std::size_t> nodeIndices =
                        edge->getElement(0)
                            ->getReferenceGeometry()
                            ->getCodim2EntityLocalIndices(
                                edge->getEdgeNumber(0));
                    firstNode = edge->getElement(0)
                                    ->getPhysicalGeometry()
                                    ->getLocalNodeCoordinates(nodeIndices[0]);
                    secondNode = edge->getElement(0)
                                     ->getPhysicalGeometry()
                                     ->getLocalNodeCoordinates(nodeIndices[1]);
                    // length is scaled in case somebody hasty decides to add
                    // smoothing at this point all edges should be squeezed a
                    // little if the algorithm is to work correctly so pretend
                    // the volume is 1.5 times as large remember to scale back
                    // from a volume measure to a length measure
                    double length =
                        (expectedLength[edge->getElement(0)
                                            ->getNode(nodeIndices[0])
                                            ->getID()] +
                         expectedLength[edge->getElement(0)
                                            ->getNode(nodeIndices[1])
                                            ->getID()]) /
                        currentLength[centerPoints.size()] *
                        std::pow(2 * totalCurrentLength / totalexpectedLength,
                                 1. / 3.);
                    centerPoints.insert({length,
                                         {(firstNode + secondNode) / 2,
                                          edge->getElement(0)
                                              ->getNode(nodeIndices[0])
                                              ->getID()}});
                }
            }
            std::vector<bool> hasTeleported(
                theMesh_.getNumberOfNodeCoordinates(), false);
            auto longEdge = centerPoints.begin();
            auto shortEdge = centerPoints.rbegin();
            for (std::size_t index : fixedPointIdxs) {
                hasTeleported[index] = true;
            }
            Geometry::PointPhysical<DIM> point;
            Geometry::PointPhysical<DIM> other;
            for (Node *node : theMesh_.getNodesList(IteratorType::GLOBAL)) {
                point = node->getElement(0)
                            ->getPhysicalGeometry()
                            ->getLocalNodeCoordinates(node->getNodeNumber(0));
                for (std::size_t i = 0; i < node->getNumberOfElements(); ++i) {
                    other =
                        node->getElement(i)
                            ->getPhysicalGeometry()
                            ->getLocalNodeCoordinates(node->getNodeNumber(i));
                    if (point != other) {
                        hasTeleported[node->getID()] = true;
                    }
                }
            }
            // remember that the size measure is inverted
            while (3 * longEdge->first < shortEdge->first) {
                if (hasTeleported[shortEdge->second.second]) {
                    shortEdge++;
                } else {
                    if (domainDescription(longEdge->second.first) < 0) {
                        maxMovement = std::max(
                            maxMovement,
                            L2Norm(
                                longEdge->second.first -
                                theMesh_.getNodeCoordinates()[shortEdge->second
                                                                  .second]));
                        // it is quite unlikely that the current triangulation
                        // suffices after randomly teleporting nodes about
                        maxShift = std::numeric_limits<double>::infinity();
                        theMesh_
                            .getNodeCoordinates()[shortEdge->second.second] =
                            longEdge->second.first;
                        hasTeleported[shortEdge->second.second] = true;
                        shortEdge++;
                    }
                    longEdge++;
                }
            }
        }
    }
    if (counter == 10000) {
        logger(WARN,
               "WARNING: Maximum iteration count reached, mesh quality may not "
               "be optimal");
    }
    // coordinate transformation may have changed, update to the current
    // situation
    for (Element *element : theMesh_.getElementsList()) {
        element->getReferenceToPhysicalMap()->reinit();
    }
}
#endif

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

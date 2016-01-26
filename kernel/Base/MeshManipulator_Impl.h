/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, University of Twente
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
//QHull uses assert internally, but the macro definition causes conflicts with the rest of hpGEM
#undef assert
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
#include "AssembleBasisFunctionSet.h"
#include "OrientedBasisFunctionSet.h"
#include "Geometry/PointPhysical.h"
#include "Geometry/ReferenceSquare.h"
#include "Geometry/GlobalNamespaceGeometry.h"
#include "Geometry/PointReference.h"
#include "ElementCacheData.h"
#include "FaceCacheData.h"
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
#include "Utilities/BasisFunctions3DH1ConformingCube.h"
#include "Utilities/BasisFunctions3DH1ConformingPrism.h"
#include "Utilities/BasisFunctions3DH1ConformingTetrahedron.h"
#include "Utilities/BasisFunctions3DNedelec.h"
#include "Utilities/BasisFunctions3DAinsworthCoyle.h"
#include "Logger.h"

#include <algorithm>
#include <iostream>
#include <unordered_set>
#include <array>
#include <vector>
#include <numeric>
#include <type_traits>
#include <typeinfo>

//(crude) fix for pre-c++14 limitations of std::hash: just use the std::hash of the underlying type
template<typename T>
class EnumHash
{
    //cause a compile error if someone uses this for non-enum types
    using onlyForEnums = typename std::enable_if<std::is_enum<T>::value, T>::type;
public:
    std::size_t operator()(const T& t) const {
        return std::hash<typename std::underlying_type<T>::type>()(static_cast<typename std::underlying_type<T>::type>(t));
    }
};

namespace Base
{
    template<std::size_t DIM>
    void MeshManipulator<DIM>::createDefaultBasisFunctions(std::size_t order)
    {
        Base::BasisFunctionSet* bFset1 = new Base::BasisFunctionSet(order);
        switch (configData_->dimension_)
        {
            case 1:
                switch (order)
                {
                    case 0:
                        Base::AssembleBasisFunctionSet_1D_Ord0_A0(*bFset1);
                        break;
                    case 1:
                        Base::AssembleBasisFunctionSet_1D_Ord1_A0(*bFset1);
                        break;
                    case 2:
                        Base::AssembleBasisFunctionSet_1D_Ord2_A0(*bFset1);
                        break;
                    case 3:
                        Base::AssembleBasisFunctionSet_1D_Ord3_A0(*bFset1);
                        break;
                    case 4:
                        Base::AssembleBasisFunctionSet_1D_Ord4_A0(*bFset1);
                        break;
                    case 5:
                        Base::AssembleBasisFunctionSet_1D_Ord5_A0(*bFset1);
                        break;
                    default:
                        logger(WARN, "WARNING: No default basisFunction sets have been defined for this polynomial order; defaulting to 2");
                        const_cast<Base::ConfigurationData*>(configData_)->polynomialOrder_ = 2;
                        delete bFset1;
                        bFset1 = new Base::BasisFunctionSet(2);
                        Base::AssembleBasisFunctionSet_1D_Ord2_A0(*bFset1);
                }
                break;
            case 2:
                switch (order)
                {
                    case 0:
                        Base::AssembleBasisFunctionSet_2D_Ord0_A0(*bFset1);
                        break;
                    case 1:
                        Base::AssembleBasisFunctionSet_2D_Ord1_A0(*bFset1);
                        break;
                    case 2:
                        Base::AssembleBasisFunctionSet_2D_Ord2_A1(*bFset1);
                        break;
                    case 3:
                        Base::AssembleBasisFunctionSet_2D_Ord3_A1(*bFset1);
                        break;
                    case 4:
                        Base::AssembleBasisFunctionSet_2D_Ord4_A1(*bFset1);
                        break;
                    case 5:
                        Base::AssembleBasisFunctionSet_2D_Ord5_A1(*bFset1);
                        break;
                    default:
                        logger(WARN, "WARNING: No default basisFunction sets have been defined for this polynomial order; defaulting to 2");
                        const_cast<Base::ConfigurationData*>(configData_)->polynomialOrder_ = 2;
                        delete bFset1;
                        bFset1 = new Base::BasisFunctionSet(2);
                        Base::AssembleBasisFunctionSet_2D_Ord2_A1(*bFset1);
                }
                break;
            case 3:
                switch (order)
                {
                    case 0:
                        Base::AssembleBasisFunctionSet_3D_Ord0_A0(*bFset1);
                        break;
                    case 1:
                        Base::AssembleBasisFunctionSet_3D_Ord1_A0(*bFset1);
                        break;
                    case 2:
                        Base::AssembleBasisFunctionSet_3D_Ord2_A1(*bFset1);
                        break;
                    case 3:
                        Base::AssembleBasisFunctionSet_3D_Ord3_A1(*bFset1);
                        break;
                    case 4:
                        Base::AssembleBasisFunctionSet_3D_Ord4_A1(*bFset1);
                        break;
                    case 5:
                        Base::AssembleBasisFunctionSet_3D_Ord5_A1(*bFset1);
                        break;
                    default:
                        logger(WARN, "WARNING: No default basisFunction sets have been defined for this polynomial order; defaulting to 2");
                        const_cast<Base::ConfigurationData*>(configData_)->polynomialOrder_ = 2;
                        delete bFset1;
                        bFset1 = new Base::BasisFunctionSet(2);
                        Base::AssembleBasisFunctionSet_3D_Ord2_A1(*bFset1);
                }
                break;
            default:
                logger(ERROR, "No basisfunctions exist in this dimension");
        }
        if (collBasisFSet_.size() == 0)
        {
            collBasisFSet_.resize(1);
        }
        collBasisFSet_[0] = std::shared_ptr<const BasisFunctionSet>(bFset1);
    }

    template<std::size_t DIM>
    void MeshManipulator<DIM>::useDefaultDGBasisFunctions()
    {
        collBasisFSet_.clear();
        std::unordered_map<Geometry::ReferenceGeometryType, std::size_t, EnumHash<Geometry::ReferenceGeometryType> > shapeToIndex;
        for(Element* element : getElementsList(IteratorType::GLOBAL))
        {
            try
            {
                element->setDefaultBasisFunctionSet(shapeToIndex.at(element->getReferenceGeometry()->getGeometryType()));
            }
            catch(std::out_of_range&)
            {
                switch(element->getReferenceGeometry()->getGeometryType())
                {
                    case Geometry::ReferenceGeometryType::LINE:
                        shapeToIndex[Geometry::ReferenceGeometryType::LINE] = collBasisFSet_.size();
                        collBasisFSet_.emplace_back(Utilities::createDGBasisFunctionSet1DH1Line(configData_->polynomialOrder_));
                        break;
                    case Geometry::ReferenceGeometryType::SQUARE:
                        shapeToIndex[Geometry::ReferenceGeometryType::SQUARE] = collBasisFSet_.size();
                        collBasisFSet_.emplace_back(Utilities::createDGBasisFunctionSet2DH1Square(configData_->polynomialOrder_));
                        break;
                    case Geometry::ReferenceGeometryType::TRIANGLE:
                        shapeToIndex[Geometry::ReferenceGeometryType::TRIANGLE] = collBasisFSet_.size();
                        collBasisFSet_.emplace_back(Utilities::createDGBasisFunctionSet2DH1Triangle(configData_->polynomialOrder_));
                        break;
                    case Geometry::ReferenceGeometryType::CUBE:
                        shapeToIndex[Geometry::ReferenceGeometryType::CUBE] = collBasisFSet_.size();
                        collBasisFSet_.emplace_back(Utilities::createDGBasisFunctionSet3DH1Cube(configData_->polynomialOrder_));
                        break;
                    case Geometry::ReferenceGeometryType::TETRAHEDRON:
                        shapeToIndex[Geometry::ReferenceGeometryType::TETRAHEDRON] = collBasisFSet_.size();
                        collBasisFSet_.emplace_back(Utilities::createDGBasisFunctionSet3DH1Tetrahedron(configData_->polynomialOrder_));
                        break;
                    case Geometry::ReferenceGeometryType::TRIANGULARPRISM:
                        shapeToIndex[Geometry::ReferenceGeometryType::TRIANGULARPRISM] = collBasisFSet_.size();
                        collBasisFSet_.emplace_back(Utilities::createDGBasisFunctionSet3DH1ConformingPrism(configData_->polynomialOrder_));
                        break;
                    case Geometry::ReferenceGeometryType::PYRAMID:
                    case Geometry::ReferenceGeometryType::HYPERCUBE:
                        logger(ERROR, "No well-conditioned basis functions have been implemented for %s", element->getReferenceGeometry()->getName());
                        break;
                    case Geometry::ReferenceGeometryType::POINT:
                        logger(ERROR, "A point is not a valid geometry for an Element!");
                        break;
                    default:
                        logger(ERROR, "A new geometry has been implemented, please add it to the cases in MeshManipulator::useDefaultDGBasisFunctions and MeshManipulator::useDefaultConformingBasisFunctions");

                }
                element->setDefaultBasisFunctionSet(shapeToIndex.at(element->getReferenceGeometry()->getGeometryType()));
            }
        }
    }
    
    template<std::size_t DIM>
    void MeshManipulator<DIM>::useNedelecDGBasisFunctions()
    {
        collBasisFSet_.clear();
        std::unordered_map<Geometry::ReferenceGeometryType, std::size_t, EnumHash<Geometry::ReferenceGeometryType> > shapeToIndex;
        for(Element* element : getElementsList(IteratorType::GLOBAL))
        {
            try
            {
                element->setDefaultBasisFunctionSet(shapeToIndex.at(element->getReferenceGeometry()->getGeometryType()));
            }
            catch(std::out_of_range&)
            {
                switch(element->getReferenceGeometry()->getGeometryType())
                {
                    case Geometry::ReferenceGeometryType::TETRAHEDRON:
                        shapeToIndex[Geometry::ReferenceGeometryType::TETRAHEDRON] = collBasisFSet_.size();
                        collBasisFSet_.emplace_back(Utilities::createDGBasisFunctionSet3DNedelec(configData_->polynomialOrder_));
                        break;
                    default:
                        logger(ERROR, "No Nedelec basis functions have been implemented for %s", element->getReferenceGeometry()->getName());

                }
                element->setDefaultBasisFunctionSet(shapeToIndex.at(element->getReferenceGeometry()->getGeometryType()));
            }
        }
    }
    
    template<std::size_t DIM>
    void MeshManipulator<DIM>::useAinsworthCoyleDGBasisFunctions()
    {
        collBasisFSet_.clear();
        std::unordered_map<Geometry::ReferenceGeometryType, std::size_t, EnumHash<Geometry::ReferenceGeometryType> > shapeToIndex;
        for(Element* element : getElementsList(IteratorType::GLOBAL))
        {
            try
            {
                element->setDefaultBasisFunctionSet(shapeToIndex.at(element->getReferenceGeometry()->getGeometryType()));
            }
            catch(std::out_of_range&)
            {
                switch(element->getReferenceGeometry()->getGeometryType())
                {
                    case Geometry::ReferenceGeometryType::TETRAHEDRON:
                        shapeToIndex[Geometry::ReferenceGeometryType::TETRAHEDRON] = collBasisFSet_.size();
                        collBasisFSet_.emplace_back(Utilities::createDGBasisFunctionSet3DAinsworthCoyle(configData_->polynomialOrder_));
                        break;
                    default:
                        logger(ERROR, "No Nedelec basis functions have been implemented for %s", element->getReferenceGeometry()->getName());

                }
                element->setDefaultBasisFunctionSet(shapeToIndex.at(element->getReferenceGeometry()->getGeometryType()));
            }
        }
    }

    template<std::size_t DIM>
    void MeshManipulator<DIM>::useDefaultConformingBasisFunctions()
    {
        logger.assert(configData_->polynomialOrder_ > 0, "Basis function may not have an empty union of supporting elements. Use a DG basis function on a single element non-periodic mesh instead");
        collBasisFSet_.clear();
        std::unordered_map<Geometry::ReferenceGeometryType, std::size_t, EnumHash<Geometry::ReferenceGeometryType> > shapeToElementIndex;
        std::unordered_map<Geometry::ReferenceGeometryType, std::size_t, EnumHash<Geometry::ReferenceGeometryType> > numberOfFaceSets;
        std::unordered_map<Geometry::ReferenceGeometryType, std::size_t, EnumHash<Geometry::ReferenceGeometryType> > numberOfEdgeSets;
        std::unordered_map<Geometry::ReferenceGeometryType, std::size_t, EnumHash<Geometry::ReferenceGeometryType> > numberOfNodeSets;
        for(Element* element : getElementsList(IteratorType::GLOBAL))
        {
            try
            {
                element->setDefaultBasisFunctionSet(shapeToElementIndex.at(element->getReferenceGeometry()->getGeometryType()));
            }
            //there is more relevant code after the huge catch block
            catch(std::out_of_range&)
            {
                auto type = element->getReferenceGeometry()->getGeometryType();
                std::vector<const Base::BasisFunctionSet*> nodeSet;
                std::vector<const Base::OrientedBasisFunctionSet*> faceSet;
                std::vector<const Base::OrientedBasisFunctionSet*> edgeSet;
                switch(type)
                {
                    case Geometry::ReferenceGeometryType::LINE:
                        shapeToElementIndex[type] = collBasisFSet_.size();
                        collBasisFSet_.emplace_back(Utilities::createInteriorBasisFunctionSet1DH1Line(configData_->polynomialOrder_));
                        faceSet = Utilities::createVertexBasisFunctionSet1DH1Line(configData_->polynomialOrder_);
                        for(const Base::BasisFunctionSet* set : faceSet)
                        {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfFaceSets[type] = collBasisFSet_.size() - shapeToElementIndex[type] - 1;
                        numberOfEdgeSets[type] = 0;
                        numberOfNodeSets[type] = 0;
                        break;
                    case Geometry::ReferenceGeometryType::SQUARE:
                        shapeToElementIndex[type] = collBasisFSet_.size();
                        collBasisFSet_.emplace_back(Utilities::createInteriorBasisFunctionSet2DH1Square(configData_->polynomialOrder_));
                        faceSet = Utilities::createFaceBasisFunctionSet2DH1Square(configData_->polynomialOrder_);
                        for(const Base::BasisFunctionSet* set : faceSet)
                        {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfFaceSets[type] = collBasisFSet_.size() - shapeToElementIndex[type] - 1;
                        numberOfEdgeSets[type] = 0;
                        nodeSet = Utilities::createVertexBasisFunctionSet2DH1Square(configData_->polynomialOrder_);
                        for(const Base::BasisFunctionSet* set : nodeSet)
                        {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfNodeSets[type] = collBasisFSet_.size() - shapeToElementIndex[type] - numberOfFaceSets[type] - numberOfEdgeSets[type] - 1;
                        break;
                    case Geometry::ReferenceGeometryType::TRIANGLE:
                        shapeToElementIndex[type] = collBasisFSet_.size();
                        collBasisFSet_.emplace_back(Utilities::createInteriorBasisFunctionSet2DH1Triangle(configData_->polynomialOrder_));
                        faceSet = Utilities::createFaceBasisFunctionSet2DH1Triangle(configData_->polynomialOrder_);
                        for(const Base::BasisFunctionSet* set : faceSet)
                        {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfFaceSets[type] = collBasisFSet_.size() - shapeToElementIndex[type] - 1;
                        numberOfEdgeSets[type] = 0;
                        nodeSet = Utilities::createVertexBasisFunctionSet2DH1Triangle(configData_->polynomialOrder_);
                        for(const Base::BasisFunctionSet* set : nodeSet)
                        {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfNodeSets[type] = collBasisFSet_.size() - shapeToElementIndex[type] - numberOfFaceSets[type] - numberOfEdgeSets[type] - 1;
                        break;
                    case Geometry::ReferenceGeometryType::CUBE:
                        shapeToElementIndex[type] = collBasisFSet_.size();
                        collBasisFSet_.emplace_back(Utilities::createInteriorBasisFunctionSet3DH1Cube(configData_->polynomialOrder_));
                        faceSet = Utilities::createFaceBasisFunctionSet3DH1Cube(configData_->polynomialOrder_);
                        for(const Base::BasisFunctionSet* set : faceSet)
                        {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfFaceSets[type] = collBasisFSet_.size() - shapeToElementIndex[type] - 1;
                        edgeSet = Utilities::createEdgeBasisFunctionSet3DH1Cube(configData_->polynomialOrder_);
                        for(const Base::BasisFunctionSet* set : edgeSet)
                        {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfEdgeSets[type] = collBasisFSet_.size() - shapeToElementIndex[type] - numberOfFaceSets[type] - 1;
                        nodeSet = Utilities::createVertexBasisFunctionSet3DH1Cube(configData_->polynomialOrder_);
                        for(const Base::BasisFunctionSet* set : nodeSet)
                        {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfNodeSets[type] = collBasisFSet_.size() - shapeToElementIndex[type] - numberOfFaceSets[type] - numberOfEdgeSets[type] - 1;
                        break;
                    case Geometry::ReferenceGeometryType::TETRAHEDRON:
                        shapeToElementIndex[type] = collBasisFSet_.size();
                        collBasisFSet_.emplace_back(Utilities::createInteriorBasisFunctionSet3DH1Tetrahedron(configData_->polynomialOrder_));
                        faceSet = Utilities::createFaceBasisFunctionSet3DH1Tetrahedron(configData_->polynomialOrder_);
                        for(const Base::BasisFunctionSet* set : faceSet)
                        {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfFaceSets[type] = collBasisFSet_.size() - shapeToElementIndex[type] - 1;
                        edgeSet = Utilities::createEdgeBasisFunctionSet3DH1Tetrahedron(configData_->polynomialOrder_);
                        for(const Base::BasisFunctionSet* set : edgeSet)
                        {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfEdgeSets[type] = collBasisFSet_.size() - shapeToElementIndex[type] - numberOfFaceSets[type] - 1;
                        nodeSet = Utilities::createVertexBasisFunctionSet3DH1Tetrahedron(configData_->polynomialOrder_);
                        for(const Base::BasisFunctionSet* set : nodeSet)
                        {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfNodeSets[type] = collBasisFSet_.size() - shapeToElementIndex[type] - numberOfFaceSets[type] - numberOfEdgeSets[type] - 1;
                        break;
                    case Geometry::ReferenceGeometryType::TRIANGULARPRISM:
                        shapeToElementIndex[type] = collBasisFSet_.size();
                        collBasisFSet_.emplace_back(Utilities::createInteriorBasisFunctionSet3DH1ConformingPrism(configData_->polynomialOrder_));
                        faceSet = Utilities::createFaceBasisFunctionSet3DH1ConformingPrism(configData_->polynomialOrder_);
                        for(const Base::BasisFunctionSet* set : faceSet)
                        {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfFaceSets[type] = collBasisFSet_.size() - shapeToElementIndex[type] - 1;
                        edgeSet = Utilities::createEdgeBasisFunctionSet3DH1ConformingPrism(configData_->polynomialOrder_);
                        for(const Base::BasisFunctionSet* set : edgeSet)
                        {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfEdgeSets[type] = collBasisFSet_.size() - shapeToElementIndex[type] - numberOfFaceSets[type] - 1;
                        nodeSet = Utilities::createVertexBasisFunctionSet3DH1ConformingPrism(configData_->polynomialOrder_);
                        for(const Base::BasisFunctionSet* set : nodeSet)
                        {
                            collBasisFSet_.emplace_back(set);
                        }
                        numberOfNodeSets[type] = collBasisFSet_.size() - shapeToElementIndex[type] - numberOfFaceSets[type] - numberOfEdgeSets[type] - 1;
                        break;
                    case Geometry::ReferenceGeometryType::PYRAMID:
                    case Geometry::ReferenceGeometryType::HYPERCUBE:
                        logger(ERROR, "No well-conditioned basis functions have been implemented for %s", element->getReferenceGeometry()->getName());
                        break;
                    case Geometry::ReferenceGeometryType::POINT:
                        logger(ERROR, "A point is not a valid geometry for an Element!");
                        break;
                    default:
                        logger(ERROR, "A new geometry has been implemented, please add it to the cases in MeshManipulator::useDefaultDGBasisFunctions and MeshManipulator::useDefaultConformingBasisFunctions");

                }
                element->setDefaultBasisFunctionSet(shapeToElementIndex.at(type));
            }
        }
        for(Face* face : getFacesList(IteratorType::GLOBAL))
        {
            std::size_t faceNumber = face->localFaceNumberLeft();
            auto type = face->getPtrElementLeft()->getReferenceGeometry()->getGeometryType();
            for (std::size_t i = shapeToElementIndex[type] + 1; i < shapeToElementIndex[type] + numberOfFaceSets[type] + 1; ++i)
            {
                logger.assert(typeid(*collBasisFSet_[i]) == typeid(const OrientedBasisFunctionSet), "This is not supposed to happen");
                if (static_cast<const OrientedBasisFunctionSet*>(collBasisFSet_[i].get())->checkOrientation(0, faceNumber))
                {
                    face->getPtrElementLeft()->setFaceBasisFunctionSet(i, faceNumber);
                    //the number of basis functions depends on the shape of the face, not on the shape of the element
                    face->setLocalNumberOfBasisFunctions(collBasisFSet_[i]->size());
                }
            }
            if (face->isInternal())
            {
                faceNumber = face->localFaceNumberRight();
                type = face->getPtrElementRight()->getReferenceGeometry()->getGeometryType();
                std::size_t orientation = face->getFaceToFaceMapIndex();
                for (std::size_t i = shapeToElementIndex[type] + 1; i < shapeToElementIndex[type] + numberOfFaceSets[type] + 1; ++i)
                {
                    logger.assert(typeid(*collBasisFSet_[i]) == typeid(const OrientedBasisFunctionSet), "This is not supposed to happen");
                    if (static_cast<const OrientedBasisFunctionSet*>(collBasisFSet_[i].get())->checkOrientation(orientation, faceNumber))
                    {
                        face->getPtrElementRight()->setFaceBasisFunctionSet(i, faceNumber);
                    }
                }
            }
        }
        for(Edge* edge : getEdgesList(IteratorType::GLOBAL))
        {
            for(std::size_t i = 0; i < edge->getNumberOfElements(); ++i)
            {
                Element* element = edge->getElement(i);
                std::size_t edgeNumber = edge->getEdgeNumber(i);
                std::size_t orientation = edge->getOrientation(i);
                auto type = element->getReferenceGeometry()->getGeometryType();
                for(std::size_t j = shapeToElementIndex[type] + numberOfFaceSets[type] + 1; j < shapeToElementIndex[type] + numberOfFaceSets[type] + numberOfEdgeSets[type] + 1; ++j)
                {
                    logger.assert(typeid(*collBasisFSet_[j]) == typeid(const OrientedBasisFunctionSet), "This is not supposed to happen");
                    if (static_cast<const OrientedBasisFunctionSet*>(collBasisFSet_[j].get())->checkOrientation(orientation, edgeNumber))
                    {
                        element->setEdgeBasisFunctionSet(j, edgeNumber);
                        edge->setLocalNumberOfBasisFunctions(collBasisFSet_[j]->size());
                    }

                }
            }
        }
        for(Node* node : getNodesList(IteratorType::GLOBAL))
        {
            if(DIM > 1)
            {
                for(std::size_t i = 0; i < node->getNumberOfElements(); ++i)
                {
                    Element* element = node->getElement(i);
                    std::size_t nodeNumber = node->getNodeNumber(i);
                    auto type = element->getReferenceGeometry()->getGeometryType();
                    element->setVertexBasisFunctionSet(shapeToElementIndex[type] + numberOfFaceSets[type] + numberOfEdgeSets[type] + 1 + nodeNumber, nodeNumber);
                    node->setLocalNumberOfBasisFunctions(collBasisFSet_[shapeToElementIndex[type] + numberOfFaceSets[type] + numberOfEdgeSets[type] + 1 + nodeNumber]->size());
                }
            }
        }
    }

    template<std::size_t DIM>
    MeshManipulator<DIM>::MeshManipulator(const ConfigurationData* config, BoundaryType xPer, BoundaryType yPer, BoundaryType zPer, std::size_t orderOfFEM, std::size_t idRangeBegin, std::size_t numberOfElementMatrices, std::size_t numberOfElementVectors, std::size_t numberOfFaceMatrices, std::size_t numberOfFaceVectors)
            : MeshManipulatorBase(config, xPer, yPer, zPer, orderOfFEM, idRangeBegin, numberOfElementMatrices, numberOfElementVectors, numberOfFaceMatrices, numberOfFaceVectors), meshMover_(nullptr)
    {
        logger.assert(DIM == config->dimension_, "Invalid configuration passed");
        logger.assert(orderOfFEM==config->polynomialOrder_, "Inconsistent redundant information passed");
        logger.assert(idRangeBegin==0, "c++ starts counting at 0");
        createDefaultBasisFunctions(orderOfFEM);
        const_cast<ConfigurationData*>(configData_)->numberOfBasisFunctions_ = collBasisFSet_[0]->size();
    }

    template<std::size_t DIM>
    MeshManipulator<DIM>::MeshManipulator(const MeshManipulator& other)
    : MeshManipulatorBase(other), theMesh_(other.theMesh_), meshMover_(other.meshMover_), collBasisFSet_(other.collBasisFSet_)
    {        
    }

    template<std::size_t DIM>
    MeshManipulator<DIM>::~MeshManipulator()
    {        
        delete meshMover_;
    }

    template<std::size_t DIM>
    void MeshManipulator<DIM>::setDefaultBasisFunctionSet(BasisFunctionSet* bFSet)
    {
        logger.assert(bFSet!=nullptr, "Invalid basis function set passed");
        collBasisFSet_[0] = std::shared_ptr<const BasisFunctionSet>(bFSet);
        const_cast<ConfigurationData*>(configData_)->numberOfBasisFunctions_ = bFSet->size();
        for (Base::Face* face : getFacesList(IteratorType::GLOBAL))
        {
            face->setLocalNumberOfBasisFunctions(0);
        }
        for (Base::Edge* edge : getEdgesList(IteratorType::GLOBAL))
        {
            edge->setLocalNumberOfBasisFunctions(0);
        }
        for (Base::Node* node : getNodesList(IteratorType::GLOBAL))
        {
            node->setLocalNumberOfBasisFunctions(0);
        }
        for (ElementIterator it = elementColBegin(IteratorType::GLOBAL); it != elementColEnd(IteratorType::GLOBAL); ++it)
        {
            (*it)->setDefaultBasisFunctionSet(0);
        }
    }

    template<std::size_t DIM>
    void MeshManipulator<DIM>::addVertexBasisFunctionSet(const std::vector<const BasisFunctionSet*>& bFsets)
    {
        std::size_t firstNewEntry = collBasisFSet_.size();
        for (const BasisFunctionSet* set : bFsets)
        {
            logger.assert(set!=nullptr, "Invalid basis function set detected");
            collBasisFSet_.emplace_back(set);
        }
        for (Node* node : getNodesList())
        {
            for (std::size_t i = 0; i < node->getNumberOfElements(); ++i)
            {
                node->getElement(i)->setVertexBasisFunctionSet(firstNewEntry + node->getNodeNumber(i), node->getNodeNumber(i));
            }
            node->setLocalNumberOfBasisFunctions(bFsets[0]->size());
        }
        const_cast<ConfigurationData*>(configData_)->numberOfBasisFunctions_ += (*elementColBegin())->getNumberOfNodes() * bFsets[0]->size();
    }

    template<std::size_t DIM>
    void MeshManipulator<DIM>::addFaceBasisFunctionSet(const std::vector<const OrientedBasisFunctionSet*>& bFsets)
    {
        std::size_t firstNewEntry = collBasisFSet_.size();
        for (const BasisFunctionSet* set : bFsets)
        {
            logger.assert(set!=nullptr, "Invalid basis function set detected");
            collBasisFSet_.emplace_back(set);
        }
        for (Face* face : getFacesList())
        {
            std::size_t faceNumber = face->localFaceNumberLeft();
            for (std::size_t i = 0; i < bFsets.size(); ++i)
            {
                if (bFsets[i]->checkOrientation(0, faceNumber))
                {
                    face->getPtrElementLeft()->setFaceBasisFunctionSet(firstNewEntry + i, faceNumber);
                }
            }
            if (face->isInternal())
            {
                faceNumber = face->localFaceNumberRight();
                std::size_t orientation = face->getFaceToFaceMapIndex();
                for (std::size_t i = 0; i < bFsets.size(); ++i)
                {
                    if (bFsets[i]->checkOrientation(orientation, faceNumber))
                    {
                        face->getPtrElementRight()->setFaceBasisFunctionSet(firstNewEntry + i, faceNumber);
                    }
                }
            }
            face->setLocalNumberOfBasisFunctions(bFsets[0]->size());
        }
        const_cast<ConfigurationData*>(configData_)->numberOfBasisFunctions_ += (*elementColBegin())->getPhysicalGeometry()->getNumberOfFaces() * bFsets[0]->size();
    }

    template<std::size_t DIM>
    void MeshManipulator<DIM>::addEdgeBasisFunctionSet(const std::vector<const OrientedBasisFunctionSet*>& bFsets)
    {
        std::size_t firstNewEntry = collBasisFSet_.size();
        for (const BasisFunctionSet* set : bFsets)
        {
            logger.assert(set!=nullptr, "Invalid basis function set detected");
            collBasisFSet_.emplace_back(set);
        }
        logger(DEBUG, "In MeshManipulator::addEdgeBasisFunctionSet: ");
        for (Edge* edge : getEdgesList())
        {
            for (std::size_t i = 0; i < edge->getNumberOfElements(); ++i)
            {
                for (std::size_t j = 0; j < bFsets.size(); ++j)
                {
                    if (bFsets[j]->checkOrientation(edge->getOrientation(i), edge->getEdgeNumber(i)))
                    {
                        edge->getElement(i)->setEdgeBasisFunctionSet(firstNewEntry + j, edge->getEdgeNumber(i));
                        logger(DEBUG, "% % %", edge->getOrientation(i), 
                               edge->getEdgeNumber(i), bFsets[j]->size());
                    }
                }
            }
            edge->setLocalNumberOfBasisFunctions(bFsets[0]->size());
        }
        const_cast<ConfigurationData*>(configData_)->numberOfBasisFunctions_ += (*elementColBegin())->getPhysicalGeometry()->getNumberOfFaces() * bFsets[0]->size();
    }

    template<std::size_t DIM>
    Base::Element*
    MeshManipulator<DIM>::addElement(const std::vector<std::size_t>& globalNodeIndexes)
    {
        logger.assert([&]()->bool{
            for(std::size_t i = 0; i < globalNodeIndexes.size(); ++i)
                for(std::size_t j = 0; j < i; ++j)
                    if(globalNodeIndexes[i] == globalNodeIndexes[j])
                        return false;
            return true;
        }(), "Trying to pass the same node twice");
        return theMesh_.addElement(globalNodeIndexes);
    }

    template<std::size_t DIM>
    void MeshManipulator<DIM>::move()
    {
        for (Geometry::PointPhysical<DIM>& p : theMesh_.getNodeCoordinates())
        {
            if (meshMover_ != nullptr)
            {
                meshMover_->movePoint(p);
            }
        }
    }

    template<std::size_t DIM>
    void MeshManipulator<DIM>::setMeshMover(const MeshMoverBase<DIM>* meshMover)
    {
        //can be set to nullptr if you don't want to move the mesh anymore
        meshMover_ = meshMover;
    }

    template<std::size_t DIM>
    bool MeshManipulator<DIM>::addFace(Element* leftElementPtr, std::size_t leftElementLocalFaceNo, Element* rightElementPtr, std::size_t rightElementLocalFaceNo, const Geometry::FaceType& faceType)
    {
        logger.assert(leftElementPtr!=nullptr, "Invalid element passed");
        //rightElementPtr may be nullptr for boundary faces
        return theMesh_.addFace(leftElementPtr, leftElementLocalFaceNo, rightElementPtr, rightElementLocalFaceNo, faceType);
    }

    template<std::size_t DIM>
    void MeshManipulator<DIM>::addEdge()
    {
        theMesh_.addEdge();
    }
    
    template<std::size_t DIM>
    void MeshManipulator<DIM>::addNode()
    {
        theMesh_.addNode();
    }

    template<std::size_t DIM>
    std::tuple<const Base::Element*, Geometry::PointReference<DIM>> MeshManipulator<DIM>::physicalToReference(Geometry::PointPhysical<DIM> pointPhysical) const
    {
        return physicalToReference_detail(pointPhysical, getElementsList(IteratorType::GLOBAL).getRootEntries());
    }

    template<std::size_t DIM>
    template<typename Iterable>
    std::tuple<const Base::Element*, Geometry::PointReference<DIM>> MeshManipulator<DIM>::physicalToReference_detail(Geometry::PointPhysical<DIM>pointPhysical, Iterable elementContainer) const
    {
        for(auto singleEntry : elementContainer)
        {
            const Base::Element* element = singleEntry->getData();
            Geometry::PointReference<DIM> pointReference = element->physicalToReference(pointPhysical);
            if(element->getReferenceGeometry()->isInternalPoint(pointReference))
            {
                if(singleEntry->hasChild())
                {
                    return physicalToReference_detail(pointPhysical, singleEntry->getChildren());
                }
                else
                {
                    return {element, pointReference};
                }
            }
        }
        logger(ERROR, "The point % lies outsize the domain", pointPhysical);
        return std::tuple<Base::Element*, Geometry::PointReference<DIM>>{nullptr, {}};
    }

    template<std::size_t DIM>
    void MeshManipulator<DIM>::addMeasurePoint(Geometry::PointPhysical<DIM> pointPhysical)
    {
        measurePoints_.push_back(physicalToReference(pointPhysical));
    }

    template<std::size_t DIM>
    const std::vector<std::tuple<const Base::Element*, Geometry::PointReference<DIM>>>& MeshManipulator<DIM>::getMeasurePoints() const
    {
        return measurePoints_;
    }
}

template<std::size_t DIM>
std::ostream& operator<<(std::ostream& os, const Base::MeshManipulator<DIM>& mesh)
{
    for (Geometry::PointPhysical<DIM> p : mesh.getNodeCoordinates())
    {
        os << "Node " << " " << p << std::endl;
    }

    for (Base::Element* element : mesh.getElementsList())
    {
        os << "Element " << element->getID() << " " << element << std::endl;
    }
    return os;
}
    
namespace Base
{

    template<std::size_t DIM>
    void MeshManipulator<DIM>::createRectangularMesh(const Geometry::PointPhysical<DIM>& bottomLeft, const Geometry::PointPhysical<DIM>& topRight, const std::vector<std::size_t>& linearNoElements)
    {
        logger.assert(bottomLeft.size()==topRight.size(), "The corners of the mesh must have the same dimension");
        logger.assert(bottomLeft.size()==configData_->dimension_, "The corners of the mesh have the wrong dimension");
        logger.assert(linearNoElements.size()==configData_->dimension_, "There are amounts of elements spicified in % dimensions, but there are % dimensions", linearNoElements.size(), configData_->dimension_);
        //set to correct value in case some other meshmanipulator changed things
        ElementFactory::instance().setCollectionOfBasisFunctionSets(&collBasisFSet_);
        ElementFactory::instance().setNumberOfMatrices(numberOfElementMatrices_);
        ElementFactory::instance().setNumberOfVectors(numberOfElementVectors_);
        ElementFactory::instance().setNumberOfTimeLevels(configData_->numberOfTimeLevels_);
        ElementFactory::instance().setNumberOfUnknowns(configData_->numberOfUnknowns_);
        FaceFactory::instance().setNumberOfFaceMatrices(numberOfFaceMatrices_);
        FaceFactory::instance().setNumberOfFaceVectors(numberOfFaceVectors_);
        
        logger.assert(linearNoElements.size() == DIM, "The number of Linear Intervals has to map the size of the problem and current it does not");
        std::vector<bool> periodicDIM;
        for (std::size_t i = 0; i < DIM; ++i)
        {
            if (i == 0)
                periodicDIM.push_back(periodicX_);
            if (i == 1)
                periodicDIM.push_back(periodicY_);
            if (i == 2)
                periodicDIM.push_back(periodicZ_);
        }
        //Stage 1 : Precompute some required values;
        ///////
        
        //This store the size length of the domain i.e. it is DIM sized vector
        Geometry::PointPhysical<DIM> delta_x;
        
        for (std::size_t i = 0; i < DIM; i++)
        {
            delta_x[i] = (topRight[i] - bottomLeft[i]) / (linearNoElements[i]);
        }
        
        //This stores the number of nodes in each coDIMension i.e. if you have 2 by 2 element it is 3 nodes 
        //nodeCoordinates mark physical location, notes mark connectivity-based location
        std::vector<std::size_t> numberOfNodeCoordinatesInEachSubspace(DIM), numberOfElementsInEachSubspace(DIM), numberOfNodesInEachSubspace(DIM);
        
        numberOfNodeCoordinatesInEachSubspace[0] = 1;
        numberOfNodesInEachSubspace[0] = 1;
        numberOfElementsInEachSubspace[0] = 1;
        
        //This will be the total number of nodes required in the problem
        std::size_t totalNumberOfNodeCoordinates, totalNumberOfNodes, totalNumberOfElements, numberOfNodesPerElement;
        
        totalNumberOfNodeCoordinates = (linearNoElements[0] + 1);
        totalNumberOfNodes = (linearNoElements[0] + (periodicDIM[0] ? 0 : 1));
        
        totalNumberOfElements = (linearNoElements[0]);
        
        numberOfNodesPerElement = 2;
        std::size_t powerOf2;
        
        for (std::size_t iDIM = 1; iDIM < DIM; ++iDIM)
        {
            totalNumberOfNodeCoordinates *= (linearNoElements[iDIM] + 1);
            totalNumberOfNodes *= (linearNoElements[iDIM] + (periodicDIM[iDIM] ? 0 : 1));
            totalNumberOfElements *= (linearNoElements[iDIM]);
            numberOfNodesPerElement *= 2;
            
            numberOfElementsInEachSubspace[iDIM] = numberOfElementsInEachSubspace[iDIM - 1] * (linearNoElements[iDIM - 1]);
            numberOfNodeCoordinatesInEachSubspace[iDIM] = numberOfNodeCoordinatesInEachSubspace[iDIM - 1] * (linearNoElements[iDIM - 1] + 1);
            numberOfNodesInEachSubspace[iDIM] = numberOfNodesInEachSubspace[iDIM - 1] * (linearNoElements[iDIM - 1] + (periodicDIM[iDIM - 1] ? 0 : 1));
        }
        
        //temp point for storing the node locations
        Geometry::PointPhysical<DIM> x;
        
        //Stage 2 : Create the nodes
        //Now loop over all the nodes and calculate the coordinates for reach DIMension (this makes the algorithm independent of DIMension
        for (std::size_t nodeIndex = 0; nodeIndex < totalNumberOfNodeCoordinates; ++nodeIndex)
        {
            std::size_t nodeIndexRemain = nodeIndex;
            
            for (int iDIM = DIM - 1; iDIM > -1; --iDIM)
            {
                x[iDIM] = bottomLeft[iDIM] + (nodeIndexRemain / numberOfNodeCoordinatesInEachSubspace[iDIM] * delta_x[iDIM]);
                nodeIndexRemain %= numberOfNodeCoordinatesInEachSubspace[iDIM];
            }
            
            //actually add the point
            theMesh_.addNodeCoordinate(x);
            
        }
        
        //stage 2.5 : create the vertices
        
        for (std::size_t nodeID = 0; nodeID < totalNumberOfNodes; ++nodeID)
        {
            theMesh_.addNode();
        }
        
        //Stage 3 : Create the elements
        
        std::vector<std::size_t> elementNdId(DIM), nodeCoordinateNdId(DIM), nodeNdId(DIM), globalNodeCoordinateID(numberOfNodesPerElement), globalNodeID(numberOfNodesPerElement);
        
        auto& nodesList = getNodesList(IteratorType::GLOBAL);
        
        //elementNdId is DIM coordinate of the bottom left node i.e. in two (0,0), (1,0) ,(2,0) ... etc are the first three (if at least three elements in x)
        for (std::size_t elementIndex = 0; elementIndex < totalNumberOfElements; ++elementIndex)
        {
            std::size_t numberOfElementsRemaining = elementIndex;
            
            for (int iDIM = DIM - 1; iDIM > -1; --iDIM)
            {
                elementNdId[iDIM] = numberOfElementsRemaining / numberOfElementsInEachSubspace[iDIM];
                numberOfElementsRemaining %= numberOfElementsInEachSubspace[iDIM];
            }
            
            // nodeNdId are the DIM coordinate of each node in the element with nodeNdId[0] being the bottom left
            for (std::size_t i = 0; i < numberOfNodesPerElement; ++i)
            {
                powerOf2 = 1;
                for (std::size_t iDIM = 0; iDIM < DIM; ++iDIM)
                {
                    nodeCoordinateNdId[iDIM] = elementNdId[iDIM] + ((i & powerOf2) != 0);
                    nodeNdId[iDIM] = elementNdId[iDIM] + ((i & powerOf2) != 0);
                    if ((nodeNdId[iDIM] >= linearNoElements[iDIM]) && (periodicDIM[iDIM] == true))
                    {
                        nodeNdId[iDIM] = 0;
                    }
                    powerOf2 *= 2;
                }
                globalNodeCoordinateID[i] = nodeCoordinateNdId[0];
                globalNodeID[i] = nodeNdId[0];
                
                //Now map to the one DIMensional global ID
                for (std::size_t iDIM = 1; iDIM < DIM; ++iDIM)
                {
                    globalNodeCoordinateID[i] += nodeCoordinateNdId[iDIM] * numberOfNodeCoordinatesInEachSubspace[iDIM];
                    globalNodeID[i] += nodeNdId[iDIM] * numberOfNodesInEachSubspace[iDIM];
                }
            }
            Element* newElement = addElement(globalNodeCoordinateID);
            for (std::size_t i = 0; i < globalNodeID.size(); ++i)
            {
                nodesList[globalNodeID[i]]->addElement(newElement, i);
            }
        }
        
        faceFactory();
        edgeFactory();
        
    }
    
    //createTrianglularMesh follows the same structure as createRectangularMesh. 
    //Where createRectangularMesh makes rectangular elements, createTrianglularMesh
    //splits the elements into a partition of triangles.
    template<std::size_t DIM>
    void MeshManipulator<DIM>::createTriangularMesh(Geometry::PointPhysical<DIM> bottomLeft, Geometry::PointPhysical<DIM> topRight, const std::vector<std::size_t>& linearNoElements)
    {
        logger.assert(bottomLeft.size()==topRight.size(), "The corners of the mesh must have the same dimension");
        logger.assert(bottomLeft.size()==configData_->dimension_, "The corners of the mesh have the wrong dimension");
        logger.assert(linearNoElements.size()==configData_->dimension_, "There are amounts of elements specified in % dimensions, but there are % dimensions", linearNoElements.size(), configData_->dimension_);
        //set to correct value in case some other meshmanipulator changed things
        ElementFactory::instance().setCollectionOfBasisFunctionSets(&collBasisFSet_);
        ElementFactory::instance().setNumberOfMatrices(numberOfElementMatrices_);
        ElementFactory::instance().setNumberOfVectors(numberOfElementVectors_);
        ElementFactory::instance().setNumberOfTimeLevels(configData_->numberOfTimeLevels_);
        ElementFactory::instance().setNumberOfUnknowns(configData_->numberOfUnknowns_);
        FaceFactory::instance().setNumberOfFaceMatrices(numberOfFaceMatrices_);
        FaceFactory::instance().setNumberOfFaceVectors(numberOfFaceVectors_);
        
        //Stage 0 : Check for required requirements
        logger.assert(linearNoElements.size() == DIM, "The number of Linear Intervals has to map the size of the problem and current it does not");
        
        logger.assert(!(DIM == 3 && periodicX_ && linearNoElements[0] % 2 == 1), "The 3D triangular grid generator can't handle an odd amount of elements in the periodic dimension X");
        logger.assert(!(DIM == 3 && periodicY_ && linearNoElements[1] % 2 == 1), "The 3D triangular grid generator can't handle an odd amount of elements in the periodic dimension Y");
        logger.assert(!(DIM == 3 && periodicZ_ && linearNoElements[2] % 2 == 1), "The 3D triangular grid generator can't handle an odd amount of elements in the periodic dimension Z");
        
        //place the boundary conditions together in a vector.
        std::vector<bool> periodicDIM;
        periodicDIM.push_back(periodicX_);
        if (DIM > 1)
        {
            periodicDIM.push_back(periodicY_);
        }
        if (DIM > 2)
        {
            periodicDIM.push_back(periodicZ_);
        }
        
        //Stage 1 : Precompute some required values
        
        Geometry::PointPhysical<DIM> delta_x;
        
        for (std::size_t i = 0; i < DIM; ++i)
        {
            delta_x[i] = (topRight[i] - bottomLeft[i]) / (linearNoElements[i]);
        }
        
        std::vector<std::size_t> numberOfNodeCoordinatesInEachSubspace(DIM), numberOfNodesInEachSubspace(DIM), numberOfElementsInEachSubspace(DIM);
        
        numberOfNodeCoordinatesInEachSubspace[0] = 1;
        numberOfNodesInEachSubspace[0] = 1;
        numberOfElementsInEachSubspace[0] = 1;
        
        std::size_t totalNumberOfNodeCoordinates, totalNumberOfElements;
        std::size_t numberOfNodesPerElement, numberOfNodesPerGroup;
        std::size_t numberOfTrianglesPerRectangle, totalNumberOfNodes;
        
        totalNumberOfNodeCoordinates = (linearNoElements[0] + 1);
        totalNumberOfNodes = (linearNoElements[0] + (periodicDIM[0] ? 0 : 1));
        //'elements' in this counter denote groups of trianglesPerRectangle elements
        totalNumberOfElements = (linearNoElements[0]);
        numberOfNodesPerElement = 2;
        numberOfNodesPerGroup = 2;
        numberOfTrianglesPerRectangle = 1;
        std::size_t powerOf2;
        
        //start with 1 because you want to ask for the entry at idim - 1
        for (std::size_t idim = 1; idim < DIM; ++idim)
        {
            totalNumberOfNodeCoordinates *= (linearNoElements[idim] + 1);
            totalNumberOfNodes *= (linearNoElements[idim] + (periodicDIM[idim] ? 0 : 1));
            totalNumberOfElements *= (linearNoElements[idim]);
            numberOfNodesPerElement += 1;
            numberOfNodesPerGroup *= 2;
            numberOfTrianglesPerRectangle += (2 * idim - 1);
            numberOfElementsInEachSubspace[idim] = numberOfElementsInEachSubspace[idim - 1] * (linearNoElements[idim - 1]);
            numberOfNodeCoordinatesInEachSubspace[idim] = numberOfNodeCoordinatesInEachSubspace[idim - 1] * (linearNoElements[idim - 1] + 1);
            numberOfNodesInEachSubspace[idim] = numberOfNodesInEachSubspace[idim - 1] * (linearNoElements[idim - 1] + (periodicDIM[idim - 1] ? 0 : 1));
        }
        
        Geometry::PointPhysical<DIM> x;
        
        //Stage 2 : Create the nodes
        
        for (std::size_t nodeIndex = 0; nodeIndex < totalNumberOfNodeCoordinates; ++nodeIndex)
        {
            std::size_t nodeIndexRemain = nodeIndex;
            for (int idim = DIM - 1; idim > -1; --idim)
            {
                x[idim] = bottomLeft[idim] + (nodeIndexRemain / numberOfNodeCoordinatesInEachSubspace[idim] * delta_x[idim]);
                nodeIndexRemain %= numberOfNodeCoordinatesInEachSubspace[idim];
            }
            theMesh_.addNodeCoordinate(x);
        }
        
        for (std::size_t nodeIndex = 0; nodeIndex < totalNumberOfNodes; ++nodeIndex)
        {
            theMesh_.addNode();
        }
        
        auto& nodes = getNodesList(IteratorType::GLOBAL);
        
        //Stage 3 : Create the elements
        
        std::vector<std::size_t> elementNdId(DIM), nodeNdId(DIM), nodeCoordinateNdId(DIM);
        std::vector<std::vector<std::size_t> > globalNodeID(numberOfTrianglesPerRectangle);
        std::vector<std::vector<std::size_t> > globalNodeCoordinateID(numberOfTrianglesPerRectangle);
        
        for (std::size_t elementGroupIndex = 0; elementGroupIndex < totalNumberOfElements; ++elementGroupIndex)
        {
            //first generate node indexes as if we are a cube
            //indicates if the element has to be rotated to make connecting faces; rotates the element 90 degrees along the y-axis if needed
            std::size_t rotate = 0;
            
            for (std::size_t i = 0; i < numberOfTrianglesPerRectangle; ++i)
            {
                globalNodeCoordinateID[i].clear();
                globalNodeID[i].clear();
            }
            int elementIndexRemainder = elementGroupIndex;
            
            for (int idim = DIM - 1; idim > -1; --idim)
            {
                elementNdId[idim] = elementIndexRemainder / numberOfElementsInEachSubspace[idim];
                elementIndexRemainder %= numberOfElementsInEachSubspace[idim];
                rotate = (elementNdId[idim] + rotate) % 2;
            }
            
            for (std::size_t i = 0; i < numberOfNodesPerGroup; ++i)
            {
                if (rotate == 0)
                {
                    powerOf2 = 1;
                    for (std::size_t idim = 0; idim < DIM; ++idim)
                    {
                        nodeCoordinateNdId[idim] = elementNdId[idim] + ((i & powerOf2) != 0);
                        nodeNdId[idim] = elementNdId[idim] + ((i & powerOf2) != 0);
                        if (nodeNdId[idim] >= linearNoElements[idim] && periodicDIM[idim])
                            nodeNdId[idim] = 0;
                        powerOf2 *= 2;
                    }
                }
                else
                {
                    powerOf2 = numberOfNodesPerGroup;
                    for (std::size_t idim = 0; idim < DIM; ++idim)
                    {
                        powerOf2 /= 2;
                        nodeCoordinateNdId[idim] = elementNdId[idim] + (((i ^ rotate) & powerOf2) != 0);
                        nodeNdId[idim] = elementNdId[idim] + (((i ^ rotate) & powerOf2) != 0);
                        if (nodeNdId[idim] >= linearNoElements[idim] && periodicDIM[idim])
                            nodeNdId[idim] = 0;
                    }
                }
                
                std::size_t nodeCoordinateIndex = nodeCoordinateNdId[0];
                std::size_t nodeIndex = nodeNdId[0];
                for (std::size_t idim = 1; idim < DIM; ++idim)
                {
                    nodeCoordinateIndex += nodeCoordinateNdId[idim] * numberOfNodeCoordinatesInEachSubspace[idim];
                    nodeIndex += nodeNdId[idim] * numberOfNodesInEachSubspace[idim];
                }
                
                //then cherrypick the element(s) these vertices should connect to (probably not the cleanest implementation; \bug doesn't work if DIM>3)
                switch (i)
                {
                    case 0:
                        globalNodeID[0].push_back(nodeIndex);
                        globalNodeCoordinateID[0].push_back(nodeCoordinateIndex);
                        break;
                    case 3:
                        //insert in the second place because the ordering of the vertices will work out better
                        globalNodeID[1].insert( ++(globalNodeID[1].begin()), nodeIndex);
                        globalNodeCoordinateID[1].insert( ++(globalNodeCoordinateID[1].begin()), nodeCoordinateIndex);
                        break;
                    case 5:
                        globalNodeID[2].push_back(nodeIndex);
                        globalNodeCoordinateID[2].push_back(nodeCoordinateIndex);
                        break;
                    case 6:
                        //insert in the second place because the ordering of the vertices will work out better
                        globalNodeID[3].insert( ++(globalNodeID[3].begin()), nodeIndex);
                        globalNodeCoordinateID[3].insert( ++(globalNodeCoordinateID[3].begin()), nodeCoordinateIndex);
                        break;
                    case 1:
                        for (std::size_t j = 0; j < numberOfTrianglesPerRectangle; ++j)
                        {
                            if (j != 3)
                            {
                                globalNodeID[j].push_back(nodeIndex);
                                globalNodeCoordinateID[j].push_back(nodeCoordinateIndex);
                            }
                        }
                        break;
                    case 2:
                        for (std::size_t j = 0; j < numberOfTrianglesPerRectangle; ++j)
                        {
                            if (j != 2)
                            {
                                globalNodeID[j].push_back(nodeIndex);
                                globalNodeCoordinateID[j].push_back(nodeCoordinateIndex);
                            }
                        }
                        break;
                    case 4:
                        for (std::size_t j = 0; j < numberOfTrianglesPerRectangle; ++j)
                        {
                            if (j != 1)
                            {
                                globalNodeID[j].push_back(nodeIndex);
                                globalNodeCoordinateID[j].push_back(nodeCoordinateIndex);
                            }
                        }
                        break;
                    case 7:
                        for (std::size_t j = 0; j < numberOfTrianglesPerRectangle; ++j)
                        {
                            if (j != 0)
                            {
                                globalNodeID[j].push_back(nodeIndex);
                                globalNodeCoordinateID[j].push_back(nodeCoordinateIndex);
                            }
                        }
                        break;
                } //switch
            } //for all vertices of the rectangle
            
            for (std::size_t i = 0; i < numberOfTrianglesPerRectangle; ++i)
            {
                Element* newElement = addElement(globalNodeCoordinateID[i]);
                for (std::size_t j = 0; j < globalNodeID[i].size(); ++j)
                {
                    logger.assert(i < globalNodeID.size(), "Requested node %, while there are only %", i, globalNodeID.size());
                    logger.assert(j < globalNodeID[i].size(), "Requested element %, but this node only has %.", j, globalNodeID[i].size());
                    logger.assert(globalNodeID[i][j] < totalNumberOfNodes, "Requested node %, while there are only %", globalNodeID[i][j], totalNumberOfNodes);
                    logger.assert(nodes.size() == totalNumberOfNodes, "Number of vertices is wrong.");
                    nodes[globalNodeID[i][j]]->addElement(newElement, j);
                }
            }
        } //for all rectangles
        
        //Stage 4 : Create the faces
        faceFactory();
        edgeFactory();
    }

    template<std::size_t DIM>
    void MeshManipulator<DIM>::readCentaurMesh(const std::string& filename)
    {
        //set to correct value in case some other meshManipulator changed things
        ElementFactory::instance().setCollectionOfBasisFunctionSets(&collBasisFSet_);
        ElementFactory::instance().setNumberOfMatrices(numberOfElementMatrices_);
        ElementFactory::instance().setNumberOfVectors(numberOfElementVectors_);
        ElementFactory::instance().setNumberOfTimeLevels(configData_->numberOfTimeLevels_);
        ElementFactory::instance().setNumberOfUnknowns(configData_->numberOfUnknowns_);
        FaceFactory::instance().setNumberOfFaceMatrices(numberOfFaceMatrices_);
        FaceFactory::instance().setNumberOfFaceVectors(numberOfFaceVectors_);
        
        //First open the file
        std::ifstream centaurFile;
        
        centaurFile.open(filename.c_str(), std::ios::binary);
        logger.assert_always(centaurFile.is_open(), "Cannot open Centaur meshfile.");
        logger.assert_always(centaurFile.good(), "Something is not so good about this mesh");
        
        switch (configData_->dimension_)
        {
            case 2:
                readCentaurMesh2D(centaurFile);
                break;
            case 3:
                readCentaurMesh3D(centaurFile);
                break;
            default:
                logger(ERROR, "Centaur mesh reader has not been implemented in this DIMension (%)", configData_->dimension_);
        }
        
        //Finally close the file
        centaurFile.close();
    }

    template<std::size_t DIM>
    void MeshManipulator<DIM>::readCentaurMesh2D(std::ifstream& centaurFile)
    {
        auto& elementslist = theMesh_.getElementsList(IteratorType::GLOBAL);
        
        //These are used to check the length of the read lines to check for read errors
        std::uint32_t sizeOfLine;
        
        //This first value in the centaur file is the size of each line in the file;
        centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
        
        // Version number of the Centaur mesh file 0.1 if using Tito's matlab generator
        float version;
        centaurFile.read(reinterpret_cast<char*>(&version), sizeof(version));
        logger(INFO, "This read mesh is in Centaur version % format", version);
        
        // Centaur File Type <0 is two DIMensional and >0 is three DIMensional
        int32_t centaurFileType;
        centaurFile.read(reinterpret_cast<char*>(&centaurFileType), sizeof(centaurFileType));
        
        logger.assert_always(centaurFileType < 0, "Incorrect mesh file. This mesh appears to contain three DIMensional data");
            
        logger(INFO, "Reading a two DIMensional centaur mesh");

        //The rest of the first line is junk
        centaurFile.ignore(sizeOfLine - sizeof(version) - sizeof(centaurFileType));
        std::uint32_t checkInt;

        //Check the first line was read correctly : each line in centaur start and end with the line size as a check
        centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
        logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");

        //Start the second line
        centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));

        //Next read the total number of nodes
        std::uint32_t numberOfNodes;
        centaurFile.read(reinterpret_cast<char*>(&numberOfNodes), sizeof(numberOfNodes));
        logger(INFO, "File contains % nodes", numberOfNodes);

        //Check the second line was read correctly : each line in centaur start and end with the line size as a check
        centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
        logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");

        //placeholder for vertices until we know where the periodic boundaries are
        std::vector<std::vector<std::size_t> > listOfElementsForEachNode(numberOfNodes);

        //Now we will read in all the nodes
        centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
        double nodeCoord[2];
        Geometry::PointPhysical<DIM> nodeCoordPointFormat;
        for (std::size_t i = 0; i < numberOfNodes; i++)
        {
            // Reads the x and y coordinates of each node.
            centaurFile.read(reinterpret_cast<char*>(nodeCoord), sizeof(nodeCoord));
            // pass the node to the nodelist.
            
            //Covert from *double to hpGEM PointPhysical format
            nodeCoordPointFormat[0] = nodeCoord[0];
            nodeCoordPointFormat[1] = nodeCoord[1];
            theMesh_.addNodeCoordinate(nodeCoordPointFormat);
            
        }
        //Now check the node line was read correctly : each line in centaur start and end with the line size as a check
        centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
        logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");

        //Now check how many triangle in the file
        centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
        // Number of triangular elements
        std::uint32_t numberOfTriangles;
        centaurFile.read(reinterpret_cast<char*>(&numberOfTriangles), sizeof(numberOfTriangles));
        logger(INFO, "File contains % triangle(s)", numberOfTriangles);

        //Check the line was read correctly : each line in centaur start and end with the line size as a check
        centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
        logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");

        centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
        if (numberOfTriangles > 0)
        {
            std::vector<std::uint32_t> globalNodeIndexes(3);
            std::vector<std::size_t> globalNodeIndexesSizeT(3);
            
            for (std::size_t i = 0; i < numberOfTriangles; i++)
            {
                
                centaurFile.read(reinterpret_cast<char*>(globalNodeIndexes.data()), sizeof(std::uint32_t) * globalNodeIndexes.size());
                for (std::size_t j = 0; j < 3; j++)
                {
                    globalNodeIndexes[j] = globalNodeIndexes[j] - 1;
                    globalNodeIndexesSizeT[j] = static_cast<std::size_t>(globalNodeIndexes[j]);
                }

                std::size_t id = addElement(globalNodeIndexesSizeT)->getID();
                
                for (std::uint32_t j : globalNodeIndexes)
                {
                    listOfElementsForEachNode[j].push_back(id);
                }
            }
            
        }
        centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
        logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");

        //Now check the number of quaduratiles in the file
        centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
        std::uint32_t numberOfQuads;
        centaurFile.read(reinterpret_cast<char*>(&numberOfQuads), sizeof(numberOfQuads));
        logger(INFO, "File contains % quaduratile(s)", numberOfQuads);

        centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
        logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");

        //Now read the quaduritles in
        centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
        if (numberOfQuads > 0)
        {
            std::vector<std::uint32_t> globalNodeIndexes(4);
            std::vector<std::size_t> globalNodeIndexesSizeT(4);
            for (std::size_t i = 0; i < numberOfQuads; i++)
            {
                //Reading the node indices of each quadrilateral.
                centaurFile.read(reinterpret_cast<char*>(globalNodeIndexes.data()), sizeof(std::uint32_t) * globalNodeIndexes.size());

                // renumbering of the vertices to match the ordering assumed by
                // hpGem:
                std::swap(globalNodeIndexes[2], globalNodeIndexes[3]);
                
                // renumber them from 1..N to 0..N-1.
                for (std::size_t j = 0; j < 4; j++)
                {
                    globalNodeIndexes[j] = globalNodeIndexes[j] - 1;
                    globalNodeIndexesSizeT[j] = static_cast<std::size_t>(globalNodeIndexes[j]);
                }
                
                std::size_t id = addElement(globalNodeIndexesSizeT)->getID();

                for (std::uint32_t j : globalNodeIndexes)
                {
                    listOfElementsForEachNode[j].push_back(id);
                }
            }
            
        }
        //Check the line was read correctly : each line in centaur start and end with the line size as a check
        centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
        logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");

        //now read boundary data
        //nodes first

        centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));

        std::uint32_t numberOfBoundaryNodes;
        centaurFile.read(reinterpret_cast<char*>(&numberOfBoundaryNodes), sizeof(numberOfBoundaryNodes));

        centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
        logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
        logger(INFO, "File contains % boundary nodes", numberOfBoundaryNodes);

        centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));

        std::vector<std::vector<std::uint32_t> > boundaryNodesForEachGroup;
        std::uint32_t boundaryNodeInformation[2];

        for (uint_fast32_t i = 0; i < numberOfBoundaryNodes; ++i)
        {
            centaurFile.read(reinterpret_cast<char*>(boundaryNodeInformation), sizeof(boundaryNodeInformation));
            if (boundaryNodesForEachGroup.size() < boundaryNodeInformation[1] + 1)
            {
                boundaryNodesForEachGroup.resize(boundaryNodeInformation[1] + 1);
            }
            boundaryNodesForEachGroup[boundaryNodeInformation[1]].push_back(boundaryNodeInformation[0] - 1);
        }

        centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
        logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");

        //then faces

        centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));

        std::uint32_t numberOfBoundaryFaces;
        centaurFile.read(reinterpret_cast<char*>(&numberOfBoundaryFaces), sizeof(numberOfBoundaryFaces));
        logger(INFO, "File contains % boundary faces", numberOfBoundaryFaces);

        centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
        logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");

        centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));

        //read the half-faces and store them until we get to the boundary information
        std::vector<HalfFaceDescription> boundaryFaces(numberOfBoundaryFaces);
        std::uint32_t readBoundaryIndex[2];

        for (std::uint_fast32_t i = 0; i < numberOfBoundaryFaces; ++i)
        {
            boundaryFaces[i].nodesList.resize(2);
            centaurFile.read(reinterpret_cast<char*>(readBoundaryIndex), sizeof(std::uint32_t) * 2);
            
            boundaryFaces[i].nodesList[0] = readBoundaryIndex[0] - 1;
            boundaryFaces[i].nodesList[1] = readBoundaryIndex[1] - 1;
            std::vector<std::size_t> candidateElements;
            std::vector<std::size_t>& leftNodeElements = listOfElementsForEachNode[boundaryFaces[i].nodesList[0]];
            std::vector<std::size_t>& rightNodeElements = listOfElementsForEachNode[boundaryFaces[i].nodesList[1]];
            
            std::set_intersection(leftNodeElements.begin(), leftNodeElements.end(), rightNodeElements.begin(), rightNodeElements.end(), std::back_inserter(candidateElements));
            
            //boundary faces *should* border only one element
            logger.assert_always(candidateElements.size() == 1, "candidate boundary face lies at two or more elements");
            boundaryFaces[i].elementNumber = candidateElements[0];
            
            Element* current = elementslist[candidateElements[0]];
            std::vector<std::size_t> faceNodes(2);
            
            for (std::size_t j = 0; j < current->getNumberOfFaces(); ++j)
            {
                faceNodes = current->getPhysicalGeometry()->getGlobalFaceNodeIndices(j);
                if ((faceNodes[0] == boundaryFaces[i].nodesList[0] || faceNodes[0] == boundaryFaces[i].nodesList[1]) && (faceNodes[1] == boundaryFaces[i].nodesList[0] || faceNodes[1] == boundaryFaces[i].nodesList[1]))
                {
                    boundaryFaces[i].localFaceIndex = j;
                }
            }
            
            candidateElements.clear();
        }

        centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
        logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");

        //now read a list of boundary segments (that will link the faces to boundary conditions later on)
        centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));

        std::vector<std::uint32_t> faceToSegment(numberOfBoundaryFaces);
        centaurFile.read(reinterpret_cast<char*>(faceToSegment.data()), sizeof(std::uint32_t) * numberOfBoundaryFaces);

        centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
        logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");

        //now couple the segments to boundary groups
        centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));

        std::uint32_t numberOfSegments;
        centaurFile.read(reinterpret_cast<char*>(&numberOfSegments), sizeof(numberOfSegments));
        logger(INFO, "Boundary is segmented into % types", numberOfSegments);

        centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
        logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");

        centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));

        std::vector<std::uint32_t> segmentToGroup(numberOfSegments);
        centaurFile.read(reinterpret_cast<char*>(segmentToGroup.data()), sizeof(std::uint32_t) * numberOfSegments);

        centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
        logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");

        //now we get the actual boundary information; centaur distinguishes the following
//            1   -1000  - Viscous Wall
//            1001-2000  - Inviscid Wall
//            2001-3000  - Symmetry
//            3001-4000  - Inlet
//            4001-5000  - Outlet
//            5001-6000  - Farfield
//            6001-7000  - Periodic
//            7001-8000  - Shadow
//            8001-8500  - Interface
//            8501-9000  - Wake Surfaces
//            9001-10000 - Moving Walls
//           10001-      - Other
        centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));

        std::uint32_t numberOfGroups;
        centaurFile.read(reinterpret_cast<char*>(&numberOfGroups), sizeof(std::uint32_t));

        centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
        logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");

        centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));

        std::vector<std::uint32_t> groupBCType(numberOfGroups);
        centaurFile.read(reinterpret_cast<char*>(groupBCType.data()), sizeof(std::uint32_t) * numberOfGroups);

        centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
        logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");

        for (std::uint_fast32_t i = 0; i < numberOfBoundaryFaces; ++i)
        {
            logger.assert_always(faceToSegment[i]-1 < numberOfSegments, "tries to find segment %, but there are only % segments", faceToSegment[i], numberOfSegments);
            logger.assert_always(segmentToGroup[faceToSegment[i]-1]-1 < numberOfGroups, "the segment to group mapping contains corrupt data");
            std::uint32_t boundaryType = groupBCType[segmentToGroup[faceToSegment[i]-1]-1];
            if (boundaryType < 1001)
            {
                logger(INFO, "Viscous Wall boundary for face % assigned as WALL_BC", i);
                addFace(elementslist[boundaryFaces[i].elementNumber], boundaryFaces[i].localFaceIndex, nullptr, 0, Geometry::FaceType::WALL_BC);
            }
            else if (boundaryType < 2001)
            {
                logger(INFO, "Inviscid Wall boundary for face % assigned as WALL_BC", i);
                addFace(elementslist[boundaryFaces[i].elementNumber], boundaryFaces[i].localFaceIndex, nullptr, 0, Geometry::FaceType::WALL_BC);
            }
            else if (boundaryType < 3001)
            {
                logger(INFO,  "symmetry plane boundary for face % assigned as WALL_BC", i);
                addFace(elementslist[boundaryFaces[i].elementNumber], boundaryFaces[i].localFaceIndex, nullptr, 0, Geometry::FaceType::WALL_BC);
            }
            else if (boundaryType < 4001)
            {
                logger(INFO, "inlet pipe boundary for face % assigned as OPEN_BC", i);
                addFace(elementslist[boundaryFaces[i].elementNumber], boundaryFaces[i].localFaceIndex, nullptr, 0, Geometry::FaceType::OPEN_BC);
            }
            else if (boundaryType < 5001)
            {
                logger(INFO, "outlet pipe boundary for face % assigned as OPEN_BC", i);
                addFace(elementslist[boundaryFaces[i].elementNumber], boundaryFaces[i].localFaceIndex, nullptr, 0, Geometry::FaceType::OPEN_BC);
            }
            else if (boundaryType < 6001)
            {
                logger(INFO, "farfield boundary for face % assigned as OPEN_BC", i);
                addFace(elementslist[boundaryFaces[i].elementNumber], boundaryFaces[i].localFaceIndex, nullptr, 0, Geometry::FaceType::OPEN_BC);
            }
            else if (boundaryType < 7001)
            {
                logger(INFO,  "periodic boundary for face % ignored for being internal; node connections will be assigned later", i);
            }
            else if (boundaryType < 8001)
            {
                logger(INFO, "shadow boundary for face % ignored for being internal; node connections will be assigned later", i);
            }
            else if (boundaryType < 8501)
            {
                logger(INFO, "interface boundary for face % ignored for being internal", i);
            }
            else if (boundaryType < 9001)
            {
                logger(INFO, "wake boundary for face % assigned as OPEN_BC", i);
                addFace(elementslist[boundaryFaces[i].elementNumber], boundaryFaces[i].localFaceIndex, nullptr, 0, Geometry::FaceType::OPEN_BC);
            }
            else if (boundaryType < 10001)
            {
                logger(INFO, "moving wall boundary for face % assigned as WALL_BC", i);
                addFace(elementslist[boundaryFaces[i].elementNumber], boundaryFaces[i].localFaceIndex, nullptr, 0, Geometry::FaceType::WALL_BC);
            }
            else
            {
                logger(INFO, "alternative boundary condition for face % assigned as WALL_BC", i);
                addFace(elementslist[boundaryFaces[i].elementNumber], boundaryFaces[i].localFaceIndex, nullptr, 0, Geometry::FaceType::WALL_BC);
            }
        }
        //I don't care about the names of the groups, just skip them
        centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));

        centaurFile.ignore(80 * numberOfGroups);

        centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
        logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");

        //now read periodic information and link corresponding vertices
        centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));

        std::uint32_t numberOfTransforms;
        centaurFile.read(reinterpret_cast<char*>(&numberOfTransforms), sizeof(numberOfTransforms));

        centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
        logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");

        for (std::uint_fast32_t i = 0; i < numberOfTransforms; ++i)
        {
            //I dont care about the actual coordinate transformations, just skip them
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            
            centaurFile.ignore(9 * sizeof(std::uint32_t));
            
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            
            centaurFile.ignore(9 * sizeof(std::uint32_t));
            
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
            
            std::uint32_t numberOfNodePairs;
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            
            centaurFile.read(reinterpret_cast<char*>(&numberOfNodePairs), sizeof(numberOfNodePairs));
            
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
            
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            
            std::uint32_t nodePair[2];
            std::vector<std::size_t> combine;
            for (std::uint_fast32_t j = 0; j < numberOfNodePairs; ++j)
            {
                centaurFile.read(reinterpret_cast<char*>(nodePair), 2 * sizeof(std::uint32_t));
                auto& firstList = listOfElementsForEachNode[nodePair[0] - 1];
                auto& secondList = listOfElementsForEachNode[nodePair[1] - 1];
                std::set_union(firstList.begin(), firstList.end(), secondList.begin(), secondList.end(), std::back_inserter(combine));
                firstList = combine;
                secondList = combine;
                combine.clear();
            }
            
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
            
        }

        //new centaur files have zones, but I'm not interested in them

        //now we know periodicity information, construct the vertices
        //remember that listOfElementsForEachNode will in general contain duplicates

        auto& listOfNodes = theMesh_.getNodesList(IteratorType::GLOBAL);

        bool addedNewNode(false);

        for (std::size_t i = 0; i < listOfElementsForEachNode.size(); ++i)
        {
            for (std::size_t j = 0; j < listOfElementsForEachNode[i].size(); ++j)
            {
                Element* current = elementslist[listOfElementsForEachNode[i][j]];
                for (std::size_t k = 0; k < current->getNumberOfNodes(); ++k)
                {
                    //if we did not jet deal with this node and it is the correct one
                    if (current->getNode(k) == nullptr && current->getPhysicalGeometry()->getNodeIndex(k) == i)
                    {
                        if (!addedNewNode)
                        {
                            addNode();
                            addedNewNode = true;
                        }
                        listOfNodes.back()->addElement(current, k);
                    }
                }
            }
            addedNewNode = false;
        }

        //construct the rest of the faces
        faceFactory();
        getElementsList();
    }

    template<std::size_t DIM>
    void MeshManipulator<DIM>::readCentaurMesh3D(std::ifstream& centaurFile)
    {
        
        auto& elementsList = theMesh_.getElementsList(IteratorType::GLOBAL);
        
        //These are used to check the length of the read lines to check for read errors
        std::uint32_t sizeOfLine;
        
        //This first value in the centaur file is the size of each line in the file;
        centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
        
        // Version number of the Centaur mesh file 0.1 if using Tito's matlab generator
        float version;
        centaurFile.read(reinterpret_cast<char*>(&version), sizeof(version));
        logger(INFO, "This read mesh is in Centaur version % format", version);
        
        // Centaur File Type <0 is two DIMensional and >0 is three DIMensional
        std::uint32_t centaurFileType;
        centaurFile.read(reinterpret_cast<char*>(&centaurFileType), sizeof(centaurFileType));
        
        logger.assert_always(centaurFileType > 0, "Incorrect mesh file. This mesh appears to contain two DIMensional data");
            logger(INFO, "Reading a three DIMensional centaur mesh");
            
            //The rest of the first line is junk
            char junk[1024];
            
            std::uint32_t checkInt;
            centaurFile.read(&junk[0], sizeOfLine - sizeof(version) - sizeof(centaurFileType));
            
            //Check the first line was read correctly : each line in centaur start and end with the line size as a check
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
            
            //Start the second line
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            
            //Next read the total number of nodes
            std::uint32_t numberOfNodes;
            centaurFile.read(reinterpret_cast<char*>(&numberOfNodes), sizeof(numberOfNodes));
            logger(INFO, "File contains % nodes", numberOfNodes);
            
            //new centaur versions support splitting this list over multiple lines
            std::uint32_t numberOfNodesPerLine(numberOfNodes);
            if (centaurFileType > 4)
            {
                centaurFile.read(reinterpret_cast<char*>(&numberOfNodesPerLine), sizeof(numberOfNodesPerLine));
                logger(INFO, "One line in the file contains at most % nodes.", numberOfNodesPerLine);
            }
            
            //Check the second line was read correctly : each line in centaur start and end with the line size as a check
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
            
            //Now we will read in all the nodes
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            double nodeCoord[3];
            Geometry::PointPhysical<3> nodeCoordPointFormat;
            for (std::size_t i = 0; i < numberOfNodes; i++)
            {
                if (i > 0 && i % numberOfNodesPerLine == 0)
                {
                    //If all the nodes on a line are read end the line and start a new one
                    centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
                    logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
                    centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
                }
                
                // Reads the x, y and z coordinates of each node.
                centaurFile.read(reinterpret_cast<char*>(nodeCoord), sizeof(nodeCoord));
                // pass the node to the nodelist.
                
                //Covert from *double to hpGEM PointPhysical format
                nodeCoordPointFormat[0] = nodeCoord[0];
                nodeCoordPointFormat[1] = nodeCoord[1];
                nodeCoordPointFormat[2] = nodeCoord[2];
                logger(DEBUG, "In MeshManipulator::readCentaurMesh3D, "
                        "nodeCoordPointFormat = %", nodeCoordPointFormat);
                theMesh_.addNodeCoordinate(nodeCoordPointFormat);
                
            }
            //Now check the node line was read correctly : each line in centaur start and end with the line size as a check
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
            
            //Keep track of the node->Element connectivity to ease face creation later on
            std::vector<std::vector<std::size_t> > listOfElementsForEachNode(numberOfNodes);
            std::vector<Element*> tempElementVector;
            
            //file version 1 has no lines about hexahedra
            if (centaurFileType > 1)
            {
                //Now check how many hexahedra in the file
                centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
                // Number of hexahedral elements
                std::uint32_t numberOfHexahedra;
                centaurFile.read(reinterpret_cast<char*>(&numberOfHexahedra), sizeof(numberOfHexahedra));
                logger(INFO,  "File contains  hexahedron(s)", numberOfHexahedra);
                
                //new centaur versions support splitting this list over multiple lines
                std::uint32_t numberOfHexahedraPerLine(numberOfHexahedra);
                if (centaurFileType > 4)
                {
                    centaurFile.read(reinterpret_cast<char*>(&numberOfHexahedraPerLine), sizeof(numberOfHexahedraPerLine));
                    logger(INFO, "One line in the file contains at most  hexahedra", numberOfHexahedraPerLine);
                }
                
                centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
                logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
                
                std::uint32_t temp;
                std::vector<std::uint32_t> globalNodeIndexes(8);
                std::vector<std::size_t> globalNodeIndexesSizeT(8);
                centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
                for (std::size_t i = 0; i < numberOfHexahedra; i++)
                {
                    if (i > 0 && i % numberOfHexahedraPerLine == 0)
                    {
                        //If all the nodes on a line are read end the line and start a new one
                        centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
                        logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
                        centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
                    }
                    
                    //Reading the node indices of each hexahedron.
                    centaurFile.read(reinterpret_cast<char*>(globalNodeIndexes.data()), sizeof(std::uint32_t) * globalNodeIndexes.size());
                    
                    // renumbering of the vertices to match the ordering assumed by
                    // hpGem: (based on the numbering in hpGEM 1)
                    
                    temp = globalNodeIndexes[2];
                    globalNodeIndexes[2] = globalNodeIndexes[3];
                    globalNodeIndexes[3] = temp;
                    
                    temp = globalNodeIndexes[6];
                    globalNodeIndexes[6] = globalNodeIndexes[7];
                    globalNodeIndexes[7] = temp;
                    
                    // renumber them from 1..N to 0..N-1.
                    for (std::size_t j = 0; j < 8; j++)
                    {
                        globalNodeIndexes[j] = globalNodeIndexes[j] - 1;
                        globalNodeIndexesSizeT[j] = static_cast<std::size_t>(globalNodeIndexes[j]);
                    }
                    
                    Base::Element* newElement = addElement(globalNodeIndexesSizeT);
                    tempElementVector.push_back(newElement);
                    
                    for (std::size_t j = 0; j < 8; ++j)
                    {
                        listOfElementsForEachNode[globalNodeIndexes[j]].push_back(newElement->getID());
                    }
                }
                centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
                logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
            }
            
            //Now check how many triangular prisms in the file
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            // Number of prismatic elements
            std::uint32_t numberOfPrisms;
            centaurFile.read(reinterpret_cast<char*>(&numberOfPrisms), sizeof(numberOfPrisms));
            logger(INFO,  "File contains % triangular prism(s)", numberOfPrisms);
            
            //new centaur versions support splitting this list over multiple lines
            std::uint32_t numberOfPrismsPerLine(numberOfPrisms);
            if (centaurFileType > 4)
            {
                centaurFile.read(reinterpret_cast<char*>(&numberOfPrismsPerLine), sizeof(numberOfPrismsPerLine));
                logger(INFO, "One line in the file contains at most % prisms", numberOfPrismsPerLine);
            }
            
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
            
            std::vector<std::uint32_t> globalNodeIndexes(6);
            std::vector<std::size_t> globalNodeIndexesSizeT(6);
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            for (std::size_t i = 0; i < numberOfPrisms; i++)
            {
                if (i > 0 && i % numberOfPrismsPerLine == 0)
                {
                    //If all the nodes on a line are read end the line and start a new one
                    centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
                    logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
                    centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
                }
                
                //Reading the node indices of each hexahedron.
                centaurFile.read(reinterpret_cast<char*>(globalNodeIndexes.data()), sizeof(std::uint32_t) * globalNodeIndexes.size());
                
                // renumber them from 1..N to 0..N-1.
                for (std::size_t j = 0; j < 6; j++)
                {
                    globalNodeIndexes[j] = globalNodeIndexes[j] - 1;
                    globalNodeIndexesSizeT[j] = static_cast<std::size_t>(globalNodeIndexes[j]);
                }
                
                Base::Element* newElement = addElement(globalNodeIndexesSizeT);
                tempElementVector.push_back(newElement);
                
                for (std::size_t j = 0; j < 6; ++j)
                {
                    listOfElementsForEachNode[globalNodeIndexes[j]].push_back(newElement->getID());
                }
            }
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
            
            //file version 1 has no lines about pyramids
            if (centaurFileType > 1)
            {
                //Now check how many pyramids in the file
                centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
                // Number of pyramids elements
                std::uint32_t numberOfPyraminds;
                centaurFile.read(reinterpret_cast<char*>(&numberOfPyraminds), sizeof(numberOfPyraminds));
                logger(INFO, "File contains % pyramid(s)", numberOfPyraminds);
                
                //new centaur versions support splitting this list over multiple lines
                std::uint32_t numberOfPyramindsPerLine(numberOfPyraminds);
                if (centaurFileType > 4)
                {
                    centaurFile.read(reinterpret_cast<char*>(&numberOfPyramindsPerLine), sizeof(numberOfPyramindsPerLine));
                    logger(INFO, "One line in the file contains at most % pyramids", numberOfPyramindsPerLine);
                }
                
                centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
                logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
                
                std::uint32_t temp;
                globalNodeIndexes.resize(5);
                globalNodeIndexesSizeT.resize(5);
                centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
                for (std::size_t i = 0; i < numberOfPyraminds; i++)
                {
                    if (i > 0 && i % numberOfPyramindsPerLine == 0)
                    {
                        //If all the nodes on a line are read end the line and start a new one
                        centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
                        logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
                        centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
                    }
                    
                    //Reading the node indices of each pyramid.
                    centaurFile.read(reinterpret_cast<char*>(globalNodeIndexes.data()), sizeof(std::uint32_t) * globalNodeIndexes.size());
                    
                    //and then the renumbering fun begins
                    //first make sure we always use the same numbering sceme even when reading from old centaur files
                    //     For centaurfiletypes <= 3, the pyramids will have the
                    //     orientation:
                    //     1-2-3-4 points away from 5
                    //     
                    //     For centaurfiletypes >  3, the pyramids will have the 
                    //     orientation:
                    //     1-2-3-4 points towards 5
                    //     
                    //     mirror the pyramid in the centaurFileType <4 case
                    //     to always get the orientation: 1-2-3-4 points towards 5
                    
                    if (centaurFileType < 4)
                    {
                        temp = globalNodeIndexes[0];
                        globalNodeIndexes[0] = globalNodeIndexes[1];
                        globalNodeIndexes[1] = temp;
                        
                        temp = globalNodeIndexes[2];
                        globalNodeIndexes[2] = globalNodeIndexes[3];
                        globalNodeIndexes[3] = temp;
                    }
                    
                    //now renumber the ordered vertices to the expected numbering in hpGEM
                    //for the moment we have the following corresponcence:
                    // Centaur | hpGEM 1
                    // -----------------
                    //  0      | 1
                    //  1      | 2
                    //  2      | 4
                    //  3      | 3
                    //  4      | 0
                    
                    temp = globalNodeIndexes[0];
                    globalNodeIndexes[0] = globalNodeIndexes[4];
                    globalNodeIndexes[4] = globalNodeIndexes[2];
                    globalNodeIndexes[2] = globalNodeIndexes[1];
                    globalNodeIndexes[1] = temp;
                    
                    // renumber them from 1..N to 0..N-1.
                    for (std::size_t j = 0; j < 5; j++)
                    {
                        globalNodeIndexes[j] = globalNodeIndexes[j] - 1;
                        globalNodeIndexesSizeT[j] = static_cast<std::size_t>(globalNodeIndexes[j]);
                    }
                    
                    Base::Element* newElement = addElement(globalNodeIndexesSizeT);
                    tempElementVector.push_back(newElement);
                    
                    for (std::size_t j = 0; j < 5; ++j)
                    {
                        listOfElementsForEachNode[globalNodeIndexes[j]].push_back(newElement->getID());
                    }
                }
                centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
                logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
            }
            
            //Now check how many tetrahedra in the file
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            // Number of tetrahedral elements
            std::uint32_t numberOfTetrahedra;
            centaurFile.read(reinterpret_cast<char*>(&numberOfTetrahedra), sizeof(numberOfTetrahedra));
            logger(INFO,  "File contains % tetrahedron(s)", numberOfTetrahedra);
            
            //new centaur versions support splitting this list over multiple lines
            std::uint32_t numberOfTetrahedraPerLine(numberOfTetrahedra);
            if (centaurFileType > 4)
            {
                centaurFile.read(reinterpret_cast<char*>(&numberOfTetrahedraPerLine), sizeof(numberOfTetrahedraPerLine));
                logger(INFO, "One line in the file contains at most % tetrahedra", numberOfTetrahedraPerLine);
            }
            
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
            
            globalNodeIndexes.resize(4);
            globalNodeIndexesSizeT.resize(4);
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            for (std::size_t i = 0; i < numberOfTetrahedra; i++)
            {
                if (i > 0 && i % numberOfTetrahedraPerLine == 0)
                {
                    //If all the tetrahedra on a line are read end the line and start a new one
                    centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
                    logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
                    centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
                }
                
                //Reading the node indices of each tetrahedron.
                centaurFile.read(reinterpret_cast<char*>(globalNodeIndexes.data()), sizeof(std::uint32_t) * globalNodeIndexes.size());
                
                // renumber them from 1..N to 0..N-1.
                for (std::size_t j = 0; j < 4; j++)
                {
                    globalNodeIndexes[j] = globalNodeIndexes[j] - 1;
                    globalNodeIndexesSizeT[j] = static_cast<std::size_t>(globalNodeIndexes[j]);
                }
                
                Base::Element* newElement = addElement(globalNodeIndexesSizeT);
                tempElementVector.push_back(newElement);
                
                for (std::size_t j = 0; j < 4; ++j)
                {
                    listOfElementsForEachNode[globalNodeIndexes[j]].push_back(newElement->getID());
                }
            }
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
            
            //    Read number of boundary faces and face to node
            //    
            //    A triangular boundary face can belong to a tetrahedron, prism 
            //    or a pyramid
            //    A quadrilateral boundary face belongs to a hexahedron,
            //    prism or a pyramid
            //    
            //    The storage in the ibfnt array is as follows:
            //                            ibfnt:  1  2  3  4  5  6  7  8
            //    quadrilateral hexahedral face   x  x  x  x  x  x  x  x
            //    triangular prism face           x  x  x  0  x  x  x  0
            //    quadrilateral prism face        x  x  x  x  x  x  0  0
            //    triangular pyramidal face       x  x  x  0  x  x  0  0
            //    quadrilateral pyramidal face    x  x  x  x  x  0  0  0
            //    tetrahedral face                x  x  x  0  x  0  0  0 
            //    
            //    In each case, the node numbers in the first 4 slots are the those
            //    that belong to the face and the node numbers after those are the
            //    rest of the nodes belonging to the cell associated with the face.
            //    For hybfiletypes 3 and before, the quadrilateral were given
            //    with the orientation 1-2-4-3, that is 3 above 1 and 4 above 2.
            //    For hybfiletypes 4 and after, the quadrilateral faces have been
            //    changed to have a more standard 1-2-3-4 ordering. The code below
            //    will convert the old ordering to the new ordering.
            //    
            //    For hybfiletypes 3 and before, for each cell type, the 8th slot
            //    in the ibfnt array is reserved for the panel number (pan)
            //    to which the face belongs. The code below will convert this old
            //    scheme into the new scheme described next.
            //
            //    For hybfiletypes 4 and after, there is another array, ibfacpan,
            //    which stores the panel number allowing the eighth spot of the
            //    ibfnt array to only be used for hexahedral quadrilateral boundary
            //    faces.
            //       
            //    This panel number is then related to the boundary condition.
            
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            
            //number of boundary faces
            std::uint32_t numberOfBoundaryFaces;
            centaurFile.read(reinterpret_cast<char*>(&numberOfBoundaryFaces), sizeof(numberOfBoundaryFaces));
            logger(INFO, "File contains % boundaryFace(s)", numberOfBoundaryFaces);
            
            std::uint32_t boundaryFacesPerLine(numberOfBoundaryFaces);
            if (centaurFileType > 4)
            {
                centaurFile.read(reinterpret_cast<char*>(&boundaryFacesPerLine), sizeof(boundaryFacesPerLine));
                logger(INFO, "One line in the file contains at most % tetrahedra", numberOfTetrahedraPerLine);
            }
            
            HalfFaceDescription *boundarFaces = new HalfFaceDescription[numberOfBoundaryFaces];
            
            std::vector<std::vector<std::uint32_t> > facesForEachCentaurPanel(0);
            //old centaur files use 7 entries to describe the face
            std::uint32_t nodalDescriptionOfTheFace[centaurFileType > 3 ? 8 : 7];
            std::uint32_t panelNumber, ijunk;
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
            
            //first read the information about the faces
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            for (std::size_t i = 0; i < numberOfBoundaryFaces; ++i)
            {
                
                if (i > 0 && i % numberOfBoundaryFaces == 0)
                {
                    //If all the tetrahedra on a line are read end the line and start a new one
                    centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
                    logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
                    centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
                }
                
                centaurFile.read(reinterpret_cast<char*>(&nodalDescriptionOfTheFace[0]), sizeof(nodalDescriptionOfTheFace));
                
                boundarFaces[i].nodesList.resize(3);
                boundarFaces[i].nodesList[0] = nodalDescriptionOfTheFace[0];
                boundarFaces[i].nodesList[1] = nodalDescriptionOfTheFace[1];
                boundarFaces[i].nodesList[2] = nodalDescriptionOfTheFace[2];
                
                //three nodes will uniquely determine the face
                auto& firstNodeList = listOfElementsForEachNode[nodalDescriptionOfTheFace[0]];
                auto& secondNodeList = listOfElementsForEachNode[nodalDescriptionOfTheFace[1]];
                auto& thirdNodeList = listOfElementsForEachNode[nodalDescriptionOfTheFace[2]];
                std::vector<std::size_t> temp, candidates, nodes(boundarFaces->nodesList), intersect;
                std::sort(nodes.begin(), nodes.end());
                std::set_intersection(firstNodeList.begin(), firstNodeList.end(), secondNodeList.begin(), secondNodeList.end(), std::back_inserter(temp));
                std::set_intersection(temp.begin(), temp.end(), thirdNodeList.begin(), thirdNodeList.end(), std::back_inserter(candidates));
                
                logger.assert_always(candidates.size() < 2, "candidate boundary face lies at two or more elements");
                
                boundarFaces[i].elementNumber = candidates[0];
                
                Element* current = elementsList[candidates[0]];
                
                for (std::size_t j = 0; j < current->getNumberOfFaces(); ++j)
                {
                    temp = current->getPhysicalGeometry()->getGlobalFaceNodeIndices(j);
                    std::sort(temp.begin(), temp.end());
                    std::set_intersection(temp.begin(), temp.end(), nodes.begin(), nodes.end(), std::back_inserter(intersect));
                    if (intersect.size() == 3)
                    {
                        boundarFaces->localFaceIndex = j;
                        if (typeid(current->getReferenceGeometry()->getCodim1ReferenceGeometry(j)) == typeid(Geometry::ReferenceSquare))
                        {
                            if (centaurFileType > 3)
                            {
                                boundarFaces[i].nodesList.push_back(boundarFaces[i].nodesList[2]);
                                boundarFaces[i].nodesList[2] = nodalDescriptionOfTheFace[3];
                            }
                            else
                            {
                                boundarFaces[i].nodesList.push_back(nodalDescriptionOfTheFace[3]);
                            }
                        }
                    }
                    intersect.clear();
                }
                
                //then read the information about the panel numbers
                if (centaurFileType < 4)
                {
                    centaurFile.read(reinterpret_cast<char*>(&panelNumber), sizeof(panelNumber));
                    if (centaurFileType > 1)
                    {
                        centaurFile.read(reinterpret_cast<char*>(&ijunk), sizeof(ijunk));
                    }
                    if (panelNumber > facesForEachCentaurPanel.size())
                    {
                        facesForEachCentaurPanel.resize(panelNumber);
                    }
                    facesForEachCentaurPanel[panelNumber].push_back(i);
                }
            }
            
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
            
            //modern centaur file version store panel information separately
            if (centaurFileType > 3)
            {
                centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
                for (std::size_t i = 0; i < numberOfBoundaryFaces; ++i)
                {
                    if (i > 0 && i % numberOfBoundaryFaces == 0)
                    {
                        //If all the tetrahedra on a line are read end the line and start a new one
                        centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
                        logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
                        centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
                    }
                    centaurFile.read(reinterpret_cast<char*>(&panelNumber), sizeof(panelNumber));
                    if (panelNumber > facesForEachCentaurPanel.size())
                    {
                        facesForEachCentaurPanel.resize(panelNumber);
                    }
                    facesForEachCentaurPanel[panelNumber - 1].push_back(i);
                }
                centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
                logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
            }
            
            //put the centaur panels in their boudary group
            std::vector<std::vector<std::size_t> > facesForEachBoundaryGroup(0);
            std::uint32_t groupOfPanelNumber;
            
            //this bit of information is a little late
            std::uint32_t numberOfPanels;
            
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            centaurFile.read(reinterpret_cast<char*>(&numberOfPanels), sizeof(numberOfPanels));
            logger.assert(numberOfPanels == facesForEachCentaurPanel.size(), "Not enough faces in centaur file");
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
            
            //then read the panel to group connections
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            for (std::size_t i = 0; i < numberOfPanels; ++i)
            {
                centaurFile.read(reinterpret_cast<char*>(&groupOfPanelNumber), sizeof(groupOfPanelNumber));
                if (groupOfPanelNumber > facesForEachBoundaryGroup.size())
                {
                    facesForEachBoundaryGroup.resize(groupOfPanelNumber);
                }
                for (std::size_t j = 0; j < facesForEachCentaurPanel[i].size(); ++j)
                {
                    facesForEachBoundaryGroup[groupOfPanelNumber - 1].push_back(facesForEachCentaurPanel[i][j]);
                }
            }
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
            
            std::uint32_t centaurBCType;
            char nameOfBoundaryCondition[80];
            
            //this bit of information is again a little late
            std::uint32_t numberOfBoundaryGroups;
            
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            centaurFile.read(reinterpret_cast<char*>(&numberOfBoundaryGroups), sizeof(numberOfBoundaryGroups));
            logger.assert(numberOfBoundaryGroups == facesForEachBoundaryGroup.size(), "Not enough boundary groups in centaur file");
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
            
            //now set the boundary conditions for each group
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            for (std::size_t i = 0; i < numberOfBoundaryGroups; ++i)
            {
                centaurFile.read(reinterpret_cast<char*>(&centaurBCType), sizeof(centaurBCType));
                if (centaurBCType < 1001)
                {
                    logger(INFO, "Viscous Wall boundary for group % assigned as WALL_BC", i);
                    for (std::size_t j = 0; j < facesForEachBoundaryGroup[i].size(); ++j)
                    {
                        //big assumption on the nature of elementIDs here...
                        addFace(tempElementVector[boundarFaces[facesForEachBoundaryGroup[i][j]].elementNumber], boundarFaces[facesForEachBoundaryGroup[i][j]].localFaceIndex, nullptr, 0, Geometry::FaceType::WALL_BC);
                    }
                }
                else if (centaurBCType < 2001)
                {
                    logger(INFO, "Inviscid Wall boundary for group % assigned as WALL_BC", i);
                    for (std::size_t j = 0; j < facesForEachBoundaryGroup[i].size(); ++j)
                    {
                        //big assumption on the nature of elementIDs here...
                        addFace(tempElementVector[boundarFaces[facesForEachBoundaryGroup[i][j]].elementNumber], boundarFaces[facesForEachBoundaryGroup[i][j]].localFaceIndex, nullptr, 0, Geometry::FaceType::WALL_BC);
                    }
                }
                else if (centaurBCType < 3001)
                {
                    logger(INFO, "symmetry plane boundary for group % assigned as WALL_BC", i);
                    for (std::size_t j = 0; j < facesForEachBoundaryGroup[i].size(); ++j)
                    {
                        //big assumption on the nature of elementIDs here...
                        addFace(tempElementVector[boundarFaces[facesForEachBoundaryGroup[i][j]].elementNumber], boundarFaces[facesForEachBoundaryGroup[i][j]].localFaceIndex, nullptr, 0, Geometry::FaceType::WALL_BC);
                    }
                }
                else if (centaurBCType < 4001)
                {
                    logger(INFO, "inlet pipe boundary for group % assigned as OPEN_BC", i);
                    for (std::size_t j = 0; j < facesForEachBoundaryGroup[i].size(); ++j)
                    {
                        //big assumption on the nature of elementIDs here...
                        addFace(tempElementVector[boundarFaces[facesForEachBoundaryGroup[i][j]].elementNumber], boundarFaces[facesForEachBoundaryGroup[i][j]].localFaceIndex, nullptr, 0, Geometry::FaceType::OPEN_BC);
                    }
                }
                else if (centaurBCType < 5001)
                {
                    logger(INFO, "outlet pipe boundary for group % assigned as OPEN_BC", i);
                    for (std::size_t j = 0; j < facesForEachBoundaryGroup[i].size(); ++j)
                    {
                        //big assumption on the nature of elementIDs here...
                        addFace(tempElementVector[boundarFaces[facesForEachBoundaryGroup[i][j]].elementNumber], boundarFaces[facesForEachBoundaryGroup[i][j]].localFaceIndex, nullptr, 0, Geometry::FaceType::OPEN_BC);
                    }
                }
                else if (centaurBCType < 6001)
                {
                    logger(INFO, "farfield boundary for group % assigned as OPEN_BC", i);
                    for (std::size_t j = 0; j < facesForEachBoundaryGroup[i].size(); ++j)
                    {
                        //big assumption on the nature of elementIDs here...
                        addFace(tempElementVector[boundarFaces[facesForEachBoundaryGroup[i][j]].elementNumber], boundarFaces[facesForEachBoundaryGroup[i][j]].localFaceIndex, nullptr, 0, Geometry::FaceType::OPEN_BC);
                    }
                }
                else if (centaurBCType < 7001)
                {
                    logger(INFO, "periodic boundary for group % ignored for being internal; node connections will be assigned later", i);
                }
                else if (centaurBCType < 8001)
                {
                    logger(INFO, "shadow boundary for group % ignored for being internal; node connections will be assigned later", i);
                }
                else if (centaurBCType < 8501)
                {
                    logger(INFO, "interface boundary for group % ignored for being internal", i);
                }
                else if (centaurBCType < 9001)
                {
                    logger(INFO,  "wake boundary for group % assigned as OPEN_BC", i);
                    for (std::size_t j = 0; j < facesForEachBoundaryGroup[i].size(); ++j)
                    {
                        addFace(tempElementVector[boundarFaces[facesForEachBoundaryGroup[i][j]].elementNumber], boundarFaces[facesForEachBoundaryGroup[i][j]].localFaceIndex, nullptr, 0, Geometry::FaceType::OPEN_BC);
                    }
                }
                else if (centaurBCType < 10001)
                {
                    logger(INFO, "moving wall boundary for group % assigned as WALL_BC", i);
                    for (std::size_t j = 0; j < facesForEachBoundaryGroup[i].size(); ++j)
                    {
                        addFace(tempElementVector[boundarFaces[facesForEachBoundaryGroup[i][j]].elementNumber], boundarFaces[facesForEachBoundaryGroup[i][j]].localFaceIndex, nullptr, 0, Geometry::FaceType::WALL_BC);
                    }
                }
                else
                {
                    logger(INFO, "alternative boundary condition for group % assigned as WALL_BC", i);
                    for (std::size_t j = 0; j < facesForEachBoundaryGroup[i].size(); ++j)
                    {
                        addFace(tempElementVector[boundarFaces[facesForEachBoundaryGroup[i][j]].elementNumber], boundarFaces[facesForEachBoundaryGroup[i][j]].localFaceIndex, nullptr, 0, Geometry::FaceType::WALL_BC);
                    }
                }
                logger(INFO, "total number of boundary faces: %", getFacesList().size());
            }
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
            
            //This is where centaur tells the names of all the boundary groups
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            for (std::size_t i = 0; i < numberOfBoundaryGroups; ++i)
            {
                centaurFile.read(reinterpret_cast<char*>(&nameOfBoundaryCondition[0]), sizeof(nameOfBoundaryCondition));
                logger(INFO, "boundary condition % is called %", i, nameOfBoundaryCondition);
            }
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
            
            //Then comes periodic boundary information
            //file versions 3 and greater store some extra information that hpGEM will be constructing itself
            //this extra information mangles the reading of the usefull information a bit
            double transformationData[16];
            std::uint32_t matchingNodes[2];
            std::uint32_t numberOfPeriodicNodes;
            
            if (centaurFileType > 3)
            {
                std::uint32_t numberOfPeriodicTransformations;
                centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
                centaurFile.read(reinterpret_cast<char*>(&numberOfPeriodicTransformations), sizeof(numberOfPeriodicTransformations));
                logger(INFO, "There are % periodic boundary -> shadow boundary transformation(s)", numberOfPeriodicTransformations);
                centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
                logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
                for (std::size_t i = 0; i < numberOfPeriodicTransformations; ++i)
                {
                    //information on how to do the transformation can be computed later so just throw it away now
                    centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
                    centaurFile.read(reinterpret_cast<char*>(&ijunk), sizeof(ijunk));
                    centaurFile.read(reinterpret_cast<char*>(&transformationData[0]), sizeof(transformationData));
                    centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
                    logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
                    centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
                    centaurFile.read(reinterpret_cast<char*>(&ijunk), sizeof(ijunk));
                    centaurFile.read(reinterpret_cast<char*>(&transformationData[0]), sizeof(transformationData));
                    centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
                    logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
                    
                    //now read the amount of periodic nodes
                    centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
                    centaurFile.read(reinterpret_cast<char*>(&numberOfPeriodicNodes), sizeof(numberOfPeriodicNodes));
                    logger(INFO, "transformation group % contains % node->node matching(s)", i, numberOfPeriodicNodes);
                    centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
                    logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
                    
                    //and the actual pairing information
                    centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
                    for (std::size_t j = 0; j < numberOfPeriodicNodes; j++)
                    {
                        centaurFile.read(reinterpret_cast<char*>(&matchingNodes[0]), sizeof(matchingNodes));
                        
                        ///\bug for EXTREMELY coarse meshes this will destroy the distinction between faces on the boundary of the domain. Workaround: use at least 3 nodes per direction on each face.
                        auto& target = listOfElementsForEachNode[matchingNodes[0]];
                        auto first = std::move(target);
                        auto& second = listOfElementsForEachNode[matchingNodes[1]];
                        
                        //We just std::move()d target, put it back in a defined state
                        target.clear();
                        target.reserve(first.size() + second.size());
                        std::set_union(first.begin(), first.end(), second.begin(), second.end(), target.begin());
                        
                        listOfElementsForEachNode[matchingNodes[1]] = listOfElementsForEachNode[matchingNodes[0]];
                    }
                    centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
                    logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
                }
            }
            else
            {
                //now read the amount of periodic nodes
                centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
                centaurFile.read(reinterpret_cast<char*>(&numberOfPeriodicNodes), sizeof(numberOfPeriodicNodes));
                logger(INFO, "the transformation group contains % node -> node matching(s)", numberOfPeriodicNodes);
                centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
                logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
                
                //and the actual pairing information
                centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
                for (std::size_t j = 0; j < numberOfPeriodicNodes; j++)
                {
                    centaurFile.read(reinterpret_cast<char*>(&matchingNodes[0]), sizeof(matchingNodes));
                    
                    auto& target = listOfElementsForEachNode[matchingNodes[0]];
                    auto first = std::move(target);
                    auto& second = listOfElementsForEachNode[matchingNodes[1]];
                    
                    //We just std::move()d target, put it back in a defined state
                    target.clear();
                    target.reserve(first.size() + second.size());
                    std::set_union(first.begin(), first.end(), second.begin(), second.end(), target.begin());
                    
                    listOfElementsForEachNode[matchingNodes[1]] = listOfElementsForEachNode[matchingNodes[0]];
                }
                centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
                logger.assert_always(checkInt == sizeOfLine && centaurFile.good(), "Error in centaur file.");
            }
            
            logger(INFO, "begin constructing internal faces and internal \"boundaries\"");
            
            //now we know periodicity information, construct the vertices
            //remember that listOfElementsForEachNode will in general contain duplicates
            
            auto& listOfNodes = theMesh_.getNodesList(IteratorType::GLOBAL);
            
            bool addedNewNode(false);
            
            for (std::size_t i = 0; i < listOfElementsForEachNode.size(); ++i)
            {
                for (std::size_t j = 0; j < listOfElementsForEachNode[i].size(); ++j)
                {
                    Element* current = elementsList[listOfElementsForEachNode[i][j]];
                    for (std::size_t k = 0; k < current->getNumberOfNodes(); ++k)
                    {
                        //if we did not jet deal with this node and it is the correct one
                        if (current->getNode(k) == nullptr && current->getPhysicalGeometry()->getNodeIndex(k) == i)
                        {
                            if (!addedNewNode)
                            {
                                addNode();
                                addedNewNode = true;
                            }
                            listOfNodes.back()->addElement(current, k);
                        }
                    }
                }
                addedNewNode = false;
            }
            
            faceFactory();
            edgeFactory();
            
            delete[] boundarFaces;
        
    }
    
#ifdef HPGEM_USE_QHULL
    template<std::size_t DIM>
    void MeshManipulator<DIM>::createUnstructuredMesh(Geometry::PointPhysical<DIM> BottomLeft, Geometry::PointPhysical<DIM> TopRight, std::size_t TotalNoNodes, std::function<double(Geometry::PointPhysical<DIM>)> domainDescription, std::vector<Geometry::PointPhysical<DIM>> fixedPoints, std::function<double(Geometry::PointPhysical<DIM>)> relativeEdgeLength, double growFactor)
    {
        //impossible to create a mesh with more fixed nodes that total nodes
        //note that when equality is met, this will only do a delaunay triangulation
        logger.assert(fixedPoints.size() <= TotalNoNodes, "Cannot create a mesh with more fixed nodes than total nodes");
        
        //set to correct value in case some other meshmanipulator changed things
        ElementFactory::instance().setCollectionOfBasisFunctionSets(&collBasisFSet_);
        ElementFactory::instance().setNumberOfMatrices(numberOfElementMatrices_);
        ElementFactory::instance().setNumberOfVectors(numberOfElementVectors_);
        ElementFactory::instance().setNumberOfTimeLevels(configData_->numberOfTimeLevels_);
        ElementFactory::instance().setNumberOfUnknowns(configData_->numberOfUnknowns_);
        FaceFactory::instance().setNumberOfFaceMatrices(numberOfFaceMatrices_);
        FaceFactory::instance().setNumberOfFaceVectors(numberOfFaceVectors_);
        
        //periodic unstructured mesh generation not yet implemented
        logger.assert(!(periodicX_ || periodicY_ || periodicZ_), "Unstructured mesh generator does not support periodic boundaries");
        
        //guess the required distance between two nodes
        double dist = std::pow(double(TotalNoNodes), -1. / double(dimension()));
        
        for (std::size_t i = 0; i < dimension(); ++i)
        {
            dist *= std::pow(TopRight[i] - BottomLeft[i], 1. / double(dimension()));
        }
        
        std::vector<Geometry::PointPhysical<DIM> > hpGEMCoordinates = fixedPoints;
        
        //seed approximately N points inside the bounding box (total amount will be tweaked later)
        Geometry::PointPhysical<DIM> nextPoint = BottomLeft;
        for (std::size_t i = 0; i < fixedPoints.size(); ++i)
        {
            logger.assert(domainDescription(fixedPoints[i]) < 1e-10, "One of the fixed points is outside of the domain");
            theMesh_.addNode();
            theMesh_.addNodeCoordinate(fixedPoints[i]);
        }
        //cant do nested for loops for generic dimension
        while (nextPoint[DIM - 1] < TopRight[DIM - 1] - 1e-10)
        {
            std::size_t incrementalDimension = 0;
            //if the point is already far enough to the right, reset and continue with the next dimension
            for (; nextPoint[incrementalDimension] > TopRight[incrementalDimension] + 1e-10; ++incrementalDimension)
            {
                nextPoint[incrementalDimension] = BottomLeft[incrementalDimension];
            }
            nextPoint[incrementalDimension] += dist;
            if (domainDescription(nextPoint) < 0)
            {
                hpGEMCoordinates.push_back(nextPoint);
                theMesh_.addNode();
                theMesh_.addNodeCoordinate(nextPoint);
            }
        }
        std::size_t nFixedPoints = fixedPoints.size();
        //there are not enough points to do a triangulation
        logger.assert(DIM < nFixedPoints, "Could not construct enough points for the initial triangulation");
        //there is inherent rounding down in the gridding and some nodes are outside the domain (so they are discarded)
        logger.assert(hpGEMCoordinates.size() <= TotalNoNodes, "Constructed too many nodes");
        
        while (hpGEMCoordinates.size() < TotalNoNodes)
        {
            //start of QHull magic to create a triangulation
            orgQhull::RboxPoints qHullCoordinates;
            qHullCoordinates.setDimension(DIM);
            qHullCoordinates.reserveCoordinates(DIM * hpGEMCoordinates.size());
            
            for (Geometry::PointPhysical<DIM> point : hpGEMCoordinates)
            {
                qHullCoordinates.append(DIM, point.data());
            }
            
            //create the triangulation, pass "d" for delaunay
            //"QJ" because there are likely to be groups of more that (d+1) cocircular nodes in a regular grid, so joggle them up a bit
            orgQhull::Qhull triangulation;
            triangulation.runQhull(qHullCoordinates, "d Qbb Qx Qc Qt");
            
            for (orgQhull::QhullFacet triangle : triangulation.facetList())
            {
                if (triangle.isGood() && !triangle.isUpperDelaunay())
                {
                    Geometry::PointPhysical<DIM> center;
                    std::vector<std::size_t> pointIndices;
                    for (auto vertex : triangle.vertices())
                    {
                        logger.assert(vertex.point().id() >= 0, "QHull breaks our assumptions on indexes");
                        center += hpGEMCoordinates[vertex.point().id()];
                        pointIndices.push_back(vertex.point().id());
                    }
                    center = center / pointIndices.size();
                    if (domainDescription(center) < 0)
                    {
                        auto newElement = addElement(pointIndices);
                        for (std::size_t i = 0; i < pointIndices.size(); ++i)
                        {
                            theMesh_.getNodesList(IteratorType::GLOBAL)[pointIndices[i]]->addElement(newElement, i);
                        }
                    }
                }
            }
            //end of QHull magic
            
            //extract connectivity information
            faceFactory();
            edgeFactory();
            
            //compute current and expected (relative) edge length
            std::vector<double> expectedLength;
            std::multimap<double, std::size_t> knownLengths;
            std::vector<double> currentLength;
            //for proper scaling
            double totalcurrentLength = 0;
            double totalexpectedLength = 0;
            bool needsExpansion = false;
            
            //compute the expected relative length at the coordinates
            for (std::size_t i = 0; i < hpGEMCoordinates.size(); ++i)
            {
                double newLength = relativeEdgeLength(hpGEMCoordinates[i]);
                expectedLength.push_back(newLength);
                if (std::isnan(newLength) || std::isinf(newLength))
                {
                    needsExpansion |= true;
                }
                else
                {
                    //cannot deliberately construct tangled meshes
                    logger.assert(newLength > 0, "Found an edge that is supposed to have a negative length");
                    knownLengths.insert( {newLength, i});
                }
            }
            
            //if the desired relative edge length is not known everywhere, slowly make them larger
            //because apparently the user is not interested in controlling edge lengths for this part
            //but sudden enlargement leads to a bad mesh quality
            if (needsExpansion)
            {
                //iterate over all nodes, sorted by edge lengths
                for (std::pair<double, std::size_t> entry : knownLengths)
                {
                    Node* current = theMesh_.getNodesList(IteratorType::GLOBAL)[entry.second];
                    for (Element* element : current->getElements())
                    {
                        for (std::size_t i = 0; i < element->getNumberOfNodes(); ++i)
                        {
                            if (std::isnan(expectedLength[element->getNode(i)->getID()]) || std::isinf(expectedLength[element->getNode(i)->getID()]))
                            {
                                expectedLength[element->getNode(i)->getID()] = growFactor * entry.first;
                                
                                //inserting does not invalidate the iterators;
                                //new node has a larger edge length, so it is guaranteed to be visited later on
                                knownLengths.insert(knownLengths.end(), {growFactor * entry.first, element->getNode(i)->getID()});
                            }
                        }
                    }
                }
            }
            //iterate over all edges to compute total length and scaling factor
            //the volume scales with (total edge length)^dimension
            //the total volume filled by the edges should be constant
            //so scale appropriately
            if (DIM == 1)
            {
                //the algorithm is mostly dimension independent, but the data type it operates on is not
                for (Element* element : theMesh_.getElementsList(IteratorType::GLOBAL))
                {
                    Geometry::PointPhysical<DIM> firstNode = element->getPhysicalGeometry()->getLocalNodeCoordinates(0);
                    Geometry::PointPhysical<DIM> secondNode = element->getPhysicalGeometry()->getLocalNodeCoordinates(1);
                    currentLength.push_back(L2Norm(firstNode - secondNode));
                    totalcurrentLength += currentLength.back();
                    totalexpectedLength += expectedLength[element->getNode(0)->getID()] / 2;
                    totalexpectedLength += expectedLength[element->getNode(1)->getID()] / 2;
                }
            }
            else if (DIM == 2)
            {
                for (Face* face : theMesh_.getFacesList(IteratorType::GLOBAL))
                {
                    Geometry::PointPhysical<DIM> firstNode, secondNode;
                    std::vector<std::size_t> nodeIndices = face->getPtrElementLeft()->getReferenceGeometry()->getCodim1EntityLocalIndices(face->localFaceNumberLeft());
                    firstNode = face->getPtrElementLeft()->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[0]);
                    secondNode = face->getPtrElementLeft()->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[1]);
                    currentLength.push_back(L2Norm(firstNode - secondNode));
                    totalcurrentLength += currentLength.back() * currentLength.back();
                    totalexpectedLength += std::pow(expectedLength[face->getPtrElementLeft()->getNode(nodeIndices[0])->getID()] + expectedLength[face->getPtrElementLeft()->getNode(nodeIndices[1])->getID()], 2.) / 4.;
                }
            }
            else
            {
                for (Edge* edge : theMesh_.getEdgesList(IteratorType::GLOBAL))
                {
                    Geometry::PointPhysical<DIM> firstNode, secondNode;
                    std::vector<std::size_t> nodeIndices = edge->getElement(0)->getReferenceGeometry()->getCodim2EntityLocalIndices(edge->getEdgeNumber(0));
                    firstNode = edge->getElement(0)->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[0]);
                    secondNode = edge->getElement(0)->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[1]);
                    currentLength.push_back(L2Norm(firstNode - secondNode));
                    totalcurrentLength += currentLength.back() * currentLength.back() * currentLength.back();
                    totalexpectedLength += std::pow(expectedLength[edge->getElement(0)->getNode(nodeIndices[0])->getID()] + expectedLength[edge->getElement(0)->getNode(nodeIndices[1])->getID()], 3.) / 8.;
                }
            }
            
            //all regions of the domain where elements are allowed to be as large as possible must be connected to regions where relativeEdgeLength provides a limitation
            logger.assert(!std::isnan(totalexpectedLength) && !std::isinf(totalexpectedLength), "could not infer edge sizes for the entirety of the domain");
            
            //sort the centers of the edges such that the centers of the large edges are indexed first
            //note that in this case the inverse measure is computed, because that will result in a more natural force computation later on
            std::multimap<double, Geometry::PointPhysical<DIM> > centerPoints;
            if (DIM == 1)
            {
                //the algorithm is mostly dimension independent, but the data type it operates on is not
                for (Element* element : theMesh_.getElementsList(IteratorType::GLOBAL))
                {
                    Geometry::PointPhysical<DIM> firstNode, secondNode;
                    firstNode = element->getPhysicalGeometry()->getLocalNodeCoordinates(0);
                    secondNode = element->getPhysicalGeometry()->getLocalNodeCoordinates(1);
                    //length is scaled in case somebody hasty decides to add smoothing at this point
                    //all edges should be squeezed a little if the algorithm is to work correctly so pretend the volume is 1.5 times as large
                    //remember to scale back from a volume measure to a length measure
                    double length = (expectedLength[element->getNode(0)->getID()] + expectedLength[element->getNode(1)->getID()]) / currentLength[centerPoints.size()] * 2 * totalcurrentLength / totalexpectedLength;
                    centerPoints.insert( {length, (firstNode + secondNode) / 2.});
                }
            }
            else if (DIM == 2)
            {
                for (Face* face : theMesh_.getFacesList(IteratorType::GLOBAL))
                {
                    Geometry::PointPhysical<DIM> firstNode, secondNode;
                    std::vector<std::size_t> nodeIndices = face->getPtrElementLeft()->getReferenceGeometry()->getCodim1EntityLocalIndices(face->localFaceNumberLeft());
                    firstNode = face->getPtrElementLeft()->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[0]);
                    secondNode = face->getPtrElementLeft()->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[1]);
                    //length is scaled in case somebody hasty decides to add smoothing at this point
                    //all edges should be squeezed a little if the algorithm is to work correctly so pretend the volume is 1.5 times as large
                    //remember to scale back from a volume measure to a length measure
                    double length = (expectedLength[face->getPtrElementLeft()->getNode(nodeIndices[0])->getID()] + expectedLength[face->getPtrElementLeft()->getNode(nodeIndices[1])->getID()]) / currentLength[centerPoints.size()] * std::pow(2 * totalcurrentLength / totalexpectedLength, 1. / 2.);
                    centerPoints.insert( {length, (firstNode + secondNode) / 2});
                }
            }
            else
            {
                for (Edge* edge : theMesh_.getEdgesList(IteratorType::GLOBAL))
                {
                    Geometry::PointPhysical<DIM> firstNode, secondNode;
                    std::vector<std::size_t> nodeIndices = edge->getElement(0)->getReferenceGeometry()->getCodim2EntityLocalIndices(edge->getEdgeNumber(0));
                    firstNode = edge->getElement(0)->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[0]);
                    secondNode = edge->getElement(0)->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[1]);
                    //length is scaled in case somebody hasty decides to add smoothing at this point
                    //all edges should be squeezed a little if the algorithm is to work correctly so pretend the volume is 1.5 times as large
                    //remember to scale back from a volume measure to a length measure
                    double length = (expectedLength[edge->getElement(0)->getNode(nodeIndices[0])->getID()] + expectedLength[edge->getElement(0)->getNode(nodeIndices[1])->getID()]) / currentLength[centerPoints.size()] * std::pow(2 * totalcurrentLength / totalexpectedLength, 1. / 3.);
                    centerPoints.insert( {length, (firstNode + secondNode) / 2});
                }
            }
            //insert nodes in the longest edges only, in an attempt to make them all equally long
            //cannot use range based loop because the centerpoint is not marked
            auto it = centerPoints.begin();
            std::size_t nNewNodes = std::min(centerPoints.size() / 2, TotalNoNodes - hpGEMCoordinates.size());
            for (std::size_t i = 0; i < nNewNodes && it != centerPoints.end(); ++i, ++it)
            {
                if (domainDescription(it->second) < 0)
                {
                    hpGEMCoordinates.push_back(it->second);
                }
                else
                {
                    --i;
                }
            }
            
            //cleanest solution, but not the fastest
            theMesh_.clear();
            for (Geometry::PointPhysical<DIM> point : hpGEMCoordinates)
            {
                theMesh_.addNode();
                theMesh_.addNodeCoordinate(point);
            }
        }
        //start of QHull magic to create a triangulation
        orgQhull::RboxPoints qHullCoordinates;
        qHullCoordinates.setDimension(DIM);
        qHullCoordinates.reserveCoordinates(DIM * hpGEMCoordinates.size());
        
        for (Geometry::PointPhysical<DIM> point : hpGEMCoordinates)
        {
            qHullCoordinates.append(DIM, point.data());
        }
        
        //create the triangulation, pass "d" for delaunay
        //no "QJ" because periodic boundary nodes must keep the current exact distaince (including direction)
        orgQhull::Qhull triangulation(qHullCoordinates, "d QbB Qx Qc Qt");
        
        for (orgQhull::QhullFacet triangle : triangulation.facetList())
        {
            if (triangle.isGood() && !triangle.isUpperDelaunay())
            {
                Geometry::PointPhysical<DIM> center;
                std::vector<std::size_t> pointIndices;
                for (auto vertexIt1 = triangle.vertices().begin(); vertexIt1 != triangle.vertices().end(); ++vertexIt1)
                {
                    center += hpGEMCoordinates[(*vertexIt1).point().id()];
                    pointIndices.push_back((*vertexIt1).point().id());
                }
                center = center / pointIndices.size();
                if (domainDescription(center) < 0)
                {
                    auto newElement = addElement(pointIndices);
                    for (std::size_t i = 0; i < pointIndices.size(); ++i)
                    {
                        theMesh_.getNodesList(IteratorType::GLOBAL)[pointIndices[i]]->addElement(newElement, i);
                    }
                }
            }
        }
        //end of QHull magic
        
        edgeFactory();
        faceFactory();
        
        std::vector<std::size_t> fixedPointIdxs;
        for (std::size_t i = 0; i < nFixedPoints; ++i)
        {
            fixedPointIdxs.push_back(i);
        }
        
        updateMesh(domainDescription, fixedPointIdxs, relativeEdgeLength, growFactor);
    }
    
    ///\bug Assumes a DG basis function set is used. (Workaround: set the basis function set again after calling this routine if you are using something conforming)
    template<std::size_t DIM>
    void MeshManipulator<DIM>::updateMesh(std::function<double(Geometry::PointPhysical<DIM>)> domainDescription, std::vector<std::size_t> fixedPointIdxs, std::function<double(Geometry::PointPhysical<DIM>)> relativeEdgeLength, double growFactor, std::function<bool(Geometry::PointPhysical<DIM>)> isOnPeriodic, std::function<Geometry::PointPhysical<DIM>(Geometry::PointPhysical<DIM>)> duplicatePeriodic, std::function<bool(Geometry::PointPhysical<DIM>)> isOnOtherPeriodic, std::function<Geometry::PointPhysical<DIM>(Geometry::PointPhysical<DIM>)> safeSpot, std::vector<std::size_t> dontConnect)
    {
        std::sort(fixedPointIdxs.begin(), fixedPointIdxs.end());
        bool needsExpansion = false;
        double totalCurrentLength = 0;
        double oldQuality = 0;
        double worstQuality = 0.5;
        
        std::set<std::pair<std::size_t, std::size_t> > periodicPairing {};

        //set to correct value in case some other meshmanipulator changed things
        ElementFactory::instance().setCollectionOfBasisFunctionSets(&collBasisFSet_);
        ElementFactory::instance().setNumberOfMatrices(numberOfElementMatrices_);
        ElementFactory::instance().setNumberOfVectors(numberOfElementVectors_);
        ElementFactory::instance().setNumberOfTimeLevels(configData_->numberOfTimeLevels_);
        ElementFactory::instance().setNumberOfUnknowns(configData_->numberOfUnknowns_);
        FaceFactory::instance().setNumberOfFaceMatrices(numberOfFaceMatrices_);
        FaceFactory::instance().setNumberOfFaceVectors(numberOfFaceVectors_);
        
        for (Node* node : theMesh_.getNodesList(IteratorType::GLOBAL))
        {
            Geometry::PointPhysical<DIM> point;
            Geometry::PointPhysical<DIM> compare;
            point = node->getElement(0)->getPhysicalGeometry()->getLocalNodeCoordinates(node->getNodeNumber(0));
            std::set<std::size_t> equivalentIndices {};
            equivalentIndices.insert(node->getElement(0)->getPhysicalGeometry()->getNodeIndex(node->getNodeNumber(0)));
            for (std::size_t i = 1; i < node->getNumberOfElements(); ++i)
            {
                compare = node->getElement(i)->getPhysicalGeometry()->getLocalNodeCoordinates(node->getNodeNumber(i));
                if (compare != point)
                {
                    equivalentIndices.insert(node->getElement(i)->getPhysicalGeometry()->getNodeIndex(node->getNodeNumber(i)));
                }
            }
            auto equivalentIterator = equivalentIndices.begin();
            ++equivalentIterator;
            for (; equivalentIterator != equivalentIndices.end(); ++equivalentIterator)
            {
                periodicPairing.insert( {*equivalentIndices.begin(), *equivalentIterator});
            }
        }
        
        //compute the lengths of the edges and how far the nodes have moved, to see if the nodes have moved so far that a retriangulation is in order
        double maxShift = 0;
        //except don't bother if a retriangulation is in order anyway
        if (oldNodeLocations_.size() == theMesh_.getNodeCoordinates().size())
        {
            std::vector<double> unscaledShift {};
            unscaledShift.reserve(theMesh_.getNumberOfNodes(IteratorType::GLOBAL));
            //compute current and expected (relative) edge length
            std::vector<double> expectedLength {};
            expectedLength.reserve(theMesh_.getNumberOfNodes(IteratorType::GLOBAL));
            std::multimap<double, std::size_t> knownLengths;
            std::vector<double> currentLength {};
            //for proper scaling
            double totalexpectedLength = 0;
            for (Node* node : theMesh_.getNodesList(IteratorType::GLOBAL))
            {
                Geometry::PointPhysical<DIM> point = node->getElement(0)->getPhysicalGeometry()->getLocalNodeCoordinates(node->getNodeNumber(0));
                unscaledShift.push_back(L2Norm(oldNodeLocations_[expectedLength.size()] - point));
                expectedLength.push_back(relativeEdgeLength(point));
                if (isnan(expectedLength.back()) || isinf(expectedLength.back()))
                {
                    needsExpansion |= true;
                }
                else
                {
                    knownLengths.insert( {expectedLength.back(), expectedLength.size() - 1});
                }
            }
            if (needsExpansion)
            {
                //iterate over all nodes, sorted by edge lengths
                for (std::pair<double, std::size_t> entry : knownLengths)
                {
                    Node* current = theMesh_.getNodesList(IteratorType::GLOBAL)[entry.second];
                    for (Element* element : current->getElements())
                    {
                        for (std::size_t i = 0; i < element->getNumberOfNodes(); ++i)
                        {
                            if (std::isnan(expectedLength[element->getNode(i)->getID()]) || std::isinf(expectedLength[element->getNode(i)->getID()]))
                            {
                                expectedLength[element->getNode(i)->getID()] = growFactor * entry.first;
                                
                                //inserting does not invalidate the iterators;
                                //new node has a larger edge length, so it is guaranteed to be visited later on
                                
                                knownLengths.insert(knownLengths.end(), {growFactor * entry.first, element->getNode(i)->getID()});
                            }
                        }
                    }
                }
            }
            //iterate over all edges to compute total length and scaling factor
            //the volume scales with (total edge length)^dimension
            //the total volume filled by the edges should be constant
            //so scale appropriately
            if (DIM == 1)
            {
                currentLength.reserve(theMesh_.getNumberOfElements(IteratorType::GLOBAL));
                //the algorithm is mostly dimension independent, but the data type it operates on is not
                for (Element* element : theMesh_.getElementsList(IteratorType::GLOBAL))
                {
                    Geometry::PointPhysical<DIM> firstNode, secondNode;
                    firstNode = element->getPhysicalGeometry()->getLocalNodeCoordinates(0);
                    secondNode = element->getPhysicalGeometry()->getLocalNodeCoordinates(1);
                    currentLength.push_back(L2Norm(firstNode - secondNode));
                    totalCurrentLength += currentLength.back();
                    totalexpectedLength += expectedLength[element->getNode(0)->getID()] / 2;
                    totalexpectedLength += expectedLength[element->getNode(1)->getID()] / 2;
                    maxShift = std::max(maxShift, L2Norm(firstNode - oldNodeLocations_[element->getPhysicalGeometry()->getNodeIndex(0)]) / currentLength.back());
                    maxShift = std::max(maxShift, L2Norm(secondNode - oldNodeLocations_[element->getPhysicalGeometry()->getNodeIndex(1)]) / currentLength.back());
                }
            }
            else if (DIM == 2)
            {
                currentLength.reserve(theMesh_.getNumberOfFaces(IteratorType::GLOBAL));
                for (Face* face : theMesh_.getFacesList(IteratorType::GLOBAL))
                {
                    Geometry::PointPhysical<DIM> firstNode, secondNode;
                    std::vector<std::size_t> nodeIndices = face->getPtrElementLeft()->getReferenceGeometry()->getCodim1EntityLocalIndices(face->localFaceNumberLeft());
                    firstNode = face->getPtrElementLeft()->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[0]);
                    secondNode = face->getPtrElementLeft()->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[1]);
                    currentLength.push_back(L2Norm(firstNode - secondNode));
                    totalCurrentLength += currentLength.back() * currentLength.back();
                    totalexpectedLength += std::pow(expectedLength[face->getPtrElementLeft()->getNode(nodeIndices[0])->getID()] + expectedLength[face->getPtrElementLeft()->getNode(nodeIndices[1])->getID()], 2.) / 4.;
                    maxShift = std::max(maxShift, L2Norm(firstNode - oldNodeLocations_[face->getPtrElementLeft()->getPhysicalGeometry()->getNodeIndex(nodeIndices[0])]) / currentLength.back());
                    maxShift = std::max(maxShift, L2Norm(secondNode - oldNodeLocations_[face->getPtrElementLeft()->getPhysicalGeometry()->getNodeIndex(nodeIndices[1])]) / currentLength.back());
                }
                worstQuality = 1;
                for (Element* element : theMesh_.getElementsList(IteratorType::GLOBAL))
                {
                    std::array<double, 3> edgeLengths;
                    for (std::size_t i = 0; i < 3; ++i)
                    {
                        edgeLengths[i] = currentLength[element->getFace(i)->getID()];
                    }
                    worstQuality = std::min(worstQuality, (edgeLengths[0] + edgeLengths[1] - edgeLengths[2]) * (edgeLengths[1] + edgeLengths[2] - edgeLengths[0]) * (edgeLengths[2] + edgeLengths[0] - edgeLengths[1]) / edgeLengths[0] / edgeLengths[1] / edgeLengths[2]);
                }
            }
            else
            {
                currentLength.reserve(theMesh_.getNumberOfEdges(IteratorType::GLOBAL));
                for (Edge* edge : theMesh_.getEdgesList(IteratorType::GLOBAL))
                {
                    Geometry::PointPhysical<DIM> firstNode, secondNode;
                    std::vector<std::size_t> nodeIndices = edge->getElement(0)->getReferenceGeometry()->getCodim2EntityLocalIndices(edge->getEdgeNumber(0));
                    firstNode = edge->getElement(0)->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[0]);
                    secondNode = edge->getElement(0)->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[1]);
                    currentLength.push_back(L2Norm(firstNode - secondNode));
                    totalCurrentLength += currentLength.back() * currentLength.back() * currentLength.back();
                    totalexpectedLength += std::pow(expectedLength[edge->getElement(0)->getNode(nodeIndices[0])->getID()] + expectedLength[edge->getElement(0)->getNode(nodeIndices[1])->getID()], 3.) / 8.;
                    maxShift = std::max(maxShift, L2Norm(firstNode - oldNodeLocations_[edge->getElement(0)->getPhysicalGeometry()->getNodeIndex(nodeIndices[0])]) / currentLength.back());
                    maxShift = std::max(maxShift, L2Norm(secondNode - oldNodeLocations_[edge->getElement(0)->getPhysicalGeometry()->getNodeIndex(nodeIndices[1])]) / currentLength.back());
                }
            }
            
            //all regions of the domain where elements are allowed to be as large as possible must be connected to regions where relativeEdgeLength provides a limitation
            logger.assert(!std::isnan(totalexpectedLength) && !std::isinf(totalexpectedLength), "Could not infer edge sizes for the entirety of the domain");
        }
        std::size_t counter = 0;
        double maxMovement = std::numeric_limits<double>::infinity();
        std::vector<double> currentLength {};
        std::vector<Geometry::PointPhysical<DIM>> movement(theMesh_.getNodeCoordinates().size());
        //stop after n iterations, or (when the nodes have stopped moving and the mesh is not becoming worse), or when the mesh is great, or when the mesh is decent, but worsening
        while ((counter < 10000 && (maxMovement > 1e-3 || oldQuality - worstQuality > 1e-3) && worstQuality < 0.8 && (worstQuality < 2. / 3. || oldQuality - worstQuality < 0)) || counter < 5)
        {
            counter++;
            if ((maxShift > 0.1 && (oldQuality - worstQuality) < 5e-3 * maxShift) || (oldNodeLocations_.size() != theMesh_.getNumberOfNodeCoordinates()) || worstQuality < 1e-6)
            {
                maxShift = 0;
                
                orgQhull::RboxPoints qHullCoordinates {};
                qHullCoordinates.setDimension(DIM);
                oldNodeLocations_.clear();
                oldNodeLocations_.reserve(theMesh_.getNumberOfNodeCoordinates());
                for (Geometry::PointPhysical<DIM>& point : theMesh_.getNodeCoordinates())
                {
                    oldNodeLocations_.push_back(point);
                }
                theMesh_.clear();
                auto pairingIterator = periodicPairing.begin();
                for (Geometry::PointPhysical<DIM> point : oldNodeLocations_)
                {
                    theMesh_.addNodeCoordinate(point);
                    if (pairingIterator == periodicPairing.end())
                    {
                        theMesh_.addNode();
                    }
                    else
                    {
                        //skip one insertion for each master/slave pair
                        pairingIterator++;
                    }
                }
                
                std::vector<std::size_t> vertexIndex {};
                vertexIndex.resize(theMesh_.getNumberOfNodeCoordinates(), std::numeric_limits<std::size_t>::max());
                pairingIterator = periodicPairing.begin();
                std::size_t currentNodeNumber = 0;
                for (std::size_t i = 0; i < theMesh_.getNumberOfNodeCoordinates();)
                {
                    vertexIndex[i] = currentNodeNumber;
                    //see if there are any new boundary nodes
                    if(isOnPeriodic(theMesh_.getNodeCoordinates()[i]) && !(pairingIterator != periodicPairing.end() && pairingIterator->first == i))
                    {
                        std::size_t j = pairingIterator->first;
                        periodicPairing.insert({i, theMesh_.getNumberOfNodeCoordinates()});
                        logger(DEBUG, "periodic pair: % % ", i, theMesh_.getNumberOfNodeCoordinates());
                        pairingIterator = std::find_if(periodicPairing.begin(), periodicPairing.end(), [=](const std::pair<std::size_t, std::size_t>& p)->bool{return p.first == std::min(i, j);});
                        Geometry::PointPhysical<DIM> newNodeCoordinate = duplicatePeriodic(theMesh_.getNodeCoordinates()[i]);
                        logger(DEBUG, "new periodic pair coordinates: % %", theMesh_.getNodeCoordinates()[i], newNodeCoordinate);
                        theMesh_.addNodeCoordinate(newNodeCoordinate);
                        oldNodeLocations_.push_back(newNodeCoordinate);
                        vertexIndex.resize(theMesh_.getNumberOfNodeCoordinates(), std::numeric_limits<std::size_t>::max());
                    }
                    //see if there are any non-boundary nodes that slipped into the boundary
                    if(isOnOtherPeriodic(theMesh_.getNodeCoordinates()[i]) && !(pairingIterator != periodicPairing.end() && pairingIterator->first == i))
                    {
                        theMesh_.getNodeCoordinates()[i] = safeSpot(theMesh_.getNodeCoordinates()[i]);
                    }
                    //assign boundary nodes
                    while (pairingIterator != periodicPairing.end() && pairingIterator->first == i)
                    {
                        logger(DEBUG, "periodic pair: % % ", pairingIterator->first, pairingIterator->second);
                        logger.assert(Base::L2Norm(duplicatePeriodic(theMesh_.getNodeCoordinates()[pairingIterator->first]) - theMesh_.getNodeCoordinates()[pairingIterator->second]) < 1e-9, "periodic pair is not moving simulateously");
                        vertexIndex[pairingIterator->second] = currentNodeNumber;
                        ++pairingIterator;
                    }
                    currentNodeNumber++;
                    //skip over already set boundary nodes
                    while (i < theMesh_.getNumberOfNodeCoordinates() && vertexIndex[i] < std::numeric_limits<std::size_t>::max())
                    {
                        ++i;
                    }
                }
                logger(DEBUG, "periodic pairs end");
                
                //all periodic boundary pairs are used
                logger.assert(pairingIterator == periodicPairing.end(), "Somehow missed some periodic pair");
                //the actual amount of vertices and the assigned amount of vertices match
                logger.assert(currentNodeNumber == theMesh_.getNumberOfNodes(IteratorType::GLOBAL), "Missed some node indexes");

                qHullCoordinates.reserveCoordinates(DIM * theMesh_.getNumberOfNodeCoordinates());
                for(Geometry::PointPhysical<DIM>& point : theMesh_.getNodeCoordinates())
                {
                    qHullCoordinates.append(DIM, point.data());
                }

                orgQhull::Qhull triangulation(qHullCoordinates, "d PF1e-10 QbB Qx Qc Qt");
                
                for (orgQhull::QhullFacet triangle : triangulation.facetList())
                {
                    if (triangle.isGood() && !triangle.isUpperDelaunay())
                    {
                        logger(DEBUG, "adding %", triangle);
                        Geometry::PointPhysical<DIM> center;
                        std::vector<std::size_t> pointIndices {};
                        bool shouldConnect = false;
                        for (auto vertexIt1 = triangle.vertices().begin(); vertexIt1 != triangle.vertices().end(); ++vertexIt1)
                        {
                            logger.assert((*vertexIt1).point().id() >= 0, "QHull breaks our assumptions on indexes");
                            logger(DEBUG, "% % %", shouldConnect, vertexIndex[(*vertexIt1).point().id()], *vertexIt1);
                            center += oldNodeLocations_[(*vertexIt1).point().id()];
                            pointIndices.push_back((*vertexIt1).point().id());
                            shouldConnect |= (std::find(dontConnect.begin(), dontConnect.end(), vertexIndex[(*vertexIt1).point().id()]) == dontConnect.end());
                        }
                        center = center / pointIndices.size();
                        if (domainDescription(center) < -1e-10 && shouldConnect)
                        {
                            auto newElement = addElement(pointIndices);
                            for (std::size_t i = 0; i < pointIndices.size(); ++i)
                            {
                                theMesh_.getNodesList(IteratorType::GLOBAL)[vertexIndex[pointIndices[i]]]->addElement(newElement, i);
                            }
                        }
                        else
                        {
                            logger(VERBOSE, "external element % ignored", triangle);
                        }
                    }
                    if (!triangle.isGood() && !triangle.isUpperDelaunay())
                    {
                        logger(VERBOSE, "small element % ignored", triangle);
                    }
                }
                for (Node* node : theMesh_.getNodesList(IteratorType::GLOBAL))
                {
                    //all of the nodes should be in the interior of the domain or near the boundary of the domain
                    if (node->getNumberOfElements() == 0)
                    {
                        for (std::size_t i = 0; i < vertexIndex.size(); ++i)
                        {
                            if (vertexIndex[i] == node->getID())
                            {
                                logger(DEBUG, "% % %", i, theMesh_.getNodeCoordinates()[i], domainDescription(theMesh_.getNodeCoordinates()[i]));
                            }
                        }
                    }
                    logger.assert(node->getNumberOfElements() > 0, "There is an node without any elements connected to it");
                }
                edgeFactory();
                faceFactory();
            }
            oldQuality = worstQuality;
            
            std::vector<double> expectedLength {};
            std::multimap<double, std::size_t> knownLengths {};
            expectedLength.reserve(theMesh_.getNumberOfNodes(IteratorType::GLOBAL));
            //for proper scaling
            totalCurrentLength = 0;
            double totalexpectedLength = 0.;
            for (Node* node : theMesh_.getNodesList(IteratorType::GLOBAL))
            {
                Geometry::PointPhysical<DIM> point = node->getElement(0)->getPhysicalGeometry()->getLocalNodeCoordinates(node->getNodeNumber(0));
                expectedLength.push_back(relativeEdgeLength(point));
                if (!isnan(expectedLength.back()) && !isinf(expectedLength.back()))
                {
                    knownLengths.insert( {expectedLength.back(), expectedLength.size() - 1});
                }
                else
                {
                    needsExpansion |= true;
                }
            }
            if (needsExpansion)
            {
                //iterate over all nodes, sorted by edge lengths
                for (std::pair<double, std::size_t> entry : knownLengths)
                {
                    Node* current = theMesh_.getNodesList(IteratorType::GLOBAL)[entry.second];
                    for (Element* element : current->getElements())
                    {
                        for (std::size_t i = 0; i < element->getNumberOfNodes(); ++i)
                        {
                            if (std::isnan(expectedLength[element->getNode(i)->getID()]) || std::isinf(expectedLength[element->getNode(i)->getID()]))
                            {
                                expectedLength[element->getNode(i)->getID()] = growFactor * entry.first;
                                
                                //inserting does not invalidate the iterators;
                                //new node has a larger edge length, so it is guaranteed to be visited later on
                                
                                knownLengths.insert(knownLengths.end(), {growFactor * entry.first, element->getNode(i)->getID()});
                            }
                        }
                    }
                }
            }
            
            //iterate over all edges to compute total length and scaling factor
            //the volume scales with (total edge length)^dimension
            //the total volume filled by the edges should be constant
            //so scale appropriately
            totalCurrentLength = 0;
            totalexpectedLength = 0;
            if (DIM == 1)
            {
                currentLength.resize(theMesh_.getNumberOfElements(IteratorType::GLOBAL));
                //the algorithm is mostly dimension independent, but the data type it operates on is not
                for (Element* element : theMesh_.getElementsList(IteratorType::GLOBAL))
                {
                    Geometry::PointPhysical<DIM> firstNode, secondNode;
                    firstNode = element->getPhysicalGeometry()->getLocalNodeCoordinates(0);
                    secondNode = element->getPhysicalGeometry()->getLocalNodeCoordinates(1);
                    currentLength[element->getID()] = L2Norm(firstNode - secondNode);
                    totalCurrentLength += currentLength[element->getID()];
                    totalexpectedLength += expectedLength[element->getNode(0)->getID()] / 2.;
                    totalexpectedLength += expectedLength[element->getNode(1)->getID()] / 2.;
                }
            }
            else if (DIM == 2)
            {
                currentLength.resize(theMesh_.getNumberOfFaces(IteratorType::GLOBAL));
                for (Face* face : theMesh_.getFacesList(IteratorType::GLOBAL))
                {
                    Geometry::PointPhysical<DIM> firstNode, secondNode;
                    std::vector<std::size_t> nodeIndices = face->getPtrElementLeft()->getReferenceGeometry()->getCodim1EntityLocalIndices(face->localFaceNumberLeft());
                    firstNode = face->getPtrElementLeft()->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[0]);
                    secondNode = face->getPtrElementLeft()->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[1]);
                    currentLength[face->getID()] = L2Norm(firstNode - secondNode);
                    totalCurrentLength += currentLength[face->getID()] * currentLength[face->getID()];
                    totalexpectedLength += std::pow(expectedLength[face->getPtrElementLeft()->getNode(nodeIndices[0])->getID()] + expectedLength[face->getPtrElementLeft()->getNode(nodeIndices[1])->getID()], 2.) / 4.;
                }
            }
            else
            {
                currentLength.resize(theMesh_.getNumberOfEdges(IteratorType::GLOBAL));
                for (Edge* edge : theMesh_.getEdgesList(IteratorType::GLOBAL))
                {
                    Geometry::PointPhysical<DIM> firstNode, secondNode;
                    std::vector<std::size_t> nodeIndices = edge->getElement(0)->getReferenceGeometry()->getCodim2EntityLocalIndices(edge->getEdgeNumber(0));
                    firstNode = edge->getElement(0)->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[0]);
                    secondNode = edge->getElement(0)->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[1]);
                    currentLength[edge->getID()] = L2Norm(firstNode - secondNode);
                    totalCurrentLength += currentLength[edge->getID()] * currentLength[edge->getID()] * currentLength[edge->getID()];
                    totalexpectedLength += std::pow(expectedLength[edge->getElement(0)->getNode(nodeIndices[0])->getID()] + expectedLength[edge->getElement(0)->getNode(nodeIndices[1])->getID()], 3.) / 8.;
                }
            }
            
            for (Geometry::PointPhysical<DIM>& point : movement)
            {
                point *= 0;
            }
            
            if (DIM == 1)
            {
                for (Element* element : theMesh_.getElementsList(IteratorType::GLOBAL))
                {
                    Geometry::PointPhysical<DIM> firstNode, secondNode;
                    firstNode = element->getPhysicalGeometry()->getLocalNodeCoordinates(0);
                    secondNode = element->getPhysicalGeometry()->getLocalNodeCoordinates(1);
                    //it is impossible to detect if a node inside the domain should be clipped to the edge
                    //instead make sure that the nodes that DO belong dont get pulled into the interior
                    //roundoff error should make sure that nodes move away from the boundary if there are too many
                    //all edges should be squeezed a little if the algorithm is to work correctly so pretend the volume is 1.4 times as large
                    //remember to scale back from a volume measure to a length measure
                    //the non-linearity makes everything slightly more robust
                    double length = (expectedLength[element->getNode(0)->getID()] + expectedLength[element->getNode(1)->getID()]) / currentLength[element->getID()] * 1.4 * totalCurrentLength / totalexpectedLength / 2.;
                    movement[element->getNode(0)->getID()] += std::max(length - 1., 0.) * (firstNode - secondNode) * (length + 1.) * 0.5;
                    movement[element->getNode(1)->getID()] += std::max(length - 1., 0.) * (secondNode - firstNode) * (length + 1.) * 0.5;
                }
            }
            else if (DIM == 2)
            {
                for (Face* face : theMesh_.getFacesList(IteratorType::GLOBAL))
                {
                    Geometry::PointPhysical<DIM> firstNode, secondNode;
                    std::vector<std::size_t> nodeIndices = face->getPtrElementLeft()->getReferenceGeometry()->getCodim1EntityLocalIndices(face->localFaceNumberLeft());
                    firstNode = face->getPtrElementLeft()->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[0]);
                    secondNode = face->getPtrElementLeft()->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[1]);
                    double length = (expectedLength[face->getPtrElementLeft()->getNode(nodeIndices[0])->getID()] + expectedLength[face->getPtrElementLeft()->getNode(nodeIndices[1])->getID()]) / currentLength[face->getID()] * std::pow(1.4 * totalCurrentLength / totalexpectedLength, 1. / 2.) / 2.;
                    movement[face->getPtrElementLeft()->getNode(nodeIndices[0])->getID()] += std::max(length - 1., 0.) * (firstNode - secondNode) * (length + 1.) * 0.5;
                    movement[face->getPtrElementLeft()->getNode(nodeIndices[1])->getID()] += std::max(length - 1., 0.) * (secondNode - firstNode) * (length + 1.) * 0.5;
                }
            }
            else if (DIM == 3)
            {
                for (Edge* edge : theMesh_.getEdgesList(IteratorType::GLOBAL))
                {
                    Geometry::PointPhysical<DIM> firstNode, secondNode;
                    std::vector<std::size_t> nodeIndices = edge->getElement(0)->getReferenceGeometry()->getCodim2EntityLocalIndices(edge->getEdgeNumber(0));
                    firstNode = edge->getElement(0)->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[0]);
                    secondNode = edge->getElement(0)->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[1]);
                    double length = (expectedLength[edge->getElement(0)->getNode(nodeIndices[0])->getID()] + expectedLength[edge->getElement(0)->getNode(nodeIndices[1])->getID()]) / currentLength[edge->getID()] * std::pow(1.4 * totalCurrentLength / totalexpectedLength, 1. / 3.) / 2.;
                    movement[edge->getElement(0)->getNode(nodeIndices[0])->getID()] += std::max(length - 1., 0.) * (firstNode - secondNode) * (length + 1.) * 0.5;
                    movement[edge->getElement(0)->getNode(nodeIndices[1])->getID()] += std::max(length - 1., 0.) * (secondNode - firstNode) * (length + 1.) * 0.5;
                }
            }
            
            //forward Euler discretisation of an optimally damped mass-spring system, with time step 0.02
            //this time step could be 0.1, but there is a stability issue where springs aligned along the periodic boundary are applied twice
            maxMovement = 0;
            auto moveIterator = movement.begin();
            auto fixIterator = fixedPointIdxs.begin();
            for (std::size_t i = 0; i < theMesh_.getNumberOfNodes(IteratorType::GLOBAL); ++moveIterator, ++i)
            {
                if (fixIterator != fixedPointIdxs.end() && i == *fixIterator)
                {
                    ++fixIterator;
                    *moveIterator *= 0;
                }
                else
                {
                    Node* node = theMesh_.getNodesList(IteratorType::GLOBAL)[i];
                    Geometry::PointPhysical<DIM>& point = theMesh_.getNodeCoordinates()[node->getElement(0)->getPhysicalGeometry()->getNodeIndex(node->getNodeNumber(0))];
                    point += 0.1 * (*moveIterator);
                    logger.assert(!(std::isnan(point[0])), "%", i);
                    bool isPeriodic = false;
                    std::map<std::size_t, bool> hasMoved {};
                    hasMoved[node->getElement(0)->getPhysicalGeometry()->getNodeIndex(node->getNodeNumber(0))] = true;
                    for (std::size_t j = 1; j < node->getNumberOfElements(); ++j)
                    {
                        if (!hasMoved[node->getElement(j)->getPhysicalGeometry()->getNodeIndex(node->getNodeNumber(j))])
                        {
                            Geometry::PointPhysical<DIM>& other = theMesh_.getNodeCoordinates()[node->getElement(j)->getPhysicalGeometry()->getNodeIndex(node->getNodeNumber(j))];
                            other += 0.1 * (*moveIterator);
                            hasMoved[node->getElement(j)->getPhysicalGeometry()->getNodeIndex(node->getNodeNumber(j))] = true;
                            isPeriodic = true;
                        }
                    }
                    if (domainDescription(point) > 0 && !isPeriodic)
                    {
                        //the point is outside of the domain, move it back inside
                        double currentValue = domainDescription(point);
                        LinearAlgebra::SmallVector<DIM> gradient;
                        LinearAlgebra::SmallVector<DIM> offset;
                        //one-sided numerical derivative
                        for (std::size_t j = 0; j < DIM; ++j)
                        {
                            offset[j] = 1e-7;
                            gradient[j] = (currentValue - domainDescription(point + offset)) * 1e7;
                            offset[j] = 0;
                        }
                        point += currentValue * gradient / L2Norm(gradient);
                        *moveIterator += 10 * currentValue * gradient / L2Norm(gradient);
                        currentValue = domainDescription(point);
                        //second step for robustness and accuracy if needed
                        if (currentValue > 0)
                        {
                            for (std::size_t j = 0; j < DIM; ++j)
                            {
                                offset[j] = 1e-7;
                                gradient[j] = (currentValue - domainDescription(point + offset)) * 1e7;
                                offset[j] = 0;
                            }
                            point += currentValue * gradient / L2Norm(gradient);
                            *moveIterator += 10 * currentValue * gradient / L2Norm(gradient);
                            //if two steps are not enough, more are also not likely to help
                            currentValue = domainDescription(point);
                            if (currentValue > 1e-10)
                            {
                                logger(WARN, "NOTE: Failed to move point % (%) back into the domain."
                                        "\n Distance from boundary is %. Algorithm may crash.\n Consider fixing "
                                        "points at corners to remedy this issue.", i, point, currentValue);
                            }
                        }
                    }
                    if (isPeriodic)
                    {
                        //do a total of four newton iteration before giving up
                        Geometry::PointPhysical<DIM> testPoint;
                        for (std::size_t j = 0; j < 4; ++j)
                        {
                            //make sure the node stays on the periodic boundary, to prevent faces with 3 or more elements connected to them
                            for (std::size_t k = 0; k < node->getNumberOfElements(); ++k)
                            {
                                testPoint = node->getElement(k)->getPhysicalGeometry()->getLocalNodeCoordinates(node->getNodeNumber(k));
                                double currentValue = domainDescription(testPoint);
                                if (currentValue > 0)
                                {
                                    LinearAlgebra::SmallVector<DIM> gradient;
                                    LinearAlgebra::SmallVector<DIM> offset;
                                    for (std::size_t l = 0; l < DIM; ++l)
                                    {
                                        offset[l] = 1e-7;
                                        gradient[l] = (currentValue - domainDescription(testPoint + offset)) * 1e7;
                                        offset[l] = 0;
                                    }
                                    hasMoved.clear();
                                    for (std::size_t l = 0; l < node->getNumberOfElements(); ++l)
                                    {
                                        if (!hasMoved[node->getElement(l)->getPhysicalGeometry()->getNodeIndex(node->getNodeNumber(l))])
                                        {
                                            Geometry::PointPhysical<DIM>& other = theMesh_.getNodeCoordinates()[node->getElement(l)->getPhysicalGeometry()->getNodeIndex(node->getNodeNumber(l))];
                                            other += currentValue * gradient / L2Norm(gradient);
                                            hasMoved[node->getElement(l)->getPhysicalGeometry()->getNodeIndex(node->getNodeNumber(l))] = true;
                                        }
                                    }
                                    *moveIterator += 10 * currentValue * gradient / L2Norm(gradient);
                                }
                            }
                        }
                        for (std::size_t j = 0; j < node->getNumberOfElements(); ++j)
                        {
                            testPoint = node->getElement(j)->getPhysicalGeometry()->getLocalNodeCoordinates(node->getNodeNumber(j));
                            if (domainDescription(testPoint) > 1e-10)
                            {
                                logger(WARN, "NOTE: Failed to move periodic testPoint % (%) back to the periodic boundary.\n "
                                        "Distance from boundary is %. Algorithm may crash.\n "
                                        "Consider fixing points at corners to remedy this issue.", i, testPoint, domainDescription(testPoint));
                            }
                        };
                    }
                }
            }
            
            worstQuality = 1;
            if (DIM == 1)
            {
                //quality measure is not an issue in 1D just create a mesh with the proper lengths
                worstQuality = 0.5;
                //the algorithm is mostly dimension independent, but the data type it operates on is not
                for (Element* element : theMesh_.getElementsList(IteratorType::GLOBAL))
                {
                    Geometry::PointPhysical<DIM> firstNode, secondNode;
                    firstNode = element->getPhysicalGeometry()->getLocalNodeCoordinates(0);
                    secondNode = element->getPhysicalGeometry()->getLocalNodeCoordinates(1);
                    maxMovement = std::max(maxMovement, L2Norm(movement[element->getNode(0)->getID()]) / 10 / currentLength[element->getID()]);
                    maxMovement = std::max(maxMovement, L2Norm(movement[element->getNode(1)->getID()]) / 10 / currentLength[element->getID()]);
                    maxShift = std::max(maxShift, L2Norm(firstNode - oldNodeLocations_[element->getPhysicalGeometry()->getNodeIndex(0)]) / currentLength[element->getID()]);
                    maxShift = std::max(maxShift, L2Norm(secondNode - oldNodeLocations_[element->getPhysicalGeometry()->getNodeIndex(1)]) / currentLength[element->getID()]);
                }
            }
            else if (DIM == 2)
            {
                //ratio between incircle and circumcircle (scaled so equilateral is quality 1 and reference is quality ~.8)
                for (Element* element : theMesh_.getElementsList(IteratorType::GLOBAL))
                {
                    std::array<double, 3> edgeLengths {};
                    for (std::size_t i = 0; i < 3; ++i)
                    {
                        edgeLengths[i] = currentLength[element->getFace(i)->getID()];
                    }
                    worstQuality = std::min(worstQuality, (edgeLengths[0] + edgeLengths[1] - edgeLengths[2]) * (edgeLengths[1] + edgeLengths[2] - edgeLengths[0]) * (edgeLengths[2] + edgeLengths[0] - edgeLengths[1]) / edgeLengths[0] / edgeLengths[1] / edgeLengths[2]);
                }
                for (Face* face : theMesh_.getFacesList(IteratorType::GLOBAL))
                {
                    Geometry::PointPhysical<DIM> firstNode, secondNode;
                    std::vector<std::size_t> nodeIndices = face->getPtrElementLeft()->getReferenceGeometry()->getCodim1EntityLocalIndices(face->localFaceNumberLeft());
                    firstNode = face->getPtrElementLeft()->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[0]);
                    secondNode = face->getPtrElementLeft()->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[1]);
                    maxMovement = std::max(maxMovement, L2Norm(movement[face->getPtrElementLeft()->getNode(nodeIndices[0])->getID()]) / 10 / currentLength[face->getID()]);
                    maxMovement = std::max(maxMovement, L2Norm(movement[face->getPtrElementLeft()->getNode(nodeIndices[1])->getID()]) / 10 / currentLength[face->getID()]);
                    maxShift = std::max(maxShift, L2Norm(firstNode - oldNodeLocations_[face->getPtrElementLeft()->getPhysicalGeometry()->getNodeIndex(nodeIndices[0])]) / currentLength[face->getID()]);
                    maxShift = std::max(maxShift, L2Norm(secondNode - oldNodeLocations_[face->getPtrElementLeft()->getPhysicalGeometry()->getNodeIndex(nodeIndices[1])]) / currentLength[face->getID()]);
                }
            }
            else
            {
                //ratio between volume and cubed average edge length (scaled so equilateral is quality 1 and reference is quality ~.8)
                for (Element* element : theMesh_.getElementsList(IteratorType::GLOBAL))
                {
                    std::array<double, 6> edgeLengths {};
                    for (std::size_t i = 0; i < 6; ++i)
                    {
                        edgeLengths[i] = currentLength[element->getEdge(i)->getID()];
                    }
                    double average = std::accumulate(edgeLengths.begin(), edgeLengths.end(), 0) / 6;
                    const Geometry::PointReference<DIM>& center = element->getReferenceGeometry()->getCenter();
                    Geometry::Jacobian<DIM, DIM> jac = element->calcJacobian(center);
                    worstQuality = std::min(worstQuality, jac.determinant() / average * std::sqrt(2));
                }
                for (Edge* edge : theMesh_.getEdgesList(IteratorType::GLOBAL))
                {
                    Geometry::PointPhysical<DIM> firstNode, secondNode;
                    std::vector<std::size_t> nodeIndices = edge->getElement(0)->getReferenceGeometry()->getCodim2EntityLocalIndices(edge->getEdgeNumber(0));
                    firstNode = edge->getElement(0)->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[0]);
                    secondNode = edge->getElement(0)->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[1]);
                    maxMovement = std::max(maxMovement, L2Norm(movement[edge->getElement(0)->getNode(nodeIndices[0])->getID()]) / 10 / currentLength[edge->getID()]);
                    maxMovement = std::max(maxMovement, L2Norm(movement[edge->getElement(0)->getNode(nodeIndices[1])->getID()]) / 10 / currentLength[edge->getID()]);
                    maxShift = std::max(maxShift, L2Norm(firstNode - oldNodeLocations_[edge->getElement(0)->getPhysicalGeometry()->getNodeIndex(nodeIndices[0])]) / currentLength[edge->getID()]);
                    maxShift = std::max(maxShift, L2Norm(secondNode - oldNodeLocations_[edge->getElement(0)->getPhysicalGeometry()->getNodeIndex(nodeIndices[1])]) / currentLength[edge->getID()]);
                }
            }
            
            //no teleporting nodes in the final iteration
            ///\todo temporarily toggled off for debug reasons
            if (counter % 50 == 1 && false)
            {
                //the actual sorting is more expensive than computing the lengths and this does not happen very often
                std::multimap<double, std::pair<Geometry::PointPhysical<DIM>, std::size_t > > centerPoints {};
                if (DIM == 1)
                {
                    //the algorithm is mostly dimension independent, but the data type it operates on is not
                    for (Element* element : theMesh_.getElementsList(IteratorType::GLOBAL))
                    {
                        Geometry::PointPhysical<DIM> firstNode, secondNode;
                        firstNode = element->getPhysicalGeometry()->getLocalNodeCoordinates(0);
                        secondNode = element->getPhysicalGeometry()->getLocalNodeCoordinates(1);
                        //length is scaled in case somebody hasty decides to add smoothing at this point
                        //all edges should be squeezed a little if the algorithm is to work correctly so pretend the volume is 1.5 times as large
                        //remember to scale back from a volume measure to a length measure
                        double length = (expectedLength[element->getNode(0)->getID()] + expectedLength[element->getNode(1)->getID()]) / currentLength[centerPoints.size()] * 2 * totalCurrentLength / totalexpectedLength;
                        centerPoints.insert( {length, {(firstNode + secondNode) / 2, element->getNode(0)->getID()}});
                    }
                }
                else if (DIM == 2)
                {
                    for (Face* face : theMesh_.getFacesList(IteratorType::GLOBAL))
                    {
                        Geometry::PointPhysical<DIM> firstNode, secondNode;
                        std::vector<std::size_t> nodeIndices = face->getPtrElementLeft()->getReferenceGeometry()->getCodim1EntityLocalIndices(face->localFaceNumberLeft());
                        firstNode = face->getPtrElementLeft()->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[0]);
                        secondNode = face->getPtrElementLeft()->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[1]);
                        //length is scaled in case somebody hasty decides to add smoothing at this point
                        //all edges should be squeezed a little if the algorithm is to work correctly so pretend the volume is 1.5 times as large
                        //remember to scale back from a volume measure to a length measure
                        double length = (expectedLength[face->getPtrElementLeft()->getNode(nodeIndices[0])->getID()] + expectedLength[face->getPtrElementLeft()->getNode(nodeIndices[1])->getID()]) / currentLength[centerPoints.size()] * std::pow(2 * totalCurrentLength / totalexpectedLength, 1. / 2.);
                        centerPoints.insert( {length, {(firstNode + secondNode) / 2, face->getPtrElementLeft()->getNode(nodeIndices[0])->getID()}});
                    }
                }
                else
                {
                    for (Edge* edge : theMesh_.getEdgesList(IteratorType::GLOBAL))
                    {
                        Geometry::PointPhysical<DIM> firstNode, secondNode;
                        std::vector<std::size_t> nodeIndices = edge->getElement(0)->getReferenceGeometry()->getCodim2EntityLocalIndices(edge->getEdgeNumber(0));
                        firstNode = edge->getElement(0)->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[0]);
                        secondNode = edge->getElement(0)->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[1]);
                        //length is scaled in case somebody hasty decides to add smoothing at this point
                        //all edges should be squeezed a little if the algorithm is to work correctly so pretend the volume is 1.5 times as large
                        //remember to scale back from a volume measure to a length measure
                        double length = (expectedLength[edge->getElement(0)->getNode(nodeIndices[0])->getID()] + expectedLength[edge->getElement(0)->getNode(nodeIndices[1])->getID()]) / currentLength[centerPoints.size()] * std::pow(2 * totalCurrentLength / totalexpectedLength, 1. / 3.);
                        centerPoints.insert( {length, {(firstNode + secondNode) / 2, edge->getElement(0)->getNode(nodeIndices[0])->getID()}});
                    }
                }
                std::vector<bool> hasTeleported(theMesh_.getNumberOfNodeCoordinates(), false);
                auto longEdge = centerPoints.begin();
                auto shortEdge = centerPoints.rbegin();
                for (std::size_t index : fixedPointIdxs)
                {
                    hasTeleported[index] = true;
                }
                Geometry::PointPhysical<DIM> point;
                Geometry::PointPhysical<DIM> other;
                for (Node* node : theMesh_.getNodesList(IteratorType::GLOBAL))
                {
                    point = node->getElement(0)->getPhysicalGeometry()->getLocalNodeCoordinates(node->getNodeNumber(0));
                    for (std::size_t i = 0; i < node->getNumberOfElements(); ++i)
                    {
                        other = node->getElement(i)->getPhysicalGeometry()->getLocalNodeCoordinates(node->getNodeNumber(i));
                        if (point != other)
                        {
                            hasTeleported[node->getID()] = true;
                        }
                    }
                }
                //remember that the size measure is inverted
                while (3 * longEdge->first < shortEdge->first)
                {
                    if (hasTeleported[shortEdge->second.second])
                    {
                        shortEdge++;
                    }
                    else
                    {
                        if (domainDescription(longEdge->second.first) < 0)
                        {
                            maxMovement = std::max(maxMovement, L2Norm(longEdge->second.first - theMesh_.getNodeCoordinates()[shortEdge->second.second]));
                            //it is quite unlikely that the current triangulation suffices after randomly teleporting nodes about
                            maxShift = std::numeric_limits<double>::infinity();
                            theMesh_.getNodeCoordinates()[shortEdge->second.second] = longEdge->second.first;
                            hasTeleported[shortEdge->second.second] = true;
                            shortEdge++;
                        }
                        longEdge++;
                    }
                }
            }
        }
        if (counter == 10000)
        {
            logger(WARN, "WARNING: Maximum iteration count reached, mesh quality may not be optimal");
        }
        //coordinate transformation may have changed, update to the current situation
        for (Element* element : theMesh_.getElementsList())
        {
            element->getReferenceToPhysicalMap()->reinit();
        }
    }
#endif
    
    /// \bug does not do the bc flags yet
    template<std::size_t DIM>
    void MeshManipulator<DIM>::faceFactory()
    {   
        std::vector<std::size_t> nodeIndices;
        std::vector<Element*> candidates;
        
        for (Element* element : theMesh_.getElementsList(IteratorType::GLOBAL))
        {
            for (std::size_t i = 0; i < element->getNumberOfFaces(); ++i)
            {
                std::vector<const Node*> localNodes;
                //if this face is not there yet
                if (element->getFace(i) == nullptr)
                {
                    localNodes.clear();
                    candidates.clear();
                    nodeIndices = element->getReferenceGeometry()->getCodim1EntityLocalIndices(i);
                    
                    candidates = element->getNode(nodeIndices[0])->getElements();
                    localNodes.push_back(element->getNode(nodeIndices[0]));
                    std::sort(candidates.begin(), candidates.end(), [](Element* left, Element* right)
                    {   return left->getID()<right->getID();});

                    logger(DEBUG, "Candidates: ");
                    for(Element* coutElement : candidates)
                    {
                        logger(DEBUG, "Element %: %", coutElement->getID(), *coutElement);
                    }
                    for (std::size_t j = 1; j < nodeIndices.size(); ++j)
                    {
                        localNodes.push_back(element->getNode(nodeIndices[j]));
                        std::vector<Element*> temp, nextIndices;
                        nextIndices = element->getNode(nodeIndices[j])->getElements();
                        std::sort(nextIndices.begin(), nextIndices.end(), [](Element* left, Element* right)
                        {   return left->getID()<right->getID();});
                        std::set_intersection(candidates.begin(), candidates.end(), nextIndices.begin(), nextIndices.end(), std::back_inserter(temp), [](Element* left, Element* right)
                        {   return left->getID()<right->getID();});
                        candidates = std::move(temp);
                        logger(DEBUG, "Candidates: ");
                        for(Element* coutElement : candidates)
                        {
                            logger(DEBUG, "Element %: %", coutElement->getID(), *coutElement);
                        }
                    }
                    
                    //the current element does not bound the face or more than two elements bound the face
                    logger.assert_always(candidates.size() == 1 || candidates.size() == 2, 
                                         "Detected % bounding elements for face %, which is impossible", candidates.size(),
                                         theMesh_.getFacesList(IteratorType::GLOBAL).size() + 1);
                    //boundary face
                    if (candidates.size() == 1)
                    {
                        logger.assert(candidates[0] == element, "dropped the original element");
                        addFace(element, i, nullptr, 0, Geometry::FaceType::WALL_BC);
                    }
                    if (candidates.size() == 2)
                    {
                        Element* other;
                        if (candidates[0] == element)
                        {
                            other = candidates[1];
                        }
                        else
                        {
                            other = candidates[0];
                        }
                        bool matchFound = false;
                        std::vector<std::size_t> otherNodeIndices;
                        for (std::size_t j = 0; j < other->getNumberOfFaces(); ++j)
                        {
                            otherNodeIndices = other->getReferenceGeometry()->getCodim1EntityLocalIndices(j);
                            bool match = true;
                            for (std::size_t k : otherNodeIndices)
                            {
                                if (std::find(localNodes.begin(), localNodes.end(), other->getNode(k)) == localNodes.end())
                                {
                                    match = false;
                                }
                            }
                            if (match)
                            {
                                logger.assert(!matchFound,"Found two opposing faces for face % " 
                                     " of element % in opposing element %." 
                                    , i, element->getID(),other->getID());
                                addFace(element, i, other, j);
                                matchFound = true;
                            }
                        }
                        logger.assert(matchFound, "Could not find matching face for face % "
                                "of element % in opposing element %." , i, element->getID(), other->getID());
                    }
                }
            }
        }
        
        logger(VERBOSE, "Total number of Faces: %", getFacesList(IteratorType::GLOBAL).size());
    }
    
    //the algorithm for the edge factory is based on that of the face factory
    //with some minor adaptation to account for the fact that there may be
    //more than two elements per edge
    ///\bug does not do 4D yet
    template<std::size_t DIM>
    void MeshManipulator<DIM>::edgeFactory()
    {
        //'edges' in DIM 2 are actually nodes
        if (DIM != 2)
        {
            std::vector<std::size_t> nodeList, otherNodeList;
            
            const Node* nodes[2];
            
            for (Element* element : theMesh_.getElementsList(IteratorType::GLOBAL))
            {
                for (std::size_t i = 0; i < element->getNumberOfEdges(); ++i)
                {
                    if (element->getEdge(i) == nullptr)
                    {
                        nodeList = element->getReferenceGeometry()->getCodim2EntityLocalIndices(i);
                        std::vector<Element*> candidates(0);
                        auto& leftElements = element->getNode(nodeList[0])->getElements();
                        auto& rightElements = element->getNode(nodeList[1])->getElements();
                        std::set_intersection(leftElements.begin(), leftElements.end(), rightElements.begin(), rightElements.end(), std::back_inserter(candidates), [](Element* a, Element* b)
                        {   return a->getID()<b->getID();});
                        logger.assert(candidates.size() > 0, "current element is not adjacent to its own edges"); 
                        addEdge();
                        nodes[0] = element->getNode(nodeList[0]);
                        nodes[1] = element->getNode(nodeList[1]);
                        Edge* newEdge = *(--theMesh_.getEdgesList(IteratorType::GLOBAL).end());
                        newEdge->addElement(element, i);
                        for (std::size_t j = 1; j < candidates.size(); ++j)
                        {
                            Element* other = candidates[j];
                            for (std::size_t k = 0; k < other->getNumberOfEdges(); ++k)
                            {
                                otherNodeList = other->getReferenceGeometry()->getCodim2EntityLocalIndices(k);
                                if ((other->getNode(otherNodeList[0]) == nodes[0] || other->getNode(otherNodeList[0]) == nodes[1]) && (other->getNode(otherNodeList[1]) == nodes[0] || other->getNode(otherNodeList[1]) == nodes[1]))
                                {
                                    newEdge->addElement(other, k);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    template<std::size_t DIM>
    Mesh<DIM>& MeshManipulator<DIM>::getMesh()
    {
        return theMesh_;
    }

    template<std::size_t DIM>
    const Mesh<DIM>& MeshManipulator<DIM>::getMesh() const
    {
        return theMesh_;
    }
}

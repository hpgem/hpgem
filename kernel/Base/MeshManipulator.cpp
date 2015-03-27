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
#include "Logger.h"

#include <algorithm>
#include <iostream>
#include <unordered_set>
#include <array>
#include <numeric>

namespace Base
{
    
    void MeshManipulator::createDefaultBasisFunctions(std::size_t order)
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
        collBasisFSet_[0] = bFset1;
    }
    
    MeshManipulator::MeshManipulator(const ConfigurationData* config, BoundaryType xPer, BoundaryType yPer, BoundaryType zPer, std::size_t orderOfFEM, std::size_t idRangeBegin, std::size_t nrOfElementMatrixes, std::size_t nrOfElementVectors, std::size_t nrOfFaceMatrtixes, std::size_t nrOfFaceVectors)
            : configData_(config),
            meshMover_(nullptr), numberOfElementMatrixes_(nrOfElementMatrixes), numberOfFaceMatrixes_(nrOfFaceMatrtixes), numberOfElementVectors_(nrOfElementVectors), numberOfFaceVectors_(nrOfFaceVectors)
    {
        logger.assert(config!=nullptr, "Invalid configuration passed");
        logger.assert(orderOfFEM==config->polynomialOrder_, "Inconsistent redundant information passed");
        logger.assert(idRangeBegin==0, "c++ starts counting at 0");
        logger(INFO, "******Mesh creation started!**************");
        std::size_t DIM = configData_->dimension_;
        periodicX_ = (xPer == BoundaryType::PERIODIC);
        periodicY_ = (yPer == BoundaryType::PERIODIC);
        periodicZ_ = (zPer == BoundaryType::PERIODIC);
        for (std::size_t i = 0; i < DIM; ++i)
        {
            if (i == 0)
                logger(INFO, "Boundries: % in X direction", (periodicX_ ? "Periodic  " : "Solid Wall"));
            if (i == 1)
                logger(INFO, "Boundries: % in Y direction", (periodicY_ ? "Periodic  " : "Solid Wall"));
            if (i == 2)
                logger(INFO, "Boundries: % in Z direction", (periodicZ_ ? "Periodic  " : "Solid Wall"));
        }
        createDefaultBasisFunctions(orderOfFEM);
        const_cast<ConfigurationData*>(configData_)->numberOfBasisFunctions_ = collBasisFSet_[0]->size();
        
        logger(INFO, "******Mesh creation is finished!**********");
        logger(VERBOSE, "nrOfElementVector = %", nrOfElementVectors);
        logger(VERBOSE, "nrOfFaceVector = %", nrOfFaceVectors);
    }
    
    MeshManipulator::MeshManipulator(const MeshManipulator& other)
            : theMesh_(other.theMesh_), configData_(other.configData_), periodicX_(other.periodicX_), periodicY_(other.periodicY_), periodicZ_(other.periodicZ_), meshMover_(other.meshMover_),
            collBasisFSet_(other.collBasisFSet_),
            numberOfElementMatrixes_(other.numberOfElementMatrixes_), numberOfFaceMatrixes_(other.numberOfFaceMatrixes_), numberOfElementVectors_(other.numberOfElementVectors_), numberOfFaceVectors_(other.numberOfFaceVectors_)
    {        
    }
    
    MeshManipulator::~MeshManipulator()
    {
        
        for (typename CollectionOfBasisFunctionSets::iterator bit = collBasisFSet_.begin(); bit != collBasisFSet_.end(); ++bit)
        {
            const BasisFunctionSetT* bf = *bit;
            ///\bug segfaults when using two meshes with the same sets of basisfunctions
            delete bf; 
        }
        
        delete meshMover_;
        
    }
    
    void MeshManipulator::setDefaultBasisFunctionSet(BasisFunctionSetT* bFSet)
    {
        logger.assert(bFSet!=nullptr, "Invalid basis function set passed");
        delete collBasisFSet_[0];
        collBasisFSet_[0] = bFSet;
        const_cast<ConfigurationData*>(configData_)->numberOfBasisFunctions_ = bFSet->size();
        for (Base::Face* face : getFacesList(IteratorType::GLOBAL))
        {
            face->setLocalNrOfBasisFunctions(0);
        }
        for (Base::Edge* edge : getEdgesList(IteratorType::GLOBAL))
        {
            edge->setLocalNrOfBasisFunctions(0);
        }
        for (Base::Node* node : getVerticesList(IteratorType::GLOBAL))
        {
            node->setLocalNrOfBasisFunctions(0);
        }
        for (ElementIterator it = elementColBegin(IteratorType::GLOBAL); it != elementColEnd(IteratorType::GLOBAL); ++it)
        {
            (*it)->setDefaultBasisFunctionSet(0);
        }
    }
    
    void MeshManipulator::addVertexBasisFunctionSet(const CollectionOfBasisFunctionSets& bFsets)
    {
        std::size_t firstNewEntry = collBasisFSet_.size();
        for (const BasisFunctionSet* set : bFsets)
        {
            logger.assert(set!=nullptr, "Invalid basis function set detected");
            collBasisFSet_.push_back(set);
        }
        for (Node* node : getVerticesList())
        {
            for (std::size_t i = 0; i < node->getNrOfElements(); ++i)
            {
                node->getElement(i)->setVertexBasisFunctionSet(firstNewEntry + node->getVertexNr(i), node->getVertexNr(i));
            }
            node->setLocalNrOfBasisFunctions(bFsets[0]->size());
        }
        const_cast<ConfigurationData*>(configData_)->numberOfBasisFunctions_ += (*elementColBegin())->getNrOfNodes() * bFsets[0]->size();
    }
    
    void MeshManipulator::addFaceBasisFunctionSet(const std::vector<const OrientedBasisFunctionSet*>& bFsets)
    {
        std::size_t firstNewEntry = collBasisFSet_.size();
        for (const BasisFunctionSet* set : bFsets)
        {
            logger.assert(set!=nullptr, "Invalid basis function set detected");
            collBasisFSet_.push_back(set);
        }
        for (Face* face : getFacesList())
        {
            std::size_t faceNr = face->localFaceNumberLeft();
            for (std::size_t i = 0; i < bFsets.size(); ++i)
            {
                if (bFsets[i]->checkOrientation(0, faceNr))
                {
                    face->getPtrElementLeft()->setFaceBasisFunctionSet(firstNewEntry + i, faceNr);
                }
            }
            if (face->isInternal())
            {
                faceNr = face->localFaceNumberRight();
                int orientation = face->getFaceToFaceMapIndex();
                for (std::size_t i = 0; i < bFsets.size(); ++i)
                {
                    if (bFsets[i]->checkOrientation(orientation, faceNr))
                    {
                        face->getPtrElementRight()->setFaceBasisFunctionSet(firstNewEntry + i, faceNr);
                    }
                }
            }
            face->setLocalNrOfBasisFunctions(bFsets[0]->size());
        }
        const_cast<ConfigurationData*>(configData_)->numberOfBasisFunctions_ += (*elementColBegin())->getPhysicalGeometry()->getNrOfFaces() * bFsets[0]->size();
    }
    
    void MeshManipulator::addEdgeBasisFunctionSet(const std::vector<const OrientedBasisFunctionSet*>& bFsets)
    {
        std::size_t firstNewEntry = collBasisFSet_.size();
        for (const BasisFunctionSet* set : bFsets)
        {
            logger.assert(set!=nullptr, "Invalid basis function set detected");
            collBasisFSet_.push_back(set);
        }
        logger(DEBUG, "In MeshManipulator::addEdgeBasisFunctionSet: ");
        for (Edge* edge : getEdgesList())
        {
            for (std::size_t i = 0; i < edge->getNrOfElements(); ++i)
            {
                for (std::size_t j = 0; j < bFsets.size(); ++j)
                {
                    if (bFsets[j]->checkOrientation(edge->getOrientation(i), edge->getEdgeNr(i)))
                    {
                        edge->getElement(i)->setEdgeBasisFunctionSet(firstNewEntry + j, edge->getEdgeNr(i));
                        logger(DEBUG, "% % %", edge->getOrientation(i), 
                               edge->getEdgeNr(i), bFsets[j]->size());
                    }
                }
            }
            edge->setLocalNrOfBasisFunctions(bFsets[0]->size());
        }
        const_cast<ConfigurationData*>(configData_)->numberOfBasisFunctions_ += (*elementColBegin())->getPhysicalGeometry()->getNrOfFaces() * bFsets[0]->size();
    }
    
    Base::Element*
    MeshManipulator::addElement(const VectorOfPointIndicesT& globalNodeIndexes)
    {
        return theMesh_.addElement(globalNodeIndexes);
    }
    
    void MeshManipulator::move()
    {
        for (Geometry::PointPhysical& p : theMesh_.getNodes())
        {
            if (meshMover_ != nullptr)
            {
                meshMover_->movePoint(p);
            }
        }
    }
    
    void MeshManipulator::setMeshMover(const MeshMoverBase* meshMover)
    {
        //can be set to nullptr if you dont want to move the mesh anymore
        meshMover_ = meshMover;
    }
    
    bool MeshManipulator::addFace(ElementT* leftElementPtr, std::size_t leftElementLocalFaceNo, ElementT* rightElementPtr, std::size_t rightElementLocalFaceNo, const Geometry::FaceType& faceType)
    {
        logger.assert(leftElementPtr!=nullptr, "Invalid element passed");
        //rightElementPtr may be nullptr for boundary faces
        return theMesh_.addFace(leftElementPtr, leftElementLocalFaceNo, rightElementPtr, rightElementLocalFaceNo, faceType);
    }
    
    void MeshManipulator::addEdge()
    {
        theMesh_.addEdge();
    }
    
    void MeshManipulator::addVertex()
    {
        theMesh_.addVertex();
    }
}

std::ostream& operator<<(std::ostream& os, const Base::MeshManipulator& mesh)
{
    for (Geometry::PointPhysical p : mesh.getNodes())
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
    
    void MeshManipulator::createRectangularMesh(const PointPhysicalT& bottomLeft, const PointPhysicalT& topRight, const VectorOfPointIndicesT& linearNoElements)
    {
        logger.assert(bottomLeft.size()==topRight.size(), "The corners of the mesh must have the same dimension");
        logger.assert(bottomLeft.size()==configData_->dimension_, "The corners of the mesh have the wrong dimension");
        logger.assert(linearNoElements.size()==configData_->dimension_, "There are amounts of elements spicified in % dimensions, but there are % dimensions", linearNoElements.size(), configData_->dimension_);
        //set to correct value in case some other meshmanipulator changed things
        ElementFactory::instance().setCollectionOfBasisFunctionSets(&collBasisFSet_);
        ElementFactory::instance().setNumberOfMatrices(numberOfElementMatrixes_);
        ElementFactory::instance().setNumberOfVectors(numberOfFaceVectors_);
        ElementFactory::instance().setNumberOfTimeLevels(configData_->numberOfTimeLevels_);
        ElementFactory::instance().setNumberOfUnknowns(configData_->numberOfUnknowns_);
        FaceFactory::instance().setNumberOfFaceMatrices(numberOfFaceMatrixes_);
        FaceFactory::instance().setNumberOfFaceVectors(numberOfFaceVectors_);
        
        std::size_t DIM = configData_->dimension_;
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
        PointPhysicalT delta_x(DIM);
        
        for (std::size_t i = 0; i < DIM; i++)
        {
            delta_x[i] = (topRight[i] - bottomLeft[i]) / (linearNoElements[i]);
        }
        
        //This stores the number of nodes in each coDIMension i.e. if you have 2 by 2 element it is 3 nodes 
        //nodes mark physical location, vertices mark connectivity-based location
        std::vector<std::size_t> numOfNodesInEachSubspace(DIM), numOfElementsInEachSubspace(DIM), numOfVerticesInEachSubspace(DIM);
        
        numOfNodesInEachSubspace[0] = 1;
        numOfVerticesInEachSubspace[0] = 1;
        numOfElementsInEachSubspace[0] = 1;
        
        //This will be the total number of nodes required in the problem
        std::size_t totalNumOfNodes, totalNumOfVertices, totalNumOfElements, verticesPerElement;
        
        totalNumOfNodes = (linearNoElements[0] + 1);
        totalNumOfVertices = (linearNoElements[0] + (periodicDIM[0] ? 0 : 1));
        
        totalNumOfElements = (linearNoElements[0]);
        
        verticesPerElement = 2;
        std::size_t powerOf2;
        
        for (std::size_t iDIM = 1; iDIM < DIM; ++iDIM)
        {
            totalNumOfNodes *= (linearNoElements[iDIM] + 1);
            totalNumOfVertices *= (linearNoElements[iDIM] + (periodicDIM[iDIM] ? 0 : 1));
            totalNumOfElements *= (linearNoElements[iDIM]);
            verticesPerElement *= 2;
            
            numOfElementsInEachSubspace[iDIM] = numOfElementsInEachSubspace[iDIM - 1] * (linearNoElements[iDIM - 1]);
            numOfNodesInEachSubspace[iDIM] = numOfNodesInEachSubspace[iDIM - 1] * (linearNoElements[iDIM - 1] + 1);
            numOfVerticesInEachSubspace[iDIM] = numOfVerticesInEachSubspace[iDIM - 1] * (linearNoElements[iDIM - 1] + (periodicDIM[iDIM - 1] ? 0 : 1));
        }
        
        //temp point for storing the node locations
        PointPhysicalT x(DIM);
        
        //Stage 2 : Create the nodes
        //Now loop over all the nodes and calculate the coordinates for reach DIMension (this makes the algorithm independent of DIMension
        for (std::size_t nodeIndex = 0; nodeIndex < totalNumOfNodes; ++nodeIndex)
        {
            std::size_t nodeIndexRemain = nodeIndex;
            
            for (int iDIM = DIM - 1; iDIM > -1; --iDIM)
            {
                x[iDIM] = bottomLeft[iDIM] + (nodeIndexRemain / numOfNodesInEachSubspace[iDIM] * delta_x[iDIM]);
                nodeIndexRemain %= numOfNodesInEachSubspace[iDIM];
            }
            
            //actually add the point
            theMesh_.addNode(x);
            
        }
        
        //stage 2.5 : create the vertices
        
        for (std::size_t vertexID = 0; vertexID < totalNumOfVertices; ++vertexID)
        {
            theMesh_.addVertex();
        }
        
        //Stage 3 : Create the elements
        
        std::vector<std::size_t> elementNdId(DIM), nodeNdId(DIM), vertexNdId(DIM), globalNodeID(verticesPerElement), globalVertexID(verticesPerElement);
        
        auto& vertexlist = getVerticesList(IteratorType::GLOBAL);
        
        //elementNdId is DIM coordinate of the bottom left node i.e. in two (0,0), (1,0) ,(2,0) ... etc are the first three (if at least three elements in x)
        for (std::size_t elementIndex = 0; elementIndex < totalNumOfElements; ++elementIndex)
        {
            std::size_t numElementsRemaining = elementIndex;
            
            for (int iDIM = DIM - 1; iDIM > -1; --iDIM)
            {
                elementNdId[iDIM] = numElementsRemaining / numOfElementsInEachSubspace[iDIM];
                numElementsRemaining %= numOfElementsInEachSubspace[iDIM];
            }
            
            // vertexNdId are the DIM coordinate of each vertex in the element with vertexNdId[0] being the bottom left
            for (std::size_t i = 0; i < verticesPerElement; ++i)
            {
                powerOf2 = 1;
                for (std::size_t iDIM = 0; iDIM < DIM; ++iDIM)
                {
                    nodeNdId[iDIM] = elementNdId[iDIM] + ((i & powerOf2) != 0);
                    vertexNdId[iDIM] = elementNdId[iDIM] + ((i & powerOf2) != 0);
                    if ((vertexNdId[iDIM] >= linearNoElements[iDIM]) && (periodicDIM[iDIM] == true))
                    {
                        vertexNdId[iDIM] = 0;
                    }
                    powerOf2 *= 2;
                }
                globalNodeID[i] = nodeNdId[0];
                globalVertexID[i] = vertexNdId[0];
                
                //Now map to the one DIMensional global ID
                for (std::size_t iDIM = 1; iDIM < DIM; ++iDIM)
                {
                    globalNodeID[i] += nodeNdId[iDIM] * numOfNodesInEachSubspace[iDIM];
                    globalVertexID[i] += vertexNdId[iDIM] * numOfVerticesInEachSubspace[iDIM];
                }
            }
            Element* newElement = addElement(globalNodeID);
            for (std::size_t i = 0; i < globalVertexID.size(); ++i)
            {
                vertexlist[globalVertexID[i]]->addElement(newElement, i);
            }
        }
        
        faceFactory();
        edgeFactory();
        
    }
    
    //createTrianglularMesh follows the same structure as createRectangularMesh. 
    //Where createRectangularMesh makes rectangular elements, createTrianglularMesh
    //splits the elements into a partition of triangles.    
    void MeshManipulator::createTriangularMesh(PointPhysicalT bottomLeft, PointPhysicalT topRight, const VectorOfPointIndicesT& linearNoElements)
    {
        logger.assert(bottomLeft.size()==topRight.size(), "The corners of the mesh must have the same dimension");
        logger.assert(bottomLeft.size()==configData_->dimension_, "The corners of the mesh have the wrong dimension");
        logger.assert(linearNoElements.size()==configData_->dimension_, "There are amounts of elements spicified in % dimensions, but there are % dimensions", linearNoElements.size(), configData_->dimension_);
        //set to correct value in case some other meshmanipulator changed things
        ElementFactory::instance().setCollectionOfBasisFunctionSets(&collBasisFSet_);
        ElementFactory::instance().setNumberOfMatrices(numberOfElementMatrixes_);
        ElementFactory::instance().setNumberOfVectors(numberOfFaceVectors_);
        ElementFactory::instance().setNumberOfTimeLevels(configData_->numberOfTimeLevels_);
        ElementFactory::instance().setNumberOfUnknowns(configData_->numberOfUnknowns_);
        FaceFactory::instance().setNumberOfFaceMatrices(numberOfFaceMatrixes_);
        FaceFactory::instance().setNumberOfFaceVectors(numberOfFaceVectors_);
        
        //Stage 0 : Check for required requirements
        std::size_t DIM = configData_->dimension_;
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
        
        PointPhysicalT delta_x(DIM);
        
        for (std::size_t i = 0; i < DIM; ++i)
        {
            delta_x[i] = (topRight[i] - bottomLeft[i]) / (linearNoElements[i]);
        }
        
        std::vector<std::size_t> numOfNodesInEachSubspace(DIM), numOfVerticesInEachSubspace(DIM), numOfElementsInEachSubspace(DIM);
        
        numOfNodesInEachSubspace[0] = 1;
        numOfVerticesInEachSubspace[0] = 1;
        numOfElementsInEachSubspace[0] = 1;
        
        std::size_t totalNumOfNodes, totalNumOfElements, verticesPerElement, verticesPerGroup, trianglesPerRectangle, totalNumOfVertices;
        
        totalNumOfNodes = (linearNoElements[0] + 1);
        totalNumOfVertices = (linearNoElements[0] + (periodicDIM[0] ? 0 : 1));
        //'elements' in this counter denote groups of trianglesPerRectangle elements
        totalNumOfElements = (linearNoElements[0]);
        verticesPerElement = 2;
        verticesPerGroup = 2;
        trianglesPerRectangle = 1;
        std::size_t powerOf2;
        
        //start with 1 because you want to ask for the entry at idim - 1
        for (std::size_t idim = 1; idim < DIM; ++idim)
        {
            totalNumOfNodes *= (linearNoElements[idim] + 1);
            totalNumOfVertices *= (linearNoElements[idim] + (periodicDIM[idim] ? 0 : 1));
            totalNumOfElements *= (linearNoElements[idim]);
            verticesPerElement += 1;
            verticesPerGroup *= 2;
            trianglesPerRectangle += (2 * idim - 1);
            numOfElementsInEachSubspace[idim] = numOfElementsInEachSubspace[idim - 1] * (linearNoElements[idim - 1]);
            numOfNodesInEachSubspace[idim] = numOfNodesInEachSubspace[idim - 1] * (linearNoElements[idim - 1] + 1);
            numOfVerticesInEachSubspace[idim] = numOfVerticesInEachSubspace[idim - 1] * (linearNoElements[idim - 1] + (periodicDIM[idim - 1] ? 0 : 1));
        }
        
        PointPhysicalT x(DIM);
        
        //Stage 2 : Create the nodes
        
        for (std::size_t nodeIndex = 0; nodeIndex < totalNumOfNodes; ++nodeIndex)
        {
            std::size_t nodeIndexRemain = nodeIndex;
            for (int idim = DIM - 1; idim > -1; --idim)
            {
                x[idim] = bottomLeft[idim] + (nodeIndexRemain / numOfNodesInEachSubspace[idim] * delta_x[idim]);
                nodeIndexRemain %= numOfNodesInEachSubspace[idim];
            }
            theMesh_.addNode(x);
        }
        
        for (std::size_t nodeIndex = 0; nodeIndex < totalNumOfVertices; ++nodeIndex)
        {
            theMesh_.addVertex();
        }
        
        auto& vertices = getVerticesList(IteratorType::GLOBAL);
        
        //Stage 3 : Create the elements
        
        std::vector<std::size_t> elementNdId(DIM), vertexNdId(DIM), nodeNdId(DIM);
        std::vector<std::vector<std::size_t> > globalVertexID(trianglesPerRectangle);
        std::vector<std::vector<std::size_t> > globalNodeID(trianglesPerRectangle);
        
        for (std::size_t elementGroupIndex = 0; elementGroupIndex < totalNumOfElements; ++elementGroupIndex)
        {
            //first generate node indexes as if we are a cube
            //indicates if the element has to be rotated to make connecting faces; rotates the element 90 degrees along the y-axis if needed
            std::size_t rotate = 0;
            
            for (std::size_t i = 0; i < trianglesPerRectangle; ++i)
            {
                globalNodeID[i].clear();
                globalVertexID[i].clear();
            }
            int elementIndexRemainder = elementGroupIndex;
            
            for (int idim = DIM - 1; idim > -1; --idim)
            {
                elementNdId[idim] = elementIndexRemainder / numOfElementsInEachSubspace[idim];
                elementIndexRemainder %= numOfElementsInEachSubspace[idim];
                rotate = (elementNdId[idim] + rotate) % 2;
            }
            
            for (std::size_t i = 0; i < verticesPerGroup; ++i)
            {
                if (rotate == 0)
                {
                    powerOf2 = 1;
                    for (std::size_t idim = 0; idim < DIM; ++idim)
                    {
                        nodeNdId[idim] = elementNdId[idim] + ((i & powerOf2) != 0);
                        vertexNdId[idim] = elementNdId[idim] + ((i & powerOf2) != 0);
                        if (vertexNdId[idim] >= linearNoElements[idim] && periodicDIM[idim])
                            vertexNdId[idim] = 0;
                        powerOf2 *= 2;
                    }
                }
                else
                {
                    powerOf2 = verticesPerGroup;
                    for (std::size_t idim = 0; idim < DIM; ++idim)
                    {
                        powerOf2 /= 2;
                        nodeNdId[idim] = elementNdId[idim] + (((i ^ rotate) & powerOf2) != 0);
                        vertexNdId[idim] = elementNdId[idim] + (((i ^ rotate) & powerOf2) != 0);
                        if (vertexNdId[idim] >= linearNoElements[idim] && periodicDIM[idim])
                            vertexNdId[idim] = 0;
                    }
                }
                
                std::size_t nodeIndex = nodeNdId[0];
                std::size_t vertexIndex = vertexNdId[0];
                for (std::size_t idim = 1; idim < DIM; ++idim)
                {
                    nodeIndex += nodeNdId[idim] * numOfNodesInEachSubspace[idim];
                    vertexIndex += vertexNdId[idim] * numOfVerticesInEachSubspace[idim];
                }
                
                //then cherrypick the element(s) these vertices should connect to (probably not the cleanest implementation; \bug doesn't work if DIM>3)
                switch (i)
                {
                    case 0:
                        globalVertexID[0].push_back(vertexIndex);
                        globalNodeID[0].push_back(nodeIndex);
                        break;
                    case 3:
                        //insert in the second place because the ordering of the vertices will work out better
                        globalVertexID[1].insert( ++(globalVertexID[1].begin()), vertexIndex);
                        globalNodeID[1].insert( ++(globalNodeID[1].begin()), nodeIndex);
                        break;
                    case 5:
                        globalVertexID[2].push_back(vertexIndex);
                        globalNodeID[2].push_back(nodeIndex);
                        break;
                    case 6:
                        //insert in the second place because the ordering of the vertices will work out better
                        globalVertexID[3].insert( ++(globalVertexID[3].begin()), vertexIndex);
                        globalNodeID[3].insert( ++(globalNodeID[3].begin()), nodeIndex);
                        break;
                    case 1:
                        for (std::size_t i = 0; i < trianglesPerRectangle; ++i)
                        {
                            if (i != 3)
                            {
                                globalVertexID[i].push_back(vertexIndex);
                                globalNodeID[i].push_back(nodeIndex);
                            }
                        }
                        break;
                    case 2:
                        for (std::size_t i = 0; i < trianglesPerRectangle; ++i)
                        {
                            if (i != 2)
                            {
                                globalVertexID[i].push_back(vertexIndex);
                                globalNodeID[i].push_back(nodeIndex);
                            }
                        }
                        break;
                    case 4:
                        for (std::size_t i = 0; i < trianglesPerRectangle; ++i)
                        {
                            if (i != 1)
                            {
                                globalVertexID[i].push_back(vertexIndex);
                                globalNodeID[i].push_back(nodeIndex);
                            }
                        }
                        break;
                    case 7:
                        for (std::size_t i = 0; i < trianglesPerRectangle; ++i)
                        {
                            if (i != 0)
                            {
                                globalVertexID[i].push_back(vertexIndex);
                                globalNodeID[i].push_back(nodeIndex);
                            }
                        }
                        break;
                } //switch
            } //for all vertices of the rectangle
            
            for (std::size_t i = 0; i < trianglesPerRectangle; ++i)
            {
                Element* newElement = addElement(globalNodeID[i]);
                for (std::size_t j = 0; j < globalVertexID[i].size(); ++j)
                {
                    logger.assert(i < globalVertexID.size(), "Requested vertex %, while there are only %", i, globalVertexID.size());
                    logger.assert(j < globalVertexID[i].size(), "Requested element %, but this vertex only has %.", j, globalVertexID[i].size());
                    logger.assert(globalVertexID[i][j] < totalNumOfVertices, "Requested vertex %, while there are only %", globalVertexID[i][j], totalNumOfVertices);
                    logger.assert(vertices.size() == totalNumOfVertices, "Number of vertices is wrong.");
                    vertices[globalVertexID[i][j]]->addElement(newElement, j);
                }
            }
        } //for all rectangles
        
        //Stage 4 : Create the faces
        faceFactory();
        edgeFactory();
    }
    
    void MeshManipulator::readCentaurMesh(const std::string& filename)
    {
        //set to correct value in case some other meshManipulator changed things
        ElementFactory::instance().setCollectionOfBasisFunctionSets(&collBasisFSet_);
        ElementFactory::instance().setNumberOfMatrices(numberOfElementMatrixes_);
        ElementFactory::instance().setNumberOfVectors(numberOfFaceVectors_);
        ElementFactory::instance().setNumberOfTimeLevels(configData_->numberOfTimeLevels_);
        ElementFactory::instance().setNumberOfUnknowns(configData_->numberOfUnknowns_);
        FaceFactory::instance().setNumberOfFaceMatrices(numberOfFaceMatrixes_);
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
    
    void MeshManipulator::readCentaurMesh2D(std::ifstream& centaurFile)
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
            char junk[1024];
            
            std::uint32_t checkInt;
            centaurFile.read(&junk[0], sizeOfLine - sizeof(version) - sizeof(centaurFileType));
            
            //Check the first line was read correctly : each line in centaur start and end with the line size as a check
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
            
            //Start the second line
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            
            //Next read the total number of nodes
            std::uint32_t numberOfNodes;
            centaurFile.read(reinterpret_cast<char*>(&numberOfNodes), sizeof(numberOfNodes));
            logger(INFO, "File contains % nodes", numberOfNodes);
            
            //Check the second line was read correctly : each line in centaur start and end with the line size as a check
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
            
            //placeholder for vertices until we know where the periodic boundaries are
            std::vector<std::vector<std::size_t> > listOfElementsForEachNode(numberOfNodes);
            
            //Now we will read in all the nodes
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            double nodeCoord[2];
            PointPhysicalT nodeCoordPointFormat(2);
            for (std::size_t i = 0; i < numberOfNodes; i++)
            {
                // Reads the x and y coordinates of each node.
                centaurFile.read(reinterpret_cast<char*>(nodeCoord), sizeof(nodeCoord));
                // pass the node to the nodelist.
                
                //Covert from *double to hpGEM PointPhysical format
                nodeCoordPointFormat[0] = nodeCoord[0];
                nodeCoordPointFormat[1] = nodeCoord[1];
                theMesh_.addNode(nodeCoordPointFormat);
                
            }
            //Now check the node line was read correctly : each line in centaur start and end with the line size as a check
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
            
            //Now check how many triangle in the file
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            // Number of triangular elements
            std::uint32_t numberOfTriangles;
            centaurFile.read(reinterpret_cast<char*>(&numberOfTriangles), sizeof(numberOfTriangles));
            logger(INFO, "File contains % triangle(s)", numberOfTriangles);
            
            //Check the line was read correctly : each line in centaur start and end with the line size as a check
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
            
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
            logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
            
            //Now check the number of quaduratiles in the file
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            std::uint32_t numberOfQuads;
            centaurFile.read(reinterpret_cast<char*>(&numberOfQuads), sizeof(numberOfQuads));
            logger(INFO, "File contains % quaduratile(s)", numberOfQuads);
            
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
            
            //Now read the quaduritles in
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            if (numberOfQuads > 0)
            {
                std::uint32_t temp;
                std::vector<std::uint32_t> globalNodeIndexes(4);
                std::vector<std::size_t> globalNodeIndexesSizeT(4);
                for (std::size_t i = 0; i < numberOfQuads; i++)
                {
                    //Reading the vertex indices of each quadrilateral.
                    centaurFile.read(reinterpret_cast<char*>(globalNodeIndexes.data()), sizeof(std::uint32_t) * globalNodeIndexes.size());
                    
                    // renumbering of the vertices to match the ordering assumed by
                    // hpGem:
                    
                    temp = globalNodeIndexes[2];
                    globalNodeIndexes[2] = globalNodeIndexes[3];
                    globalNodeIndexes[3] = temp;
                    
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
            logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
            
            //now read boundary data
            //nodes first
            
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            
            std::uint32_t numberOfBoundaryNodes;
            centaurFile.read(reinterpret_cast<char*>(&numberOfBoundaryNodes), sizeof(numberOfBoundaryNodes));
            
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
            
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
            logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
            
            //then faces
            
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            
            std::uint32_t numberOfBoundaryFaces;
            centaurFile.read(reinterpret_cast<char*>(&numberOfBoundaryFaces), sizeof(numberOfBoundaryFaces));
            
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
            
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            
            //read the half-faces and store them until we get to the boundary information
            std::vector<HalfFaceDescription> boundaryFaces(numberOfBoundaryFaces);
            
            for (uint_fast32_t i = 0; i < numberOfBoundaryFaces; ++i)
            {
                boundaryFaces[i].nodeList.resize(2);
                centaurFile.read(reinterpret_cast<char*>(boundaryFaces[i].nodeList.data()), sizeof(std::uint32_t) * 2);
                
                boundaryFaces[i].nodeList[0] -= 1;
                boundaryFaces[i].nodeList[1] -= 1;
                
                std::vector<std::size_t> candidateElements;
                std::vector<std::size_t>& leftNodeElements = listOfElementsForEachNode[boundaryFaces[i].nodeList[0]];
                std::vector<std::size_t>& rightNodeElements = listOfElementsForEachNode[boundaryFaces[i].nodeList[1]];
                
                std::set_intersection(leftNodeElements.begin(), leftNodeElements.end(), rightNodeElements.begin(), rightNodeElements.end(), std::back_inserter(candidateElements));
                
                //boundary faces *should* border only one element
                logger.assert_always(candidateElements.size() < 2, "candidate boundary face lies at two or more elements"); 
                boundaryFaces[i].elementNum = candidateElements[0];
                
                Element* current = elementslist[candidateElements[0]];
                std::vector<std::size_t> faceVertices(2);
                
                for (std::size_t j = 0; j < current->getNrOfFaces(); ++j)
                {
                    faceVertices = current->getPhysicalGeometry()->getGlobalFaceNodeIndices(j);
                    if ((faceVertices[0] == boundaryFaces[i].nodeList[0] || faceVertices[0] == boundaryFaces[i].nodeList[1]) && (faceVertices[1] == boundaryFaces[i].nodeList[0] || faceVertices[1] == boundaryFaces[i].nodeList[1]))
                    {
                        boundaryFaces[i].localFaceIndex = j;
                    }
                }
                
                candidateElements.clear();
            }
            
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
            
            //now read a list of boundary segments (that will link the faces to boundary conditions later on)
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            
            std::vector<std::uint32_t> faceToSegment(numberOfBoundaryFaces);
            centaurFile.read(reinterpret_cast<char*>(faceToSegment.data()), sizeof(std::uint32_t) * numberOfBoundaryFaces);
            
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
            
            //now couple the segments to boundary groups
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            
            std::uint32_t numberOfSegments;
            centaurFile.read(reinterpret_cast<char*>(&numberOfSegments), sizeof(numberOfSegments));
            
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
            
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            
            std::vector<std::uint32_t> segmentToGroup(numberOfSegments);
            centaurFile.read(reinterpret_cast<char*>(segmentToGroup.data()), sizeof(std::uint32_t) * numberOfSegments);
            
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
            
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
            logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
            
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            
            std::vector<std::uint32_t> groupBCType(numberOfGroups);
            centaurFile.read(reinterpret_cast<char*>(groupBCType.data()), sizeof(std::uint32_t) * numberOfGroups);
            
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
            
            for (uint_fast32_t i = 0; i < numberOfBoundaryFaces; ++i)
            {
                std::uint32_t boundaryType = groupBCType[segmentToGroup[faceToSegment[i]]];
                if (boundaryType < 1001)
                {
                    logger(INFO, "Viscous Wall boundary for face % assigned as WALL_BC", i);
                    addFace(elementslist[boundaryFaces[i].elementNum], boundaryFaces[i].localFaceIndex, nullptr, 0, Geometry::FaceType::WALL_BC);
                }
                else if (boundaryType < 2001)
                {
                    logger(INFO, "Inviscid Wall boundary for face % assigned as WALL_BC", i);
                    addFace(elementslist[boundaryFaces[i].elementNum], boundaryFaces[i].localFaceIndex, nullptr, 0, Geometry::FaceType::WALL_BC);
                }
                else if (boundaryType < 3001)
                {
                    logger(INFO,  "symmetry plane boundary for face % assigned as WALL_BC", i);
                    addFace(elementslist[boundaryFaces[i].elementNum], boundaryFaces[i].localFaceIndex, nullptr, 0, Geometry::FaceType::WALL_BC);
                }
                else if (boundaryType < 4001)
                {
                    logger(INFO, "inlet pipe boundary for face % assigned as OPEN_BC", i);
                    addFace(elementslist[boundaryFaces[i].elementNum], boundaryFaces[i].localFaceIndex, nullptr, 0, Geometry::FaceType::OPEN_BC);
                }
                else if (boundaryType < 5001)
                {
                    logger(INFO, "outlet pipe boundary for face % assigned as OPEN_BC", i);
                    addFace(elementslist[boundaryFaces[i].elementNum], boundaryFaces[i].localFaceIndex, nullptr, 0, Geometry::FaceType::OPEN_BC);
                }
                else if (boundaryType < 6001)
                {
                    logger(INFO, "farfield boundary for face % assigned as OPEN_BC", i);
                    addFace(elementslist[boundaryFaces[i].elementNum], boundaryFaces[i].localFaceIndex, nullptr, 0, Geometry::FaceType::OPEN_BC);
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
                    addFace(elementslist[boundaryFaces[i].elementNum], boundaryFaces[i].localFaceIndex, nullptr, 0, Geometry::FaceType::OPEN_BC);
                }
                else if (boundaryType < 10001)
                {
                    logger(INFO, "moving wall boundary for face % assigned as WALL_BC", i);
                    addFace(elementslist[boundaryFaces[i].elementNum], boundaryFaces[i].localFaceIndex, nullptr, 0, Geometry::FaceType::WALL_BC);
                }
                else
                {
                    logger(INFO, "alternative boundary condition for face % assigned as WALL_BC", i);
                    addFace(elementslist[boundaryFaces[i].elementNum], boundaryFaces[i].localFaceIndex, nullptr, 0, Geometry::FaceType::WALL_BC);
                }
            }
            
            //I don't care about the names of the groups, just skip them
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            
            centaurFile.ignore(80 * numberOfGroups);
            
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
            
            //now read periodic information and link corresponding vertices
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            
            std::uint32_t numberOfTransforms;
            centaurFile.read(reinterpret_cast<char*>(&numberOfTransforms), sizeof(numberOfTransforms));
            
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
            
            for (uint_fast32_t i = 0; i < numberOfTransforms; ++i)
            {
                //I dont care about the actual coordinate transformations, just skip them
                centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
                
                centaurFile.ignore(9 * sizeof(std::uint32_t));
                
                centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
                logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
                centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
                
                centaurFile.ignore(9 * sizeof(uint));
                
                centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
                logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
                
                std::uint32_t numberOfNodePairs;
                centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
                
                centaurFile.read(reinterpret_cast<char*>(&numberOfNodePairs), sizeof(numberOfNodePairs));
                
                centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
                logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
                
                centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
                
                std::uint32_t nodePair[2];
                std::vector<std::size_t> combine;
                for (uint_fast32_t j = 0; j < numberOfNodePairs; ++j)
                {
                    centaurFile.read(reinterpret_cast<char*>(nodePair), sizeof(nodePair));
                    auto& firstList = listOfElementsForEachNode[nodePair[0]];
                    auto& secondList = listOfElementsForEachNode[nodePair[1]];
                    std::set_union(firstList.begin(), firstList.end(), secondList.begin(), secondList.end(), std::back_inserter(combine));
                    firstList = combine;
                    secondList = combine;
                    combine.clear();
                }
                
                centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
                logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
                
            }
            
            //new centaur files have zones, but I'm not interested in them
            
            //now we know periodicity information, construct the vertices
            //remember that listOfElementsForEachNode will in general contain duplicates
            
            auto& listOfVertices = theMesh_.getVerticesList(IteratorType::GLOBAL);
            
            bool addedNewVertex(false);
            
            for (std::size_t i = 0; i < listOfElementsForEachNode.size(); ++i)
            {
                for (std::size_t j = 0; j < listOfElementsForEachNode[i].size(); ++j)
                {
                    Element* current = elementslist[listOfElementsForEachNode[i][j]];
                    for (std::size_t k = 0; k < current->getNrOfNodes(); ++k)
                    {
                        //if we did not jet deal with this node and it is the correct one
                        if (current->getNode(k) == nullptr && current->getPhysicalGeometry()->getNodeIndex(k) == i)
                        {
                            if (!addedNewVertex)
                            {
                                addVertex();
                                addedNewVertex = true;
                            }
                            listOfVertices.back()->addElement(current, k);
                        }
                    }
                }
                addedNewVertex = false;
            }
            
            //construct the rest of the faces
            faceFactory();
            
    }
    
    void MeshManipulator::readCentaurMesh3D(std::ifstream& centaurFile)
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
            logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
            
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
            logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
            
            //Now we will read in all the nodes
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            double nodeCoord[3];
            PointPhysicalT nodeCoordPointFormat(3);
            for (std::size_t i = 0; i < numberOfNodes; i++)
            {
                if (i > 0 && i % numberOfNodesPerLine == 0)
                {
                    //If all the nodes on a line are read end the line and start a new one
                    centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
                    logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
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
                theMesh_.addNode(nodeCoordPointFormat);
                
            }
            //Now check the node line was read correctly : each line in centaur start and end with the line size as a check
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
            
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
                logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
                
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
                        logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
                        centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
                    }
                    
                    //Reading the vertex indices of each hexahedron.
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
                logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
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
            logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
            
            std::vector<std::uint32_t> globalNodeIndexes(6);
            std::vector<std::size_t> globalNodeIndexesSizeT(6);
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            for (std::size_t i = 0; i < numberOfPrisms; i++)
            {
                if (i > 0 && i % numberOfPrismsPerLine == 0)
                {
                    //If all the nodes on a line are read end the line and start a new one
                    centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
                    logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
                    centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
                }
                
                //Reading the vertex indices of each hexahedron.
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
            logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
            
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
                logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
                
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
                        logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
                        centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
                    }
                    
                    //Reading the vertex indices of each pyramid.
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
                logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
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
            logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
            
            globalNodeIndexes.resize(4);
            globalNodeIndexesSizeT.resize(4);
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            for (std::size_t i = 0; i < numberOfTetrahedra; i++)
            {
                if (i > 0 && i % numberOfTetrahedraPerLine == 0)
                {
                    //If all the tetrahedra on a line are read end the line and start a new one
                    centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
                    logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
                    centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
                }
                
                //Reading the vertex indices of each tetrahedron.
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
            logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
            
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
            logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
            
            //first read the information about the faces
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            for (std::size_t i = 0; i < numberOfBoundaryFaces; ++i)
            {
                
                if (i > 0 && i % numberOfBoundaryFaces == 0)
                {
                    //If all the tetrahedra on a line are read end the line and start a new one
                    centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
                    logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
                    centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
                }
                
                centaurFile.read(reinterpret_cast<char*>(&nodalDescriptionOfTheFace[0]), sizeof(nodalDescriptionOfTheFace));
                
                boundarFaces[i].nodeList.resize(3);
                boundarFaces[i].nodeList[0] = nodalDescriptionOfTheFace[0];
                boundarFaces[i].nodeList[1] = nodalDescriptionOfTheFace[1];
                boundarFaces[i].nodeList[2] = nodalDescriptionOfTheFace[2];
                
                //three nodes will uniquely determine the face
                auto& firstNodeList = listOfElementsForEachNode[nodalDescriptionOfTheFace[0]];
                auto& secondNodeList = listOfElementsForEachNode[nodalDescriptionOfTheFace[1]];
                auto& thirdNodeList = listOfElementsForEachNode[nodalDescriptionOfTheFace[2]];
                std::vector<std::size_t> temp, candidates, nodes(boundarFaces->nodeList), intersect;
                std::sort(nodes.begin(), nodes.end());
                std::set_intersection(firstNodeList.begin(), firstNodeList.end(), secondNodeList.begin(), secondNodeList.end(), std::back_inserter(temp));
                std::set_intersection(temp.begin(), temp.end(), thirdNodeList.begin(), thirdNodeList.end(), std::back_inserter(candidates));
                
                logger.assert_always(candidates.size() < 2, "candidate boundary face lies at two or more elements");
                
                boundarFaces[i].elementNum = candidates[0];
                
                Element* current = elementsList[candidates[0]];
                
                for (std::size_t j = 0; j < current->getNrOfFaces(); ++j)
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
                                boundarFaces[i].nodeList.push_back(boundarFaces[i].nodeList[2]);
                                boundarFaces[i].nodeList[2] = nodalDescriptionOfTheFace[3];
                            }
                            else
                            {
                                boundarFaces[i].nodeList.push_back(nodalDescriptionOfTheFace[3]);
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
            logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
            
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
                        logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
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
                logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
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
            logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
            
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
            logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
            
            std::uint32_t centaurBCType;
            char nameOfBoundaryCondition[80];
            
            //this bit of information is again a little late
            std::uint32_t numberOfBoundaryGroups;
            
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            centaurFile.read(reinterpret_cast<char*>(&numberOfBoundaryGroups), sizeof(numberOfBoundaryGroups));
            logger.assert(numberOfBoundaryGroups == facesForEachBoundaryGroup.size(), "Not enough boundary groups in centaur file");
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
            
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
                        addFace(tempElementVector[boundarFaces[facesForEachBoundaryGroup[i][j]].elementNum], boundarFaces[facesForEachBoundaryGroup[i][j]].localFaceIndex, nullptr, 0, Geometry::FaceType::WALL_BC);
                    }
                }
                else if (centaurBCType < 2001)
                {
                    logger(INFO, "Inviscid Wall boundary for group % assigned as WALL_BC", i);
                    for (std::size_t j = 0; j < facesForEachBoundaryGroup[i].size(); ++j)
                    {
                        //big assumption on the nature of elementIDs here...
                        addFace(tempElementVector[boundarFaces[facesForEachBoundaryGroup[i][j]].elementNum], boundarFaces[facesForEachBoundaryGroup[i][j]].localFaceIndex, nullptr, 0, Geometry::FaceType::WALL_BC);
                    }
                }
                else if (centaurBCType < 3001)
                {
                    logger(INFO, "symmetry plane boundary for group % assigned as WALL_BC", i);
                    for (std::size_t j = 0; j < facesForEachBoundaryGroup[i].size(); ++j)
                    {
                        //big assumption on the nature of elementIDs here...
                        addFace(tempElementVector[boundarFaces[facesForEachBoundaryGroup[i][j]].elementNum], boundarFaces[facesForEachBoundaryGroup[i][j]].localFaceIndex, nullptr, 0, Geometry::FaceType::WALL_BC);
                    }
                }
                else if (centaurBCType < 4001)
                {
                    logger(INFO, "inlet pipe boundary for group % assigned as OPEN_BC", i);
                    for (std::size_t j = 0; j < facesForEachBoundaryGroup[i].size(); ++j)
                    {
                        //big assumption on the nature of elementIDs here...
                        addFace(tempElementVector[boundarFaces[facesForEachBoundaryGroup[i][j]].elementNum], boundarFaces[facesForEachBoundaryGroup[i][j]].localFaceIndex, nullptr, 0, Geometry::FaceType::OPEN_BC);
                    }
                }
                else if (centaurBCType < 5001)
                {
                    logger(INFO, "outlet pipe boundary for group % assigned as OPEN_BC", i);
                    for (std::size_t j = 0; j < facesForEachBoundaryGroup[i].size(); ++j)
                    {
                        //big assumption on the nature of elementIDs here...
                        addFace(tempElementVector[boundarFaces[facesForEachBoundaryGroup[i][j]].elementNum], boundarFaces[facesForEachBoundaryGroup[i][j]].localFaceIndex, nullptr, 0, Geometry::FaceType::OPEN_BC);
                    }
                }
                else if (centaurBCType < 6001)
                {
                    logger(INFO, "farfield boundary for group % assigned as OPEN_BC", i);
                    for (std::size_t j = 0; j < facesForEachBoundaryGroup[i].size(); ++j)
                    {
                        //big assumption on the nature of elementIDs here...
                        addFace(tempElementVector[boundarFaces[facesForEachBoundaryGroup[i][j]].elementNum], boundarFaces[facesForEachBoundaryGroup[i][j]].localFaceIndex, nullptr, 0, Geometry::FaceType::OPEN_BC);
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
                        addFace(tempElementVector[boundarFaces[facesForEachBoundaryGroup[i][j]].elementNum], boundarFaces[facesForEachBoundaryGroup[i][j]].localFaceIndex, nullptr, 0, Geometry::FaceType::OPEN_BC);
                    }
                }
                else if (centaurBCType < 10001)
                {
                    logger(INFO, "moving wall boundary for group % assigned as WALL_BC", i);
                    for (std::size_t j = 0; j < facesForEachBoundaryGroup[i].size(); ++j)
                    {
                        addFace(tempElementVector[boundarFaces[facesForEachBoundaryGroup[i][j]].elementNum], boundarFaces[facesForEachBoundaryGroup[i][j]].localFaceIndex, nullptr, 0, Geometry::FaceType::WALL_BC);
                    }
                }
                else
                {
                    logger(INFO, "alternative boundary condition for group % assigned as WALL_BC", i);
                    for (std::size_t j = 0; j < facesForEachBoundaryGroup[i].size(); ++j)
                    {
                        addFace(tempElementVector[boundarFaces[facesForEachBoundaryGroup[i][j]].elementNum], boundarFaces[facesForEachBoundaryGroup[i][j]].localFaceIndex, nullptr, 0, Geometry::FaceType::WALL_BC);
                    }
                }
                logger(INFO, "total number of boundary faces: %", getFacesList().size());
            }
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
            
            //This is where centaur tells the names of all the boundary groups
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            for (std::size_t i = 0; i < numberOfBoundaryGroups; ++i)
            {
                centaurFile.read(reinterpret_cast<char*>(&nameOfBoundaryCondition[0]), sizeof(nameOfBoundaryCondition));
                logger(INFO, "boundary condition % is called %", i, nameOfBoundaryCondition);
            }
            centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
            logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
            
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
                logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
                for (std::size_t i = 0; i < numberOfPeriodicTransformations; ++i)
                {
                    //information on how to do the transformation can be computed later so just throw it away now
                    centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
                    centaurFile.read(reinterpret_cast<char*>(&ijunk), sizeof(ijunk));
                    centaurFile.read(reinterpret_cast<char*>(&transformationData[0]), sizeof(transformationData));
                    centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
                    logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
                    centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
                    centaurFile.read(reinterpret_cast<char*>(&ijunk), sizeof(ijunk));
                    centaurFile.read(reinterpret_cast<char*>(&transformationData[0]), sizeof(transformationData));
                    centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
                    logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
                    
                    //now read the amount of periodic nodes
                    centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
                    centaurFile.read(reinterpret_cast<char*>(&numberOfPeriodicNodes), sizeof(numberOfPeriodicNodes));
                    logger(INFO, "transformation group % contains % node->node matching(s)", i, numberOfPeriodicNodes);
                    centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
                    logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
                    
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
                    logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
                }
            }
            else
            {
                //now read the amount of periodic nodes
                centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
                centaurFile.read(reinterpret_cast<char*>(&numberOfPeriodicNodes), sizeof(numberOfPeriodicNodes));
                logger(INFO, "the transformation group contains % node -> node matching(s)", numberOfPeriodicNodes);
                centaurFile.read(reinterpret_cast<char*>(&checkInt), sizeof(checkInt));
                logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
                
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
                logger.assert_always(checkInt == sizeOfLine, "Error in centaur file.");
            }
            
            logger(INFO, "begin constructing internal faces and internal \"boundaries\"");
            
            //now we know periodicity information, construct the vertices
            //remember that listOfElementsForEachNode will in general contain duplicates
            
            auto& listOfVertices = theMesh_.getVerticesList(IteratorType::GLOBAL);
            
            bool addedNewVertex(false);
            
            for (std::size_t i = 0; i < listOfElementsForEachNode.size(); ++i)
            {
                for (std::size_t j = 0; j < listOfElementsForEachNode[i].size(); ++j)
                {
                    Element* current = elementsList[listOfElementsForEachNode[i][j]];
                    for (std::size_t k = 0; k < current->getNrOfNodes(); ++k)
                    {
                        //if we did not jet deal with this node and it is the correct one
                        if (current->getNode(k) == nullptr && current->getPhysicalGeometry()->getNodeIndex(k) == i)
                        {
                            if (!addedNewVertex)
                            {
                                addVertex();
                                addedNewVertex = true;
                            }
                            listOfVertices.back()->addElement(current, k);
                        }
                    }
                }
                addedNewVertex = false;
            }
            
            faceFactory();
            edgeFactory();
            
            delete[] boundarFaces;
        
    }
    
#ifdef HPGEM_USE_QHULL
    void MeshManipulator::createUnstructuredMesh(PointPhysicalT BottomLeft, PointPhysicalT TopRight, std::size_t TotalNoNodes, std::function<double(PointPhysicalT)> domainDescription, std::vector<PointPhysicalT> fixedPoints, std::function<double(PointPhysicalT)> relativeEdgeLength, double growFactor)
    {
        //impossible to create a mesh with more fixed nodes that total nodes
        //note that when equality is met, this will only do a delaunay triangulation
        logger.assert(fixedPoints.size() <= TotalNoNodes, "Cannot create a mesh with more fixed nodes than total nodes");
        
        //set to correct value in case some other meshmanipulator changed things
        ElementFactory::instance().setCollectionOfBasisFunctionSets(&collBasisFSet_);
        ElementFactory::instance().setNumberOfMatrices(numberOfElementMatrixes_);
        ElementFactory::instance().setNumberOfVectors(numberOfFaceVectors_);
        ElementFactory::instance().setNumberOfTimeLevels(configData_->numberOfTimeLevels_);
        ElementFactory::instance().setNumberOfUnknowns(configData_->numberOfUnknowns_);
        FaceFactory::instance().setNumberOfFaceMatrices(numberOfFaceMatrixes_);
        FaceFactory::instance().setNumberOfFaceVectors(numberOfFaceVectors_);
        
        //periodic unstructured mesh generation not yet implemented
        logger.assert(!(periodicX_ || periodicY_ || periodicZ_), "Unstructured mesh generator does not support periodic boundaries");
        
        //guess the required distance between two nodes
        double dist = std::pow(double(TotalNoNodes), -1. / double(dimension()));
        
        for (std::size_t i = 0; i < dimension(); ++i)
        {
            dist *= std::pow(TopRight[i] - BottomLeft[i], 1. / double(dimension()));
        }
        
        std::vector<PointPhysicalT> hpGEMCoordinates = fixedPoints;
        
        //seed approximately N points inside the bounding box (total amount will be tweaked later)
        std::size_t DIM = dimension();
        PointPhysicalT nextPoint = BottomLeft;
        for (std::size_t i = 0; i < fixedPoints.size(); ++i)
        {
            logger.assert(domainDescription(fixedPoints[i]) < 1e-10, "One of the fixed points is outside of the domain");
            theMesh_.addVertex();
            theMesh_.addNode(fixedPoints[i]);
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
                theMesh_.addVertex();
                theMesh_.addNode(nextPoint);
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
            
            for (PointPhysicalT point : hpGEMCoordinates)
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
                    PointPhysicalT center(DIM);
                    std::vector<std::size_t> pointIndices;
                    for (auto vertex : triangle.vertices())
                    {
                        center += hpGEMCoordinates[vertex.point().id()];
                        pointIndices.push_back(vertex.point().id());
                    }
                    center = center / pointIndices.size();
                    if (domainDescription(center) < 0)
                    {
                        auto newElement = addElement(pointIndices);
                        for (std::size_t i = 0; i < pointIndices.size(); ++i)
                        {
                            theMesh_.getVerticesList(IteratorType::GLOBAL)[pointIndices[i]]->addElement(newElement, i);
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
                    Node* current = theMesh_.getVerticesList(IteratorType::GLOBAL)[entry.second];
                    for (Element* element : current->getElements())
                    {
                        for (std::size_t i = 0; i < element->getNrOfNodes(); ++i)
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
                    Geometry::PointPhysical firstNode = element->getPhysicalGeometry()->getLocalNodeCoordinates(0);
                    Geometry::PointPhysical secondNode = element->getPhysicalGeometry()->getLocalNodeCoordinates(1);
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
                    PointPhysicalT firstNode(DIM), secondNode(DIM);
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
                    PointPhysicalT firstNode(DIM), secondNode(DIM);
                    std::vector<std::size_t> nodeIndices = edge->getElement(0)->getReferenceGeometry()->getCodim2EntityLocalIndices(edge->getEdgeNr(0));
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
            std::multimap<double, PointPhysicalT> centerPoints;
            if (DIM == 1)
            {
                //the algorithm is mostly dimension independent, but the data type it operates on is not
                for (Element* element : theMesh_.getElementsList(IteratorType::GLOBAL))
                {
                    PointPhysicalT firstNode(DIM), secondNode(DIM);
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
                    PointPhysicalT firstNode(DIM), secondNode(DIM);
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
                    PointPhysicalT firstNode(DIM), secondNode(DIM);
                    std::vector<std::size_t> nodeIndices = edge->getElement(0)->getReferenceGeometry()->getCodim2EntityLocalIndices(edge->getEdgeNr(0));
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
            for (PointPhysicalT point : hpGEMCoordinates)
            {
                theMesh_.addVertex();
                theMesh_.addNode(point);
            }
        }
        //start of QHull magic to create a triangulation
        orgQhull::RboxPoints qHullCoordinates;
        qHullCoordinates.setDimension(DIM);
        qHullCoordinates.reserveCoordinates(DIM * hpGEMCoordinates.size());
        
        for (PointPhysicalT point : hpGEMCoordinates)
        {
            qHullCoordinates.append(2, point.data());
        }
        
        //create the triangulation, pass "d" for delaunay
        //"QJ" because there are likely to be groups of more that (d+1) cocircular nodes in a regular grid, so joggle them up a bit
        orgQhull::Qhull triangulation(qHullCoordinates, "d QbB Qx Qc Qt");
        
        for (orgQhull::QhullFacet triangle : triangulation.facetList())
        {
            if (triangle.isGood() && !triangle.isUpperDelaunay())
            {
                PointPhysicalT center {DIM};
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
                        theMesh_.getVerticesList(IteratorType::GLOBAL)[pointIndices[i]]->addElement(newElement, i);
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
    
    void MeshManipulator::updateMesh(std::function<double(PointPhysicalT)> domainDescription, std::vector<std::size_t> fixedPointIdxs, std::function<double(PointPhysicalT)> relativeEdgeLength, double growFactor)
    {
        std::sort(fixedPointIdxs.begin(), fixedPointIdxs.end());
        std::size_t DIM = dimension();
        bool needsExpansion = false;
        double totalCurrentLength = 0;
        double oldQuality = 0;
        double worstQuality = 0.5;
        
        std::set<std::pair<std::size_t, std::size_t> > periodicPairing {};
        
        for (Node* node : theMesh_.getVerticesList(IteratorType::GLOBAL))
        {
            PointPhysicalT point {DIM};
            PointPhysicalT compare {DIM};
            point = node->getElement(0)->getPhysicalGeometry()->getLocalNodeCoordinates(node->getVertexNr(0));
            std::set<std::size_t> equivalentIndices {};
            equivalentIndices.insert(node->getElement(0)->getPhysicalGeometry()->getNodeIndex(node->getVertexNr(0)));
            for (std::size_t i = 1; i < node->getNrOfElements(); ++i)
            {
                compare = node->getElement(i)->getPhysicalGeometry()->getLocalNodeCoordinates(node->getVertexNr(i));
                if (compare != point)
                {
                    equivalentIndices.insert(node->getElement(i)->getPhysicalGeometry()->getNodeIndex(node->getVertexNr(i)));
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
        //except dont bother if a retriangulation is in order anyway
        if (oldNodeLocations_.size() == theMesh_.getNodes().size())
        {
            std::vector<double> unscaledShift {};
            unscaledShift.reserve(theMesh_.getNumberOfNodes());
            //compute current and expected (relative) edge length
            std::vector<double> expectedLength {};
            expectedLength.reserve(theMesh_.getNumberOfNodes());
            std::multimap<double, std::size_t> knownLengths;
            std::vector<double> currentLength {};
            //for proper scaling
            double totalexpectedLength = 0;
            for (Node* node : theMesh_.getVerticesList(IteratorType::GLOBAL))
            {
                PointPhysicalT point = node->getElement(0)->getPhysicalGeometry()->getLocalNodeCoordinates(node->getVertexNr(0));
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
                    Node* current = theMesh_.getVerticesList(IteratorType::GLOBAL)[entry.second];
                    for (Element* element : current->getElements())
                    {
                        for (std::size_t i = 0; i < element->getNrOfNodes(); ++i)
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
                    PointPhysicalT firstNode(DIM), secondNode(DIM);
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
                    PointPhysicalT firstNode(DIM), secondNode(DIM);
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
                    PointPhysicalT firstNode(DIM), secondNode(DIM);
                    std::vector<std::size_t> nodeIndices = edge->getElement(0)->getReferenceGeometry()->getCodim2EntityLocalIndices(edge->getEdgeNr(0));
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
        std::vector<PointPhysicalT> movement(theMesh_.getNodes().size(), DIM);
        //stop after n iterations, or (when the nodes have stopped moving and the mesh is not becoming worse), or when the mesh is great, or when the mesh is decent, but worsening
        while ((counter < 10000 && (maxMovement > 1e-3 || oldQuality - worstQuality > 1e-3) && worstQuality < 0.8 && (worstQuality < 2. / 3. || oldQuality - worstQuality < 0)) || counter < 5)
        {
            counter++;
            if ((maxShift > 0.1 && (oldQuality - worstQuality) < 5e-3 * maxShift) || (oldNodeLocations_.size() != theMesh_.getNumberOfNodes()) || worstQuality < 1e-6)
            {
                maxShift = 0;
                
                orgQhull::RboxPoints qHullCoordinates {};
                qHullCoordinates.setDimension(DIM);
                qHullCoordinates.reserveCoordinates(DIM * theMesh_.getNumberOfNodes());
                oldNodeLocations_.clear();
                oldNodeLocations_.reserve(theMesh_.getNumberOfNodes());
                for (PointPhysicalT point : theMesh_.getNodes())
                {
                    qHullCoordinates.append(DIM, point.data());
                    oldNodeLocations_.push_back(point);
                }
                theMesh_.clear();
                auto pairingIterator = periodicPairing.begin();
                for (PointPhysicalT point : oldNodeLocations_)
                {
                    theMesh_.addNode(point);
                    if (pairingIterator == periodicPairing.end())
                    {
                        theMesh_.addVertex();
                    }
                    else
                    {
                        //skip one insertion for each master/slave pair
                        pairingIterator++;
                    }
                }
                
                std::vector<std::size_t> vertexIndex {};
                vertexIndex.resize(theMesh_.getNumberOfNodes(), std::numeric_limits<std::size_t>::max());
                pairingIterator = periodicPairing.begin();
                std::size_t currentVertexNumber = 0;
                for (std::size_t i = 0; i < theMesh_.getNumberOfNodes();)
                {
                    vertexIndex[i] = currentVertexNumber;
                    //assign boundary nodes
                    while (pairingIterator != periodicPairing.end() && pairingIterator->first == i)
                    {
                        vertexIndex[pairingIterator->second] = currentVertexNumber;
                        ++pairingIterator;
                    }
                    currentVertexNumber++;
                    //skip over already set boundary nodes
                    while (vertexIndex[i] < std::numeric_limits<std::size_t>::max() && i < theMesh_.getNumberOfNodes())
                    {
                        ++i;
                    }
                }
                
                //all periodic boundary pairs are used
                logger.assert(pairingIterator == periodicPairing.end(), "Somehow missed some periodic");
                //the actual amount of vertices and the assigned amount of vertices match
                logger.assert(currentVertexNumber == theMesh_.getNumberOfVertices(IteratorType::GLOBAL), "Missed some node indexes");
                
                orgQhull::Qhull triangulation(qHullCoordinates, "d PF1e-10 QbB Qx Qc Qt");
                
                for (orgQhull::QhullFacet triangle : triangulation.facetList())
                {
                    if (triangle.isGood() && !triangle.isUpperDelaunay())
                    {
                        PointPhysicalT center {DIM};
                        std::vector<std::size_t> pointIndices {};
                        for (auto vertexIt1 = triangle.vertices().begin(); vertexIt1 != triangle.vertices().end(); ++vertexIt1)
                        {
                            center += oldNodeLocations_[(*vertexIt1).point().id()];
                            pointIndices.push_back((*vertexIt1).point().id());
                        }
                        center = center / pointIndices.size();
                        if (domainDescription(center) < -1e-10)
                        {
                            auto newElement = addElement(pointIndices);
                            for (std::size_t i = 0; i < pointIndices.size(); ++i)
                            {
                                theMesh_.getVerticesList(IteratorType::GLOBAL)[vertexIndex[pointIndices[i]]]->addElement(newElement, i);
                            }
                        }
                    }
                    if (!triangle.isGood() && !triangle.isUpperDelaunay())
                    {
                        logger(INFO, "small element % ignored", triangle);
                    }
                }
                for (Node* node : theMesh_.getVerticesList(IteratorType::GLOBAL))
                {
                    //all of the nodes should be in the interior of the domain or near the boundary of the domain
                    if (node->getNrOfElements() == 0)
                    {
                        for (std::size_t i = 0; i < vertexIndex.size(); ++i)
                        {
                            if (vertexIndex[i] == node->getID())
                            {
                                logger(DEBUG, "% % %", i, theMesh_.getNodes()[i], domainDescription(theMesh_.getNodes()[i]));
                            }
                        }
                    }
                    logger.assert(node->getNrOfElements() > 0, "There is an node without any elements connected to it");
                }
                edgeFactory();
                faceFactory();
            }
            oldQuality = worstQuality;
            
            std::vector<double> expectedLength {};
            std::multimap<double, std::size_t> knownLengths {};
            expectedLength.reserve(theMesh_.getNumberOfNodes());
            //for proper scaling
            totalCurrentLength = 0;
            double totalexpectedLength = 0.;
            for (Node* node : theMesh_.getVerticesList(IteratorType::GLOBAL))
            {
                PointPhysicalT point = node->getElement(0)->getPhysicalGeometry()->getLocalNodeCoordinates(node->getVertexNr(0));
                expectedLength.push_back(relativeEdgeLength(point));
                if (!isnan(expectedLength.back()) && !isinf(expectedLength.back()))
                {
                    knownLengths.insert( {expectedLength.back(), expectedLength.size() - 1});
                    needsExpansion |= true;
                }
            }
            if (needsExpansion)
            {
                //iterate over all nodes, sorted by edge lengths
                for (std::pair<double, std::size_t> entry : knownLengths)
                {
                    Node* current = theMesh_.getVerticesList(IteratorType::GLOBAL)[entry.second];
                    for (Element* element : current->getElements())
                    {
                        for (std::size_t i = 0; i < element->getNrOfNodes(); ++i)
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
                    PointPhysicalT firstNode(DIM), secondNode(DIM);
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
                    PointPhysicalT firstNode(DIM), secondNode(DIM);
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
                    PointPhysicalT firstNode(DIM), secondNode(DIM);
                    std::vector<std::size_t> nodeIndices = edge->getElement(0)->getReferenceGeometry()->getCodim2EntityLocalIndices(edge->getEdgeNr(0));
                    firstNode = edge->getElement(0)->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[0]);
                    secondNode = edge->getElement(0)->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[1]);
                    currentLength[edge->getID()] = L2Norm(firstNode - secondNode);
                    totalCurrentLength += currentLength[edge->getID()] * currentLength[edge->getID()] * currentLength[edge->getID()];
                    totalexpectedLength += std::pow(expectedLength[edge->getElement(0)->getNode(nodeIndices[0])->getID()] + expectedLength[edge->getElement(0)->getNode(nodeIndices[1])->getID()], 3.) / 8.;
                }
            }
            
            for (PointPhysicalT& point : movement)
            {
                point *= 0;
            }
            
            if (DIM == 1)
            {
                for (Element* element : theMesh_.getElementsList(IteratorType::GLOBAL))
                {
                    PointPhysicalT firstNode(DIM), secondNode(DIM);
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
                    PointPhysicalT firstNode(DIM), secondNode(DIM);
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
                    PointPhysicalT firstNode(DIM), secondNode(DIM);
                    std::vector<std::size_t> nodeIndices = edge->getElement(0)->getReferenceGeometry()->getCodim2EntityLocalIndices(edge->getEdgeNr(0));
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
            for (std::size_t i = 0; i < theMesh_.getNumberOfVertices(IteratorType::GLOBAL); ++moveIterator, ++i)
            {
                if (fixIterator != fixedPointIdxs.end() && i == *fixIterator)
                {
                    ++fixIterator;
                    *moveIterator *= 0;
                }
                else
                {
                    Node* node = theMesh_.getVerticesList(IteratorType::GLOBAL)[i];
                    PointPhysicalT& point = theMesh_.getNodes()[node->getElement(0)->getPhysicalGeometry()->getNodeIndex(node->getVertexNr(0))];
                    point += 0.1 * (*moveIterator);
                    logger.assert(!(std::isnan(point[0])), "%", i);
                    bool isPeriodic = false;
                    std::map<std::size_t, bool> hasMoved {};
                    hasMoved[node->getElement(0)->getPhysicalGeometry()->getNodeIndex(node->getVertexNr(0))] = true;
                    for (std::size_t j = 1; j < node->getNrOfElements(); ++j)
                    {
                        if (!hasMoved[node->getElement(j)->getPhysicalGeometry()->getNodeIndex(node->getVertexNr(j))])
                        {
                            PointPhysicalT& other = theMesh_.getNodes()[node->getElement(j)->getPhysicalGeometry()->getNodeIndex(node->getVertexNr(j))];
                            other += 0.1 * (*moveIterator);
                            hasMoved[node->getElement(j)->getPhysicalGeometry()->getNodeIndex(node->getVertexNr(j))] = true;
                            isPeriodic = true;
                        }
                    }
                    if (domainDescription(point) > 0 && !isPeriodic)
                    {
                        //the point is outside of the domain, move it back inside
                        double currentValue = domainDescription(point);
                        LinearAlgebra::NumericalVector gradient (DIM);
                        LinearAlgebra::NumericalVector offset (DIM);
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
                        //do a total of tree newton iteration before giving up
                        PointPhysicalT point {DIM};
                        for (std::size_t j = 0; j < 4; ++j)
                        {
                            //make sure the node stays on the periodic boundary, to prevent faces with 3 or more elements connected to them
                            for (std::size_t k = 0; k < node->getNrOfElements(); ++k)
                            {
                                point = node->getElement(k)->getPhysicalGeometry()->getLocalNodeCoordinates(node->getVertexNr(k));
                                double currentValue = domainDescription(point);
                                if (currentValue > 0)
                                {
                                    LinearAlgebra::NumericalVector gradient (DIM);
                                    LinearAlgebra::NumericalVector offset (DIM);
                                    for (std::size_t l = 0; l < DIM; ++l)
                                    {
                                        offset[l] = 1e-7;
                                        gradient[l] = (currentValue - domainDescription(point + offset)) * 1e7;
                                        offset[l] = 0;
                                    }
                                    std::map<std::size_t, bool> hasMoved;
                                    for (std::size_t l = 0; l < node->getNrOfElements(); ++l)
                                    {
                                        if (!hasMoved[node->getElement(l)->getPhysicalGeometry()->getNodeIndex(node->getVertexNr(l))])
                                        {
                                            PointPhysicalT& other = theMesh_.getNodes()[node->getElement(l)->getPhysicalGeometry()->getNodeIndex(node->getVertexNr(l))];
                                            other += currentValue * gradient / L2Norm(gradient);
                                            hasMoved[node->getElement(l)->getPhysicalGeometry()->getNodeIndex(node->getVertexNr(l))] = true;
                                        }
                                    }
                                    *moveIterator += 10 * currentValue * gradient / L2Norm(gradient);
                                }
                            }
                        }
                        for (std::size_t j = 0; j < node->getNrOfElements(); ++j)
                        {
                            point = node->getElement(j)->getPhysicalGeometry()->getLocalNodeCoordinates(node->getVertexNr(j));
                            if (domainDescription(point) > 1e-10)
                            {
                                logger(WARN, "NOTE: Failed to move periodic point % (%) back to the periodic boundary.\n "
                                        "Distance from boundary is %. Algorithm may crash.\n "
                                        "Consider fixing points at corners to remedy this issue.", i, point, domainDescription(point));
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
                    PointPhysicalT firstNode(DIM), secondNode(DIM);
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
                    PointPhysicalT firstNode(DIM), secondNode(DIM);
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
                    Geometry::PointReference center = element->getReferenceGeometry()->getCenter();
                    Geometry::Jacobian jac = element->calcJacobian(center);
                    worstQuality = std::min(worstQuality, jac.determinant() / average * std::sqrt(2));
                }
                for (Edge* edge : theMesh_.getEdgesList(IteratorType::GLOBAL))
                {
                    PointPhysicalT firstNode(DIM), secondNode(DIM);
                    std::vector<std::size_t> nodeIndices = edge->getElement(0)->getReferenceGeometry()->getCodim2EntityLocalIndices(edge->getEdgeNr(0));
                    firstNode = edge->getElement(0)->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[0]);
                    secondNode = edge->getElement(0)->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[1]);
                    maxMovement = std::max(maxMovement, L2Norm(movement[edge->getElement(0)->getNode(nodeIndices[0])->getID()]) / 10 / currentLength[edge->getID()]);
                    maxMovement = std::max(maxMovement, L2Norm(movement[edge->getElement(0)->getNode(nodeIndices[1])->getID()]) / 10 / currentLength[edge->getID()]);
                    maxShift = std::max(maxShift, L2Norm(firstNode - oldNodeLocations_[edge->getElement(0)->getPhysicalGeometry()->getNodeIndex(nodeIndices[0])]) / currentLength[edge->getID()]);
                    maxShift = std::max(maxShift, L2Norm(secondNode - oldNodeLocations_[edge->getElement(0)->getPhysicalGeometry()->getNodeIndex(nodeIndices[1])]) / currentLength[edge->getID()]);
                }
            }
            
            //no teleporting nodes in the final iteration
            if (counter % 50 == 1 && false)
            {
                //the actual sorting is more expensive than computing the lengths and this does not happen very often
                std::multimap<double, std::pair<PointPhysicalT, PointIndexT> > centerPoints {};
                if (DIM == 1)
                {
                    //the algorithm is mostly dimension independent, but the data type it operates on is not
                    for (Element* element : theMesh_.getElementsList(IteratorType::GLOBAL))
                    {
                        PointPhysicalT firstNode(DIM), secondNode(DIM);
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
                        PointPhysicalT firstNode(DIM), secondNode(DIM);
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
                        PointPhysicalT firstNode(DIM), secondNode(DIM);
                        std::vector<std::size_t> nodeIndices = edge->getElement(0)->getReferenceGeometry()->getCodim2EntityLocalIndices(edge->getEdgeNr(0));
                        firstNode = edge->getElement(0)->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[0]);
                        secondNode = edge->getElement(0)->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIndices[1]);
                        //length is scaled in case somebody hasty decides to add smoothing at this point
                        //all edges should be squeezed a little if the algorithm is to work correctly so pretend the volume is 1.5 times as large
                        //remember to scale back from a volume measure to a length measure
                        double length = (expectedLength[edge->getElement(0)->getNode(nodeIndices[0])->getID()] + expectedLength[edge->getElement(0)->getNode(nodeIndices[1])->getID()]) / currentLength[centerPoints.size()] * std::pow(2 * totalCurrentLength / totalexpectedLength, 1. / 3.);
                        centerPoints.insert( {length, {(firstNode + secondNode) / 2, edge->getElement(0)->getNode(nodeIndices[0])->getID()}});
                    }
                }
                std::vector<bool> hasTeleported(theMesh_.getNumberOfNodes(), false);
                auto longEdge = centerPoints.begin();
                auto shortEdge = centerPoints.rbegin();
                for (std::size_t index : fixedPointIdxs)
                {
                    hasTeleported[index] = true;
                }
                PointPhysicalT point {DIM};
                PointPhysicalT other {DIM};
                for (Node* node : theMesh_.getVerticesList(IteratorType::GLOBAL))
                {
                    point = node->getElement(0)->getPhysicalGeometry()->getLocalNodeCoordinates(node->getVertexNr(0));
                    for (std::size_t i = 0; i < node->getNrOfElements(); ++i)
                    {
                        other = node->getElement(i)->getPhysicalGeometry()->getLocalNodeCoordinates(node->getVertexNr(i));
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
                            maxMovement = std::max(maxMovement, L2Norm(longEdge->second.first - theMesh_.getNodes()[shortEdge->second.second]));
                            //it is quite unlikely that the current triangulation suffices after randomly teleporting nodes about
                            maxShift = std::numeric_limits<double>::infinity();
                            theMesh_.getNodes()[shortEdge->second.second] = longEdge->second.first;
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
            element->getReferenceToPhysicalMap()->reinit(element->getPhysicalGeometry());
        }
    }
#endif
    
    /// \bug does not do the bc flags yet
    void MeshManipulator::faceFactory()
    {   
        std::vector<std::size_t> nodeIndices;
        std::vector<Element*> candidates;
        
        for (Element* element : theMesh_.getElementsList(IteratorType::GLOBAL))
        {
            for (std::size_t i = 0; i < element->getNrOfFaces(); ++i)
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
                    }
                    
                    //the current element does not bound the face or more than two elements bound the face
                    logger.assert_always(candidates.size() == 1 || candidates.size() == 2, 
                                         "Invalid number of bounding elements detected for face %", 
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
                        for (std::size_t j = 0; j < other->getNrOfFaces(); ++j)
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
        
        logger(INFO, "Total number of Faces: %", getFacesList(IteratorType::GLOBAL).size());
    }
    
    //the algorithm for the edge factory is based on that of the face factory
    //with some minor adaptation to account for the fact that there may be
    //more than two elements per edge
    ///\bug does not do 4D yet
    void MeshManipulator::edgeFactory()
    {
        std::size_t DIM(configData_->dimension_);
        //'edges' in DIM 2 are actually nodes
        if (DIM != 2)
        {
            std::vector<std::size_t> nodeList, otherNodeList;
            
            const Node* nodes[2];
            
            for (Element* element : theMesh_.getElementsList(IteratorType::GLOBAL))
            {
                for (std::size_t i = 0; i < element->getNrOfEdges(); ++i)
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
                        Edge* newEdge = theMesh_.getEdgesList(IteratorType::GLOBAL).back();
                        newEdge->addElement(element, i);
                        for (std::size_t j = 1; j < candidates.size(); ++j)
                        {
                            Element* other = candidates[j];
                            for (std::size_t k = 0; k < other->getNrOfEdges(); ++k)
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
    
    Mesh& MeshManipulator::getMesh()
    {
        return theMesh_;
    }
    
    const Mesh& MeshManipulator::getMesh() const
    {
        return theMesh_;
    }
    
    std::size_t MeshManipulator::dimension() const
    {
        return configData_->dimension_;
    }
}

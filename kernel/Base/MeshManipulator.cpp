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

#ifndef meshmanipulator_impl_h_
#define meshmanipulator_impl_h_

#include "MeshManipulator.hpp"

#include "Geometry/PhysicalGeometry.hpp"
#include "Geometry/ReferenceTriangle.hpp"
#include "Edge.hpp"
#include "Base/BasisFunctionSet.hpp"
#include "ConfigurationData.hpp"
#include "Element.hpp"
#include "Face.hpp"
#include "MeshMoverBase.hpp"
#include "AssembleBasisFunctionSet.hpp"
#include "OrientedBasisFunctionSet.hpp"
#include "Geometry/PointPhysical.hpp"
#include "Geometry/ReferenceSquare.hpp"
#include "Geometry/GlobalNamespaceGeometry.hpp"
#include "Geometry/PointReference.hpp"
#include "ElementCacheData.hpp"
#include "FaceCacheData.hpp"
#include "BaseBasisFunction.hpp"
#include "Geometry/Mappings/MappingReferenceToPhysical.hpp"
#include "ElementFactory.hpp"
#include "FaceFactory.hpp"

#include <cassert>
#include <algorithm>
#include <iostream>

//
//  MeshManipulator.cpp
//
//
//  Created by Shavarsh Nurijanyan on 9/4/13.
//
//


namespace Base {
    /// \brief This is the function that checks if two halfFaces nned to swapped or not.
    /// \details

    /*!
     *  Implements a named version of operator<() that sorts first on first index, then on second index and so on
     !*/
    bool compareHalfFace(HalfFaceDescription first, HalfFaceDescription second) {
        throw "please send your stack trace to freekjan";
        assert(first.nodeList.size() == second.nodeList.size());
        for (int i = 0; i < first.nodeList.size(); ++i) {
            if (first.nodeList[i] > second.nodeList[i]) {
                return true;
            } else if (second.nodeList[i] > first.nodeList[i]) {
                return false;
            }
        }
        //give a consistent return value
        return first.elementNum > second.elementNum;
    }

    //deleted a failed attempt to implement compareHalfFace that was already commented out see revision <325 for details -FB 
    //***********************************************************************************************************************
    //***********************************************************************************************************************
    //***********************************************************************************************************************

    void
    MeshManipulator::createDefaultBasisFunctions(unsigned int order) {
        Base::BasisFunctionSet* bFset1 = new Base::BasisFunctionSet(order);
        switch (configData_->dimension_) {
            case 1:
                switch (order) {
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
                        std::cout << "WARNING: No default basisFunction sets have been defined for this polynomial order; defaulting to 2" << std::endl;
                        const_cast<Base::ConfigurationData*> (configData_)->polynomialOrder_ = 2;
                        delete bFset1;
                        bFset1 = new Base::BasisFunctionSet(2);
                        Base::AssembleBasisFunctionSet_1D_Ord2_A0(*bFset1);
                }
                break;
            case 2:
                switch (order) {
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
                        std::cout << "WARNING: No default basisFunction sets have been defined for this polynomial order; defaulting to 2" << std::endl;
                        const_cast<Base::ConfigurationData*> (configData_)->polynomialOrder_ = 2;
                        delete bFset1;
                        bFset1 = new Base::BasisFunctionSet(2);
                        Base::AssembleBasisFunctionSet_2D_Ord2_A1(*bFset1);
                }
                break;
            case 3:
                switch (order) {
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
                        std::cout << "WARNING: No default basisFunction sets have been defined for this polynomial order; defaulting to 2" << std::endl;
                        const_cast<Base::ConfigurationData*> (configData_)->polynomialOrder_ = 2;
                        delete bFset1;
                        bFset1 = new Base::BasisFunctionSet(2);
                        Base::AssembleBasisFunctionSet_3D_Ord2_A1(*bFset1);
                }
                break;
            default:
                throw "No basisfunctions exist in this dimension";
        }
        if (collBasisFSet_.size() == 0) {
            collBasisFSet_.resize(1);
        }
        collBasisFSet_[0] = bFset1;
    }

    MeshManipulator::MeshManipulator(const ConfigurationData* config, bool xPer, bool yPer, bool zPer, unsigned int orderOfFEM, unsigned int idRangeBegin, int nrOfElementMatrixes, int nrOfElementVectors, int nrOfFaceMatrtixes, int nrOfFaceVectors) :
    configData_(config),
    periodicX_(xPer),
    periodicY_(yPer),
    periodicZ_(zPer),
    //activeMeshTree_(0),
    //numMeshTree_(0),
    numberOfElementMatrixes_(nrOfElementMatrixes),
    numberOfElementVectors_(nrOfElementVectors),
    numberOfFaceMatrixes_(nrOfFaceMatrtixes),
    numberOfFaceVectors_(nrOfFaceVectors),
    meshMover_(nullptr) {

        std::cout << "******Mesh creation started!**************" << std::endl;
        unsigned int DIM = configData_->dimension_;
        for (unsigned int i = 0; i < DIM; ++i) {
            if (i == 0)
                std::cout << "Boundries: " << (periodicX_ ? "Periodic  " : "Solid Wall") << " in X direction" << std::endl;
            if (i == 1)
                std::cout << "Boundries: " << (periodicY_ ? "Periodic " : "Solid Wall") << " in Y direction" << std::endl;
            if (i == 2)
                std::cout << "Boundries: " << (periodicZ_ ? "Periodic " : "Solid Wall") << " in Z direction" << std::endl;
        }
        createDefaultBasisFunctions(orderOfFEM);
        const_cast<ConfigurationData*> (configData_)->numberOfBasisFunctions_ = collBasisFSet_[0]->size();

        //createNewMeshTree();
        std::cout << "******Mesh creation is finished!**********" << std::endl;
    }

    MeshManipulator::MeshManipulator(const MeshManipulator& other) :
    configData_(other.configData_),
    theMesh_(other.theMesh_),
    periodicX_(other.periodicX_),
    periodicY_(other.periodicY_),
    periodicZ_(other.periodicZ_),
    meshMover_(other.meshMover_),
    //defaultSetOfBasisFunctions_(other.defaultSetOfBasisFunctions_),
    collBasisFSet_(other.collBasisFSet_),
    //activeMeshTree_(other.activeMeshTree_),
    //numMeshTree_(other.numMeshTree_),
    //vecOfElementTree_(other.vecOfElementTree_),
    //vecOfFaceTree_(other.vecOfFaceTree_),
    numberOfElementMatrixes_(other.numberOfElementMatrixes_),
    numberOfElementVectors_(other.numberOfElementVectors_),
    numberOfFaceMatrixes_(other.numberOfFaceMatrixes_),
    numberOfFaceVectors_(other.numberOfFaceVectors_) {

    }

    MeshManipulator::~MeshManipulator() {

        for (typename CollectionOfBasisFunctionSets::iterator bit = collBasisFSet_.begin(); bit != collBasisFSet_.end(); ++bit) {
            const BasisFunctionSetT* bf = *bit;
            delete bf; ///\bug segfaults when using two meshes with the same sets of basisfunctions
        }

        delete meshMover_;
        //delete defaultSetOfBasisFunctions_;


        // Kill all elements in all mesh-tree
        //while (!vecOfElementTree_.empty()) {
            //delete vecOfElementTree_.back();
            //vecOfElementTree_.pop_back();
        //}
    }

    void
    MeshManipulator::setDefaultBasisFunctionSet(BasisFunctionSetT* bFSet) {
        delete collBasisFSet_[0];
        collBasisFSet_[0] = bFSet;
        const_cast<ConfigurationData*> (configData_)->numberOfBasisFunctions_ = bFSet->size();
        for (Base::Face* face : getFacesList(IteratorType::GLOBAL)) {
            face->setLocalNrOfBasisFunctions(0);
        }
        for (Base::Edge* edge : getEdgesList(IteratorType::GLOBAL)) {
            edge->setLocalNrOfBasisFunctions(0);
        }
        for (Base::Node* node : getVerticesList(IteratorType::GLOBAL)) {
            node->setLocalNrOfBasisFunctions(0);
        }
        for (ElementIterator it = elementColBegin(IteratorType::GLOBAL); it != elementColEnd(IteratorType::GLOBAL); ++it) {
            (*it)->setDefaultBasisFunctionSet(0);
        }
    }

    void
    MeshManipulator::addVertexBasisFunctionSet(CollectionOfBasisFunctionSets& bFsets) {
        int firstNewEntry = collBasisFSet_.size();
        for (const BasisFunctionSet* it : bFsets) {
            collBasisFSet_.push_back(it);
        }
        for (Node* node : getVerticesList()) {
            for (int i = 0; i < node->getNrOfElements(); ++i) {
                node->getElement(i)->setVertexBasisFunctionSet(firstNewEntry + node->getVertexNr(i), node->getVertexNr(i));
            }
            node->setLocalNrOfBasisFunctions(bFsets[0]->size());
        }
        const_cast<ConfigurationData*> (configData_)->numberOfBasisFunctions_ += (*elementColBegin())->getNrOfNodes() * bFsets[0]->size();
    }

    void
    MeshManipulator::addFaceBasisFunctionSet(std::vector<const OrientedBasisFunctionSet*>& bFsets) {
        int firstNewEntry = collBasisFSet_.size();
        for (const BasisFunctionSet* it : bFsets) {
            collBasisFSet_.push_back(it);
        }
        for (Face* face : getFacesList()) {
            int faceNr = face->localFaceNumberLeft();
            for (int i = 0; i < bFsets.size(); ++i) {
                if (bFsets[i]->checkOrientation(0, faceNr)) {
                    face->getPtrElementLeft()->setFaceBasisFunctionSet(firstNewEntry + i, faceNr);
                }
            }
            if (face->isInternal()) {
                faceNr = face->localFaceNumberRight();
                int orientation = face->getFaceToFaceMapIndex();
                for (int i = 0; i < bFsets.size(); ++i) {
                    if (bFsets[i]->checkOrientation(orientation, faceNr)) {
                        face->getPtrElementRight()->setFaceBasisFunctionSet(firstNewEntry + i, faceNr);
                    }
                }
            }
            face->setLocalNrOfBasisFunctions(bFsets[0]->size());
        }
        const_cast<ConfigurationData*> (configData_)->numberOfBasisFunctions_ += (*elementColBegin())->getPhysicalGeometry()->getNrOfFaces() * bFsets[0]->size();
    }

    void
    MeshManipulator::addEdgeBasisFunctionSet(std::vector<const OrientedBasisFunctionSet*>& bFsets) {
        int firstNewEntry = collBasisFSet_.size();
        for (const BasisFunctionSet* it : bFsets) {
            collBasisFSet_.push_back(it);
        }
        for (Edge* edge : getEdgesList()) {
            for (int i = 0; i < edge->getNrOfElements(); ++i) {
                for (int j = 0; j < bFsets.size(); ++j) {
                    if (bFsets[j]->checkOrientation(edge->getOrientation(i), edge->getEdgeNr(i))) {
                        edge->getElement(i)->setEdgeBasisFunctionSet(firstNewEntry + j, edge->getEdgeNr(i));
                        //cout<<edge->getOrientation(i)<<edge->getEdgeNr(i)<<bFsets[j]->size()<<endl;
                    }
                }
            }
            edge->setLocalNrOfBasisFunctions(bFsets[0]->size());
        }
        const_cast<ConfigurationData*> (configData_)->numberOfBasisFunctions_ += (*elementColBegin())->getPhysicalGeometry()->getNrOfFaces() * bFsets[0]->size();
    }

    Base::Element*
    MeshManipulator::addElement(const VectorOfPointIndicesT& globalNodeIndexes) {
        return theMesh_.addElement(globalNodeIndexes);
    }

    void
    MeshManipulator::move() {
        for (Geometry::PointPhysical& p : theMesh_.getNodes()) {
            if (meshMover_ != NULL) {
                meshMover_->movePoint(p);
            }
        }
        /*for(Base::Element* element:getElementsList()){
            const_cast<Geometry::MappingReferenceToPhysical*>(element->getReferenceToPhysicalMap())->reinit(element->getPhysicalGeometry());
        }*/
    }

    void
    MeshManipulator::setMeshMover(const MeshMoverBase* meshMover) {
        meshMover_ = meshMover;
    }

    bool
    MeshManipulator::addFace(ElementT* leftElementPtr, unsigned int leftElementLocalFaceNo, ElementT* rightElementPtr, unsigned int rightElementLocalFaceNo, const Geometry::FaceType& faceType) {
        return theMesh_.addFace(leftElementPtr, leftElementLocalFaceNo, rightElementPtr, rightElementLocalFaceNo, faceType);
    }

    //void
    //MeshManipulator::addEdge(std::vector< Element*> elements, std::vector<unsigned int> localEdgeNrs) {
    //    theMesh_.addEdge(elements, localEdgeNrs);
    //}
    
    void MeshManipulator::addEdge()
    {
        theMesh_.addEdge();
    }

    void
    MeshManipulator::addVertex() {
        theMesh_.addVertex();
    }

    void
    MeshManipulator::outputMesh(std::ostream& os)const {
        for (Geometry::PointPhysical p : getNodes()) {
            os << "Node " << " " << p << std::endl;
        }

        int elementNum = 0;

        for (Element* element : getElementsList()) {
            os << "Element " << element->getID() << " " << element << std::endl;
            elementNum++;
        }

        /*int faceNum=0;
        for (typename ListOfFacesT::const_iterator cit=faces_.begin(); cit !=faces_.end(); ++cit)
        {
            /// \bug need at output routine for the face to test this.
            // os << "Face " <<faceNum <<" " <<(*cit)<<endl;
            faceNum++;
        }*/

    }

    void
    MeshManipulator::createRectangularMesh(const PointPhysicalT& BottomLeft, const PointPhysicalT& TopRight, const VectorOfPointIndicesT& linearNoElements) {
        //set to correct value in case some other meshmanipulator changed things
        ElementFactory::instance().setCollectionOfBasisFunctionSets(&collBasisFSet_);
        ElementFactory::instance().setNumberOfMatrices(numberOfElementMatrixes_);
        ElementFactory::instance().setNumberOfVectors(numberOfFaceVectors_);
        ElementFactory::instance().setNumberOfTimeLevels(configData_->numberOfTimeLevels_);
        ElementFactory::instance().setNumberOfUnknowns(configData_->numberOfUnknowns_);
        FaceFactory::instance().setNumberOfFaceMatrices(numberOfFaceMatrixes_);
        FaceFactory::instance().setNumberOfFaceVectors(numberOfFaceVectors_);

        unsigned int DIM = configData_->dimension_;
        if (linearNoElements.size() != DIM) {
            std::cout << "The number of Linear Intervals has to map the size of the problem and current it does not" << std::endl;
            throw (10);
        }
        std::vector<bool> periodicDIM;
        for(size_t i=0;i<DIM;++i)
        {
            if(i==0) periodicDIM.push_back(periodicX_);
            if(i==1) periodicDIM.push_back(periodicY_);
            if(i==2) periodicDIM.push_back(periodicZ_);
        }
        //Stage 1 : Precompute some required values;
        ///////

        //This store the size length of the domain i.e. it is DIM sized vector
        PointPhysicalT delta_x(DIM);

        for (int i = 0; i < DIM; i++) {
            delta_x[i] = (TopRight[i] - BottomLeft[i]) / (linearNoElements[i]);
        }

        //This stores the number of nodes in each coDIMension i.e. if you have 2 by 2 element it is 3 nodes 
        //nodes mark physical location, vertices mark connectivity-based location
        std::vector<unsigned int> numOfNodesInEachSubspace(DIM), numOfElementsInEachSubspace(DIM), numOfVerticesInEachSubspace(DIM);

        numOfNodesInEachSubspace[0] = 1;
        numOfVerticesInEachSubspace[0] = 1;
        numOfElementsInEachSubspace[0] = 1;

        //This will be the total number of nodes required in the problem
        unsigned int totalNumOfNodes, totalNumOfVertices, totalNumOfElements, verticesPerElement;

        totalNumOfNodes = (linearNoElements[0] + 1);
        totalNumOfVertices = (linearNoElements[0] + (periodicDIM[0]?0:1));

        totalNumOfElements = (linearNoElements[0]);

        verticesPerElement = 2;
        int powerOf2;

        for (unsigned int iDIM = 1; iDIM < DIM; ++iDIM) {
            totalNumOfNodes *= (linearNoElements[iDIM] + 1);
            totalNumOfVertices *= (linearNoElements[iDIM] + (periodicDIM[iDIM]?0:1));
            totalNumOfElements *= (linearNoElements[iDIM]);
            verticesPerElement *= 2;

            numOfElementsInEachSubspace[iDIM] = numOfElementsInEachSubspace[iDIM - 1]*(linearNoElements[iDIM - 1]);
            numOfNodesInEachSubspace[iDIM] = numOfNodesInEachSubspace[iDIM - 1]*(linearNoElements[iDIM - 1] + 1);
            numOfVerticesInEachSubspace[iDIM] = numOfVerticesInEachSubspace[iDIM - 1]*(linearNoElements[iDIM - 1] + (periodicDIM[iDIM-1]?0:1));
        }

        //temp point for storing the node locations
        PointPhysicalT x(DIM);


        //Stage 2 : Create the nodes
        //Now loop over all the nodes and calculate the coordinates for rach DIMension (this makes the algorithm independent of DIMension
        for (unsigned int nodeIndex = 0; nodeIndex < totalNumOfNodes; ++nodeIndex) {
            unsigned int nodeIndexRemain = nodeIndex;



            for (int iDIM = DIM - 1; iDIM>-1; --iDIM) {
                x[iDIM] = BottomLeft[iDIM] + (nodeIndexRemain / numOfNodesInEachSubspace[iDIM] * delta_x[iDIM]);
                nodeIndexRemain %= numOfNodesInEachSubspace[iDIM];
            }

            //actally add the point
            theMesh_.addNode(x);

        }

        //stage 2.5 : create the vertices
        
        for(size_t vertexID = 0; vertexID < totalNumOfVertices; ++vertexID)
        {
            theMesh_.addVertex();
        }

        //Stage 3 : Create the elements

        std::vector<unsigned int> elementNdId(DIM), nodeNdId(DIM), vertexNdId(DIM), globalNodeID(verticesPerElement), globalVertexID(verticesPerElement);

        //This will store a two array, which goes elementNumber and the n-DIMensional coordinate of that element. This is only used for the face creation (step 4).
        //std::vector<std::vector<unsigned int> > allElementNdId;
        //allElementNdId.resize(totalNumOfElements);

        auto& vertexlist=getVerticesList(IteratorType::GLOBAL);
        auto& elementlist=getElementsList(IteratorType::GLOBAL);

        //elementNdId is DIM coordinate of the bottom left node i.e. in two (0,0), (1,0) ,(2,0) ... etc are the first three (if at least three elements in x)
        for (unsigned int elementIndex = 0; elementIndex < totalNumOfElements; ++elementIndex) {
            unsigned int numElementsRemaining = elementIndex;

            for (int iDIM = DIM - 1; iDIM>-1; --iDIM) {
                elementNdId[iDIM] = numElementsRemaining / numOfElementsInEachSubspace[iDIM];
                numElementsRemaining %= numOfElementsInEachSubspace[iDIM];
            }


            // vertexNdId are the DIM coordinate of each vertex in the element with vertexNdId[0] being the bottom left
            for (unsigned int i = 0; i < verticesPerElement; ++i) {
                powerOf2 = 1;
                for (unsigned int iDIM = 0; iDIM < DIM; ++iDIM) {
                    nodeNdId[iDIM] = elementNdId[iDIM] + ((i & powerOf2) != 0);
                    vertexNdId[iDIM] = elementNdId[iDIM] + ((i & powerOf2) != 0);
                    if((vertexNdId[iDIM]>=linearNoElements[iDIM]) && (periodicDIM[iDIM]==true))
                    {
                        vertexNdId[iDIM]=0;
                    }
                    powerOf2 *= 2;
                }
                globalNodeID[i] = nodeNdId[0];
                globalVertexID[i] = vertexNdId[0];

                //Now map to the one DIMensional global ID
                for (unsigned int iDIM = 1; iDIM < DIM; ++iDIM) {
                    globalNodeID[i] += nodeNdId[iDIM] * numOfNodesInEachSubspace[iDIM];
                    globalVertexID[i] += vertexNdId[iDIM] * numOfVerticesInEachSubspace[iDIM];
                }
            }
            Element* newElement = addElement(globalNodeID);
            for(size_t i=0;i<globalVertexID.size();++i)
            {
                vertexlist[globalVertexID[i]]->addElement(newElement,i);
            }
        }

        //Stage 4 : Create the faces

        faceFactory();
        edgeFactory();

        /*switch (DIM) {
            case 1:
                rectangularCreateFaces1D(tempElementVector, linearNoElements);
                break;

            case 2:
                rectangularCreateFaces2D(tempElementVector, linearNoElements);
                break;

            case 3:
                rectangularCreateFaces3D(tempElementVector, linearNoElements);
                break;

            default:
                throw ("Face generator not implemented in this DIMension");

        }//end case*/


    }

    /*void MeshManipulator::rectangularCreateFaces1D(VectorOfElementPtrT& tempElementVector, const VectorOfPointIndicesT& linearNoElements) {

        //First the internal faces
        for (int i = 0; i < linearNoElements[0] - 1; i++) {
            addFace(tempElementVector[i], 1, tempElementVector[i + 1], 0);
        }

        if (periodicX_ == 1)
 {
            addFace(tempElementVector[0], 0, tempElementVector[linearNoElements[0] - 1], 1);
        }
        else {
            addFace(tempElementVector[0], 0, NULL, 0);
            addFace(tempElementVector[linearNoElements[0] - 1], 1, NULL, 0);
        }

    }

    void MeshManipulator::rectangularCreateFaces2D(VectorOfElementPtrT& tempElementVector, const VectorOfPointIndicesT& linearNoElements) {

        //first do the faces with normals pointing in the x-direction
        Geometry::FaceType xFace = Geometry::WALL_BC;
        unsigned int index;

        for (int j = 0; j < linearNoElements[1]; j++) {
            for (int i = 0; i < linearNoElements[0] - 1; i++) {

                index = i + j * linearNoElements[0];
                addFace(tempElementVector[index], 2, tempElementVector[index + 1], 1);

            }


            index = j * linearNoElements[0];
            if (periodicX_ == 1) {
                addFace(tempElementVector[index], 1, tempElementVector[index + linearNoElements[0] - 1], 2);
            } else {
                //Left bounday face
                //cout<<"index1="<<index<<endl;
                addFace(tempElementVector[index], 1, NULL, 0);

                //Right boundary face
                index = index + linearNoElements[0] - 1;
                //cout<<"index2="<<index<<endl;
                addFace(tempElementVector[index], 2, NULL, 0);
            }
        }

        //now do the faces with normals pointing in the y-direction
        Geometry::FaceType yFace = Geometry::WALL_BC;
        for (int j = 0; j < linearNoElements[0]; j++) {
            for (int i = 0; i < linearNoElements[1] - 1; i++) {
                index = j + i * linearNoElements[0];
                //cout<<"index3="<<index<<endl;
                //cout<<"index3.5="<<index+linearNoElements[0]<<endl;

                addFace(tempElementVector[index], 3, tempElementVector[index + linearNoElements[0]], 0);

            }

            index = j;

            if (periodicY_ == 1) {

                addFace(tempElementVector[index], 0, tempElementVector[index + (linearNoElements[1] - 1) * linearNoElements[0]], 3);

            } else {

                //Bottom boundary face

                // cout<<"index4="<<index<<endl;
                addFace(tempElementVector[index], 0, NULL, 0);

                //Top boundary face
                index = index + (linearNoElements[1] - 1) * linearNoElements[0];

                //cout<<"index5="<<index<<endl;
                addFace(tempElementVector[index], 3, NULL, 0);
            }

        }


    }

    void MeshManipulator::rectangularCreateFaces3D(VectorOfElementPtrT& tempElementVector, const VectorOfPointIndicesT& linearNoElements) {
        unsigned int index;
        //first do the faces in x-direction
        //counter in y
        for (int j = 0; j < linearNoElements[1]; j++) {
            //counter in z
            for (int k = 0; k < linearNoElements[2]; k++) {
                //counter in x
                for (int i = 0; i < linearNoElements[0] - 1; i++) {

                    index = i + j * linearNoElements[0] + k * linearNoElements[0] * linearNoElements[1];
                    addFace(tempElementVector[index], 3, tempElementVector[index + 1], 2);

                }//end loop over x


                index = j * linearNoElements[0] + k * linearNoElements[0] * linearNoElements[1];
                if (periodicX_ == 1) {
                    addFace(tempElementVector[index], 2, tempElementVector[index + linearNoElements[0] - 1], 3);
                } else {

                    //Left boundary

                    // cout<<"index1="<<index<<endl;
                    addFace(tempElementVector[index], 2, NULL, 0);

                    //Right boundary face
                    index = index + linearNoElements[0] - 1;
                    // cout<<"index2="<<index<<endl;
                    addFace(tempElementVector[index], 3, NULL, 0);
                } // end boundary cases

            } // end loop over z
        } //end loop over y


        //Now do the faces pointing in the y-direction
        //count in z
        for (int k = 0; k < linearNoElements[2]; k++) {
            //counter in x
            for (int j = 0; j < linearNoElements[0]; j++) {
                //counter in y
                for (int i = 0; i < linearNoElements[1] - 1; i++) {

                    index = j + i * linearNoElements[0] + k * linearNoElements[0] * linearNoElements[1];
                    addFace(tempElementVector[index], 4, tempElementVector[index + linearNoElements[0]], 1);

                }//end loop over y


                index = j + k * linearNoElements[0] * linearNoElements[1];
                if (periodicY_ == 1) {
                    //   cout << index << " , " << index+(linearNoElements[1]-1)*linearNoElements[0] << endl;
                    addFace(tempElementVector[index], 1, tempElementVector[(index + (linearNoElements[1] - 1) * linearNoElements[0])], 4);
                } else {

                    //front bounday face

                    // cout<<"index1="<<index<<endl;
                    addFace(tempElementVector[index], 1, NULL, 0);

                    //Back boundary face
                    index = index + (linearNoElements[1] - 1) * linearNoElements[0];
                    //   cout<<"index2="<<index<<endl;
                    addFace(tempElementVector[index], 4, NULL, 0);
                } // end boundary cases

            } // end loop over x
        } //end loop over z

        //Now finally do the face in the z-direction

        //count in x direction 
        for (int k = 0; k < linearNoElements[0]; k++) {
            //counter in y
            for (int j = 0; j < linearNoElements[1]; j++) {
                //counter in z
                for (int i = 0; i < linearNoElements[2] - 1; i++) {

                    index = k + j * linearNoElements[0] + i * linearNoElements[0] * linearNoElements[1];
                    addFace(tempElementVector[index], 5, tempElementVector[index + linearNoElements[0] * linearNoElements[1]], 0);

                }//end loop over z

                index = k + j * linearNoElements[0];
                if (periodicZ_ == 1) {
                    addFace(tempElementVector[index], 0, tempElementVector[index + (linearNoElements[2] - 1) * linearNoElements[1] * linearNoElements[0]], 5);
                } else {

                    //bottom boundary
                    //cout<<"index1="<<index<<endl;
                    addFace(tempElementVector[index], 0, NULL, 0);

                    //Top boundary face
                    index = index + (linearNoElements[2] - 1) * linearNoElements[1] * linearNoElements[0];
                    //cout<<"index2="<<index<<endl;
                    addFace(tempElementVector[index], 5, NULL, 0);
                } // end boundary cases

            } // end loop over y
        } //end loop over x
        edgeFactory(tempElementVector);
    }*/

    //createTrianglularMesh follows the same structure as createRectangularMesh. Where createRectangularMesh makes rectangular elements, createTrianglularMesh splits the elements into a partition of triangles.

    void MeshManipulator::createTriangularMesh(PointPhysicalT BottomLeft, PointPhysicalT TopRight, const VectorOfPointIndicesT& linearNoElements) {
        //set to correct value in case some other meshmanipulator changed things
        ElementFactory::instance().setCollectionOfBasisFunctionSets(&collBasisFSet_);
        ElementFactory::instance().setNumberOfMatrices(numberOfElementMatrixes_);
        ElementFactory::instance().setNumberOfVectors(numberOfFaceVectors_);
        ElementFactory::instance().setNumberOfTimeLevels(configData_->numberOfTimeLevels_);
        ElementFactory::instance().setNumberOfUnknowns(configData_->numberOfUnknowns_);
        FaceFactory::instance().setNumberOfFaceMatrices(numberOfFaceMatrixes_);
        FaceFactory::instance().setNumberOfFaceVectors(numberOfFaceVectors_);

        //Stage 0 : Check for required requirements
        unsigned int DIM = configData_->dimension_;
        if (linearNoElements.size() != DIM) {
            std::cout << "The number of Linear Intervals has to map the size of the problem and current it does not" << std::endl;
            throw (10);
        }

        if (DIM == 3 && periodicX_ == 1 && linearNoElements[0] % 2 == 1) {
            throw "The 3D triangular grid generator can't handle an odd amount of elements in the periodic dimension X";
        }

        if (DIM == 3 && periodicY_ == 1 && linearNoElements[1] % 2 == 1) {
            throw "The 3D triangular grid generator can't handle an odd amount of elements in the periodic dimension Y";
        }

        if (DIM == 3 && periodicZ_ == 1 && linearNoElements[2] % 2 == 1) {
            throw "The 3D triangular grid generator can't handle an odd amount of elements in the periodic dimension Z";
        }

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

        for (int i = 0; i < DIM; ++i) {
            delta_x[i] = (TopRight[i] - BottomLeft[i]) / (linearNoElements[i]);
        }

        std::vector<size_t> numOfNodesInEachSubspace(DIM), numOfVerticesInEachSubspace(DIM), numOfElementsInEachSubspace(DIM);

        numOfNodesInEachSubspace[0] = 1;
        numOfVerticesInEachSubspace[0] = 1;
        numOfElementsInEachSubspace[0] = 1;

        unsigned int totalNumOfNodes, totalNumOfElements, verticesPerElement,
                verticesPerGroup, trianglesPerRectangle, totalNumOfVertices;

        totalNumOfNodes = (linearNoElements[0] + 1);
        totalNumOfVertices = (linearNoElements[0] + (periodicDIM[0]?0:1));
        //'elements' in this counter denote groups of trianglesPerRectangle elements
        totalNumOfElements = (linearNoElements[0]);
        verticesPerElement = 2;
        verticesPerGroup = 2;
        trianglesPerRectangle = 1;
        int powerOf2;

        //start with 1 because you want to ask for the entry at idim - 1
        for (int idim = 1; idim < DIM; ++idim) {
            totalNumOfNodes *= (linearNoElements[idim] + 1);
            totalNumOfVertices *= (linearNoElements[idim] + (periodicDIM[idim]?0:1));
            totalNumOfElements *= (linearNoElements[idim]);
            verticesPerElement += 1;
            verticesPerGroup *= 2;
            trianglesPerRectangle += (2 * idim - 1);
            numOfElementsInEachSubspace[idim] = numOfElementsInEachSubspace[idim - 1] * (linearNoElements[idim - 1]);
            numOfNodesInEachSubspace[idim] = numOfNodesInEachSubspace[idim - 1] * (linearNoElements[idim - 1] + 1);
            numOfVerticesInEachSubspace[idim] = numOfNodesInEachSubspace[idim - 1] * (linearNoElements[idim - 1] + (periodicDIM[idim-1]?0:1));
        }

        //'elements' in this counter denote groups of trianglesPerRectangle elements
        //totalNumOfElements*=(trianglesPerRectangle);
        //for(int idim=0;idim<DIM;++idim)
        //{
        // numOfElementsInEachSubspace[idim]*=(trianglesPerRectangle);
        //}

        PointPhysicalT x(DIM);

        //Stage 2 : Create the nodes

        for (int nodeIndex = 0; nodeIndex < totalNumOfNodes; ++nodeIndex) {
            int nodeIndexRemain = nodeIndex;
            for (int idim = DIM - 1; idim > -1; --idim) {
                x[idim] = BottomLeft[idim] + (nodeIndexRemain / numOfNodesInEachSubspace[idim] * delta_x[idim]);
                nodeIndexRemain %= numOfNodesInEachSubspace[idim];
            }
            theMesh_.addNode(x);
        }
        
        for(size_t nodeIndex=0; nodeIndex < totalNumOfVertices; ++nodeIndex)
        {
            theMesh_.addVertex();
        }

        auto& vertices = getVerticesList(IteratorType::GLOBAL);
        auto& elements = getElementsList(IteratorType::GLOBAL);
        
        //Stage 3 : Create the elements

        std::vector<unsigned int> elementNdId(DIM), vertexNdId(DIM), nodeNdId(DIM);
        std::vector < std::vector<unsigned int> > globalVertexID(trianglesPerRectangle);
        std::vector < std::vector<unsigned int> > globalNodeID(trianglesPerRectangle);

        for (int elementGroupIndex = 0; elementGroupIndex < totalNumOfElements; ++elementGroupIndex) {
            //first generate node indexes as if we are a cube
            //indicates if the element has to be rotated to make connecting faces; rotates the element 90 degrees along the y-axis if needed
            int rotate = 0;

            for (int i = 0; i < trianglesPerRectangle; ++i) {
                globalNodeID[i].clear();
                globalVertexID[i].clear();
            }
            int elementIndexRemainder = elementGroupIndex;

            for (int idim = DIM - 1; idim > -1; --idim) {
                elementNdId[idim] = elementIndexRemainder / numOfElementsInEachSubspace[idim];
                elementIndexRemainder %= numOfElementsInEachSubspace[idim];
                rotate = (elementNdId[idim] + rotate) % 2;
            }

            for (int i = 0; i < verticesPerGroup; ++i) {
                if (rotate == 0) {
                    powerOf2 = 1;
                    for (int idim = 0; idim < DIM; ++idim) {
                        nodeNdId[idim] = elementNdId[idim] + ((i & powerOf2) != 0);
                        vertexNdId[idim] = elementNdId[idim] + ((i & powerOf2) != 0);
                        if(vertexNdId[idim]>=linearNoElements[idim] && periodicDIM[idim]) vertexNdId[idim]=0;
                        powerOf2 *= 2;
                    }
                } else {
                    powerOf2 = verticesPerGroup;
                    for (int idim = 0; idim < DIM; ++idim) {
                        powerOf2 /= 2;
                        nodeNdId[idim] = elementNdId[idim] + (((i ^ rotate) & powerOf2) != 0);
                        vertexNdId[idim] = elementNdId[idim] + (((i ^ rotate) & powerOf2) != 0);
                        if(vertexNdId[idim]>=linearNoElements[idim] && periodicDIM[idim]) vertexNdId[idim]=0;
                    }
                }
                
                int nodeIndex = nodeNdId[0];
                size_t vertexIndex = vertexNdId[0];
                for (int idim = 1; idim < DIM; ++idim) {
                    nodeIndex += nodeNdId[idim] * numOfNodesInEachSubspace[idim];
                    vertexIndex += vertexNdId[idim] * numOfVerticesInEachSubspace[idim];
                }
                
                //then cherrypick the element(s) these vertices should connect to (probably not the cleanest implementation; \bug doesn't work if DIM>3)
                switch (i) {
                    case 0:
                        globalVertexID[0].push_back(vertexIndex);
                        globalNodeID[0].push_back(nodeIndex);
                        break;
                    case 3:
                        //insert in the second place because the ordering of the vertices will work out better
                        globalVertexID[1].insert(++(globalVertexID[1].begin()), vertexIndex);
                        globalNodeID[1].insert(++(globalNodeID[1].begin()), nodeIndex);
                        break;
                    case 5:
                        globalVertexID[2].push_back(vertexIndex);
                        globalNodeID[2].push_back(nodeIndex);
                        break;
                    case 6:
                        //insert in the second place because the ordering of the vertices will work out better
                        globalVertexID[3].insert(++(globalVertexID[3].begin()), vertexIndex);
                        globalNodeID[3].insert(++(globalNodeID[3].begin()), nodeIndex);
                        break;
                    case 1:
                        for (int i = 0; i < trianglesPerRectangle; ++i) {
                            if (i != 3) {
                                globalVertexID[i].push_back(vertexIndex);
                                globalNodeID[i].push_back(nodeIndex);
                            }
                        }
                        break;
                    case 2:
                        for (int i = 0; i < trianglesPerRectangle; ++i) {
                            if (i != 2) {
                                globalVertexID[i].push_back(vertexIndex);
                                globalNodeID[i].push_back(nodeIndex);
                            }
                        }
                        break;
                    case 4:
                        for (int i = 0; i < trianglesPerRectangle; ++i) {
                            if (i != 1) {
                                globalVertexID[i].push_back(vertexIndex);
                                globalNodeID[i].push_back(nodeIndex);
                            }
                        }
                        break;
                    case 7:
                        for (int i = 0; i < trianglesPerRectangle; ++i) {
                            if (i != 0) {
                                globalVertexID[i].push_back(vertexIndex);
                                globalNodeID[i].push_back(nodeIndex);
                            }
                        }
                        break;
                } //switch
            } //for all vertices of the rectangle

            for (int i = 0; i < trianglesPerRectangle; ++i) {
                Element* newElement = addElement(globalNodeID[i]);
                for(size_t j=0;j<globalVertexID[i].size();++j)
                {
                    assert(i < globalVertexID.size());
                    assert(j < globalVertexID[i].size());
                    assert(globalVertexID[i][j]<totalNumOfVertices);
                    assert(globalVertexID[i][j] >= 0);
                    assert(vertices.size() == totalNumOfVertices);
                    vertices[globalVertexID[i][j]]->addElement(newElement,j);
                }
            }
        } //for all rectangles

        //Stage 4 : Create the faces

        /*switch (DIM) {
            case 1:
                //only provided for completeness, it is always better to just use rectangularMeshGenerator for line elements
                triangularCreateFaces1D(tempElementVector, linearNoElements);
                break;
            case 2:
                triangularCreateFaces2D(tempElementVector, linearNoElements);
                break;
            case 3:
                triangularCreateFaces3D(tempElementVector, linearNoElements);
                break;
        }*/
        
        faceFactory();
        edgeFactory();
    }

    /*void MeshManipulator::triangularCreateFaces1D(std::vector<Base::Element*>& tempElementVector, const std::vector<unsigned int>& linearNoElements) {

        //First the internal faces
        for (int i = 0; i + 1 < linearNoElements[0]; i++) {
            addFace(tempElementVector[i], 1 - i % 2, tempElementVector[i + 1], 1 - i % 2);
        }

        if (periodicX_ == 1) {
            addFace(tempElementVector[0], 0, tempElementVector[linearNoElements[0] - 1], linearNoElements[0] % 2);
        } else {
            addFace(tempElementVector[0], 0, NULL, 0, WALL_BC);
            addFace(tempElementVector[linearNoElements[0] - 1], linearNoElements[0] % 2, NULL, 0, WALL_BC);
        }

    }

    void MeshManipulator::triangularCreateFaces2D(std::vector<Base::Element*>& tempElementVector, const std::vector<unsigned int>& linearNoElements) {
        //this face creator assumes the grid has been divided in rectangles. The rectangles are then sub-divided into two triangles, which are alternatively rotated and not rotated, in a checkerboard pattern. The rectangle at DIM index (0,0,0) is not rotated. There is an ordering of the vertexes assumed that looks like the following:
        // 2 21 02 2
        // 01 0 1 01
        //not rotated rotated
        // 0 1 0 1 (element numbers within the rectangle)

        //first do the faces with normals pointing in the x-direction
        Geometry::FaceType xFace = Geometry::WALL_BC;
        unsigned int index;
        unsigned int rotate;
        for (int j = 0; j < linearNoElements[1]; ++j) {
            for (int i = 0; i + 1 < linearNoElements[0]; ++i) {
                //hard-coded number of triangles per rectangle, being able to change a constant is not going to make this easier to change
                index = 2 * (i + j * linearNoElements[0]);
                rotate = (i + j) % 2;
                addFace(tempElementVector[index + 1], 2 * rotate, tempElementVector[index + 2], rotate);
            }

            index = 2 * j * linearNoElements[0];
            rotate = (j + linearNoElements[0] - 1) % 2;
            if (periodicX_ == 1) {
                addFace(tempElementVector[index], 1 - j % 2,
                        tempElementVector[index + 2 * linearNoElements[0] - 1],
                        2 * rotate);
            } else {
                //Left bounday face

                // 	cout << "index1=" << index << endl;
                addFace(tempElementVector[index], 1 - j % 2, NULL, 0, WALL_BC);

                //Right boundary face
                index = index + 2 * linearNoElements[0] - 1;
                // 	cout << "index2=" << index << endl;
                addFace(tempElementVector[index], 2 * rotate, NULL, 0, WALL_BC);
            }
        }

        //now do the faces with normals pointing in the y-direction
        //changing loop ordering is risky+makes the code a bit less efficient (?)
        Geometry::FaceType yFace = Geometry::WALL_BC;
        for (int j = 0; j < linearNoElements[0]; ++j) {
            for (int i = 0; i + 1 < linearNoElements[1]; ++i) {
                index = 2 * (j + i * linearNoElements[0]);
                rotate = (i + j) % 2;
                // 	    cout << "index3=" << index << endl;
                // 	    cout << "index3.5=" << index + 2 * linearNoElements[0] << endl;
                addFace(tempElementVector[index + 1 - rotate], 2 - rotate, tempElementVector[index + 2 * linearNoElements[0] + 1 - rotate], 0);
            }

            index = 2 * j;
            rotate = (j + linearNoElements[1] - 1) % 2;
            if (periodicY_ == 1) {
                addFace(tempElementVector[index + j % 2], 0, tempElementVector[index + 2 * (linearNoElements[1] - 1) * linearNoElements[0] + 1 - rotate], 2 - rotate);
            } else {
                //Bottom boundary face	  
                // 	    cout << "index4=" << index + j % 2 << endl;
                addFace(tempElementVector[index + j % 2], 0, NULL, 0, WALL_BC);
                //Top boundary face
                index = index + 2 * (linearNoElements[1] - 1) * linearNoElements[0] + 1 - rotate;
                // 	    cout << "index5=" << index << endl;
                addFace(tempElementVector[index], 2 - rotate, NULL, 0, WALL_BC);
            }
        }

        //now do the diagonal faces
        for (int j = 0; j < linearNoElements[1]; ++j) {
            for (int i = 0; i < linearNoElements[0]; ++i) {
                index = 2 * (i + j * linearNoElements[0]);
                addFace(tempElementVector[index], 2, tempElementVector[index + 1], 1);
                //diagonal faces are never boudary faces
            }
        }
    }

    void MeshManipulator::triangularCreateFaces3D(std::vector<Base::Element*>& tempElementVector, const std::vector<unsigned int>& linearNoElements) {
        //my ascii-art skill are sadly not good enough to show how the rectangular element has been devided
        unsigned int index;
        unsigned int rotate;
        //counter in z
        for (int k = 0; k < linearNoElements[2]; k++) {
            //counter in y
            for (int j = 0; j < linearNoElements[1]; j++) {
                //counter in x
                for (int i = 0; i + 1 < linearNoElements[0]; i++) {
                    // 	        cout<<"("<<i<<","<<j<<","<<k<<")"<<endl;
                    //hard-coded number of triangles per rectangle, being able to change a constant is not going to make this easier to change
                    index = 5 * (i + j * linearNoElements[0] + k * linearNoElements[0] * linearNoElements[1]);
                    rotate = (i + j + k) % 2;
                    addFace(tempElementVector[index + 2], 3 * rotate,
                            tempElementVector[index + 5], 2 - 2 * rotate);
                    addFace(tempElementVector[index + 1 + 2 * rotate], 1 + 2 * rotate,
                            tempElementVector[index + 6 + 2 * rotate], 2);
                } //end loop over x

                index = 5 * (j * linearNoElements[0] + k * linearNoElements[0] * linearNoElements[1]);
                rotate = (j + k + linearNoElements[0] - 1) % 2;
                if (periodicX_ == 1) {
                    addFace(tempElementVector[index], 2 * ((k + j) % 2),
                            tempElementVector[index + 5 * linearNoElements[0] - 3], 3 * rotate);
                    addFace(tempElementVector[index + 3 - 2 * ((k + j) % 2)], 2,
                            tempElementVector[index + 5 * linearNoElements[0] - 4 + 2 * rotate], 1 + 2 * rotate);
                } else {
                    //Left boundary
                    // 	        cout << "index1=" << index << endl;
                    addFace(tempElementVector[index], 2 * ((k + j) % 2), NULL, 0, WALL_BC);
                    addFace(tempElementVector[index + 3 - 2 * ((k + j) % 2)], 2, NULL, 0, WALL_BC);
                    //Right boundary face
                    index = index + 5 * linearNoElements[0] - 5;
                    //              cout << "index2=" << index << endl;
                    addFace(tempElementVector[index + 2], 3 * rotate, NULL, 0, WALL_BC);
                    addFace(tempElementVector[index + 1 + 2 * rotate], 1 + 2 * rotate, NULL, 0, WALL_BC);
                } // end boundary cases
            } // end loop over y
        } //end loop over z

        //Now do the faces pointing in the y-direction
        //count in z
        for (int k = 0; k < linearNoElements[2]; k++) {
            //counter in x
            for (int j = 0; j < linearNoElements[0]; j++) {
                //counter in y
                for (int i = 0; i + 1 < linearNoElements[1]; i++) {
                    // 		cout<<"("<<i<<","<<j<<","<<k<<")"<<endl;
                    index = 5 * (j + i * linearNoElements[0] + k * linearNoElements[0] * linearNoElements[1]);
                    rotate = (i + j + k) % 2;
                    addFace(tempElementVector[index + 1], 3,
                            tempElementVector[index + 5 * linearNoElements[0] + 2 - 2 * rotate], 2 - rotate);
                    addFace(tempElementVector[index + 3], 1,
                            tempElementVector[index + 5 * linearNoElements[0] + 2 * rotate], 1 + rotate);
                } //end loop over x

                index = 5 * (j + k * linearNoElements[0] * linearNoElements[1]);
                rotate = (j + k) % 2;
                if (periodicY_ == 1) {
                    //cout << "Hmmm this is the problem " << endl;
                    // 		cout << index << " , " << index + 5 * (linearNoElements[1] - 1) * linearNoElements[0] << endl;
                    addFace(tempElementVector[index + 2 * rotate], 1 + rotate,
                            tempElementVector[(index + 5 * (linearNoElements[1] - 1) * linearNoElements[0]) + 1], 3);
                    addFace(tempElementVector[index + 2 - 2 * rotate], 2 - rotate,
                            tempElementVector[(index + 5 * (linearNoElements[1] - 1) * linearNoElements[0]) + 3], 1);
                } else {
                    //front bounday face
                    // 		cout << "index1=" << index << endl;
                    addFace(tempElementVector[index + 2 - 2 * rotate], 2 - rotate, NULL, 0, WALL_BC);
                    addFace(tempElementVector[index + 2 * rotate], 1 + rotate, NULL, 0, WALL_BC);
                    //Back boundary face
                    index = index + 5 * (linearNoElements[1] - 1) * linearNoElements[0];
                    // 		cout << "index2=" << index << endl;
                    addFace(tempElementVector[index + 1], 3, NULL, 0, WALL_BC);
                    addFace(tempElementVector[index + 3], 1, NULL, 0, WALL_BC);
                } // end boundary cases
            } // end loop over y
        } //end loop over z

        //Now do the face in the z-direction
        //count in x direction
        for (int k = 0; k < linearNoElements[0]; k++) {
            //counter in y
            for (int j = 0; j < linearNoElements[1]; j++) {
                //counter in z
                for (int i = 0; i + 1 < linearNoElements[2]; i++) {
                    // 		cout<<"("<<i<<","<<j<<","<<k<<")"<<endl;
                    index = 5 * (k + j * linearNoElements[0] + i * linearNoElements[0] * linearNoElements[1]);
                    rotate = (i + j + k) % 2;
                    addFace(tempElementVector[index + 2 - 2 * rotate], 3 - 3 * rotate,
                            tempElementVector[index + 5 * linearNoElements[0] * linearNoElements[1] + 2 - 2 * rotate], 2 * rotate);
                    addFace(tempElementVector[index + 3], 3 - rotate,
                            tempElementVector[index + 5 * linearNoElements[0] * linearNoElements[1] + 1], 1 + rotate);
                } //end loop over x

                index = 5 * (k + j * linearNoElements[0]);
                rotate = (k + j + linearNoElements[2] - 1) % 2;
                if (periodicZ_ == 1) {
                    // 		cout << "This is the final problem " << endl;
                    // 		cout << index << " , ";
                    // 		cout << index + 5 * (linearNoElements[2] - 1) * linearNoElements[1] * linearNoElements[0] << endl;
                    addFace(tempElementVector[index + 1 + ((j + k) % 2)], 2 - 2 * ((j + k) % 2),
                            tempElementVector[index + 5 * (linearNoElements[2] - 1) * linearNoElements[1] * linearNoElements[0] + 2 + rotate], 3 - rotate);
                    addFace(tempElementVector[index + (j + k) % 2], 2 - (j + k) % 2,
                            tempElementVector[index + 5 * (linearNoElements[2] - 1) * linearNoElements[1] * linearNoElements[0] + 3 - 3 * rotate], 3 - 3 * rotate);
                } else {
                    //bottom boundary
                    // 		cout << "index1=" << index << endl;
                    addFace(tempElementVector[index + 2 * ((j + k) % 2)], 2 - 2 * ((j + k) % 2), NULL, 0, WALL_BC);
                    addFace(tempElementVector[index + 1 ], 2 - ((j + k) % 2), NULL, 0, WALL_BC);
                    //Top boundary face
                    index = index + 5 * (linearNoElements[2] - 1) * linearNoElements[1] * linearNoElements[0];
                    // 		cout << "index2=" << index << endl;
                    addFace(tempElementVector[index + 2 - 2 * rotate], 3 - 3 * rotate, NULL, 0, WALL_BC);
                    addFace(tempElementVector[index + 3], 3 - rotate, NULL, 0, WALL_BC);
                } // end boundary cases
            } // end loop over y
        } //end loop over z
        //Now do the diagonal faces
        //counter in z
        for (int k = 0; k < linearNoElements[2]; k++) {
            //counter in y
            for (int j = 0; j < linearNoElements[1]; j++) {
                //counter in x
                for (int i = 0; i < linearNoElements[0]; i++) {
                    // 		cout<<"("<<i<<","<<j<<","<<k<<")"<<endl;	      
                    index = 5 * (i + j * linearNoElements[0] + k * linearNoElements[0] * linearNoElements[1]);
                    addFace(tempElementVector[index + 4], 0, tempElementVector[index + 2], 1);
                    addFace(tempElementVector[index + 4], 1, tempElementVector[index + 1], 0);
                    addFace(tempElementVector[index + 4], 2, tempElementVector[index], 3);
                    addFace(tempElementVector[index + 4], 3, tempElementVector[index + 3], 0);
                } //end loop over x
            } //end loop over y
        } //end loop over z  
        edgeFactory(tempElementVector);
    }*/

    void
    MeshManipulator::readCentaurMesh(const std::string& filename) {
        //set to correct value in case some other meshmanipulator changed things
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
        if (!centaurFile.is_open()) {
            throw ("Cannot open Centaur meshfile.");
        }

        switch (configData_->dimension_) {
            case 2:
                readCentaurMesh2D(centaurFile);
                break;
            case 3:
                readCentaurMesh3D(centaurFile);
                break;
            default:
                std::cerr << "Centaur mesh reader has not been implemented in this DIMension" << std::endl;
        }

        //Finally close the file
        centaurFile.close();
    }

    void MeshManipulator::readCentaurMesh2D(std::ifstream& centaurFile) {

        auto& elementslist = theMesh_.getElementsList(IteratorType::GLOBAL);
        

        //These are used to check the length of the read lines to check for read errors
        uint32_t sizeOfLine;

        //This first value in the centaur file is the size of each line in the file;
        centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));

        // Version number of the Centaur mesh file 0.1 if using Tito's matlab generator
        float version;
        centaurFile.read(reinterpret_cast<char*> (&version), sizeof (version));
        std::cout << "This read mesh is in Centaur version " << version << " format" << std::endl;

        // Centaur File Type <0 is two DIMensional and >0 is three DIMensional
        int32_t centaurFileType;
        centaurFile.read(reinterpret_cast<char*> (&centaurFileType), sizeof (centaurFileType));



        if (centaurFileType < 0) {
            std::cout << "Reading a two DIMensional centaur mesh" << std::endl;

            //The rest of the first line is junk
            char junk[1024];

            uint32_t checkInt;
            centaurFile.read(&junk[0], sizeOfLine - sizeof (version) - sizeof (centaurFileType));

            //Check the first line was read correctly : each line in centaur start and end with the line size as a check
            centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
            if (checkInt != sizeOfLine) {
                std::cerr << "Error in centaur file " << std::endl;
                return;
            }

            //Start the second line
            centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));

            //Next read the total number of nodes
            uint32_t numberOfNodes;
            centaurFile.read(reinterpret_cast<char*> (&numberOfNodes), sizeof (numberOfNodes));
            std::cout << "File contains " << numberOfNodes << " nodes" << std::endl;

            //Check the second line was read correctly : each line in centaur start and end with the line size as a check
            centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
            if (checkInt != sizeOfLine) {
                std::cerr << "Error in centaur file " << std::endl;
                return;
            }

            //placeholder for vertices until we know where the periodic boundaries are
            std::vector<std::vector<int> > listOfElementsForEachNode(numberOfNodes);

            //Now we will read in all the nodes
            centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
            double nodeCoord[2];
            PointPhysicalT nodeCoordPointFormat(2);
            for (int i = 0; i < numberOfNodes; i++) {
                // Reads the x and y coordinates of each node.
                centaurFile.read(reinterpret_cast<char*> (nodeCoord), sizeof(nodeCoord));
                // pass the node to the nodelist.

                //Covert from *double to hpGEM PointPhysical format
                nodeCoordPointFormat[0] = nodeCoord[0];
                nodeCoordPointFormat[1] = nodeCoord[1];
                theMesh_.addNode(nodeCoordPointFormat);

            }
            //Now check the node line was read correctly : each line in centaur start and end with the line size as a check
            centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
            if (checkInt != sizeOfLine) {
                std::cerr << "Error in centaur file " << std::endl;
                return;
            }



            //Now check how many triangle in the file
            centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
            // Number of triangular elements
            uint32_t numberOfTriangles;
            centaurFile.read(reinterpret_cast<char*> (&numberOfTriangles), sizeof (numberOfTriangles));
            std::cout << "File contains " << numberOfTriangles << " triangle(s)" << std::endl;

            //Check the line was read correctly : each line in centaur start and end with the line size as a check
            centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
            if (checkInt != sizeOfLine) {
                std::cerr << "Error in centaur file " << std::endl;
                return;
            }

            centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
            if (numberOfTriangles > 0) {
                //unsigned int triGlobalNodeIndexes[3];
                std::vector<uint32_t> globalNodeIndexes(3);

                for (int i = 0; i < numberOfTriangles; i++) {

                    centaurFile.read(reinterpret_cast<char*> (globalNodeIndexes.data()), sizeof(uint32_t)*globalNodeIndexes.size());
                    for (int j = 0; j < 3; j++) {
                        globalNodeIndexes[j] = globalNodeIndexes[j] - 1;
                    }

                    
                    int id = addElement(globalNodeIndexes)->getID();

                    for(uint32_t j:globalNodeIndexes)
                    {
                        listOfElementsForEachNode[j].push_back(id);
                    }
                }


            }
            centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
            if (checkInt != sizeOfLine) {
                std::cerr << "Error in centaur file " << std::endl;
                return;
            }

            //Now check the number of quaduratiles in the file
            centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
            uint32_t numberOfQuads;
            centaurFile.read(reinterpret_cast<char*> (&numberOfQuads), sizeof (numberOfQuads));
            std::cout << "File contains " << numberOfQuads << " quaduratile(s)" << std::endl;

            centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
            if (checkInt != sizeOfLine) {
                std::cerr << "Error in centaur file " << std::endl;
                return;
            }

            //Now read the quaduritles in
            centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
            if (numberOfQuads > 0) {
                //unsigned int quadGlobalNodeIndexes[4];
                uint32_t temp;
                std::vector<uint32_t> globalNodeIndexes(4);
                for (int i = 0; i < numberOfQuads; i++) {
                    //Reading the vertex indices of each quadrilateral.
                    centaurFile.read(reinterpret_cast<char*> (globalNodeIndexes.data()), sizeof (uint32_t) * globalNodeIndexes.size());

                    // renumbering of the vertices to match the ordering assumed by
                    // hpGem:

                    temp = globalNodeIndexes[2];
                    globalNodeIndexes[2] = globalNodeIndexes[3];
                    globalNodeIndexes[3] = temp;

                    
                    // renumber them from 1..N to 0..N-1.
                    for (int j = 0; j < 4; j++)
                        globalNodeIndexes[j] = globalNodeIndexes[j] - 1;

                    int id = addElement(globalNodeIndexes)->getID();
                    
                    for(uint32_t j:globalNodeIndexes)
                    {
                        listOfElementsForEachNode[j].push_back(id);
                    }
                }

            }
            //Check the line was read correctly : each line in centaur start and end with the line size as a check
            centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
            if (checkInt != sizeOfLine) {
                std::cerr << "Error in centaur file " << std::endl;
                return;
            }

            //now read boundary data
            //nodes first
            
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            
            uint32_t numberOfBoundaryNodes;
            centaurFile.read(reinterpret_cast<char*>(&numberOfBoundaryNodes), sizeof(numberOfBoundaryNodes));
            
            centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
            if (checkInt != sizeOfLine) {
                std::cerr << "Error in centaur file " << std::endl;
                return;
            }
            
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            
            std::vector<std::vector<uint32_t> > boundaryNodesForEachGroup;
            uint32_t boundaryNodeInformation[2];
            
            for(uint_fast32_t i=0; i<numberOfBoundaryNodes; ++i)
            {
                centaurFile.read(reinterpret_cast<char*> (boundaryNodeInformation), sizeof(boundaryNodeInformation));
                if(boundaryNodesForEachGroup.size()<boundaryNodeInformation[1]+1)
                {
                    boundaryNodesForEachGroup.resize(boundaryNodeInformation[1]+1);
                }
                boundaryNodesForEachGroup[boundaryNodeInformation[1]].push_back(boundaryNodeInformation[0]-1);
            }
            
            centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
            if (checkInt != sizeOfLine) {
                std::cerr << "Error in centaur file " << std::endl;
                return;
            }
            
            //then faces
            
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            
            uint32_t numberOfBoundaryFaces;
            centaurFile.read(reinterpret_cast<char*>(&numberOfBoundaryFaces), sizeof(numberOfBoundaryFaces));
            
            centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
            if (checkInt != sizeOfLine) {
                std::cerr << "Error in centaur file " << std::endl;
                return;
            }
            
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            
            //read the half-faces and store them until we get to the boundary information
            //uint32_t nodesForEachFace[2];
            std::vector<HalfFaceDescription> boundaryFaces(numberOfBoundaryFaces);
            
            for(uint_fast32_t i=0;i<numberOfBoundaryFaces;++i)
            {
                boundaryFaces[i].nodeList.resize(2);
                centaurFile.read(reinterpret_cast<char*>(boundaryFaces[i].nodeList.data()), sizeof(uint32_t)*2);
                
                boundaryFaces[i].nodeList[0]-=1;
                boundaryFaces[i].nodeList[1]-=1;
                
                std::vector<int> candidateElements;
                std::vector<int>& leftNodeElements = listOfElementsForEachNode[boundaryFaces[i].nodeList[0]];
                std::vector<int>& rightNodeElements =  listOfElementsForEachNode[boundaryFaces[i].nodeList[1]];
                
                std::set_intersection(leftNodeElements.begin(),leftNodeElements.end(),rightNodeElements.begin(),rightNodeElements.end(),std::back_inserter(candidateElements));
                
                //boundary faces *should* border only one element
                if(candidateElements.size()>1)
                {
                    std::cerr << "candidate boundary face lies at two or more elements" << std::endl;
                    return;
                }
                boundaryFaces[i].elementNum = candidateElements[0];
                
                Element* current = elementslist[candidateElements[0]];
                std::vector<unsigned int> faceVertices(2);
                
                for(size_t j=0;j<current->getNrOfFaces();++j)
                {
                    current->getPhysicalGeometry()->getGlobalFaceNodeIndices(j,faceVertices);
                    if((faceVertices[0]==boundaryFaces[i].nodeList[0] || faceVertices[0]==boundaryFaces[i].nodeList[1]) &&
                       (faceVertices[1]==boundaryFaces[i].nodeList[0] || faceVertices[1]==boundaryFaces[i].nodeList[1]))
                    {
                        boundaryFaces[i].localFaceIndex=j;
                    }
                }
                
                candidateElements.clear();
            }
            
            centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
            if (checkInt != sizeOfLine) {
                std::cerr << "Error in centaur file " << std::endl;
                return;
            }
            
            //now read a list of boundary segments (that will link the faces to boundary conditions later on)
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            
            std::vector<uint32_t> faceToSegment(numberOfBoundaryFaces);
            centaurFile.read(reinterpret_cast<char*> (faceToSegment.data()), sizeof (uint32_t) * numberOfBoundaryFaces);
            
            centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
            if (checkInt != sizeOfLine) {
                std::cerr << "Error in centaur file " << std::endl;
                return;
            }
            
            //now couple the segments to boundary groups
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            
            uint32_t numberOfSegments;
            centaurFile.read(reinterpret_cast<char*>(&numberOfSegments), sizeof(numberOfSegments));
            
            centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
            if (checkInt != sizeOfLine) {
                std::cerr << "Error in centaur file " << std::endl;
                return;
            }
            
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            
            std::vector<uint32_t> segmentToGroup(numberOfSegments);
            centaurFile.read(reinterpret_cast<char*>(segmentToGroup.data()), sizeof(uint32_t) * numberOfSegments);
            
            centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
            if (checkInt != sizeOfLine) {
                std::cerr << "Error in centaur file " << std::endl;
                return;
            }
            
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
            
            uint32_t numberOfGroups;
            centaurFile.read(reinterpret_cast<char*>(numberOfGroups), sizeof(uint32_t));
            
            centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
            if (checkInt != sizeOfLine) {
                std::cerr << "Error in centaur file " << std::endl;
                return;
            }
            
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            
            std::vector<uint32_t> groupBCType(numberOfGroups);
            centaurFile.read(reinterpret_cast<char*>(groupBCType.data()), sizeof(uint32_t) * numberOfGroups);
            
            centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
            if (checkInt != sizeOfLine) {
                std::cerr << "Error in centaur file " << std::endl;
                return;
            }
            
            for(uint_fast32_t i=0;i<numberOfBoundaryFaces; ++i)
            {
                uint32_t boundaryType = groupBCType[segmentToGroup[faceToSegment[i]]];
                if(boundaryType < 1001) {
                    std::cout << "Viscous Wall boundary for face " << i << " assigned as WALL_BC" << std::endl;
                    addFace(elementslist[boundaryFaces[i].elementNum],boundaryFaces[i].localFaceIndex,NULL,0,Geometry::FaceType::WALL_BC);
                } else if (boundaryType < 2001) {
                    std::cout << "Inviscid Wall boundary for face " << i << " assigned as WALL_BC" << std::endl;
                    addFace(elementslist[boundaryFaces[i].elementNum],boundaryFaces[i].localFaceIndex,NULL,0,Geometry::FaceType::WALL_BC);
                } else if (boundaryType < 3001) {
                    std::cout << "symmetry plane boundary for face " << i << " assigned as WALL_BC" << std::endl;
                    addFace(elementslist[boundaryFaces[i].elementNum],boundaryFaces[i].localFaceIndex,NULL,0,Geometry::FaceType::WALL_BC);
                } else if (boundaryType < 4001) {
                    std::cout << "inlet pipe boundary for face " << i << " assigned as OPEN_BC" << std::endl;
                    addFace(elementslist[boundaryFaces[i].elementNum],boundaryFaces[i].localFaceIndex,NULL,0,Geometry::FaceType::OPEN_BC);
                } else if (boundaryType < 5001) {
                    std::cout << "outlet pipe boundary for face " << i << " assigned as OPEN_BC" << std::endl;
                    addFace(elementslist[boundaryFaces[i].elementNum],boundaryFaces[i].localFaceIndex,NULL,0,Geometry::FaceType::OPEN_BC);
                } else if (boundaryType < 6001) {
                    std::cout << "farfield boundary for face " << i << " assigned as OPEN_BC" << std::endl;
                    addFace(elementslist[boundaryFaces[i].elementNum],boundaryFaces[i].localFaceIndex,NULL,0,Geometry::FaceType::OPEN_BC);
                } else if (boundaryType < 7001) {
                    std::cout << "periodic boundary for face " << i << " ignored for being internal; node connections will be assigned later" << std::endl;
                } else if (boundaryType < 8001) {
                    std::cout << "shadow boundary for face " << i << " ignored for being internal; node connections will be assigned later" << std::endl;
                } else if (boundaryType < 8501) {
                    std::cout << "interface boundary for face " << i << " ignored for being internal" << std::endl;
                } else if (boundaryType < 9001) {
                    std::cout << "wake boundary for face " << i << " assigned as OPEN_BC" << std::endl;
                    addFace(elementslist[boundaryFaces[i].elementNum],boundaryFaces[i].localFaceIndex,NULL,0,Geometry::FaceType::OPEN_BC);
                } else if (boundaryType < 10001) {
                    std::cout << "moving wall boundary for face " << i << " assigned as WALL_BC" << std::endl;
                    addFace(elementslist[boundaryFaces[i].elementNum],boundaryFaces[i].localFaceIndex,NULL,0,Geometry::FaceType::WALL_BC);
                } else {
                    std::cout << "alternative boundary condition for face " << i << " assigned as WALL_BC" << std::endl;
                    addFace(elementslist[boundaryFaces[i].elementNum],boundaryFaces[i].localFaceIndex,NULL,0,Geometry::FaceType::WALL_BC);
                }
            }
            
            //I dont care about the names of the groups, just skip them
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            
            centaurFile.ignore(80*numberOfGroups);
            
            centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
            if (checkInt != sizeOfLine) {
                std::cerr << "Error in centaur file " << std::endl;
                return;
            }
            
            //now read periodic information and link corresponding vertices
            centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
            
            uint32_t numberOfTransforms;
            centaurFile.read(reinterpret_cast<char*> (&numberOfTransforms), sizeof (numberOfTransforms));
            
            centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
            if (checkInt != sizeOfLine) {
                std::cerr << "Error in centaur file " << std::endl;
                return;
            }
            
            for(uint_fast32_t i=0; i<numberOfTransforms; ++i)
            {
                //I dont care about the actual coordinate transformations, just skip them
                centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));

                centaurFile.ignore(9*sizeof(uint32_t));

                centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
                if (checkInt != sizeOfLine) {
                    std::cerr << "Error in centaur file " << std::endl;
                    return;
                }
                centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));

                centaurFile.ignore(9*sizeof(uint));

                centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
                if (checkInt != sizeOfLine) {
                    std::cerr << "Error in centaur file " << std::endl;
                    return;
                }
                
                uint32_t numberOfNodePairs;
                centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));

                centaurFile.read(reinterpret_cast<char*> (&numberOfNodePairs), sizeof (numberOfNodePairs));

                centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
                if (checkInt != sizeOfLine) {
                    std::cerr << "Error in centaur file " << std::endl;
                    return;
                }
                
                centaurFile.read(reinterpret_cast<char*>(&sizeOfLine), sizeof(sizeOfLine));
                
                uint32_t nodePair[2];
                std::vector<int> combine;
                for(uint_fast32_t j=0;j<numberOfNodePairs;++j)
                {
                    centaurFile.read(reinterpret_cast<char*> (nodePair), sizeof (nodePair));
                    auto& firstList = listOfElementsForEachNode[nodePair[0]];
                    auto& secondList = listOfElementsForEachNode[nodePair[1]];
                    std::set_union(firstList.begin(),firstList.end(),secondList.begin(),secondList.end(),std::back_inserter(combine));
                    firstList=combine;
                    secondList=combine;
                    combine.clear();
                }
                
                centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
                if (checkInt != sizeOfLine) {
                    std::cerr << "Error in centaur file " << std::endl;
                    return;
                }
                
            }
            
            //new centaur files have zones, but I'm not interested in them
            
            //now we know periodicity information, construct the vertices
            //remember that listOfElementsForEachNode will in general contain duplicates
            
            auto& listOfVertices = theMesh_.getVerticesList(IteratorType::GLOBAL);
            
            bool addedNewVertex(false);
            
            for(size_t i=0;i<listOfElementsForEachNode.size();++i)
            {
                for(size_t j=0;j<listOfElementsForEachNode[i].size();++j)
                {
                    Element* current = elementslist[listOfElementsForEachNode[i][j]];
                    for(size_t k=0;k<current->getNrOfNodes();++k)
                    {
                        //if we did not jet deal with this node and it is the correct one
                        if(current->getNode(k)==nullptr && current->getPhysicalGeometry()->getNodeIndex(k)==i)
                        {
                            if(!addedNewVertex)
                            {
                                addVertex();
                                addedNewVertex=true;
                            }
                            listOfVertices.back()->addElement(current,k);
                        }
                    }
                }
                addedNewVertex=false;
            }
            
            //construct the rest of the faces
            faceFactory();

        } else {
            std::cerr << "Incorrect mesh file. This mesh appears to contain three DIMensional data" << std::endl;

        }
    }

    void MeshManipulator::readCentaurMesh3D(std::ifstream& centaurFile) {

        auto& elementsList = theMesh_.getElementsList(IteratorType::GLOBAL);

        //These are used to check the length of the read lines to check for read errors
        uint32_t sizeOfLine;

        //This first value in the centaur file is the size of each line in the file;
        centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));

        // Version number of the Centaur mesh file 0.1 if using Tito's matlab generator
        float version;
        centaurFile.read(reinterpret_cast<char*> (&version), sizeof (version));
        std::cout << "This read mesh is in Centaur version " << version << " format" << std::endl;

        // Centaur File Type <0 is two DIMensional and >0 is three DIMensional
        uint32_t centaurFileType;
        centaurFile.read(reinterpret_cast<char*> (&centaurFileType), sizeof (centaurFileType));



        if (centaurFileType > 0) {
            std::cout << "Reading a three DIMensional centaur mesh" << std::endl;

            //The rest of the first line is junk
            char junk[1024];

            uint32_t checkInt;
            centaurFile.read(&junk[0], sizeOfLine - sizeof (version) - sizeof (centaurFileType));

            //Check the first line was read correctly : each line in centaur start and end with the line size as a check
            centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
            if (checkInt != sizeOfLine) {
                std::cerr << "Error in centaur file " << std::endl;
                return;
            }

            //Start the second line
            centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));

            //Next read the total number of nodes
            uint32_t numberOfNodes;
            centaurFile.read(reinterpret_cast<char*> (&numberOfNodes), sizeof (numberOfNodes));
            std::cout << "File contains " << numberOfNodes << " nodes" << std::endl;

            //new centaur versions support splitting this list over multiple lines
            uint32_t numberOfNodesPerLine(numberOfNodes);
            if (centaurFileType > 4) {
                centaurFile.read(reinterpret_cast<char*> (&numberOfNodesPerLine), sizeof (numberOfNodesPerLine));
                std::cout << "One line in the file contains at most " << numberOfNodesPerLine << " nodes" << std::endl;
            }

            //Check the second line was read correctly : each line in centaur start and end with the line size as a check
            centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
            if (checkInt != sizeOfLine) {
                std::cerr << "Error in centaur file " << std::endl;
                return;
            }

            //Now we will read in all the nodes
            centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
            double nodeCoord[3];
            PointPhysicalT nodeCoordPointFormat(3);
            for (int i = 0; i < numberOfNodes; i++) {
                if (i > 0 && i % numberOfNodesPerLine == 0) {
                    //If all the nodes on a line are read end the line and start a new one
                    centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
                    if (checkInt != sizeOfLine) {
                        std::cerr << "Error in centaur file " << std::endl;
                        return;
                    }
                    centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
                }

                // Reads the x, y and z coordinates of each node.
                centaurFile.read(reinterpret_cast<char*> (nodeCoord), sizeof (nodeCoord));
                // pass the node to the nodelist.

                //Covert from *double to hpGEM PointPhysical format
                nodeCoordPointFormat[0] = nodeCoord[0];
                nodeCoordPointFormat[1] = nodeCoord[1];
                nodeCoordPointFormat[2] = nodeCoord[2];
                //std::cout << nodeCoordPointFormat << std::endl;
                theMesh_.addNode(nodeCoordPointFormat);

            }
            //Now check the node line was read correctly : each line in centaur start and end with the line size as a check
            centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
            if (checkInt != sizeOfLine) {
                std::cerr << "Error in centaur file " << std::endl;
                return;
            }

            //Keep track of the node->Element connectivity to ease face creation later on
            std::vector<std::vector<int> > listOfElementsForEachNode(numberOfNodes);
            std::vector<Element* > tempElementVector;

            //file version 1 has no lines about hexahedra
            if (centaurFileType > 1) {
                //Now check how many hexahedra in the file
                centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
                // Number of hexahedral elements
                uint32_t numberOfHexahedra;
                centaurFile.read(reinterpret_cast<char*> (&numberOfHexahedra), sizeof (numberOfHexahedra));
                std::cout << "File contains " << numberOfHexahedra << " hexahedron(s)" << std::endl;

                //new centaur versions support splitting this list over multiple lines
                uint32_t numberOfHexahedraPerLine(numberOfHexahedra);
                if (centaurFileType > 4) {
                    centaurFile.read(reinterpret_cast<char*> (&numberOfHexahedraPerLine), sizeof (numberOfHexahedraPerLine));
                    std::cout << "One line in the file contains at most " << numberOfHexahedraPerLine << " hexahedra" << std::endl;
                }

                centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
                if (checkInt != sizeOfLine) {
                    std::cerr << "Error in centaur file " << std::endl;
                    return;
                }


                //uint32_t hexahedralGlobalNodeIndexes[8];
                uint32_t temp;
                std::vector<uint32_t> globalNodeIndexes(8);
                centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
                for (int i = 0; i < numberOfHexahedra; i++) {
                    if (i > 0 && i % numberOfHexahedraPerLine == 0) {
                        //If all the nodes on a line are read end the line and start a new one
                        centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
                        if (checkInt != sizeOfLine) {
                            std::cerr << "Error in centaur file " << std::endl;
                            return;
                        }
                        centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
                    }

                    //Reading the vertex indices of each hexahedron.
                    centaurFile.read(reinterpret_cast<char*> (globalNodeIndexes.data()), sizeof (uint32_t) * globalNodeIndexes.size());

                    // renumbering of the vertices to match the ordering assumed by
                    // hpGem: (based on the numbering in hpGEM 1)

                    temp = globalNodeIndexes[2];
                    globalNodeIndexes[2] = globalNodeIndexes[3];
                    globalNodeIndexes[3] = temp;

                    temp = globalNodeIndexes[6];
                    globalNodeIndexes[6] = globalNodeIndexes[7];
                    globalNodeIndexes[7] = temp;

                    // renumber them from 1..N to 0..N-1.
                    for (int j = 0; j < 8; j++)
                        globalNodeIndexes[j] = globalNodeIndexes[j] - 1;

                    Base::Element* newElement = addElement(globalNodeIndexes);
                    tempElementVector.push_back(newElement);

                    for (int j = 0; j < 8; ++j) {
                        listOfElementsForEachNode[globalNodeIndexes[j] ].push_back(newElement->getID());
                    }
                }
                centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
                if (checkInt != sizeOfLine) {
                    std::cerr << "Error in centaur file " << std::endl;
                    return;
                }
            }


            //Now check how many triangular prisms in the file
            centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
            // Number of prismatic elements
            uint32_t numberOfPrisms;
            centaurFile.read(reinterpret_cast<char*> (&numberOfPrisms), sizeof (numberOfPrisms));
            std::cout << "File contains " << numberOfPrisms << " triangular prism(s)" << std::endl;

            //new centaur versions support splitting this list over multiple lines
            uint32_t numberOfPrismsPerLine(numberOfPrisms);
            if (centaurFileType > 4) {
                centaurFile.read(reinterpret_cast<char*> (&numberOfPrismsPerLine), sizeof (numberOfPrismsPerLine));
                std::cout << "One line in the file contains at most " << numberOfPrismsPerLine << " prisms" << std::endl;
            }

            centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
            if (checkInt != sizeOfLine) {
                std::cerr << "Error in centaur file " << std::endl;
                return;
            }


            //unsigned int prismaticGlobalNodeIndexes[6];
            std::vector<uint32_t> globalNodeIndexes(6);
            centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
            for (int i = 0; i < numberOfPrisms; i++) {
                if (i > 0 && i % numberOfPrismsPerLine == 0) {
                    //If all the nodes on a line are read end the line and start a new one
                    centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
                    if (checkInt != sizeOfLine) {
                        std::cerr << "Error in centaur file " << std::endl;
                        return;
                    }
                    centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
                }

                //Reading the vertex indices of each hexahedron.
                centaurFile.read(reinterpret_cast<char*> (globalNodeIndexes.data()), sizeof (uint32_t) * globalNodeIndexes.size());

                // renumber them from 1..N to 0..N-1.
                for (int j = 0; j < 6; j++)
                    globalNodeIndexes[j] = globalNodeIndexes[j] - 1;

                Base::Element* newElement = addElement(globalNodeIndexes);
                tempElementVector.push_back(newElement);

                for (int j = 0; j < 6; ++j) {
                    listOfElementsForEachNode[globalNodeIndexes[j] ].push_back(newElement->getID());
                }
            }
            centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
            if (checkInt != sizeOfLine) {
                std::cerr << "Error in centaur file " << std::endl;
                return;
            }

            //file version 1 has no lines about pyramids
            if (centaurFileType > 1) {
                //Now check how many pyramids in the file
                centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
                // Number of pyramids elements
                uint32_t numberOfPyraminds;
                centaurFile.read(reinterpret_cast<char*> (&numberOfPyraminds), sizeof (numberOfPyraminds));
                std::cout << "File contains " << numberOfPyraminds << " pyramid(s)" << std::endl;

                //new centaur versions support splitting this list over multiple lines
                uint32_t numberOfPyramindsPerLine(numberOfPyraminds);
                if (centaurFileType > 4) {
                    centaurFile.read(reinterpret_cast<char*> (&numberOfPyramindsPerLine), sizeof (numberOfPyramindsPerLine));
                    std::cout << "One line in the file contains at most " << numberOfPyramindsPerLine << " pyramids" << std::endl;
                }

                centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
                if (checkInt != sizeOfLine) {
                    std::cerr << "Error in centaur file " << std::endl;
                    return;
                }


                //unsigned int pyramidGlobalNodeIndexes[5];
                uint32_t temp;
                globalNodeIndexes.resize(5);
                centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
                for (int i = 0; i < numberOfPyraminds; i++) {
                    if (i > 0 && i % numberOfPyramindsPerLine == 0) {
                        //If all the nodes on a line are read end the line and start a new one
                        centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
                        if (checkInt != sizeOfLine) {
                            std::cerr << "Error in centaur file " << std::endl;
                            return;
                        }
                        centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
                    }

                    //Reading the vertex indices of each pyramid.
                    centaurFile.read(reinterpret_cast<char*> (globalNodeIndexes.data()), sizeof (uint32_t) * globalNodeIndexes.size());

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


                    if (centaurFileType < 4) {
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
                    for (int j = 0; j < 5; j++)
                        globalNodeIndexes[j] = globalNodeIndexes[j] - 1;

                    Base::Element* newElement = addElement(globalNodeIndexes);
                    tempElementVector.push_back(newElement);

                    for (int j = 0; j < 5; ++j) {
                        listOfElementsForEachNode[globalNodeIndexes[j] ].push_back(newElement->getID());
                    }
                }
                centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
                if (checkInt != sizeOfLine) {
                    std::cerr << "Error in centaur file " << std::endl;
                    return;
                }
            }

            //Now check how many tetrahedra in the file
            centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
            // Number of tetrahedral elements
            uint32_t numberOfTetrahedra;
            centaurFile.read(reinterpret_cast<char*> (&numberOfTetrahedra), sizeof (numberOfTetrahedra));
            std::cout << "File contains " << numberOfTetrahedra << " tetrahedron(s)" << std::endl;

            //new centaur versions support splitting this list over multiple lines
            uint32_t numberOfTetrahedraPerLine(numberOfTetrahedra);
            if (centaurFileType > 4) {
                centaurFile.read(reinterpret_cast<char*> (&numberOfTetrahedraPerLine), sizeof (numberOfTetrahedraPerLine));
                std::cout << "One line in the file contains at most " << numberOfTetrahedraPerLine << " tetrahedra" << std::endl;
            }

            centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
            if (checkInt != sizeOfLine) {
                std::cerr << "Error in centaur file " << std::endl;
                return;
            }


            //unsigned int tetrahedralGlobalNodeIndexes[4];
            globalNodeIndexes.resize(4);
            centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
            for (int i = 0; i < numberOfTetrahedra; i++) {
                if (i > 0 && i % numberOfTetrahedraPerLine == 0) {
                    //If all the tetrahedra on a line are read end the line and start a new one
                    centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
                    if (checkInt != sizeOfLine) {
                        std::cerr << "Error in centaur file " << std::endl;
                        return;
                    }
                    centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
                }

                //Reading the vertex indices of each tetrahedron.
                centaurFile.read(reinterpret_cast<char*> (globalNodeIndexes.data()), sizeof (uint32_t) * globalNodeIndexes.size());

                // renumber them from 1..N to 0..N-1.
                for (int j = 0; j < 4; j++)
                    globalNodeIndexes[j] = globalNodeIndexes[j] - 1;

                Base::Element* newElement = addElement(globalNodeIndexes);
                tempElementVector.push_back(newElement);

                for (int j = 0; j < 4; ++j) {
                    listOfElementsForEachNode[globalNodeIndexes[j] ].push_back(newElement->getID());
                }
            }
            centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
            if (checkInt != sizeOfLine) {
                std::cerr << "Error in centaur file " << std::endl;
                return;
            }


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



            centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));

            //number of boundary faces
            uint32_t numberOfBoundaryFaces;
            centaurFile.read(reinterpret_cast<char*> (&numberOfBoundaryFaces), sizeof (numberOfBoundaryFaces));
            std::cout << "File contains " << numberOfBoundaryFaces << " boundaryFace(s)" << std::endl;

            uint32_t boundaryFacesPerLine(numberOfBoundaryFaces);
            if (centaurFileType > 4) {
                centaurFile.read(reinterpret_cast<char*> (&boundaryFacesPerLine), sizeof (boundaryFacesPerLine));
                std::cout << "One line in the file contains at most " << boundaryFacesPerLine << " tetrahedra" << std::endl;
            }

            HalfFaceDescription *boundarFaces = new HalfFaceDescription[numberOfBoundaryFaces];

            std::vector<std::vector<uint32_t> > facesForEachCentaurPanel(0);
            //old centaur files use 7 entries to describe the face
            uint32_t nodalDescriptionOfTheFace[centaurFileType > 3 ? 8 : 7];
            uint32_t panelNumber, ijunk;
            centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
            if (checkInt != sizeOfLine) {
                std::cerr << "Error in centaur file " << std::endl;
                return;
            }

            //first read the information about the faces
            centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
            for (int i = 0; i < numberOfBoundaryFaces; ++i) {

                if (i > 0 && i % numberOfBoundaryFaces == 0) {
                    //If all the tetrahedra on a line are read end the line and start a new one
                    centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
                    if (checkInt != sizeOfLine) {
                        std::cerr << "Error in centaur file " << std::endl;
                        return;
                    }
                    centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
                }

                centaurFile.read(reinterpret_cast<char*> (&nodalDescriptionOfTheFace[0]), sizeof (nodalDescriptionOfTheFace));

                boundarFaces[i].nodeList.resize(3);
                boundarFaces[i].nodeList[0] = nodalDescriptionOfTheFace[0];
                boundarFaces[i].nodeList[1] = nodalDescriptionOfTheFace[1];
                boundarFaces[i].nodeList[2] = nodalDescriptionOfTheFace[2];
                
                //three nodes will uniquely determine the face
                auto& firstNodeList = listOfElementsForEachNode[nodalDescriptionOfTheFace[0]]; 
                auto& secondNodeList = listOfElementsForEachNode[nodalDescriptionOfTheFace[1]]; 
                auto& thirdNodeList = listOfElementsForEachNode[nodalDescriptionOfTheFace[2]]; 
                std::vector<unsigned int> temp, candidates,nodes(boundarFaces->nodeList), intersect;
                std::sort(nodes.begin(),nodes.end());
                std::set_intersection(firstNodeList.begin(),firstNodeList.end(), secondNodeList.begin(), secondNodeList.end(), std::back_inserter(temp));
                std::set_intersection(temp.begin(),temp.end(), thirdNodeList.begin(), thirdNodeList.end(), std::back_inserter(candidates));
                
                if(candidates.size() > 1)
                {
                    std::cerr << "candidate boundary face lies at two or more elements" << std::endl;
                    return;
                }
                
                boundarFaces[i].elementNum = candidates[0];
                
                Element* current = elementsList[candidates[0]];
                
                for(size_t j=0;j<current->getNrOfFaces();++j)
                {
                    current->getPhysicalGeometry()->getGlobalFaceNodeIndices(j,temp);
                    std::sort(temp.begin(),temp.end());
                    std::set_intersection(temp.begin(),temp.end(),nodes.begin(),nodes.end(),std::back_inserter(intersect));
                    if(intersect.size()==3)
                    {
                        boundarFaces->localFaceIndex=j;
                        if(typeid(current->getReferenceGeometry()->getCodim1ReferenceGeometry(j))==typeid(Geometry::ReferenceSquare)) {
                            if(centaurFileType>3)
                            {
                                boundarFaces[i].nodeList.push_back(boundarFaces[i].nodeList[2]);
                                boundarFaces[i].nodeList[2]=nodalDescriptionOfTheFace[3];
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
                if (centaurFileType < 4) {
                    centaurFile.read(reinterpret_cast<char*> (&panelNumber), sizeof (panelNumber));
                    if (centaurFileType > 1) {
                        centaurFile.read(reinterpret_cast<char*> (&ijunk), sizeof (ijunk));
                    }
                    if (panelNumber > facesForEachCentaurPanel.size()) {
                        facesForEachCentaurPanel.resize(panelNumber);
                    }
                    facesForEachCentaurPanel[panelNumber].push_back(i);
                }
            }

            centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
            if (checkInt != sizeOfLine) {
                std::cerr << "Error in centaur file " << std::endl;
                return;
            }

            //modern centaur file version store panel information separately
            if (centaurFileType > 3) {
                centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
                for (int i = 0; i < numberOfBoundaryFaces; ++i) {
                    if (i > 0 && i % numberOfBoundaryFaces == 0) {
                        //If all the tetrahedra on a line are read end the line and start a new one
                        centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
                        if (checkInt != sizeOfLine) {
                            std::cerr << "Error in centaur file " << std::endl;
                            return;
                        }
                        centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
                    }
                    centaurFile.read(reinterpret_cast<char*> (&panelNumber), sizeof (panelNumber));
                    if (panelNumber > facesForEachCentaurPanel.size()) {
                        facesForEachCentaurPanel.resize(panelNumber);
                    }
                    facesForEachCentaurPanel[panelNumber - 1].push_back(i);
                }
                centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
                if (checkInt != sizeOfLine) {
                    std::cerr << "Error in centaur file " << std::endl;
                    return;
                }
            }

            //put the centaur panels in their boudary group
            std::vector<std::vector<int> > facesForEachBoundaryGroup(0);
            uint32_t groupOfPanelNumber;

            //this bit of information is a little late
            uint32_t numberOfPanels;

            centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
            centaurFile.read(reinterpret_cast<char*> (&numberOfPanels), sizeof (numberOfPanels));
            assert(numberOfPanels == facesForEachCentaurPanel.size());
            centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
            if (checkInt != sizeOfLine) {
                std::cerr << "Error in centaur file " << std::endl;
                return;
            }

            //then read the panel to group connections
            centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
            for (int i = 0; i < numberOfPanels; ++i) {
                centaurFile.read(reinterpret_cast<char*> (&groupOfPanelNumber), sizeof (groupOfPanelNumber));
                if (groupOfPanelNumber > facesForEachBoundaryGroup.size()) {
                    facesForEachBoundaryGroup.resize(groupOfPanelNumber);
                }
                for (int j = 0; j < facesForEachCentaurPanel[i].size(); ++j) {
                    facesForEachBoundaryGroup[groupOfPanelNumber - 1].push_back(facesForEachCentaurPanel[i][j]);
                }
            }
            centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
            if (checkInt != sizeOfLine) {
                std::cerr << "Error in centaur file " << std::endl;
                return;
            }

            uint32_t centaurBCType;
            char nameOfBoundaryCondition[80];

            //this bit of information is again a little late
            uint32_t numberOfBoundaryGroups;

            centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
            centaurFile.read(reinterpret_cast<char*> (&numberOfBoundaryGroups), sizeof (numberOfBoundaryGroups));
            assert(numberOfBoundaryGroups == facesForEachBoundaryGroup.size());
            centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
            if (checkInt != sizeOfLine) {
                std::cerr << "Error in centaur file " << std::endl;
                return;
            }

            //now set the boundary conditions for each group
            centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
            for (int i = 0; i < numberOfBoundaryGroups; ++i) {
                centaurFile.read(reinterpret_cast<char*> (&centaurBCType), sizeof (centaurBCType));
                if (centaurBCType < 1001) {
                    std::cout << "Viscous Wall boundary for group " << i << " assigned as WALL_BC" << std::endl;
                    for (int j = 0; j < facesForEachBoundaryGroup[i].size(); ++j) {
                        //big assumption on the nature of elementIDs here...
                        addFace(tempElementVector[boundarFaces[facesForEachBoundaryGroup[i][j]].elementNum], boundarFaces[facesForEachBoundaryGroup[i][j]].localFaceIndex, NULL, 0, Geometry::FaceType::WALL_BC);
                    }
                } else if (centaurBCType < 2001) {
                    std::cout << "Inviscid Wall boundary for group " << i << " assigned as WALL_BC" << std::endl;
                    for (int j = 0; j < facesForEachBoundaryGroup[i].size(); ++j) {
                        //big assumption on the nature of elementIDs here...
                        addFace(tempElementVector[boundarFaces[facesForEachBoundaryGroup[i][j]].elementNum], boundarFaces[facesForEachBoundaryGroup[i][j]].localFaceIndex, NULL, 0, Geometry::FaceType::WALL_BC);
                    }
                } else if (centaurBCType < 3001) {
                    std::cout << "symmetry plane boundary for group " << i << " assigned as WALL_BC" << std::endl;
                    for (int j = 0; j < facesForEachBoundaryGroup[i].size(); ++j) {
                        //big assumption on the nature of elementIDs here...
                        addFace(tempElementVector[boundarFaces[facesForEachBoundaryGroup[i][j]].elementNum], boundarFaces[facesForEachBoundaryGroup[i][j]].localFaceIndex, NULL, 0, Geometry::FaceType::WALL_BC);
                    }
                } else if (centaurBCType < 4001) {
                    std::cout << "inlet pipe boundary for group " << i << " assigned as OPEN_BC" << std::endl;
                    for (int j = 0; j < facesForEachBoundaryGroup[i].size(); ++j) {
                        //big assumption on the nature of elementIDs here...
                        addFace(tempElementVector[boundarFaces[facesForEachBoundaryGroup[i][j]].elementNum], boundarFaces[facesForEachBoundaryGroup[i][j]].localFaceIndex, NULL, 0, Geometry::FaceType::OPEN_BC);
                    }
                } else if (centaurBCType < 5001) {
                    std::cout << "outlet pipe boundary for group " << i << " assigned as OPEN_BC" << std::endl;
                    for (int j = 0; j < facesForEachBoundaryGroup[i].size(); ++j) {
                        //big assumption on the nature of elementIDs here...
                        addFace(tempElementVector[boundarFaces[facesForEachBoundaryGroup[i][j]].elementNum], boundarFaces[facesForEachBoundaryGroup[i][j]].localFaceIndex, NULL, 0, Geometry::FaceType::OPEN_BC);
                    }
                } else if (centaurBCType < 6001) {
                    std::cout << "farfield boundary for group " << i << " assigned as OPEN_BC" << std::endl;
                    for (int j = 0; j < facesForEachBoundaryGroup[i].size(); ++j) {
                        //big assumption on the nature of elementIDs here...
                        addFace(tempElementVector[boundarFaces[facesForEachBoundaryGroup[i][j]].elementNum], boundarFaces[facesForEachBoundaryGroup[i][j]].localFaceIndex, NULL, 0, Geometry::FaceType::OPEN_BC);
                    }
                } else if (centaurBCType < 7001) {
                    std::cout << "periodic boundary for group " << i << " ignored for being internal; node connections will be assigned later" << std::endl;
                } else if (centaurBCType < 8001) {
                    std::cout << "shadow boundary for group " << i << " ignored for being internal; node connections will be assigned later" << std::endl;
                } else if (centaurBCType < 8501) {
                    std::cout << "interface boundary for group " << i << " ignored for being internal" << std::endl;
                } else if (centaurBCType < 9001) {
                    std::cout << "wake boundary for group " << i << " assigned as OPEN_BC" << std::endl;
                    for (int j = 0; j < facesForEachBoundaryGroup[i].size(); ++j) {
                        //big assumption on the nature of elementIDs here...
                        addFace(tempElementVector[boundarFaces[facesForEachBoundaryGroup[i][j]].elementNum], boundarFaces[facesForEachBoundaryGroup[i][j]].localFaceIndex, NULL, 0, Geometry::FaceType::OPEN_BC);
                    }
                } else if (centaurBCType < 10001) {
                    std::cout << "moving wall boundary for group " << i << " assigned as WALL_BC" << std::endl;
                    for (int j = 0; j < facesForEachBoundaryGroup[i].size(); ++j) {
                        //big assumption on the nature of elementIDs here...
                        addFace(tempElementVector[boundarFaces[facesForEachBoundaryGroup[i][j]].elementNum], boundarFaces[facesForEachBoundaryGroup[i][j]].localFaceIndex, NULL, 0, Geometry::FaceType::WALL_BC);
                    }
                } else {
                    std::cout << "alternative boundary condition for group " << i << " assigned as WALL_BC" << std::endl;
                    for (int j = 0; j < facesForEachBoundaryGroup[i].size(); ++j) {
                        //big assumption on the nature of elementIDs here...
                        addFace(tempElementVector[boundarFaces[facesForEachBoundaryGroup[i][j]].elementNum], boundarFaces[facesForEachBoundaryGroup[i][j]].localFaceIndex, NULL, 0, Geometry::FaceType::WALL_BC);
                    }
                }
                std::cout << "total number of boundary faces: " << getFacesList().size() << std::endl;
            }
            centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
            if (checkInt != sizeOfLine) {
                std::cerr << "Error in centaur file " << std::endl;
                return;
            }

            //This is where centaur tells the names of all the boundary groups
            centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
            for (int i = 0; i < numberOfBoundaryGroups; ++i) {
                centaurFile.read(reinterpret_cast<char*> (&nameOfBoundaryCondition[0]), sizeof (nameOfBoundaryCondition));
                std::cout << "boundary condition " << i << " is called " << nameOfBoundaryCondition << std::endl;
            }
            centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
            if (checkInt != sizeOfLine) {
                std::cerr << "Error in centaur file " << std::endl;
                return;
            }

            //Then comes periodic boundary information
            //file versions 3 and greater store some extra information that hpGEM will be constructing itself
            //this extra information mangles the reading of the usefull information a bit
            double transformationData[16];
            uint32_t matchingNodes[2];
            uint32_t numberOfPeriodicNodes;

            if (centaurFileType > 3) {
                uint32_t numberOfPeriodicTransformations;
                centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
                centaurFile.read(reinterpret_cast<char*> (&numberOfPeriodicTransformations), sizeof (numberOfPeriodicTransformations));
                std::cout << "There are " << numberOfPeriodicTransformations << " periodic boundary -> shadow boundary transformation(s)" << std::endl;
                centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
                if (checkInt != sizeOfLine) {
                    std::cerr << "Error in centaur file " << std::endl;
                    return;
                }
                for (int i = 0; i < numberOfPeriodicTransformations; ++i) {
                    //information on how to do the transformation can be computed later so just throw it away now
                    centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
                    centaurFile.read(reinterpret_cast<char*> (&ijunk), sizeof (ijunk));
                    centaurFile.read(reinterpret_cast<char*> (&transformationData[0]), sizeof (transformationData));
                    centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
                    if (checkInt != sizeOfLine) {
                        std::cerr << "Error in centaur file " << std::endl;
                        return;
                    }
                    centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
                    centaurFile.read(reinterpret_cast<char*> (&ijunk), sizeof (ijunk));
                    centaurFile.read(reinterpret_cast<char*> (&transformationData[0]), sizeof (transformationData));
                    centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
                    if (checkInt != sizeOfLine) {
                        std::cerr << "Error in centaur file " << std::endl;
                        return;
                    }




                    //now read the amount of periodic nodes
                    centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
                    centaurFile.read(reinterpret_cast<char*> (&numberOfPeriodicNodes), sizeof (numberOfPeriodicNodes));
                    std::cout << "transformation group " << i << " contains " << numberOfPeriodicNodes << " node -> node matching(s)" << std::endl;
                    centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
                    if (checkInt != sizeOfLine) {
                        std::cerr << "Error in centaur file " << std::endl;
                        return;
                    }

                    //and the actual pairing information
                    centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
                    for (int j = 0; j < numberOfPeriodicNodes; j++) {
                        centaurFile.read(reinterpret_cast<char*> (&matchingNodes[0]), sizeof (matchingNodes));

                        ///\bug for EXTREMELY coarse meshes this will destroy the distinction between faces on the boundary of the domain. Workaround: use at least 3 nodes per direction on each face.
                        auto& target = listOfElementsForEachNode[matchingNodes[0]];
                        auto first = std::move(target);
                        auto& second = listOfElementsForEachNode[matchingNodes[1]];

                        //We just std::move()d target, put it back in a defined state
                        target.clear();
                        target.reserve(first.size() + second.size());
                        std::set_union(first.begin(), first.end(), second.begin(), second.end(), target.begin());

                        listOfElementsForEachNode[matchingNodes[1] ] = listOfElementsForEachNode[matchingNodes[0] ];
                    }
                    centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
                    if (checkInt != sizeOfLine) {
                        std::cerr << "Error in centaur file " << std::endl;
                        return;
                    }
                }
            } else {
                //now read the amount of periodic nodes
                centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
                centaurFile.read(reinterpret_cast<char*> (&numberOfPeriodicNodes), sizeof (numberOfPeriodicNodes));
                std::cout << "the transformation group contains " << numberOfPeriodicNodes << " node -> node matching(s)" << std::endl;
                centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
                if (checkInt != sizeOfLine) {
                    std::cerr << "Error in centaur file " << std::endl;
                    return;
                }

                //and the actual pairing information
                centaurFile.read(reinterpret_cast<char*> (&sizeOfLine), sizeof (sizeOfLine));
                for (int j = 0; j < numberOfPeriodicNodes; j++) {
                    centaurFile.read(reinterpret_cast<char*> (&matchingNodes[0]), sizeof (matchingNodes));

                    auto& target = listOfElementsForEachNode[matchingNodes[0]];
                    auto first = std::move(target);
                    auto& second = listOfElementsForEachNode[matchingNodes[1]];

                    //We just std::move()d target, put it back in a defined state
                    target.clear();
                    target.reserve(first.size() + second.size());
                    std::set_union(first.begin(), first.end(), second.begin(), second.end(), target.begin());

                    listOfElementsForEachNode[matchingNodes[1] ] = listOfElementsForEachNode[matchingNodes[0] ];
                }
                centaurFile.read(reinterpret_cast<char*> (&checkInt), sizeof (checkInt));
                if (checkInt != sizeOfLine) {
                    std::cerr << "Error in centaur file " << std::endl;
                    return;
                }
            }

            std::cout << "begin constructing internal faces and internal \"boundaries\"" << std::endl;
            

            //now we know periodicity information, construct the vertices
            //remember that listOfElementsForEachNode will in general contain duplicates
            
            auto& listOfVertices = theMesh_.getVerticesList(IteratorType::GLOBAL);
            
            bool addedNewVertex(false);
            
            for(size_t i=0;i<listOfElementsForEachNode.size();++i)
            {
                for(size_t j=0;j<listOfElementsForEachNode[i].size();++j)
                {
                    Element* current = elementsList[listOfElementsForEachNode[i][j]];
                    for(size_t k=0;k<current->getNrOfNodes();++k)
                    {
                        //if we did not jet deal with this node and it is the correct one
                        if(current->getNode(k)==nullptr && current->getPhysicalGeometry()->getNodeIndex(k)==i)
                        {
                            if(!addedNewVertex)
                            {
                                addVertex();
                                addedNewVertex=true;
                            }
                            listOfVertices.back()->addElement(current,k);
                        }
                    }
                }
                addedNewVertex=false;
            }
            
            
            
            
            faceFactory();
            edgeFactory();

            delete[] boundarFaces;
        } else {

            std::cerr << "Incorrect mesh file. This mesh appears to contain two DIMensional data" << std::endl;

        }

    }

    /*void MeshManipulator::findElementNumber(std::vector<int>& a, std::vector<int>& b, std::vector<int>& c, int aNumber, int bNumber, int cNumber, std::vector<int>& notOnFace, HalfFaceDescription& face, std::vector<Element*>& vectorOfElements) {
        //step 1: the element should be connected to all of the given nodes

        std::vector<int>::iterator aEntry(a.begin()), bEntry(b.begin()), cEntry(c.begin()), otherEntry(notOnFace.begin());
        //std::cout<<a.size()<<" "<<b.size()<<" "<<c.size()<<" "<<notOnFace.size()<<std::endl;

        //assumes the lists are already sorted
        while (!(*aEntry == *bEntry && *aEntry == *cEntry && *aEntry == *otherEntry)) {
            if (*aEntry <*bEntry) {
                aEntry++;
            }
            if (*bEntry <*cEntry) {
                bEntry++;
            }
            if (*cEntry <*otherEntry) {
                cEntry++;
            }
            if (*otherEntry<*aEntry) {
                otherEntry++;
            }
            //std::cout<<*aEntry<<" "<<*bEntry<<" "<<*cEntry<<" "<<*otherEntry<<std::endl;
            assert(aEntry != a.end());
            assert(bEntry != b.end());
            assert(cEntry != c.end());
            assert(otherEntry != notOnFace.end());
        }
        //std::cout<<*aEntry<<" "<<*bEntry<<" "<<*cEntry<<" "<<*otherEntry<<std::endl;

        face.elementNum = *aEntry;

        //step 2: find the actual face
        std::vector<unsigned int> facenumbers;
        for (int i = 0; i < vectorOfElements[*aEntry]->getPhysicalGeometry()->getNrOfFaces(); ++i) {
            bool aFound(false), bFound(false), cFound(false);
            vectorOfElements[*aEntry]->getPhysicalGeometry()->getGlobalFaceNodeIndices(i, facenumbers);
            for (int j = 0; j < facenumbers.size(); ++j) {
                aFound |= (facenumbers[j] == aNumber);
                bFound |= (facenumbers[j] == bNumber);
                cFound |= (facenumbers[j] == cNumber);
            }
            if (aFound && bFound && cFound) {
                face.localFaceIndex = i;
                return;
            }
        }
        throw "The centaur mesh file contains a face with wrong bounding nodes!";
    }

    //in principle this will also work for 2D but fixing DIM makes the code a lot easier to read

    void MeshManipulator::constructInternalFaces(std::vector<std::vector<int> >& listOfElementsForEachNode, std::vector<Element*>& vectorOfElements) {
        int numberOfFaces(0);
        for (Element* currentElement : vectorOfElements) {
            for (int currentFace = 0; currentFace < (currentElement)->getPhysicalGeometry()->getNrOfFaces(); ++currentFace) {

                //step 1: find the other element
                std::vector<int>::iterator elementListEntry[3];
                std::vector<unsigned int> faceNode;
                (currentElement)->getPhysicalGeometry()->getGlobalFaceNodeIndices(currentFace, faceNode);
                for (int i = 0; i < 3; ++i) {
                    elementListEntry[i] = listOfElementsForEachNode[faceNode[i]].begin();
                }
                while (!(*elementListEntry[0] == *elementListEntry[1] && *elementListEntry[0] == *elementListEntry[2])) {
                    if (*elementListEntry[0]<*elementListEntry[1]) {
                        elementListEntry[0]++;
                    }
                    if (*elementListEntry[1]<*elementListEntry[2]) {
                        elementListEntry[1]++;
                    }
                    if (*elementListEntry[2]<*elementListEntry[0]) {
                        elementListEntry[2]++;
                    }
                }
                if (*elementListEntry[0] != (currentElement)->getID()) {
                    //std::cout<<"found candidate matching "<<*elementListEntry[0]<<"->"<<(*currentElement)->getID()<<std::endl;

                    //step 2: find the other face
                    //the trick used for the boundary faces wont work because of periodic 'internal' faces
                    for (int otherFace = 0; otherFace < vectorOfElements[*elementListEntry[0]]->getPhysicalGeometry()->getNrOfFaces(); ++otherFace) {

                        //but we can just find the original element back
                        std::vector<int>::iterator otherElementListEntry[3];
                        std::vector<unsigned int> otherFaceNode;
                        vectorOfElements[*elementListEntry[0]]->getPhysicalGeometry()->getGlobalFaceNodeIndices(otherFace, otherFaceNode);
                        for (int i = 0; i < 3; ++i) {
                            otherElementListEntry[i] = listOfElementsForEachNode[otherFaceNode[i]].begin();
                        }
                        while (!(*otherElementListEntry[0] == *otherElementListEntry[1] && *otherElementListEntry[0] == *otherElementListEntry[2])) {
                            if (*otherElementListEntry[0]<*otherElementListEntry[1]) {
                                otherElementListEntry[0]++;
                            }
                            if (*otherElementListEntry[1]<*otherElementListEntry[2]) {
                                otherElementListEntry[1]++;
                            }
                            if (*otherElementListEntry[2]<*otherElementListEntry[0]) {
                                otherElementListEntry[2]++;
                            }
                        }
                        //std::cout<<"finding the other element again("<<*otherElementListEntry[0]<<")"<<std::endl;
                        //assert(*otherElementListEntry[0]==*elementListEntry[0]);//for other faces this element may not be the first match found
                        if (*otherElementListEntry[0] == (currentElement)->getID() && *otherElementListEntry[0] == *otherElementListEntry[1] && *otherElementListEntry[0] == *otherElementListEntry[2]) {
                            //std::cout<<"found the fist element("<<*otherElementListEntry[0]<<")"<<std::endl;
                            addFace(currentElement, currentFace, vectorOfElements[*elementListEntry[0]], otherFace);
                            numberOfFaces++;
                        } else {
                            otherElementListEntry[0]++;
                            while (!(*otherElementListEntry[0] == *otherElementListEntry[1] && *otherElementListEntry[0] == *otherElementListEntry[2]) &&
                                    otherElementListEntry[0] != listOfElementsForEachNode[otherFaceNode[0]].end() &&
                                    otherElementListEntry[1] != listOfElementsForEachNode[otherFaceNode[1]].end() &&
                                    otherElementListEntry[2] != listOfElementsForEachNode[otherFaceNode[2]].end()) {
                                if (*otherElementListEntry[0]<*otherElementListEntry[1]) {
                                    otherElementListEntry[0]++;
                                }
                                if (*otherElementListEntry[1]<*otherElementListEntry[2]) {
                                    otherElementListEntry[1]++;
                                }
                                if (*otherElementListEntry[2]<*otherElementListEntry[0]) {
                                    otherElementListEntry[2]++;
                                }
                            }
                            if (*otherElementListEntry[0] == (currentElement)->getID() && *otherElementListEntry[0] == *otherElementListEntry[1] && *otherElementListEntry[0] == *otherElementListEntry[2]) {
                                //std::cout<<"found the fist element("<<*otherElementListEntry[0]<<")"<<std::endl;
                                addFace(currentElement, currentFace, vectorOfElements[*elementListEntry[0]], otherFace);
                                numberOfFaces++;
                            }
                        }
                    }

                } else {
                    //this face was already treated or it is a boundary face
                }
            }
        }
        std::cout << "Total number of Faces: " << numberOfFaces << std::endl;
    }*/


    /// \bug does not do the bc flags yet
    void MeshManipulator::faceFactory() {
        unsigned int DIM = configData_->dimension_;

        std::vector<unsigned int> nodeIndices;
        std::vector<Element*> candidates;
        
        for(Element* element:theMesh_.getElementsList(IteratorType::GLOBAL))
        {
            for(size_t i=0;i<element->getNrOfFaces();++i)
            {
                std::vector<const Node*> localNodes;
                //if this face is not there yet
                if(element->getFace(i)==nullptr)
                {
                    localNodes.clear();
                    candidates.clear();
                    element->getReferenceGeometry()->getCodim1EntityLocalIndices(i,nodeIndices);
                    
                    candidates=element->getNode(nodeIndices[0])->getElements();
                    localNodes.push_back(element->getNode(nodeIndices[0]));
                    std::sort(candidates.begin(),candidates.end(),[](Element* left, Element* right){return left->getID()<right->getID();});
                    for(size_t j=1;j<nodeIndices.size();++j)
                    {
                        localNodes.push_back(element->getNode(nodeIndices[j]));
                        std::vector<Element*> temp,nextIndices;
                        nextIndices=element->getNode(nodeIndices[j])->getElements();
                        std::sort(nextIndices.begin(),nextIndices.end(),[](Element* left, Element* right){return left->getID()<right->getID();});
                        std::set_intersection(candidates.begin(),candidates.end(),nextIndices.begin(),nextIndices.end(),std::back_inserter(temp),[](Element* left, Element* right){return left->getID()<right->getID();});
                        candidates=std::move(temp);
                    }
                    
                    //the current element does not bound the face or more than two elements bound the face
                    if(candidates.size()==0||candidates.size()>2)
                    {
                        std::cerr << "Invalid number of bounding elements detected for face " << theMesh_.getFacesList(IteratorType::GLOBAL).size() + 1 << std::endl;
                        return;
                    }
                    //boundary face
                    if(candidates.size() == 1)
                    {
                        assert(candidates[0]==element);
                        addFace(element,i,nullptr,0,Geometry::FaceType::WALL_BC);
                    }
                    if(candidates.size() == 2)
                    {
                        Element* other;
                        if(candidates[0]==element)
                        {
                            other=candidates[1];
                        }
                        else
                        {
                            other=candidates[0];
                        }
                        std::vector<unsigned int> otherNodeIndices;
                        for(size_t j=0;j<other->getNrOfFaces();++j)
                        {
                            other->getReferenceGeometry()->getCodim1EntityLocalIndices(j,otherNodeIndices);
                            bool match =true;
                            for(unsigned int k:otherNodeIndices)
                            {
                                if(std::find(localNodes.begin(),localNodes.end(),other->getNode(k))==localNodes.end())
                                {
                                    match=false;
                                }
                            }
                            if(match)
                            {
                                addFace(element,i,other,j);
                            }
                        }
                    }
                }
            }
        }
        
        
        
        
        /* //This will store the number of faces.
        unsigned int numOfFaces;

        //Half face holds global node number of both nodes, the element number and the local face number in that element.
        HalfFaceDescription halfFace;

        //List to hold the half faces. List is used for the quick sorting of the halfFaces
        std::vector<HalfFaceDescription> halfFaceList;
        VectorOfElementPtrT tempElementVector(getElementsList().size());

        VectorOfPointIndicesT globalFaceIndexes;
        int insertposition = 0;

        //first loop over the elements reading in all the data of the half faces
        unsigned int elementID = 0;
        for (typename ListOfElementsT::iterator it = elementColBegin(); it != elementColEnd(); ++it) {
            const Geometry::PhysicalGeometry * const myPhysicalGeometry = (*it)->getPhysicalGeometry();
            const Geometry::ReferenceGeometry * const myReferenceGeometry = (*it)->getReferenceGeometry();
            tempElementVector[elementID] = &(**it);

            numOfFaces = myReferenceGeometry->getNrOfCodim1Entities();

            //Loop over the faces on the element
            for (int i = 0; i < numOfFaces; i++) {
                halfFace.nodeList.clear();
                halfFace.nodeList.resize(DIM);
                myPhysicalGeometry->getGlobalFaceNodeIndices(i, globalFaceIndexes);

                halfFace.elementNum = elementID;
                halfFace.localFaceIndex = i;
                //make sure the node numbers in the face description are sorted
                for (int j = 0; j < globalFaceIndexes.size(); ++j) {
                    while (insertposition < DIM && halfFace.nodeList[insertposition] < globalFaceIndexes[j]) {
                        halfFace.nodeList[insertposition] = halfFace.nodeList[insertposition + 1];
                        ++insertposition;
                    }
                    if (insertposition <= DIM) {
                        halfFace.nodeList[insertposition - 1] = globalFaceIndexes[j];
                    }
                    insertposition = 0;
                }
                //Make sure the firstNode number is the smaller of the two global node number, required for the sorting.
                //                 if (globalFaceIndexes[0]<globalFaceIndexes[1])
                //                 {
                //                     halfFace.firstNode=globalFaceIndexes[0];
                //                     halfFace.secondNode=globalFaceIndexes[1];
                //                 }
                //                 else 
                //                 {
                //                     halfFace.firstNode=globalFaceIndexes[1];
                //                     halfFace.secondNode=globalFaceIndexes[0];
                //                 }

                //Add the halfFace to the list
                halfFaceList.push_back(halfFace);

            }

            elementID++;
        }

        //Now sort the list on the two value (firstNode, secondNode) so it order and the two pair internal halfFaces are next to each other
        std::sort(halfFaceList.begin(), halfFaceList.end(), compareHalfFace);

        //Writing out for testing
        for (typename std::vector<HalfFaceDescription>::const_iterator cit = halfFaceList.begin(); cit != halfFaceList.end(); ++cit) {
            //std::cout << (*cit).elementNum << " " << (*cit).localFaceIndex<< " "<< (*cit).nodeList[0]<<" " << (*cit).nodeList[1]<<" "<<  std::endl;
        }

        //Now create the faces
        HalfFaceDescription current;
        HalfFaceDescription next;

        //Loop over all the half faces
        for (typename std::vector<HalfFaceDescription>::const_iterator cit = halfFaceList.begin(); cit != halfFaceList.end(); ++cit) {
            //For the first halfFace just read it in, as we cannot tell if is an internal or external yet.
            if (cit == halfFaceList.begin()) {
                current = *cit;
            } else {
                next = *cit;

                //This is an interal face
                if ((current.nodeList[0] == next.nodeList[0]) && (current.nodeList[1] == next.nodeList[1]) && ((DIM < 3) || (current.nodeList[2] == next.nodeList[2]))) {

                    addFace(tempElementVector[current.elementNum], current.localFaceIndex, tempElementVector[next.elementNum], next.localFaceIndex);
                    //jump a face
                    ++cit;
                    current = *cit;

                    //There was only one face left, so it now the end of the list
                    if (cit == halfFaceList.end()) {
                        break;
                    }
                } else //it is a boundary face 
                {
                    addFace(tempElementVector[current.elementNum], current.localFaceIndex, NULL, 0, Geometry::WALL_BC);
                    current = next;
                }

            }
        }*/
        std::cout << "Total number of Faces: " << getFacesList(IteratorType::GLOBAL).size() << std::endl;
    }

    //the algorithm for the edge factory is based on that of the face factory
    //with some minor adaptation to account for the fact that there may be
    //more than two elements per edge
    
    ///\bug does not do 4D yet
    void MeshManipulator::edgeFactory() {
        unsigned int DIM(configData_->dimension_);
        //'edges' in DIM 2 are actually nodes
        if(DIM != 2) {
            std::vector<unsigned int> nodeList, otherNodeList;

            const Node* nodes[2];

            for(Element* element:theMesh_.getElementsList(IteratorType::GLOBAL))
            {
                for(size_t i=0;i<element->getNrOfEdges();++i)
                {
                    if(element->getEdge(i)==nullptr)
                    {
                        element->getReferenceGeometry()->getCodim2EntityLocalIndices(i,nodeList);
                        std::vector<Element*> candidates(0);
                        auto& leftElements = element->getNode(nodeList[0])->getElements();
                        auto& rightElements = element->getNode(nodeList[1])->getElements();
                        std::set_intersection(leftElements.begin(),leftElements.end(),
                                rightElements.begin(),rightElements.end(),std::back_inserter(candidates),[](Element* a, Element* b){return a->getID()<b->getID();});
                        if(candidates.size()==0)
                        {
                            std::cerr << "current element is not adjacent to its own edges" << std::endl;
                            return;
                        }
                        addEdge();
                        nodes[0]=element->getNode(nodeList[0]);
                        nodes[1]=element->getNode(nodeList[1]);
                        Edge* newEdge = theMesh_.getEdgesList(IteratorType::GLOBAL).back();
                        newEdge->addElement(element,i);
                        for(size_t j=1;j<candidates.size();++j)
                        {
                            Element* other = candidates[j];
                            for(size_t k=0;k<other->getNrOfEdges();++k)
                            {
                                other->getReferenceGeometry()->getCodim2EntityLocalIndices(k,otherNodeList);
                                if((other->getNode(otherNodeList[0])==nodes[0] || other->getNode(otherNodeList[0])==nodes[1]) &&
                                   (other->getNode(otherNodeList[1])==nodes[0] || other->getNode(otherNodeList[1])==nodes[1]))
                                {
                                    newEdge->addElement(other,k);
                                }
                            }
                        }
                    }
                }
            }
        }
        
        
        
        
        /*//the halfFaceDescription is designed to store partial information for objects that still need to be linked, so it will also work for edges
        HalfFaceDescription halfEdge;

        std::vector<HalfFaceDescription> halfEdgeList;

        VectorOfPointIndicesT globalEdgeIndexes;
        int insertposition = 0;
        unsigned int temp, dummy;

        for (auto it = tempElementVector.begin(); it != tempElementVector.end(); ++it) {
            const Geometry::PhysicalGeometry * const myPhysicalGeometry = (*it)->getPhysicalGeometry();
            const Geometry::ReferenceGeometry * const myReferenceGeometry = (*it)->getReferenceGeometry();
            tempElementVector[(*it)->getID()] = &(**it);
            bool inserted(false);

            numberOfEdges = myReferenceGeometry->getNrOfCodim2Entities();

            //Loop over the edges on the element
            for (int i = 0; i < numberOfEdges; ++i) {
                halfEdge.nodeList.clear();
                myReferenceGeometry->getCodim2EntityLocalIndices(i, globalEdgeIndexes);
                for (int j = 0; j < globalEdgeIndexes.size(); ++j)
                    globalEdgeIndexes[j] = myPhysicalGeometry->getNodeIndex(globalEdgeIndexes[j]);

                halfEdge.elementNum = (*it)->getID();
                halfEdge.localFaceIndex = i;
                //make sure the node numbers in the edge description are sorted
                halfEdge.nodeList.push_back(globalEdgeIndexes[0]);
                for (int j = 1; j < globalEdgeIndexes.size(); ++j) {
                    for (auto iit = halfEdge.nodeList.begin(); iit != halfEdge.nodeList.end(); ++iit) {
                        if (!inserted && (*iit) > globalEdgeIndexes[j]) {
                            halfEdge.nodeList.insert(iit, globalEdgeIndexes[j]);
                            inserted = true;
                        }
                    }
                    if (!inserted) {
                        halfEdge.nodeList.push_back(globalEdgeIndexes[j]);
                    }
                    inserted = false;
                }
                //Add the halfFace to the list
                halfEdgeList.push_back(halfEdge);
            }
        }
        std::sort(halfEdgeList.begin(), halfEdgeList.end(), compareHalfFace);
        HalfFaceDescription current;

        std::vector<Element*> elements;
        std::vector<unsigned int> edgeNrs;

        for (typename std::vector<HalfFaceDescription>::const_iterator cit = halfEdgeList.begin(); cit != halfEdgeList.end(); ++cit) {
            //For the first halfFace just read it in, as we cannot tell if is an internal or external yet.
            if (elements.empty()) {
                elements.push_back(tempElementVector[cit->elementNum]);
                edgeNrs.push_back(cit->localFaceIndex);
                current = *cit;
            } else {
                //This edge is not complete yet
                if ((current.nodeList[0] == cit->nodeList[0]) && (current.nodeList[1] == cit->nodeList[1])) {
                    elements.push_back(tempElementVector[cit->elementNum]);
                    edgeNrs.push_back(cit->localFaceIndex);
                } else //cit point to a new edge
                {
                    addEdge(elements, edgeNrs);
                    elements.clear();
                    edgeNrs.clear();
                    elements.push_back(tempElementVector[cit->elementNum]);
                    edgeNrs.push_back(cit->localFaceIndex);
                    current = *cit;
                }
            }
        }
        //dont forget to add the final edge
        addEdge(elements, edgeNrs);*/
    }



    //---------------------------------------------------------------------

    /*int MeshManipulator::getNumberOfMeshes() const {
        return vecOfElementTree_.size();
    }

    //! Create a new (empty) mesh-tree.

    void MeshManipulator::createNewMeshTree() {
        //vecOfElementTree_.push_back(new ElementLevelTreeT);
        //vecOfFaceTree_.push_back(new FaceLevelTreeT);
        //++numMeshTree_;
        //setActiveMeshTree(numMeshTree_ - 1);
    }

    //! Get the element container of a specific mesh-tree.

    typename MeshManipulator::ElementLevelTreeT*
    MeshManipulator::ElCont(int meshTreeIdx) const {
        int onIndex;
        if ((meshTreeIdx >= 0) && (meshTreeIdx < numMeshTree_)) {
            onIndex = meshTreeIdx;
        } else if ((meshTreeIdx == -1) && (activeMeshTree_ >= 0)) {
            onIndex = activeMeshTree_;
        } else
            throw "MeshManipulator::ElCont(): invalid mesh-tree index or no active mesh-tree.";

        return vecOfElementTree_[onIndex];
    }

    //! Get the face container of a specific mesh-tree.

    typename MeshManipulator::FaceLevelTreeT*
    MeshManipulator::FaCont(int meshTreeIdx) const {
        int onIndex;
        if ((meshTreeIdx >= 0) && (meshTreeIdx < numMeshTree_)) {
            onIndex = meshTreeIdx;
        } else if ((meshTreeIdx == -1) && (activeMeshTree_ >= 0)) {
            onIndex = activeMeshTree_;
        } else
            throw "MeshManipulator::FaCont(): invalid mesh-tree index or no active mesh-tree.";

        return vecOfFaceTree_[onIndex];
    }

    //! Some mesh generator: centaur / rectangular / triangle / tetrahedra / triangular-prism.

    void
    MeshManipulator::someMeshGenerator(int meshTreeIdx) {
        //set to correct value in case some other meshmanipulator changed things
        ElementFactory::instance().setCollectionOfBasisFunctionSets(&collBasisFSet_);
        ElementFactory::instance().setNumberOfMatrices(numberOfElementMatrixes_);
        ElementFactory::instance().setNumberOfVectors(numberOfFaceVectors_);
        ElementFactory::instance().setNumberOfTimeLevels(configData_->numberOfTimeLevels_);
        ElementFactory::instance().setNumberOfUnknowns(configData_->numberOfUnknowns_);
        FaceFactory::instance().setNumberOfFaceMatrices(numberOfFaceMatrixes_);
        FaceFactory::instance().setNumberOfFaceVectors(numberOfFaceVectors_);
        int onIndex;
        if ((meshTreeIdx >= 0) && (meshTreeIdx < numMeshTree_))
        {
          onIndex = meshTreeIdx;
        }
        else if ((meshTreeIdx == -1) && (activeMeshTree_ >= 0))
        {
          onIndex = activeMeshTree_;
        }
        else
          throw "MeshManipulator::someMeshGenerator(): invalid mesh-tree index or no active mesh-tree.";
          
        const int numberOfElement = 1 + (rand() % 10);
        const int startElementId = (onIndex+1)*1000;
        for (int id=startElementId; id<startElementId+numberOfElement; ++id)
        {
          ElementT el(id);
          vecOfElementTree_[onIndex]->addEntry(el);
        }

        const int numberOfFace = 1 + (rand() % 10);
        const int startFaceId = (onIndex+1)*1000;
        for (int id=startFaceId; id<startFaceId+numberOfFace; ++id)
        {
          FaceT fa(id);
          vecOfFaceTree_[onIndex]->addEntry(fa);
        }
    }

    //! Set active mesh-tree.

    void
    MeshManipulator::setActiveMeshTree(unsigned int meshTreeIdx) {
        if (meshTreeIdx < numMeshTree_)
            activeMeshTree_ = meshTreeIdx;
        else
            throw "MeshManipulator<DIM>::setActiveMeshTree(): invalid mesh-tree index.\n";
    }

    //! Get active mesh-tree index.

    int
    MeshManipulator::getActiveMeshTree() const {
        return activeMeshTree_;
    }

    //! Reset active mesh-tree.

    void
    MeshManipulator::resetActiveMeshTree() {
        activeMeshTree_ = -1;
    }*/

    Mesh&
    MeshManipulator::getMesh() {
        return theMesh_;
    }

    const Mesh&
    MeshManipulator::getMesh() const {
        return theMesh_;
    }


    //! Get maximum h-level of a specific mesh-tree.

    /*unsigned int
    MeshManipulator::getMaxLevel(int meshTreeIdx) const {
        int onIndex;
        if ((meshTreeIdx >= 0) && (meshTreeIdx < numMeshTree_))
        {
          onIndex = meshTreeIdx;
        }
        else if ((meshTreeIdx == -1) && (activeMeshTree_ >= 0))
        {
          onIndex = activeMeshTree_;
        }
        else
          throw "MeshManipulator::getMaxLevel(): invalid mesh-tree index or no active mesh-tree.";
        
        return vecOfElementTree_[onIndex]->maxLevel();
    }

    //! Set active level of a specific mesh-tree.

    void
    MeshManipulator::setActiveLevel(unsigned int meshTreeIdx, int level) {
        int onIndex;
        if ((meshTreeIdx >= 0) && (meshTreeIdx < numMeshTree_)) {
            onIndex = meshTreeIdx;
        } else if ((meshTreeIdx == -1) && (activeMeshTree_ >= 0)) {
            onIndex = activeMeshTree_;
        } else
            throw "MeshManipulator::setActiveLevel(): invalid mesh-tree index or no active mesh-tree.";

        //if ((level >= 0) && (level <= vecOfElementTree_[onIndex]->maxLevel()))
        //{
        //  vecOfElementTree_[onIndex]->setActiveLevel(level);
        //  vecOfFaceTree_[onIndex]->setActiveLevel(level);
        //}
        //else
        throw "MeshManipulator::setActiveLevel(): invalid level.";

    }

    //! Get active level of a specific mesh-tree.

    int
    MeshManipulator::getActiveLevel(int meshTreeIdx) const {
        int onIndex;
        if ((meshTreeIdx >= 0) && (meshTreeIdx < numMeshTree_)) {
            onIndex = meshTreeIdx;
        } else if ((meshTreeIdx == -1) && (activeMeshTree_ >= 0)) {
            onIndex = activeMeshTree_;
        } else
            throw "MeshManipulator::getActiveLevel(): invalid mesh-tree index or no active mesh-tree.";

        //return vecOfElementTree_[onIndex]->getActiveLevel();
    }

    //! Reset active level of a specific mesh-tree.

    void
    MeshManipulator::resetActiveLevel(int meshTreeIdx) {
        int onIndex;
        if ((meshTreeIdx >= 0) && (meshTreeIdx < numMeshTree_)) {
            onIndex = meshTreeIdx;
        } else if ((meshTreeIdx == -1) && (activeMeshTree_ >= 0)) {
            onIndex = activeMeshTree_;
        } else
            throw "MeshManipulator::resetActiveLevel(): invalid mesh-tree index or no active mesh-tree.";

        //vecOfElementTree_[onIndex]->resetActiveLevel();
        //vecOfFaceTree_[onIndex]->resetActiveLevel();
    }

    //! Duplicate mesh contents including all refined meshes.

    void
    MeshManipulator::duplicate(unsigned int fromMeshTreeIdx, unsigned int toMeshTreeIdx, unsigned int upToLevel) {
        ///\todo implement method
    }

    //! Refine a specific mesh-tree.

    void
    MeshManipulator::doRefinement(unsigned int meshTreeIdx, int refinementType) {
        int level = getMaxLevel(meshTreeIdx);
        setActiveLevel(meshTreeIdx, level);

        if (refinementType == -1) {
            // elements should have been flagged before calling this
            std::cout << "MeshManipulator::doRefinement(Mesh(" << meshTreeIdx << "))\n";
            doElementRefinement(meshTreeIdx);
            doFaceRefinement(meshTreeIdx);
            //for (ElementIteratorT it=ElCont(meshTreeIdx)->beginLevel(level); it != ElCont(meshTreeIdx)->end(); ++it)
            {
                // reset element's marking for refinement
                //           it->unsetRefineType();
                //           it->setBeingRefinedOff();
            }
        } else {
             // apply the uniform mesh refinement
             for (ElementIteratorT it=ElCont(meshTreeIdx)->beginLevel(level); it != ElCont(meshTreeIdx)->end(); ++it)
             {
                 // mark element for refinement
                 it->setRefineType(refinementType);
             }

             std::cout << "MeshManipulator<" << DIM << ">::doRefinement(" << meshTreeIdx << "," << refinementType << ")\n";
             doElementRefinement(meshTreeIdx);
             doFaceRefinement(meshTreeIdx);
             for (ElementIteratorT it=ElCont(meshTreeIdx).beginLevel(level); it != ElCont(meshTreeIdx).end(); ++it)
             {
                 // reset element's marking for refinement
     //           it->unsetRefineType();
     //           it->setBeingRefinedOff();
             }
        }

        level = getMaxLevel(meshTreeIdx);
        setActiveLevel(meshTreeIdx, level);
    }

    //! Do refinement on the elements.

    void
    MeshManipulator::doElementRefinement(unsigned int meshTreeIdx) {
        unsigned int needDummyFaceOnLevel = 0;
        std::vector<ElementT*> vecElementsToRefined; // list of unrefined elements
        //for (ElementIteratorT el=ElCont(meshTreeIdx)->begin(); el != ElCont(meshTreeIdx)->end(); ++el)
        {
               Geometry::RefinementGeometry* RG = el->getRefinementGeometry();
            
                int refineType;
        //         refineType = el->getRefineType();   // TODO: add this to Element



                if (refineType  < 0)   // not refined, just duplicate the element
                {
                    vecElementsToRefined.push_back(&(*el));
                    needDummyFaceOnLevel = el->getLevel()+1;
                    continue;
                }

                // ********* Add new nodes to the container
                //----------------------------------------
                RG->getAllNodes(refineType, VectorOfPointPhysicalsT& nodes);
                unsigned int nrNewNodes = nrOfNewNodes(refineType);
                unsigned int nrAllNodes = nodes.size();
                unsigned int nrOldNodes = nrAllNodes - nrNewNodes;
            
                // get physical nodes of this elements
                VectorOfPhysicalPointsT nodesPhys;    // vector of physical points
                VectorOfPointIndexesT   nodesIdx;     // vector of global indices
                for (PointIndexT j=0; j<nrVertices; ++j)
                {
                    PointPhysicalT p;       // a physical node
                    PointIndexT pIdx;       // a global index
                    PG->getVertexPoint (j, p);
                    pIdx = PG->getVertexIndex (j);
                    nodesPhys.push_back(p);
                    nodesIdx.push_back(pIdx);
                }

                // get new physical nodes due to the refinement, and add them
                // up to the nodes collection of this element
                PG->NewPhysicalNodes(refineType, nodesPhys);

                // add the new physical nodes into the nodes container
                // and get their global indices
                for (LocalPointIndexT j=nrVertices; j<nodesPhys.size(); ++j)
                {
                    int pIdx = PCptr->getGlobalIndex(nodesPhys[j]);
                    if (pIdx >= 0)
                    {
                        // it's already exist
                        nodesIdx.push_back(pIdx);
                    }
                    else
                    {
                        // it's not exist yet.  Add it to the node container
                        MeshBase<DIM>::AllNodes.addRoot(nodesPhys[j]);
                        nodesIdx.push_back(PCptr->size()-1);
                    }
                }
                // Now we already have all nodes: being used by this element and
                // to be used by the new sub-elements.

                // ********* Add sub-elements to the container
                //----------------------------------------
                unsigned int nrNewElements = PG->nrOfSubElements(refineType);
                std::vector<ElementBase*> vecSubElements;
                for (unsigned int j=0; j<nrNewElements; ++j)
                {
                    VectorOfPointIndexesT localNodeIndexes;
                    PG->SubElementNodeIndexes(refineType, j, localNodeIndexes);
                    ElementDescriptor elDescriptor(localNodeIndexes.size());
                    for (VectorOfPointIndexesT k=0; k<localNodeIndexes.size(); ++k)
                    {
                        elDescriptor.addNode(nodesIdx[localNodeIndexes[k]]);
                    }

        //             ElementFactory<DIM> elFactory(*this);
        //             ElementT *elem = elFactory.makeElement(&elDescriptor);
        //             elem->setRefinementType(refineType);
        //             vecSubElements.push_back(elem);
                }
            
                // Add the sub-elements as the children of this element
                ElementIteratorT elIt = ElCont(meshTreeIdx).addChildren(el, vecSubElements);
            
                // Clear up the storage
                nodesPhys.clear(); // clear vector of physical points
                nodesIdx.clear();  // clear vector of global index

            
                // ********* Add sub-Internal Faces to the container
                //----------------------------------------
                VectorOfPointIndexesT elementIdx1;
                VectorOfPointIndexesT elementIdx2;
                VectorOfPointIndexesT localFaceIdx1;
                VectorOfPointIndexesT localFaceIdx2;
                PG->AdjacentSubElementsPairs(refineType, elementIdx1,localFaceIdx1, elementIdx2,localFaceIdx2);
                unsigned int nrOfNewFaces = elementIdx1.size();

                std::vector<FaceT*> vecSubFaces;
                for (unsigned int j=0; j<nrOfNewFaces; ++j)
                {
                    // Create the face, connecting the two elements
                    ElementT* el1 = vecSubElements[elementIdx1[j]];
                    ElementT* el2 = vecSubElements[elementIdx2[j]];
                    FaceT* iFace = new FaceT( el1, localFaceIdx1[j], el2, localFaceIdx2[j]);
                    vecSubFaces.push_back(iFace);
                }
                // Add them to the container as the children of a dummy face
                FaceIteratorT fa = FaCont(meshTreeIdx).getDummyFace(el->getLevel());

                FaCont(meshTreeIdx).addChildren(fa, vecSubFaces);
                vecSubFaces.clear();
            
                // Mark that this element was just refined
                el->setJustRefined();
                vecSubElements.clear();

            } // end of loop over elements

            if (needDummyFaceOnLevel)
            {
                FaCont(meshTreeIdx).getDummyFace(needDummyFaceOnLevel);
            }
        
            std::vector<ElementT*> vecEmpty;  // empty list of new elements
            while(!vecElementsToRefined.empty()) 
            {
                ElementT *elem = vecElementsToRefined.back();
                ElementIteratorT elIt = ElCont(meshTreeIdx).addChildren(elem->getIterator(), vecEmpty);
                vecElementsToRefined.pop_back();
        }
    }

    //! Do refinement on the faces.

    void
    MeshManipulator::doFaceRefinement(unsigned int meshTreeIdx) {
        ///\bug nothing happens
    }


    //! Check whether the two elements may be connected by a face or not.

    void
    MeshManipulator::pairingCheck(const ElementIterator elL, unsigned int locFaceNrL,
            const ElementIterator elR, unsigned int locFaceNrR,
            int& pairingValue, bool& sizeOrder)
    // pairingValue: 0=not match, 1=partial match, 2=perfect match
    // sizeOrder: true = LR,  false = RL
    {
        // get node numbers from left (and right) sides
        std::vector<PointIndexT> globNodesL;
        std::vector<PointIndexT> globNodesR;

        const Geometry::PhysicalGeometry * const leftPG = (*elL)->getPhysicalGeometry();
        const Geometry::PhysicalGeometry* rightPG;

        leftPG->getGlobalFaceNodeIndices(locFaceNrL, globNodesL);
        if (*elR != 0) {
            rightPG = (*elR)->getPhysicalGeometry();
            rightPG->getGlobalFaceNodeIndices(locFaceNrR, globNodesR);
        }

        // store them as sets
        std::set<PointIndexT> setL(globNodesL.data(), globNodesL.data() + globNodesL.size());
        std::set<PointIndexT> setR(globNodesR.data(), globNodesR.data() + globNodesR.size());

        if (setL == setR) {
            // nodes of the faces are match
            pairingValue = 2;
            sizeOrder = true;
            return;
        }

        std::set<PointIndexT> commNodes;
        std::set_intersection(setL.begin(), setL.end(), setR.begin(), setR.end(), std::inserter(commNodes, commNodes.begin()));

        unsigned int nrNodesL = globNodesL.size();
        unsigned int nrNodesR = globNodesR.size();
        if (nrNodesL != nrNodesR) {
            pairingValue = 0;
            return;
        }

        std::set<PointIndexT> diffNodesL;
        std::set_difference(setL.begin(), setL.end(), commNodes.begin(), commNodes.end(), std::inserter(diffNodesL, diffNodesL.end()));

        std::set<PointIndexT> diffNodesR;
        std::set_difference(setR.begin(), setR.end(), commNodes.begin(), commNodes.end(), std::inserter(diffNodesR, diffNodesR.end()));

        typedef std::set<PointIndexT>::iterator SetIterType;

        // do collinear test for all possible pairs
        pairingValue = 1;

        // order of face sizes: true = LR;  false = RL
        sizeOrder = true;

        unsigned int DIM = configData_->dimension_;

        // loop until empty or until the faces proved not collinear
        while (!diffNodesL.empty()) {
            SetIterType itDiffL = diffNodesL.begin();
            PointPhysicalT ppL(DIM);
            leftPG->getGlobalNodeCoordinates(*itDiffL, ppL);

            bool foundCollinear(false);
            SetIterType itDiffL2;
            SetIterType itDiffR;
            SetIterType itDiffR2;
            for (itDiffR = diffNodesR.begin(); itDiffR != diffNodesR.end(); ++itDiffR) {
                PointPhysicalT ppR(DIM);
                rightPG->getGlobalNodeCoordinates(*itDiffR, ppR);

                if (commNodes.size() > 0) {
                    for (SetIterType itCommon = commNodes.begin(); itCommon != commNodes.end(); ++itCommon) {
                        PointPhysicalT pp0(DIM);
                        leftPG->getGlobalNodeCoordinates(*itCommon, pp0);

                        //-------------
                        // perform 3-nodes collinear test
                        PointPhysicalT dL(ppL - pp0);
                        PointPhysicalT dR(ppR - pp0);

                        double ratio = 0.;
                        unsigned int d;
                        for (d = 0; d < DIM; ++d) {
                            if (std::abs(dR[d]) > Geometry::SmallerDoubleThanMinimalSizeOfTheMesh) {
                                ratio = dL[d] / dR[d];
                                break;
                            }
                        }

                        if (ratio < 0.) {
                            pairingValue = 0;
                            return;
                        }

                        foundCollinear = true;
                        for (; d < DIM; ++d) {
                            if (std::abs(dR[d]) > Geometry::SmallerDoubleThanMinimalSizeOfTheMesh) {
                                if (dL[d] / dR[d] < 0.) {
                                    foundCollinear = false;
                                    break;
                                }

                                if (std::abs(dL[d] / dR[d] - ratio) > Geometry::SmallerDoubleThanMinimalSizeOfTheMesh) {
                                    foundCollinear = false;
                                    break;
                                }
                            }
                        }

                        // order of face sizes: true = LR;  false = RL
                        sizeOrder = (sizeOrder && (ratio <= 1.0));

                        // found a collinear nodes pair
                        if (foundCollinear) break;
                    } // for itCommon
                }// if (commNodes.size() > 0)
                else
                    // no common nodes
                {
                    for (itDiffL2 = diffNodesL.begin(); itDiffL2 != diffNodesL.end(); ++itDiffL2) {
                        if ((itDiffL2 == itDiffL) || (itDiffL2 == itDiffR))
                            continue;

                        PointPhysicalT ppL2(DIM);
                        leftPG->getGlobalNodeCoordinates(*itDiffL2, ppL2);

                        for (itDiffR2 = diffNodesR.begin(); itDiffR2 != diffNodesR.end(); ++itDiffR2) {
                            if ((*itDiffR2 == *itDiffR) || (*itDiffR2 == *itDiffL) || (*itDiffR2 == *itDiffL2))
                                continue;

                            PointPhysicalT ppR2(DIM);
                            rightPG->getGlobalNodeCoordinates(*itDiffR2, ppR2);

                            //-------------
                            // perform two 3-nodes pairs collinear test 
                            PointPhysicalT dL1(ppL - ppR);
                            PointPhysicalT dR1(ppL2 - ppR);

                            PointPhysicalT dL2(ppL2 - ppR);
                            PointPhysicalT dR2(ppR2 - ppR);

                            PointPhysicalT dL3(ppL - ppR2);
                            PointPhysicalT dR3(ppR - ppR2);

                            double ratio1 = 0.;
                            double ratio2 = 0.;
                            double ratio3 = 0.;
                            unsigned int d1;
                            unsigned int d2;
                            unsigned int d3;
                            for (d1 = 0; d1 < DIM; ++d1) {
                                if (std::abs(dR1[d1]) > Geometry::SmallerDoubleThanMinimalSizeOfTheMesh) {
                                    ratio1 = dL1[d1] / dR1[d1];
                                    break;
                                }
                            }

                            for (d2 = 0; d2 < DIM; ++d2) {
                                if (std::abs(dR2[d2]) > Geometry::SmallerDoubleThanMinimalSizeOfTheMesh) {
                                    ratio2 = dL2[d2] / dR2[d2];
                                    break;
                                }
                            }

                            for (d3 = 0; d3 < DIM; ++d3) {
                                if (std::abs(dR3[d3]) > Geometry::SmallerDoubleThanMinimalSizeOfTheMesh) {
                                    ratio3 = dL3[d3] / dR3[d3];
                                    break;
                                }
                            }

                            if ((ratio1 < 0) || (ratio2 < 0) || (ratio3 < 0)) {
                                pairingValue = 0;
                                return;
                            }

                            foundCollinear = true;
                            for (; d1 < DIM; ++d1) {
                                if (std::abs(dR1[d1]) > Geometry::SmallerDoubleThanMinimalSizeOfTheMesh) {
                                    if (dL1[d1] / dR1[d1] < 0.) {
                                        pairingValue = 0;
                                        return;
                                    }

                                    if (std::abs(dL1[d1] / dR1[d1] - ratio1) > Geometry::SmallerDoubleThanMinimalSizeOfTheMesh) {
                                        foundCollinear = false;
                                        break;
                                    }
                                }
                            }

                            for (; (d2 < DIM) && foundCollinear; ++d2) {
                                if (std::abs(dR2[d2]) > Geometry::SmallerDoubleThanMinimalSizeOfTheMesh) {
                                    if (dL2[d2] / dR2[d2] < 0.) {
                                        pairingValue = 0;
                                        return;
                                    }

                                    if (std::abs(dL2[d2] / dR2[d2] - ratio2) > Geometry::SmallerDoubleThanMinimalSizeOfTheMesh) {
                                        foundCollinear = false;
                                        break;
                                    }
                                }
                            }

                            for (; (d3 < DIM) && foundCollinear; ++d3) {
                                if (std::abs(dR3[d3]) > Geometry::SmallerDoubleThanMinimalSizeOfTheMesh) {
                                    if (dL3[d3] / dR3[d3] < 0.) {
                                        pairingValue = 0;
                                        return;
                                    }

                                    if (std::abs(dL3[d3] / dR3[d3] - ratio3) > Geometry::SmallerDoubleThanMinimalSizeOfTheMesh) {
                                        foundCollinear = false;
                                        break;
                                    }
                                }
                            }

                            // found a collinear nodes pair
                            if (foundCollinear) break; // itDiffR2 will be deleted
                        }
                        // found a collinear nodes pair
                        if (foundCollinear) break; // itDiffL2 will be deleted
                    } // for itDiffL2
                } //  else, no common nodes

                // found a collinear nodes pair
                if (foundCollinear) break; // itDiffR will be deleted
            } // for itDiffR

            if (foundCollinear)
                // found a collinear nodes pair
            {
                if (commNodes.size() > 0) {
                    diffNodesR.erase(itDiffR);
                } else {
                    diffNodesL.erase(itDiffL2);
                    diffNodesR.erase(itDiffR);
                    diffNodesR.erase(itDiffR2);
                }
            } else
                // the faces are not collinear
            {
                pairingValue = 0;
                break;
            }

            diffNodesL.erase(itDiffL);
        } // end while-loop

        if (!diffNodesR.empty())
            pairingValue = 0;

    } // end pairingCheck

    //! Check whether the two elements may be connected by a face or not in periodic face case.
    void
    MeshManipulator::periodicPairingCheck(const FaceIteratorT fa,
                                               const ElementIteratorT elL, unsigned int localFaceNrL,
                                               const ElementIteratorT elR, unsigned int localFaceNrR,
                                               int& pairingValue, bool& sizeOrder)
    {
        
    }*/
    //---------------------------------------------------------------------

    int MeshManipulator //: public MeshRefiner <DIM>
    ::dimension() const {
        return configData_->dimension_;
    }

}

#endif

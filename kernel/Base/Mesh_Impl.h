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


#include <limits>

#include "MpiContainer.h"
#include "Mesh.h"
#include "Element.h"
#include "Face.h"
#include "Edge.h"
#include "ElementFactory.h"
#include "FaceFactory.h"
#include "Geometry/PointPhysical.h"
#include "Geometry/ReferenceGeometry.h"
#include "Geometry/PointReference.h"
#include "FaceCacheData.h"
#include "ElementCacheData.h"

#ifdef HPGEM_USE_METIS
#include <metis.h>
#else
//replace the metis internal size type by something reasonably sane
using idx_t = std::size_t;
#endif

namespace Base
{

    template<std::size_t DIM>
    Mesh<DIM>::Mesh()
            : hasToSplit_(false)
    {
    }

    template<std::size_t DIM>
    Mesh<DIM>::Mesh(const Mesh& orig)
            : hasToSplit_(true),
        nodeCoordinates_(orig.nodeCoordinates_)
    {
        logger(ERROR, "To be done again with the help of Boost::serialize");
        //Make elements. Note: each element gets a new unique ID
        /*for(Element* element : orig.elements_)
        {
            elements_.addRootEntry(element->copyWithoutFacesEdgesNodes());
            ++elementCounter_;
        }
        logger.assert(orig.elementCounter_ == elementCounter_, "In the copy constructor of Mesh,"
            "there is a different number of elements in the new Mesh than in the old one.");
        
        //Make nodes and couple them with elements
        for(Node* node : orig.nodes_)
        {
            nodes_.push_back(new Node(node->getID()));
            std::vector<Element*> nodeElements = node->getElements();
            logger.assert(nodeElements.size() > 0, "There are no elements at this node.");
            for(std::size_t i = 0; i < nodeElements.size(); ++i)
            {
                std::size_t id = nodeElements[i]->getID();
                nodes_.back()->addElement(elements_[id], node->getNodeNumber(i));
            }
            ++nodeCounter_;
        }
        logger.assert(orig.nodeCounter_ == nodeCounter_, "In the copy constructor of Mesh,"
            "there is a different number of nodes in the new Mesh than in the old one.");
        
        //Make faces and couple them with elements
        for(Face* face : orig.faces_)
        {
            Element* elLeft = elements_[face->getPtrElementLeft()->getID()];
            std::size_t idOnLeft = face->localFaceNumberLeft();
            Element* elRight = nullptr;
            std::size_t idOnRight = 0;
            if (face->isInternal())
            {
                elRight = elements_[face->getPtrElementRight()->getID()];
                idOnRight = face->localFaceNumberRight();
            }
            faces_.addRootEntry(new Face(*face, elLeft, idOnLeft, elRight, idOnRight));
            ++faceCounter_;
        }
        logger.assert(orig.faceCounter_ == faceCounter_, "In the copy constructor of Mesh,"
            "there is a different number of faces in the new Mesh than in the old one.");
        
        //Make nodes and couple them with elements
        for(Edge* edge : orig.edges_)
        {
            edges_.addRootEntry(new Edge(edge->getID()));
            std::vector<Element*> edgeElements = edge->getElements();
            logger.assert(edgeElements.size() > 0, "There are no elements at this edge.");
            for(std::size_t i = 0; i < edgeElements.size(); ++i)
            {
                std::size_t id = edgeElements[i]->getID();
                (*(--edges_.end()))->addElement(elements_[id], edge->getEdgeNumber(i));
            }
            ++edgeCounter_;
        }
        logger.assert(orig.edgeCounter_ == edgeCounter_, "In the copy constructor of Mesh,"
            "there is a different number of nodes in the new Mesh than in the old one.");
        
        //call split() to make get...List() const to work with local iterators
        split();*/
    }

    template<std::size_t DIM>
    Mesh<DIM>::~Mesh()
    {
        clear();        
    }

    template<std::size_t DIM>
    Element* Mesh<DIM>::addElement(const std::vector<std::size_t>& globalNodeIndexes)
    {
        Element* newElement = ElementFactory::instance().makeElement(globalNodeIndexes, nodeCoordinates_);
        elements_.addRootEntry(newElement);
        hasToSplit_ = true;
        newElement->setPositionInTree((--elements_.end()).getTreeEntry());
        return newElement;
    }

    template<std::size_t DIM>
    void Mesh<DIM>::addSubElements(Base::Element* parent, const std::vector<Base::Element*> subElements)
    {
        elements_.addChildren(parent->getPositionInTree()->getIterator(Base::TreeTraversalMethod::ALLLEVEL), subElements);
        for(std::size_t i = 0; i < parent->getPositionInTree()->getNumberOfChildren(); ++i)
        {
            subElements[i]->setPositionInTree(parent->getPositionInTree()->getChild(i));
        }
        hasToSplit_ = true;
    }

    template<std::size_t DIM>
    bool Mesh<DIM>::addFace(Element* leftElementPtr, std::size_t leftElementLocalFaceNo, Element* rightElementPtr, std::size_t rightElementLocalFaceNo, const Geometry::FaceType& faceType)
    {
        Face* newFace = nullptr;
        logger.assert(leftElementPtr!=nullptr, "Invalid element passed");
        if (rightElementPtr == nullptr)
        {
            newFace = FaceFactory::instance().makeFace(leftElementPtr, leftElementLocalFaceNo, faceType);
        }
        else
        {
            newFace = FaceFactory::instance().makeFace(leftElementPtr, leftElementLocalFaceNo, rightElementPtr, rightElementLocalFaceNo);
        }
        faces_.addRootEntry(newFace);
        newFace->setPositionInTree((--faces_.end()).getTreeEntry());
        hasToSplit_ = true;
        return true;
    }

    template<std::size_t DIM>
    void Mesh<DIM>::addSubFaces(const Base::Face* parent, const std::vector<Base::Face*> subFaces)
    {
        faces_.addChildren(parent->getPositionInTree()->getIterator(Base::TreeTraversalMethod::ALLLEVEL), subFaces);
        for(std::size_t i = 0; i < parent->getPositionInTree()->getNumberOfChildren(); ++i)
        {
            subFaces[i]->setPositionInTree(parent->getPositionInTree()->getChild(i));
        }
        hasToSplit_ = true;
    }

    template<std::size_t DIM>
    void Mesh<DIM>::addEdge()
    {
        Edge* newEdge = new Edge(GlobalUniqueIndex::instance().getEdgeIndex());
        edges_.addRootEntry(newEdge);
        newEdge->setPositionInTree((--edges_.end()).getTreeEntry());
        hasToSplit_ = true;
    }
    
    template<std::size_t DIM>
    void Mesh<DIM>::addNodeCoordinate(Geometry::PointPhysical<DIM> node)
    {
        nodeCoordinates_.push_back(node);
        //don't distribute the points here, it will confuse the elements
    }
    
    template<std::size_t DIM>
    void Mesh<DIM>::addNode()
    {
        nodes_.push_back(new Node(GlobalUniqueIndex::instance().getNodeIndex()));
        hasToSplit_ = true;
    }

    template<std::size_t DIM>
    void Mesh<DIM>::split()
    {
        //split the mesh
        int pid = 0, elementsProcessed(0);
        std::map<std::size_t, idx_t> contiguousPositionOfElement;
        elements_.setSingleLevelTraversal(0);
        for (Element* element : elements_)
        {
            contiguousPositionOfElement[element->getID()] = elementsProcessed++;
        }
#ifdef HPGEM_USE_MPI
#ifdef HPGEM_USE_METIS
        std::vector<idx_t> partition(elements_.size()); //output
        pid = MPIContainer::Instance().getProcessorID();
        idx_t nProcs = MPIContainer::Instance().getNumberOfProcessors();

        if (pid == 0 && nProcs > 1)
        {   
            logger(INFO, "start of metis");

            idx_t one = 1; //actually the number of constraints. This can be increased for example when we want to distribute an entire mesh tree in one go (while keeping each of the levels balanced) - increasing this number turns imbalance into a vector
            idx_t numberOfElements = elements_.size();

            //int mpiCommSize=4;
            float imbalance = 1.001;//explicitly put the default for later manipulation
            idx_t totalCutSize;//output
            
            std::vector<idx_t> xadj(numberOfElements + 1);//make sure not to put this data on the stack
            std::vector<idx_t> adjncy(2 * faces_.size());//if this basic connectivity structure turns out to be very slow for conforming meshes, some improvements can be made
            int connectionsUsed(0), xadjCounter(0);
            for (Element* element : elements_)
            {   
                xadj[xadjCounter] = connectionsUsed;
                xadjCounter++;
                for (int i = 0; i < element->getReferenceGeometry()->getNumberOfCodim1Entities(); ++i)
                {   
                    const Face* face = element->getFace(i);
                    if (face->isInternal())
                    {   
                        if (element == face->getPtrElementLeft())
                        {   
                            adjncy[connectionsUsed] = contiguousPositionOfElement[face->getPtrElementRight()->getID()];
                            connectionsUsed++;
                        }
                        else
                        {   
                            adjncy[connectionsUsed] = contiguousPositionOfElement[face->getPtrElementLeft()->getID()];
                            connectionsUsed++;
                        }
                    } //boundary faces don't generate connections
                }
            }
            xadj[xadjCounter] = connectionsUsed;

            idx_t metisOptions[METIS_NOPTIONS];
            METIS_SetDefaultOptions(metisOptions);

            metisOptions[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM;
            metisOptions[METIS_OPTION_RTYPE] = METIS_RTYPE_FM;

            //the empty arguments provide options for fine-tuning the weights of nodes, edges and processors, these are currently assumed to be the same
            METIS_PartGraphKway(&numberOfElements, &one, xadj.data(), adjncy.data(), NULL, NULL, NULL, &nProcs, NULL, &imbalance, metisOptions, &totalCutSize, partition.data());
            //mpiCommunicator.Bcast((void *)&partition[0],partition.size(),MPI::INT,0);//broadcast the computed partition to all the nodes
            logger(INFO, "Done splitting mesh.");

        }

        MPIContainer::Instance().broadcast(partition, 0);
#endif
#else
        std::vector<std::size_t> partition(elements_.size()); //output
#endif
        submeshes_.clear();
        auto elementIterator = elements_.begin();
        elementIterator.setSingleLevelTraversal(0);
        for (auto targetIterator = partition.begin(); targetIterator != partition.end(); ++targetIterator, ++elementIterator)
        {
            if (pid == *targetIterator)
            {
                submeshes_.add(*elementIterator);
            }
        }
        faces_.setSingleLevelTraversal(0);
        for (Base::Face* face : faces_)
        {
            //Are we part of this face?
            if (partition[contiguousPositionOfElement[face->getRootElement()->getID()]] == pid
                || (face->isInternal() &&
                    partition[contiguousPositionOfElement[face->getPtrElementRight()->getID()]] == pid))
            {
                //if we are a part of this face
                submeshes_.add(face);

                if (face->isInternal() && contiguousPositionOfElement[face->getPtrElementRight()->getID()] < elements_.size() &&
                    (partition[contiguousPositionOfElement[face->getRootElement()->getID()]] !=
                     partition[contiguousPositionOfElement[face->getPtrElementRight()->getID()]]))
                {
                    if (face->getFaceType() == Geometry::FaceType::INTERNAL)
                    {
                        face->setFaceType(Geometry::FaceType::SUBDOMAIN_BOUNDARY);
                    }
                    else if (face->getFaceType() == Geometry::FaceType::PERIODIC_BC)
                    {
                        face->setFaceType(Geometry::FaceType::PERIODIC_SUBDOMAIN_BC);
                    }
                    //else the facetype was already updated during a previous partitioning pass
                    if (partition[contiguousPositionOfElement[face->getRootElement()->getID()]] == pid)
                    {
                        //don't send to yourself, ask the element on the other side what pid to sent to
                        //use the private addPush and addPull because they are more efficient
                        submeshes_.addPush(face->getRootElement(), partition[contiguousPositionOfElement[face->getPtrElementRight()->getID()]]);
                        //if you receive, the source is the owner of the element
                        submeshes_.addPull(face->getPtrElementRight(), partition[contiguousPositionOfElement[face->getPtrElementRight()->getID()]]);
                    }
                    else
                    {
                        submeshes_.addPush(face->getPtrElementRight(), partition[contiguousPositionOfElement[face->getRootElement()->getID()]]);
                        submeshes_.addPull(face->getRootElement(), partition[contiguousPositionOfElement[face->getRootElement()->getID()]]);
                    }
                }
            }
        }
        
        //edges and nodes are only used to construct conforming basis-functions (i.e. GlobalAssembly))
        if(edges_.size() > 0) edges_.setSingleLevelTraversal(0);
        for (Base::Edge* edge : edges_)
        {
            if (partition[contiguousPositionOfElement[edge->getRootElement()->getID()]] == pid)
            {
                submeshes_.add(edge);
            }
        }
        for (Base::Node* node : nodes_)
        {
            if (partition[contiguousPositionOfElement[node->getRootElement()->getID()]] == pid)
            {
                submeshes_.add(node);
            }
        }
        elements_.setPreOrderTraversal();
        submeshes_.getElementsList().setPreOrderTraversal();
        auto subMeshElementIterator = submeshes_.getElementsList().begin();
        for(elementIterator = elements_.begin(); subMeshElementIterator != submeshes_.getElementsList().end(); ++elementIterator)
        {
            if(*elementIterator == *subMeshElementIterator){
                if (elementIterator.getTreeEntry()->hasChild())
                {
                    std::vector<Element *> children;
                    auto childIterator = elementIterator;
                    childIterator.setPreOrderTraversal();
                    (++childIterator).setAllLevelTraversal();
                    for (auto child : elementIterator.getTreeEntry()->getChildren())
                    {
                        children.push_back(*childIterator++);
                    }
                    submeshes_.getElementsList().addChildren(subMeshElementIterator, children);
                }
                ++subMeshElementIterator;
            }
        }
        faces_.setPreOrderTraversal();
        submeshes_.getFacesList().setPreOrderTraversal();
        auto faceIterator = faces_.begin();
        auto subMeshFaceIterator = submeshes_.getFacesList().begin();
        for(; subMeshFaceIterator != submeshes_.getFacesList().end(); ++faceIterator)
        {
            if(*faceIterator == *subMeshFaceIterator)
            {
                if (faceIterator.getTreeEntry()->hasChild())
                {
                    std::vector<Face *> children;
                    auto childIterator = faceIterator;
                    childIterator.setPreOrderTraversal();
                    (++childIterator).setAllLevelTraversal();
                    for (auto child : faceIterator.getTreeEntry()->getChildren())
                    {
                        children.push_back(*childIterator++);
                    }
                    submeshes_.getFacesList().addChildren(subMeshFaceIterator, children);
                }
                ++subMeshFaceIterator;
            }
        }
        edges_.setPreOrderTraversal();
        submeshes_.getEdgesList().setPreOrderTraversal();
        auto edgeIterator = edges_.begin();
        auto subMeshEdgeIterator = submeshes_.getEdgesList().begin();
        for(; subMeshEdgeIterator != submeshes_.getEdgesList().end(); ++edgeIterator)
        {
            if(*edgeIterator == *subMeshEdgeIterator)
            {
                if (edgeIterator.getTreeEntry()->hasChild())
                {
                    std::vector<Edge *> children;
                    auto childIterator = edgeIterator;
                    childIterator.setPreOrderTraversal();
                    (++childIterator).setAllLevelTraversal();
                    for (auto child : edgeIterator.getTreeEntry()->getChildren())
                    {
                        children.push_back(*childIterator++);
                    }
                    submeshes_.getEdgesList().addChildren(subMeshEdgeIterator, children);
                }
                ++subMeshEdgeIterator;
            }
        }
        hasToSplit_ = false;
    }

    template<std::size_t DIM>
    void Mesh<DIM>::clear()
    {
        for (Element* element : elements_)
        {
            delete element;
        }
        for (Face* face : faces_)
        {
            delete face;
        }
        for (Edge* edge : edges_)
        {
            delete edge;
        }
        for (Node* node : nodes_)
        {
            delete node;
        }
        elements_.clear();
        faces_.clear();
        edges_.clear();
        nodes_.clear();
        submeshes_.clear();
        nodeCoordinates_.clear();
    }

    template<std::size_t DIM>
    const LevelTree<Element*>& Mesh<DIM>::getElementsList(IteratorType part) const
    {
        if (part == IteratorType::LOCAL)
        {
            logger.assert_always(!hasToSplit_, "Please call getElementsList() on a modifiable mesh at least once"
                    "\nbefore calling getElementsList() const");
            return submeshes_.getElementsList();
        }
        else
        {
            return elements_;
        }
    }

    template<std::size_t DIM>
    LevelTree<Element*>& Mesh<DIM>::getElementsList(IteratorType part)
    {
        if (part == IteratorType::LOCAL)
        {
            if (hasToSplit_)
            {
                split();
            }
            return submeshes_.getElementsList();
        }
        else
        {
            return elements_;
        }
    }

    template<std::size_t DIM>
    const LevelTree<Face*>& Mesh<DIM>::getFacesList(IteratorType part) const
    {
        if (part == IteratorType::LOCAL)
        {
            logger.assert_always(!hasToSplit_, "Please call getFacesList() on a modifiable mesh at least once"
                    "\nbefore calling getFacesList() const");
            return submeshes_.getFacesList();
        }
        else
        {
            return faces_;
        }
    }

    template<std::size_t DIM>
    LevelTree<Face*>& Mesh<DIM>::getFacesList(IteratorType part)
    {
        if (part == IteratorType::LOCAL)
        {
            if (hasToSplit_)
            {
                split();
            }
            return submeshes_.getFacesList();
        }
        else
        {
            return faces_;
        }
    }

    template<std::size_t DIM>
    const LevelTree<Edge*>& Mesh<DIM>::getEdgesList(IteratorType part) const
    {
        if (part == IteratorType::LOCAL)
        {
            logger.assert_always(!hasToSplit_, "Please call getEdgesList() on a modifiable mesh at least once"
                    "\nbefore calling getEdgesList() const");
            return submeshes_.getEdgesList();
        }
        else
        {
            return edges_;
        }
    }

    template<std::size_t DIM>
    LevelTree<Edge*>& Mesh<DIM>::getEdgesList(IteratorType part)
    {
        if (part == IteratorType::LOCAL)
        {
            if (hasToSplit_)
            {
                split();
            }
            return submeshes_.getEdgesList();
        }
        else
        {
            return edges_;
        }
    }
    
    template<std::size_t DIM>
    const std::vector<Node*>& Mesh<DIM>::getNodesList(IteratorType part) const
    {
        if (part == IteratorType::LOCAL)
        {
            logger.assert_always(!hasToSplit_, "Please call getNodesList() on a modifiable mesh at least once"
                    "\nbefore calling getNodesList() const");
            return submeshes_.getNodesList();
        }
        else
        {
            return nodes_;
        }
    }
    
    template<std::size_t DIM>
    std::vector<Node*>& Mesh<DIM>::getNodesList(IteratorType part)
    {
        if (part == IteratorType::LOCAL)
        {
            if (hasToSplit_)
            {
                split();
            }
            return submeshes_.getNodesList();
        }
        else
        {
            return nodes_;
        }
    }

    template<std::size_t DIM>
    const std::vector<Geometry::PointPhysical<DIM> >& Mesh<DIM>::getNodeCoordinates() const
    {
        //for historic reasons points_ is referenced directly during element 
        //creation and therefore cannot be distributed
        return nodeCoordinates_;
    }

    template<std::size_t DIM>
    std::vector<Geometry::PointPhysical<DIM> >& Mesh<DIM>::getNodeCoordinates()
    {
        //for historic reasons points_ is referenced directly during element 
        //creation and therefore cannot be distributed
        return nodeCoordinates_;
    }

}

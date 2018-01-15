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
    Mesh<DIM>::Mesh() {}

    template<std::size_t DIM>
    Mesh<DIM>::Mesh(const Mesh& orig)
            : nodeCoordinates_(orig.nodeCoordinates_)
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
            "there is a different number of nodes in the new Mesh than in the old one.");*/
    }

    template<std::size_t DIM>
    Mesh<DIM>::~Mesh()
    {
        clear();        
    }

    template<std::size_t DIM>
    Element* Mesh<DIM>::addElement(const std::vector<std::size_t>& globalNodeIndexes)
    {
        //some users don't want to see their elements added locally so don't mess with the submesh here
        Element* newElement = ElementFactory::instance().makeElement(globalNodeIndexes, nodeCoordinates_);
        elements_.addRootEntry(newElement);
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
    }

    template<std::size_t DIM>
    Edge* Mesh<DIM>::addEdge()
    {
        Edge* newEdge = new Edge(GlobalUniqueIndex::instance().getEdgeIndex());
        edges_.addRootEntry(newEdge);
        newEdge->setPositionInTree((--edges_.end()).getTreeEntry());
        return newEdge;
    }
    
    template<std::size_t DIM>
    void Mesh<DIM>::addNodeCoordinate(Geometry::PointPhysical<DIM> node)
    {
        nodeCoordinates_.push_back(node);
    }
    
    template<std::size_t DIM>
    void Mesh<DIM>::addNode()
    {
        nodes_.push_back(new Node(GlobalUniqueIndex::instance().getNodeIndex()));
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
            return submeshes_.getElementsList();
        }
        else
        {
            logger(WARN, "getElementsList will only return elements this processor knows about (including shadow elements)");
            return elements_;
        }
    }

    template<std::size_t DIM>
    LevelTree<Element*>& Mesh<DIM>::getElementsList(IteratorType part)
    {
        if (part == IteratorType::LOCAL)
        {
            return submeshes_.getElementsList();
        }
        else
        {
            logger(WARN, "getElementsList will only return elements this processor knows about (including shadow elements)");
            return elements_;
        }
    }

    template<std::size_t DIM>
    const LevelTree<Face*>& Mesh<DIM>::getFacesList(IteratorType part) const
    {
        if (part == IteratorType::LOCAL)
        {
            return submeshes_.getFacesList();
        }
        else
        {
            logger(WARN, "getFacesList will only return faces this processor knows about");
            return faces_;
        }
    }

    template<std::size_t DIM>
    LevelTree<Face*>& Mesh<DIM>::getFacesList(IteratorType part)
    {
        if (part == IteratorType::LOCAL)
        {
            return submeshes_.getFacesList();
        }
        else
        {
            logger(WARN, "getFacesList will only return faces this processor knows about");
            return faces_;
        }
    }

    template<std::size_t DIM>
    const LevelTree<Edge*>& Mesh<DIM>::getEdgesList(IteratorType part) const
    {
        if (part == IteratorType::LOCAL)
        {
            return submeshes_.getEdgesList();
        }
        else
        {
            logger(WARN, "getEdgesList will only return edges this processor knows about");
            return edges_;
        }
    }

    template<std::size_t DIM>
    LevelTree<Edge*>& Mesh<DIM>::getEdgesList(IteratorType part)
    {
        if (part == IteratorType::LOCAL)
        {
            return submeshes_.getEdgesList();
        }
        else
        {
            logger(WARN, "getEdgesList will only return edges this processor knows about");
            return edges_;
        }
    }
    
    template<std::size_t DIM>
    const std::vector<Node*>& Mesh<DIM>::getNodesList(IteratorType part) const
    {
        if (part == IteratorType::LOCAL)
        {
            return submeshes_.getNodesList();
        }
        else
        {
            logger(WARN, "getNodesList will only return nodes this processor knows about");
            return nodes_;
        }
    }
    
    template<std::size_t DIM>
    std::vector<Node*>& Mesh<DIM>::getNodesList(IteratorType part)
    {
        if (part == IteratorType::LOCAL)
        {
            return submeshes_.getNodesList();
        }
        else
        {
            logger(WARN, "getNodesList will only return nodes this processor knows about");
            return nodes_;
        }
    }

    template<std::size_t DIM>
    const std::vector<Geometry::PointPhysical<DIM> >& Mesh<DIM>::getNodeCoordinates() const
    {
        return nodeCoordinates_;
    }

    template<std::size_t DIM>
    std::vector<Geometry::PointPhysical<DIM> >& Mesh<DIM>::getNodeCoordinates()
    {
        return nodeCoordinates_;
    }

}

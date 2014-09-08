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

#include "Mesh.hpp"
#include "Element.hpp"
#include "Face.hpp"
#include "Edge.hpp"
#include "ElementFactory.hpp"
#include "FaceFactory.hpp"
#include "Geometry/PointPhysical.hpp"
#include "FaceCacheData.hpp"
#include "ElementCacheData.hpp"

namespace Base{

Mesh::Mesh():
        elementcounter_(0),
        faceCounter_(0),
        edgeCounter_(0)
{
}

Mesh::Mesh(const Mesh& orig):
        elements_(orig.elements_),
        edges_(orig.edges_),
        faces_(orig.faces_),
        points_(orig.points_),
        elementcounter_(orig.elementcounter_),
        faceCounter_(orig.faceCounter_),
        edgeCounter_(orig.edgeCounter_)
{
}

Mesh::~Mesh() {
    for(Element* element:elements_){
        delete element;
    }
    for(Face* face:faces_){
        delete face;
    }
    for(Edge* edge:edges_){
        delete edge;
    }
    
}

    Element*                            Mesh::addElement(const std::vector<unsigned int>& globalNodeIndexes){
        elements_.push_back(ElementFactory::instance().makeElement(globalNodeIndexes,points_,elementcounter_));
        ++elementcounter_;
        return elements_.back();
    }

    bool                                Mesh::addFace(Element* leftElementPtr, unsigned int leftElementLocalFaceNo, 
                                                Element* rightElementPtr, unsigned int rightElementLocalFaceNo,
                                                const Geometry::FaceType& faceType){
        if(rightElementPtr==NULL){
            faces_.push_back(FaceFactory::instance().makeFace(leftElementPtr,leftElementLocalFaceNo,faceType,faceCounter_));
        }else{
            faces_.push_back(FaceFactory::instance().makeFace(leftElementPtr,leftElementLocalFaceNo,rightElementPtr,rightElementLocalFaceNo,faceCounter_));
        }
        ++faceCounter_;
        return true;
    }

    void                                Mesh::addEdge(std::vector< Element*> elements, std::vector<unsigned int> localEdgeNrs){
        edges_.push_back(new Edge(elements,localEdgeNrs,edgeCounter_));
        ++edgeCounter_;
    }
    
    void                                Mesh::addNode(Geometry::PointPhysical node){
        points_.push_back(node);
    }

}
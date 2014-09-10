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
#include "Geometry/ReferenceGeometry.hpp"
#include "Geometry/PointReference.hpp"
#include "FaceCacheData.hpp"
#include "ElementCacheData.hpp"

//#define hpGEM_INCLUDE_METIS_SUPPORT//temporarily activating this definition makes development easier on some IDEs
#ifdef hpGEM_INCLUDE_METIS_SUPPORT
#include "metis.h"
#endif

namespace Base{

Mesh::Mesh():
        elementcounter_(0),
        faceCounter_(0),
        edgeCounter_(0),
        localProcessorID_(0),
        submeshes_(1)//TODO get the last two initializers from MPI
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

void Mesh::split(){
        //split the mesh
#ifdef hpGEM_INCLUDE_METIS_SUPPORT
    std::cout<<"start of metis"<<std::endl;
    int one(1);//actually the number of constraints. This can be increased for example when we want to distribute an entire mesh tree in one go (while keeping each of the levels balanced) - increasing this number turns imbalance into a vector
    int numberOfElements(elements_.size());
    int PlaceholderForTheNumberOfMPINondes(4);
    float imbalance(1.001);//explicitly put the default for later manipulation
    int totalCutSize;//output
    std::vector<int> partition(numberOfElements);//output
    std::vector<int> xadj(numberOfElements+1);//for some reason c-style arrays break somewhere near 1e6 faces in a mesh, so use vectors
    std::vector<int> adjncy(2*faces_.size());//if this basic connectivity structure turns out to be very slow for conforming meshes, some improvements can be made
    int connectionsUsed(0),xadjCounter(0);
    for(Element* element:elements_){
        xadj[xadjCounter]=connectionsUsed;
        xadjCounter++;
        for(int i=0;i<element->getReferenceGeometry()->getNrOfCodim1Entities();++i){
            const Face* face=element->getFace(i);
            if(face->isInternal()){
                if(element==face->getPtrElementLeft()){
                    adjncy[connectionsUsed]=face->getPtrElementRight()->getID();
                    connectionsUsed++;
                }else{
                    adjncy[connectionsUsed]=face->getPtrElementLeft()->getID();
                    connectionsUsed++;
                }
            }//boundary faces dont generate connections
        }
    }
    xadj[xadjCounter]=connectionsUsed;
    
    int metisOptions[METIS_NOPTIONS];
    METIS_SetDefaultOptions(metisOptions);
    
    metisOptions[METIS_OPTION_CTYPE]=METIS_CTYPE_SHEM;
    metisOptions[METIS_OPTION_RTYPE]=METIS_RTYPE_FM;
    
    //the empty arguments provide options for fine-tuning the weights of nodes, edges and processors, these are currently assumed to be the same
    METIS_PartGraphKway(&numberOfElements,&one,&xadj[0],&adjncy[0],NULL,NULL,NULL,&PlaceholderForTheNumberOfMPINondes,NULL,&imbalance,metisOptions,&totalCutSize,&partition[0]);
    
    //temporary debug statements
    for(int i=0;i<numberOfElements;++i){
        std::cout<<partition[i]<<std::endl;
    }
#endif
        while(!elements_.empty()){
            //if this element belongs to this mesh
            submeshes_[localProcessorID_].add(elements_.front());
            elements_.pop_front();
        }
        while(!faces_.empty()){
            //if the left element OR the right element belongs to this mesh
            submeshes_[localProcessorID_].add(faces_.front());
            faces_.pop_front();
            //if the left element XOR the right element belongs to this mesh
            //do some fiddling with push and pull elements
        }
        while(!edges_.empty()){
            //if one of the elements adjacent to this edge belong to this mesh
            submeshes_[localProcessorID_].add(edges_.front());
            edges_.pop_front();
            //if above AND one of the elements adjacent to this edge belong to another mesh
            //AND there are conforming basis functions
            //do some fiddling with push and pull elements
        }
}

const std::list<Element*>&          Mesh::getElementsList() const {
    if(elements_.empty()){
        return submeshes_[localProcessorID_].getElementsList();
    }else{
        throw "Please call getElementsList() on a modifiable mesh at least once before calling getElementsList() const";
    }
}
std::list<Element*>&                Mesh::getElementsList() { 
    if(!elements_.empty()){
        split();
    }
    return submeshes_[localProcessorID_].getElementsList();
}

const std::list<Face*>&             Mesh::getFacesList() const { 
    if(faces_.empty()){
        return submeshes_[localProcessorID_].getFacesList();
    }else{
        throw "Please call getFacesList() on a modifiable mesh at least once before calling getFacesList() const";
    }
}
std::list<Face*>&                   Mesh::getFacesList() {
    if(!faces_.empty()){
        split();
    }
    return submeshes_[localProcessorID_].getFacesList();
}

const std::list<Edge*>&             Mesh::getEdgesList() const {
    if(edges_.empty()){
        return submeshes_[localProcessorID_].getEdgesList();
    }else{
        throw "Please call getEdgesList() on a modifiable mesh at least once before calling getEdgesList() const";
    }
}
std::list<Edge*>&                   Mesh::getEdgesList() {
    if(!edges_.empty()){
        split();
    }
    return submeshes_[localProcessorID_].getEdgesList();
}

const std::vector<Geometry::PointPhysical>&  Mesh::getNodes()const{
    //for historic reasons points_ is referenced directly during element creation and therefore cannot be distributed
    return points_;
}
std::vector<Geometry::PointPhysical>&        Mesh::getNodes(){
    //for historic reasons points_ is referenced directly during element creation and therefore cannot be distributed
    return points_;
}

}
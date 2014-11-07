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
#ifdef HPGEM_USE_MPI
#include <mpi.h>
#endif

#include <cassert>

#include "MpiContainer.hpp"
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

//#define HPGEM_USE_METIS//temporarily activating this definition makes development easier on some IDEs
#ifdef HPGEM_USE_METIS
#include <metis.h>
#endif

namespace Base{

Mesh::Mesh():
        elementcounter_(0),
        faceCounter_(0),
        edgeCounter_(0),
        localProcessorID_(0)
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
    std::vector<int> partition(elements_.size());//output
        //split the mesh
    int pid = 0;
#ifdef HPGEM_USE_MPI
#ifdef HPGEM_USE_METIS
    pid = MPIContainer::Instance().getProcessorID();
    int nProcs = MPIContainer::Instance().getNumProcessors();

    if(pid==0 && nProcs>1){
        std::cout<<"start of metis"<<std::endl;
        
        int one = 1;//actually the number of constraints. This can be increased for example when we want to distribute an entire mesh tree in one go (while keeping each of the levels balanced) - increasing this number turns imbalance into a vector
        int numberOfElements = elements_.size();
        
        //int mpiCommSize=4;
        float imbalance = 1.001;//explicitly put the default for later manipulation
        int totalCutSize;//output

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
                }//boundary faces don't generate connections
            }
        }
        xadj[xadjCounter]=connectionsUsed;

        int metisOptions[METIS_NOPTIONS];
        METIS_SetDefaultOptions(metisOptions);

        metisOptions[METIS_OPTION_CTYPE]=METIS_CTYPE_SHEM;
        metisOptions[METIS_OPTION_RTYPE]=METIS_RTYPE_FM;

        //the empty arguments provide options for fine-tuning the weights of nodes, edges and processors, these are currently assumed to be the same
        METIS_PartGraphKway(&numberOfElements,&one,&xadj[0],&adjncy[0],NULL,NULL,NULL,&nProcs,NULL,&imbalance,metisOptions,&totalCutSize,&partition[0]);
        //mpiCommunicator.Bcast((void *)&partition[0],partition.size(),MPI::INT,0);//broadcast the computed partition to all the nodes
        std::cout<<"done splitting mesh"<<std::endl;
        
    }
    
    MPIContainer::Instance().broadcast(partition, 0);
    
#endif
#endif
    auto elementIterator=elements_.begin();
    for(auto targetIterator=partition.begin();targetIterator!=partition.end();++targetIterator,++elementIterator){

        if(pid == *targetIterator){
            submeshes_.add(*elementIterator);
        }
    }
    for(Base::Face* face : faces_){
        //Are we part of this face?
        if(partition[face->getPtrElementLeft()->getID()] == pid ||
            (face->isInternal() && partition[face->getPtrElementRight()->getID()]==pid)){
            //yeah we are
            submeshes_.add(face);
            
            if(face->isInternal()&&(partition[face->getPtrElementLeft()->getID()]!=partition[face->getPtrElementRight()->getID()])){
                if(partition[face->getPtrElementLeft()->getID()]== pid ){
                    //dont send to yourself, ask the element on the other side what pid to sent to
                    submeshes_.addPush(face->getPtrElementLeft(),partition[face->getPtrElementRight()->getID()]);
                    //if you recieve, the source is the owner of the element
                    submeshes_.addPull(face->getPtrElementRight(),partition[face->getPtrElementRight()->getID()]);
                }else{
                    submeshes_.addPush(face->getPtrElementRight(),partition[face->getPtrElementLeft()->getID()]);
                    submeshes_.addPull(face->getPtrElementLeft(),partition[face->getPtrElementLeft()->getID()]);
                }
            }
        }
    }
    for(Base::Edge* edge:edges_){
        bool done(false);
        for(int i=0;(i<edge->getNrOfElements())&&(!done);++i){
            if(partition[edge->getElement(i)->getID()] == pid){
                submeshes_.add(edge);
                done=true;
            }
        }
    }
    std::cout << "#YOLOSWAG!" << std::endl;
    elements_.clear();
    faces_.clear();
    edges_.clear();
    
}

const std::vector<Element*>&          Mesh::getElementsList() const {
    if(elements_.empty()){
        return submeshes_.getElementsList();
    }else{
        throw "Please call getElementsList() on a modifiable mesh at least once before calling getElementsList() const";
    }
}
std::vector<Element*>&                Mesh::getElementsList() { 
    if(!elements_.empty()){
        split();
    }
    return submeshes_.getElementsList();
}

const std::vector<Face*>&             Mesh::getFacesList() const { 
    if(faces_.empty()){
        return submeshes_.getFacesList();
    }else{
        throw "Please call getFacesList() on a modifiable mesh at least once before calling getFacesList() const";
    }
}
std::vector<Face*>&                   Mesh::getFacesList() {
    if(!faces_.empty()){
        split();
    }
    return submeshes_.getFacesList();
}

const std::vector<Edge*>&             Mesh::getEdgesList() const {
    if(edges_.empty()){
        return submeshes_.getEdgesList();
    }else{
        throw "Please call getEdgesList() on a modifiable mesh at least once before calling getEdgesList() const";
    }
}
std::vector<Edge*>&                   Mesh::getEdgesList() {
    if(!edges_.empty()){
        split();
    }
    return submeshes_.getEdgesList();
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
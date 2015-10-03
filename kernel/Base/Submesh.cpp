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

#include "PhysicalElement.h"
#include "Submesh.h"
#include "Mesh.h"
#include "Element.h"
#include "Face.h"
#include "Edge.h"
#include "Geometry/PointPhysical.h" 
#include "FaceCacheData.h"
#include "ElementCacheData.h"

namespace Base
{
    void Submesh::add(Element* element)
    {
        logger.assert(element!=nullptr, "Invalid element passed");
        elements_.push_back(element);
    }
    
    void Submesh::addPush(Element* element, int processorID)
    {
        logger.assert(element!=nullptr, "Invalid element passed");
        pushElements_[processorID].push_back(element);
    }
    
    void Submesh::addPull(Element* element, int processorID)
    {
        logger.assert(element!=nullptr, "Invalid element passed");
        pullElements_[processorID].push_back(element);
    }
    
    void Submesh::add(Face* face)
    {
        logger.assert(face!=nullptr, "Invalid face passed");
        faces_.push_back(face);
    }
    
    void Submesh::add(Edge* edge)
    {
        logger.assert(edge!=nullptr, "Invalid edge passed");
        edges_.push_back(edge);
    }
    
    void Submesh::add(Node* node)
    {
        logger.assert(node!=nullptr, "Invalid node passed");
        nodes_.push_back(node);
    }
    
    void Submesh::clear()
    {
        elements_.clear();
        faces_.clear();
        edges_.clear();
        nodes_.clear();
        pullElements_.clear();
        pushElements_.clear();
    }

}

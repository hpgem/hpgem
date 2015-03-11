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

#include "ReferenceTetrahedron.h"
#include "ReferenceLine.h"
#include "ReferenceTriangle.h"
#include "Geometry/PointReference.h"
#include "Mappings/MappingToRefTriangleToTetrahedron.h"

namespace Geometry
{
    /* The ordering of the vertex and faces in a tetrahedron:
     * 3 o
     *   |\
     *   |  \2
     *   |  / \
     *   |/     \
     * 0 o--------o 1
     */
    std::size_t ReferenceTetrahedron::localNodeIndexes_[4][3] = { {0, 3, 2}, {0, 1, 3}, {0, 2, 1}, {1, 2, 3}};
    
    std::size_t ReferenceTetrahedron::localNodesOnEdge_[6][2] = { {0, 1}, {0, 2}, {0, 3}, {2, 3}, {1, 3}, {1, 2}, };
    
    ReferenceTetrahedron::ReferenceTetrahedron()
            : /// Tetrahedron has four nodes 3D + 1
            ReferenceGeometry(4, 3, TETRAHEDRON), referenceGeometryCodim1Ptr_(&ReferenceTriangle::Instance()), referenceGeometryCodim2Ptr_(&ReferenceLine::Instance())
    {
        PointReferenceT p1(3), p2(3), p3(3), p4(3);
        
        p1[0] = +0.0;
        p1[1] = +0.0;
        p1[2] = +0.0;
        p2[0] = +1.0;
        p2[1] = +0.0;
        p2[2] = +0.0;
        p3[0] = +0.0;
        p3[1] = +1.0;
        p3[2] = +0.0;
        p4[0] = +0.0;
        p4[1] = +0.0;
        p4[2] = +1.0;
        
        points_[0] = p1;
        points_[1] = p2;
        points_[2] = p3;
        points_[3] = p4;
        
        mappingsTriangleToTetrahedron_[0] = &MappingToRefTriangleToTetrahedron0::Instance();
        mappingsTriangleToTetrahedron_[1] = &MappingToRefTriangleToTetrahedron1::Instance();
        mappingsTriangleToTetrahedron_[2] = &MappingToRefTriangleToTetrahedron2::Instance();
        mappingsTriangleToTetrahedron_[3] = &MappingToRefTriangleToTetrahedron3::Instance();
        
        /// TODO: Implement MappingTetrahedronToTetrahedron.
        mappingsTetrahedronToTetrahedron_[0] = 0;
    }
    
    ReferenceTetrahedron::ReferenceTetrahedron(const ReferenceTetrahedron& copy)
            : ReferenceGeometry(copy), referenceGeometryCodim1Ptr_(copy.referenceGeometryCodim1Ptr_), referenceGeometryCodim2Ptr_(copy.referenceGeometryCodim2Ptr_)
    {
    }
    
    bool ReferenceTetrahedron::isInternalPoint(const PointReferenceT& p) const
    {
        logger.assert(p.size()==3, "The dimension of the reference point is incorrect");
        return ((p[0] >= 0.) && (p[0] <= 1.) && (p[1] >= 0.) && (p[1] <= 1. - p[0]) && (p[2] >= 0.) && (p[2] <= 1. - p[0] - p[1]));
    }
    
    PointReference ReferenceTetrahedron::getCenter() const
    {
        PointReference p(3);
        p[0] = p[1] = p[2] = 1. / 4.;
        return p;
    }
    
    const PointReference& ReferenceTetrahedron::getNode(const IndexT& i) const
    {
        logger.assert(i<getNumberOfNodes(), "Asked for node %, but there are only % nodes", i, getNumberOfNodes());
        return points_[i];
    }
    
    std::ostream& operator<<(std::ostream& os, const ReferenceTetrahedron& tetra)
    {
        os << tetra.getName() << " =( ";
        ReferenceTetrahedron::const_iterator it = tetra.points_.begin();
        ReferenceTetrahedron::const_iterator end = tetra.points_.end();
        
        for (; it != end; ++it)
        {
            os << (*it) << '\t';
        }
        os << ')' << std::endl;
        
        return os;
    }
    
    // ================================== Codimension 0 ============================================
    
    std::size_t ReferenceTetrahedron::getCodim0MappingIndex(const ListOfIndexesT& list1, const ListOfIndexesT& list2) const
    {
        /// TODO: Implement tetrahedron to tetrahedron mappings.
        throw "ERROR: Tetrahedron to tetrahedron mappings do not exist";
    }
    
    const MappingReferenceToReference*
    ReferenceTetrahedron::getCodim0MappingPtr(const IndexT i) const
    {
        /// TODO: Implement tetrahedron to tetrahedron mappings.
        throw "ERROR: Tetrahedron to tetrahedron mappings do not exist";
    }
    
    // ================================== Codimension 1 ============================================
    
    std::vector<std::size_t> ReferenceTetrahedron::getCodim1EntityLocalIndices(const IndexT faceIndex) const
    {
        if (faceIndex < 4)
        {
            return std::vector<std::size_t>(localNodeIndexes_[faceIndex], localNodeIndexes_[faceIndex] + 3);
        }
        else
        {
            throw "ERROR: Index out of range. Tetrahedron has only 3 faces.";
        }
    }
    
    const ReferenceGeometry*
    ReferenceTetrahedron::getCodim1ReferenceGeometry(const IndexT faceIndex) const
    {
        if (faceIndex < 4)
        {
            return referenceGeometryCodim1Ptr_;
        }
        else
        {
            throw "ERROR: Index out of range. Tetrahedron has only 3 faces.";
        }
    }
    
    const MappingReferenceToReference*
    ReferenceTetrahedron::getCodim1MappingPtr(const IndexT faceIndex) const
    {
        if (faceIndex < 4)
        {
            return mappingsTriangleToTetrahedron_[faceIndex];
        }
        else
        {
            throw "ERROR: Asked for a square point index larger than 3. There are only 4 nodes in a square!";
        }
    }
    
    // ================================== Codimension 2 ============================================
    
    std::vector<std::size_t> ReferenceTetrahedron::getCodim2EntityLocalIndices(const IndexT edgeIndex) const
    {
        if (edgeIndex < 6)
        {
            return std::vector<std::size_t>(localNodesOnEdge_[edgeIndex], localNodesOnEdge_[edgeIndex] + 2);
        }
        else
        {
            throw "ERROR: Index out of range. Tetrahedron has only 6 edges.";
        }
    }
    
    const ReferenceGeometry*
    ReferenceTetrahedron::getCodim2ReferenceGeometry(const IndexT edgeIndex) const
    {
        if (edgeIndex < 6)
        {
            return referenceGeometryCodim2Ptr_;
        }
        else
        {
            throw "ERROR: Index out of range. Tetrahedron has only 6 edges.";
        }
    }
    
    const MappingReferenceToReference*
    ReferenceTetrahedron::getCodim2MappingPtr(const IndexT faceIndex) const
    {
        /// TODO: Implement line to tetrahedron mappings.
        throw "ERROR: Line to tetrahedron mappings do not exist";
    }
    
    // ================================== Codimension 3 ============================================
    
    std::vector<std::size_t> ReferenceTetrahedron::getCodim3EntityLocalIndices(const IndexT nodeIndex) const
    {
        if (nodeIndex < 4)
        {
            return std::vector<std::size_t>(1, nodeIndex);
        }
        else
        {
            throw "ERROR: Index out of range. Tetrahedron has only 4 nodes.";
        }
    }

}
;


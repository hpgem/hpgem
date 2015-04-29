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
            ReferenceGeometry(4, 3, ReferenceGeometryType::TETRAHEDRON, {1./4., 1./4., 1./4.}), referenceGeometryCodim1Ptr_(&ReferenceTriangle::Instance()), referenceGeometryCodim2Ptr_(&ReferenceLine::Instance())
    {
        name = "ReferenceTetrahedron";
        
        points_[0] = PointReferenceFactory::instance()->makePoint({0., 0., 0.});
        points_[1] = PointReferenceFactory::instance()->makePoint({1., 0., 0.});
        points_[2] = PointReferenceFactory::instance()->makePoint({0., 1., 0.});
        points_[3] = PointReferenceFactory::instance()->makePoint({0., 0., 1.});
        
        mappingsTriangleToTetrahedron_[0] = &MappingToRefTriangleToTetrahedron0::Instance();
        mappingsTriangleToTetrahedron_[1] = &MappingToRefTriangleToTetrahedron1::Instance();
        mappingsTriangleToTetrahedron_[2] = &MappingToRefTriangleToTetrahedron2::Instance();
        mappingsTriangleToTetrahedron_[3] = &MappingToRefTriangleToTetrahedron3::Instance();
        
        /// TODO: Implement MappingTetrahedronToTetrahedron.
        mappingsTetrahedronToTetrahedron_[0] = 0;
    }
    
    bool ReferenceTetrahedron::isInternalPoint(const PointReference& p) const
    {
        logger.assert(p.size()==3, "The dimension of the reference point is incorrect");
        return ((p[0] >= 0.) && (p[0] <= 1.) && (p[1] >= 0.) && (p[1] <= 1. - p[0]) && (p[2] >= 0.) && (p[2] <= 1. - p[0] - p[1]));
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
        logger(FATAL, "Tetrahedron to tetrahedron mappings do not exist.\n");
        return 0;
    }
    
    const MappingReferenceToReference*
    ReferenceTetrahedron::getCodim0MappingPtr(const std::size_t i) const
    {
        /// \TODO: Implement tetrahedron to tetrahedron mappings.
        logger(FATAL, "ERROR: Tetrahedron to tetrahedron mappings do not exist.\n");
        return 0;
    }
    
    // ================================== Codimension 1 ============================================
    
    std::vector<std::size_t> ReferenceTetrahedron::getCodim1EntityLocalIndices(const std::size_t faceIndex) const
    {
        logger.assert((faceIndex < 4), "ERROR: Index out of range. Tetrahedron has only 3 faces.\n");
        return std::vector<std::size_t>(localNodeIndexes_[faceIndex], localNodeIndexes_[faceIndex] + 3);
    }
    
    const ReferenceGeometry*
    ReferenceTetrahedron::getCodim1ReferenceGeometry(const std::size_t faceIndex) const
    {
        logger.assert((faceIndex < 4), "ERROR: Index out of range. Tetrahedron has only 3 faces.\n");
        return referenceGeometryCodim1Ptr_;
    }
    
    const MappingReferenceToReference*
    ReferenceTetrahedron::getCodim1MappingPtr(const std::size_t faceIndex) const
    {
        logger.assert((faceIndex < 4), "ERROR: Asked for a square point index larger than 3. There are only 4 nodes in a square.\n");
        return mappingsTriangleToTetrahedron_[faceIndex];
    }
    
    // ================================== Codimension 2 ============================================
    
    std::vector<std::size_t> ReferenceTetrahedron::getCodim2EntityLocalIndices(const std::size_t edgeIndex) const
    {
        logger.assert((edgeIndex < 6), "ERROR: Index out of range. Tetrahedron has only 6 edges.\n");
        return std::vector<std::size_t>(localNodesOnEdge_[edgeIndex], localNodesOnEdge_[edgeIndex] + 2);
    }
    
    const ReferenceGeometry*
    ReferenceTetrahedron::getCodim2ReferenceGeometry(const std::size_t edgeIndex) const
    {
        logger.assert((edgeIndex < 6), "ERROR: Index out of range. Tetrahedron has only 6 edges.\n");
        return referenceGeometryCodim2Ptr_;
    }
    
    const MappingReferenceToReference*
    ReferenceTetrahedron::getCodim2MappingPtr(const std::size_t faceIndex) const
    {
        /// \TODO: Implement line to tetrahedron mappings.
        logger(FATAL, "ERROR: Line to tetrahedron mappings do not exist.\n");
        return 0;
    }
    
    // ================================== Codimension 3 ============================================
    
    std::vector<std::size_t> ReferenceTetrahedron::getCodim3EntityLocalIndices(const std::size_t nodeIndex) const
    {
        logger.assert((nodeIndex < 4), "ERROR: Index out of range. Tetrahedron has only 4 nodes.\n");
        return std::vector<std::size_t>(1, nodeIndex);
    }

}


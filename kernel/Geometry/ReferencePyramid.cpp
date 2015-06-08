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

#include "ReferencePyramid.h"
#include "ReferenceTriangle.h"
#include "ReferenceSquare.h"
#include "ReferenceLine.h"
#include "Geometry/PointReference.h"
#include "Mappings/MappingToRefFaceToPyramid.h"

#include <cmath>

namespace Geometry
{
    
    std::size_t ReferencePyramid::localNodeIndexes_[5][4] = { {3, 4, 1, 2}, {3, 1, 0}, {2, 4, 0}, {1, 2, 0}, {4, 3, 0}, };
    
    std::size_t ReferencePyramid::localNodesOnEdge_[8][2] = { {0, 1}, {0, 2}, {0, 3}, {0, 4}, {1, 2}, {2, 4}, {4, 3}, {3, 1}, };
    
    ReferencePyramid::ReferencePyramid()
            : /// pyramid has three nodes 3D + 2
            ReferenceGeometry(5, 3, ReferenceGeometryType::PYRAMID, {0., 0., 1./4.}), referenceGeometryCodim1SquarePtr_(&ReferenceSquare::Instance()), referenceGeometryCodim1TrianglePtr_(&ReferenceTriangle::Instance()), referenceGeometryCodim2Ptr_(&ReferenceLine::Instance()), points_(5)
    {
        name = "ReferencePyramid";
        
        points_[0] = PointReferenceFactory<3>::instance()->makePoint({ 0.,  0.,  1.});
        points_[1] = PointReferenceFactory<3>::instance()->makePoint({-1., -1.,  0.});
        points_[2] = PointReferenceFactory<3>::instance()->makePoint({ 1., -1.,  0.});
        points_[3] = PointReferenceFactory<3>::instance()->makePoint({-1.,  1.,  0.});
        points_[4] = PointReferenceFactory<3>::instance()->makePoint({ 1.,  1.,  0.});
        center_ = PointReferenceFactory<3>::instance()->makePoint({0., 0., 1./4.});
        
        mappingsFaceToPyramid_[0] = &MappingToRefFaceToPyramid0::Instance();
        mappingsFaceToPyramid_[1] = &MappingToRefFaceToPyramid1::Instance();
        mappingsFaceToPyramid_[2] = &MappingToRefFaceToPyramid2::Instance();
        mappingsFaceToPyramid_[3] = &MappingToRefFaceToPyramid3::Instance();
        mappingsFaceToPyramid_[4] = &MappingToRefFaceToPyramid4::Instance();
        
        /// There are no MappingPyramidToPyramid. TODO: Implement.
        mappingsPyramidToPyramid_[0] = 0;
    }
    
    bool ReferencePyramid::isInternalPoint(const PointReference<3>& p) const
    {
        logger.assert(p.size()==3, "The reference point has the wrong dimension");
        return ((0. <= p[2]) && (1. >= p[2]) && (std::abs(p[0]) <= (1. - p[2])) && (std::abs(p[1]) <= (1. - p[2])));
    }
    
    std::ostream& operator<<(std::ostream& os, const ReferencePyramid& pyramid)
    {
        os << pyramid.getName() << " =( ";
        auto it = pyramid.points_.begin();
        auto end = pyramid.points_.end();
        
        for (; it != end; ++it)
        {
            os << (*it) << '\t';
        }
        os << ')' << std::endl;
        
        return os;
    }
    
    // ================================== Codimension 0 ============================================
    
    std::size_t ReferencePyramid::getCodim0MappingIndex(const ListOfIndexesT& list1, const ListOfIndexesT& list2) const
    {
        logger(FATAL, "ReferencePyramid::getCodim0MappingIndex: there are no Codim0 mappings for Pyramid.\n");
        return 0;
    }
    
    const MappingReferenceToReference<0>* ReferencePyramid::getCodim0MappingPtr(const std::size_t i) const
    {
        logger(FATAL, "ReferencePyramid::getCodim0MappingIndex: there are no Codim0 mappings for Pyramid.\n");
        return 0;
    }
    
    // ================================== Codimension 1 ============================================
    
    const MappingReferenceToReference<1>*
    ReferencePyramid::getCodim1MappingPtr(const std::size_t faceIndex) const
    {
        logger.assert((faceIndex < 5), "ReferencePyramid::getCodim1MappingPtr Index out of range. Only 5 faces for pyramid.\n");
        return mappingsFaceToPyramid_[faceIndex];
    }
    
    const ReferenceGeometry*
    ReferencePyramid::getCodim1ReferenceGeometry(const std::size_t faceIndex) const
    {
        logger.assert((faceIndex < 5), "ReferencePyramid::getCodim1ReferenceGeometry Index out of range. Only 5 faces for pyramid.\n");
        if (faceIndex == 0)
                return referenceGeometryCodim1SquarePtr_;
            else
                return referenceGeometryCodim1TrianglePtr_;
    }
    
    std::vector<std::size_t> ReferencePyramid::getCodim1EntityLocalIndices(const std::size_t faceIndex) const
    {
        logger.assert((faceIndex < 5), "ReferencePyramid::getCodim1EntityLocalIndices Index out of range. Only 5 faces for pyramid.\n");
        if (faceIndex == 0)
            {
                return std::vector<std::size_t>(localNodeIndexes_[faceIndex], localNodeIndexes_[faceIndex] + 4);
            }
            else
            {
                return std::vector<std::size_t>(localNodeIndexes_[faceIndex], localNodeIndexes_[faceIndex] + 3);
            }
    }
    
    // ================================== Codimension 2 ============================================
    
    const MappingReferenceToReference<2>*
    ReferencePyramid::getCodim2MappingPtr(const std::size_t edgeIndex) const
    {
        logger(ERROR, "There are no edge to pyramid mappings. \n");
        return 0;
    }
    
    const ReferenceGeometry*
    ReferencePyramid::getCodim2ReferenceGeometry(const std::size_t edgeIndex) const
    {
        logger.assert((edgeIndex < 8), "ReferencePyramid::getCodim2ReferenceGeometry Index out of range. Only 8 edges in pyramid.\n");
       return referenceGeometryCodim2Ptr_;
    }
    
    std::vector<std::size_t> ReferencePyramid::getCodim2EntityLocalIndices(const std::size_t edgeIndex) const
    {
        logger.assert((edgeIndex < 8), "ReferencePyramid::getCodim1EntityLocalIndices Index out of range. Only 8 edges in pyramid.\n");
        return std::vector<std::size_t>(localNodesOnEdge_[edgeIndex], localNodesOnEdge_[edgeIndex] + 2);
    }
    
    // ================================== Codimension 3 ============================================
    
    std::vector<std::size_t> ReferencePyramid::getCodim3EntityLocalIndices(const std::size_t nodeIndex) const
    {
        logger.assert((nodeIndex < 5), "ReferencePyramid::getCodim3EntityLocalIndices Index out of range. Pyramid has 5 nodes.\n");
        return std::vector<std::size_t>(1, nodeIndex);
    }
}

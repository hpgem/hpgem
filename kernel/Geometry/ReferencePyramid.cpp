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

#include "ReferencePyramid.hpp"

namespace Geometry
{

    int ReferencePyramid::localNodeIndexes_[5][4] =
    {
        { 3, 4, 1, 2 },
        { 3, 1, 0 },
        { 2, 4, 0 },
        { 1, 2, 0 },
        { 4, 3, 0 },
    };

    int ReferencePyramid::localNodesOnEdge_[8][2] =
    {
        { 0 , 1 },
        { 0 , 2 },
        { 0 , 3 },
        { 0 , 4 },
        { 1 , 2 },
        { 2 , 4 },
        { 4 , 3 },
        { 3 , 1 },
    };

    ReferencePyramid::ReferencePyramid():/// pyramid has three nodes 3D + 2
        ReferenceGeometry(ThreeD+2,3, PYRAMID),
        referenceGeometryCodim1TrianglePtr_(&ReferenceTriangle::Instance()),
        referenceGeometryCodim1SquarePtr_(&ReferenceSquare::Instance()),
        referenceGeometryCodim2Ptr_(&ReferenceLine::Instance())
    {
        PointReferenceT p1(3), p2(3), p3(3), p4(3), p5(3);

        p1[0] = +0.0; p1[1] = +0.0; p1[2] = +1.0;
        p2[0] = -1.0; p2[1] = -1.0; p2[2] = +0.0;
        p3[0] = +1.0; p3[1] = -1.0; p3[2] = +0.0;
        p4[0] = -1.0; p4[1] = +1.0; p4[2] = +0.0;
        p5[0] = +1.0; p5[1] = +1.0; p5[2] = +0.0;

        points_[0] = p1;
        points_[1] = p2;
        points_[2] = p3;
        points_[3] = p4;
        points_[4] = p5;

        mappingsFaceToPyramid_[0] = &MappingToRefFaceToPyramid0::Instance();
        mappingsFaceToPyramid_[1] = &MappingToRefFaceToPyramid1::Instance();
        mappingsFaceToPyramid_[2] = &MappingToRefFaceToPyramid2::Instance();
        mappingsFaceToPyramid_[3] = &MappingToRefFaceToPyramid3::Instance();
        mappingsFaceToPyramid_[4] = &MappingToRefFaceToPyramid4::Instance();

        /// There are no MappingPyramidToPyramid. TODO: Implement.
        mappingsPyramidToPyramid_[0] = 0;
    }

    ReferencePyramid::ReferencePyramid(const ReferencePyramid& copy):
        ReferenceGeometry(copy),
        referenceGeometryCodim1TrianglePtr_(copy.referenceGeometryCodim1TrianglePtr_),
        referenceGeometryCodim1SquarePtr_(copy.referenceGeometryCodim1SquarePtr_),
        referenceGeometryCodim2Ptr_(copy.referenceGeometryCodim2Ptr_) { }
    
    bool ReferencePyramid::isInternalPoint(const PointReferenceT& p) const
    {
        return ((0. <= p[2]) && (1. >= p[2]) &&
               (std::abs(p[0]) <= (1. - p[2])) &&
               (std::abs(p[1]) <= (1. - p[2])));
    }
    
    void ReferencePyramid::getCenter(PointReferenceT& p) const
    {
        p[0] = 0.;
        p[1] = 0.;
        p[2] = 1. / 4.;
    }
    
    void ReferencePyramid::getNode(const IndexT& nodeIndex, PointReferenceT& point) const
    {
        if (nodeIndex < 5)
        {
            point = points_[nodeIndex];
        }
        else
        {
            throw "ERROR: ReferencePyramid::getNode Index out of range. Only 5 faces for pyramid.";
        }
    }

    std::ostream& operator<<(std::ostream& os, const ReferencePyramid& pyramid)
    {
        os << pyramid.getName()<<" =( ";
        ReferencePyramid::const_iterator it = pyramid.points_.begin();
        ReferencePyramid::const_iterator end = pyramid.points_.end();

        for ( ; it != end; ++it)
        {
            os << (*it) << '\t';
        }
        os <<')'<<std::endl;

        return os;
    }

    // ================================== Codimension 0 ============================================

    int ReferencePyramid::getCodim0MappingIndex(const ListOfIndexesT& list1, const ListOfIndexesT& list2) const
    {
        throw "ReferencePyramid::getCodim0MappingIndex: there are no Codim0 mappings for Pyramid.";
    }

    const MappingReferenceToReference* ReferencePyramid::getCodim0MappingPtr(const IndexT i) const
    {
        throw "ReferencePyramid::getCodim0MappingIndex: there are no Codim0 mappings for Pyramid.";
    }

    // ================================== Codimension 1 ============================================

    const MappingReferenceToReference*
    ReferencePyramid::getCodim1MappingPtr(const IndexT faceIndex) const
    {
        if (faceIndex < 5)
        {
            return mappingsFaceToPyramid_[faceIndex];
        }
        else
        {
            throw "ReferencePyramid::getCodim1MappingPtr Index out of range. Only 5 faces for pyramid.";
        }
    }

    const ReferenceGeometry*
    ReferencePyramid::getCodim1ReferenceGeometry(const IndexT faceIndex) const
    {
        if (faceIndex < 5)
        {
            if (faceIndex == 0)
                return referenceGeometryCodim1SquarePtr_;
            else
                return referenceGeometryCodim1TrianglePtr_;
        }
        else
        {
            throw "ReferencePyramid::getCodim1ReferenceGeometry Index out of range. Only 5 faces for pyramid.";
        }
    }

    void ReferencePyramid::getCodim1EntityLocalIndices(const IndexT faceIndex, ListOfIndexesT& faceNodesLocal) const
    {
        if (faceIndex < 5)
        {
            if (faceIndex == 0)
            {
                faceNodesLocal.resize(4); // 2 nodes per edge
                faceNodesLocal[0] = localNodeIndexes_[faceIndex][0];
                faceNodesLocal[1] = localNodeIndexes_[faceIndex][1];
                faceNodesLocal[2] = localNodeIndexes_[faceIndex][2];
                faceNodesLocal[3] = localNodeIndexes_[faceIndex][3];
            }
            else
            {
                faceNodesLocal.resize(3); // 2 nodes per edge
                faceNodesLocal[0] = localNodeIndexes_[faceIndex][0];
                faceNodesLocal[1] = localNodeIndexes_[faceIndex][1];
                faceNodesLocal[2] = localNodeIndexes_[faceIndex][2];
            }
        }
        else
        {
            throw "ReferencePyramid::getCodim1EntityLocalIndices Index out of range. Only 5 faces for pyramid.";
        }

    }

    // ================================== Codimension 2 ============================================

    const MappingReferenceToReference*
    ReferencePyramid::getCodim2MappingPtr(const IndexT edgeIndex) const
    {
        throw "ERROR: There are no edge to pyramid mappings.";
    }

    const ReferenceGeometry*
    ReferencePyramid::getCodim2ReferenceGeometry(const IndexT edgeIndex) const
    {
        if (edgeIndex < 8)
        {
            return referenceGeometryCodim2Ptr_;
        }
        else
        {
            throw "ReferencePyramid::getCodim2ReferenceGeometry Index out of range. Only 8 edges in pyramid.";
        }
    }

    void ReferencePyramid::getCodim2EntityLocalIndices(const IndexT edgeIndex, ListOfIndexesT& faceNodesLocal) const
    {
        if (edgeIndex < 8)
        {
            faceNodesLocal.resize(2); // 2 nodes per edge
            faceNodesLocal[0] = localNodesOnEdge_[edgeIndex][0];
            faceNodesLocal[1] = localNodesOnEdge_[edgeIndex][1];
        }
        else
        {
            throw "ReferencePyramid::getCodim1EntityLocalIndices Index out of range. Only 8 edges in pyramid.";
        }
    }

    // ================================== Codimension 3 ============================================

    void ReferencePyramid::
    getCodim3EntityLocalIndices(const IndexT nodeIndex, ListOfIndexesT& nodeNodesLocal) const
    {
        if (nodeIndex < 5)
        {
            nodeNodesLocal.resize(1);
            nodeNodesLocal[0] = nodeIndex;
        }
        else
        {
            throw "ReferencePyramid::getCodim3EntityLocalIndices Index out of range. Pyramid has 5 nodes.";
        }
    }


    // ================================== Quadrature rules =====================================

    /// Add a quadrature rule into the list of valid quadrature rules for this geometry.
    void ReferencePyramid::addGaussQuadratureRule(QuadratureRules::GaussQuadratureRule* const qr)
    {
        std::list<QuadratureRules::GaussQuadratureRule*>::iterator it = lstGaussQuadratureRules_.begin();
        while (it != lstGaussQuadratureRules_.end())
        {
          if ((*it)->order() < qr->order()) ++it;
          else break;
        }
        lstGaussQuadratureRules_.insert(it,qr);
    }

    /// Get a valid quadrature for this geometry.
    QuadratureRules::GaussQuadratureRule* const ReferencePyramid::getGaussQuadratureRule(int order) const
    {
        for (std::list<QuadratureRules::GaussQuadratureRule*>::const_iterator it = lstGaussQuadratureRules_.begin();
              it != lstGaussQuadratureRules_.end(); ++it)
          if ((*it)->order() >= order) return *it;

        return NULL;
    }
};

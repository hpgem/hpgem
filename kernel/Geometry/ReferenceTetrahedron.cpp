/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, Univesity of Twenete
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "ReferenceTetrahedron.hpp"
#include "ReferenceLine.hpp"
#include "ReferenceTriangle.hpp"

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
     int ReferenceTetrahedron::localNodeIndexes_[4][3] =
     {
         { 0, 3, 2 },
         { 0, 1, 3 },
         { 0, 2, 1 },
         { 1, 2, 3 }
     };

     int ReferenceTetrahedron::localNodesOnEdge_[6][2] =
     {
         { 0, 1 },
         { 0, 2 },
         { 0, 3 },
         { 2, 3 },
         { 1, 3 },
         { 1, 2 },
     };

    ReferenceTetrahedron::ReferenceTetrahedron():/// Tetrahedron has four nodes 3D + 1
        ReferenceGeometry(ThreeD+1,3,TETRAHEDRON),
        referenceGeometryCodim1Ptr_(&ReferenceTriangle::Instance()),
        referenceGeometryCodim2Ptr_(&ReferenceLine::Instance())
    {
        PointReferenceT p1(3), p2(3), p3(3), p4(3);
        
        p1[0] = +0.0; p1[1] = +0.0; p1[2] = +0.0;
        p2[0] = +1.0; p2[1] = +0.0; p2[2] = +0.0;
        p3[0] = +0.0; p3[1] = +1.0; p3[2] = +0.0;
        p4[0] = +0.0; p4[1] = +0.0; p4[2] = +1.0;
        
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
    
    ReferenceTetrahedron::ReferenceTetrahedron(const ReferenceTetrahedron& copy):
        ReferenceGeometry(copy),
        referenceGeometryCodim1Ptr_(copy.referenceGeometryCodim1Ptr_),
        referenceGeometryCodim2Ptr_(copy.referenceGeometryCodim2Ptr_)
    {
    }
    
    bool ReferenceTetrahedron::isInternalPoint(const PointReferenceT& p) const
    {
        return ((p[0] >= 0.) && (p[0] <= 1.) && (p[1] >= 0.) && (p[1] <= 1. - p[0]) &&
                (p[2] >= 0.) && (p[2] <= 1. - p[0] - p[1]));
    }
    
    void ReferenceTetrahedron::getCenter(PointReferenceT& p) const
    {
        p[0] = p[1] = p[2] = 1. / 4.;
    }
    
    void ReferenceTetrahedron::getNode(const IndexT& i, PointReferenceT& point) const
    {
        point = points_[i];
    }
    
    std::ostream& operator<<(std::ostream& os, const ReferenceTetrahedron& tetra)
    {
        os <<tetra.getName()<<" =( ";
        ReferenceTetrahedron::const_iterator it = tetra.points_.begin();
        ReferenceTetrahedron::const_iterator end = tetra.points_.end();

        for ( ; it != end; ++it)
        {
            os << (*it) << '\t';
        }
        os <<')'<<std::endl;

        return os;
    }

    // ================================== Codimension 0 ============================================

    int ReferenceTetrahedron::
    getCodim0MappingIndex(const ListOfIndexesT& list1, const ListOfIndexesT& list2) const
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

    void ReferenceTetrahedron::
    getCodim1EntityLocalIndices(const IndexT faceIndex, ListOfIndexesT& faceNodesLocal) const
    {
        if (faceIndex < 4)
        {
            faceNodesLocal.resize(3); // 3 nodes per face
            faceNodesLocal[0] = (IndexT) localNodeIndexes_[faceIndex][0];
            faceNodesLocal[1] = (IndexT) localNodeIndexes_[faceIndex][1];
            faceNodesLocal[2] = (IndexT) localNodeIndexes_[faceIndex][2];
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

    void ReferenceTetrahedron::
    getCodim2EntityLocalIndices(const IndexT edgeIndex, ListOfIndexesT& edgeNodesLocal) const
    {
        if (edgeIndex < 6)
        {
            edgeNodesLocal.resize(2); // 2 nodes per edge
            edgeNodesLocal[0] = (IndexT) localNodesOnEdge_[edgeIndex][0];
            edgeNodesLocal[1] = (IndexT) localNodesOnEdge_[edgeIndex][1];
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

    void ReferenceTetrahedron::
    getCodim3EntityLocalIndices(const IndexT nodeIndex, ListOfIndexesT& nodeNodesLocal) const
    {
        if (nodeIndex < 4)
        {
            nodeNodesLocal.resize(1); // 2 nodes per edge
            nodeNodesLocal[0] = nodeIndex;
        }
        else
        {
            throw "ERROR: Index out of range. Tetrahedron has only 4 nodes.";
        }
    }


    // ================================== Quadrature rules =====================================

    /// Add a quadrature rule into the list of valid quadrature rules for this geometry.
    void ReferenceTetrahedron::addGaussQuadratureRule(QuadratureRules::GaussQuadratureRule* const qr)
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
    QuadratureRules::GaussQuadratureRule* const ReferenceTetrahedron::getGaussQuadratureRule(int order) const
    {
        for (std::list<QuadratureRules::GaussQuadratureRule*>::const_iterator it = lstGaussQuadratureRules_.begin();
              it != lstGaussQuadratureRules_.end(); ++it)
          if ((*it)->order() >= order) return *it;

        return NULL;
    }
            
};


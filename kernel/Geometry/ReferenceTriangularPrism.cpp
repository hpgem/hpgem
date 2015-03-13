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
#include "ReferenceTriangularPrism.h"
#include "ReferenceTriangle.h"
#include "ReferenceSquare.h"
#include "ReferenceLine.h"
#include "Geometry/PointReference.h"
#include "Mappings/MappingToRefFaceToTriangularPrism.h"
#include "LinearAlgebra/Matrix.h"

namespace Geometry
{
    std::size_t ReferenceTriangularPrism::localNodeIndexes_[5][4] = { {0, 2, 1}, {3, 4, 5}, {2, 0, 5, 3}, {0, 1, 3, 4}, {1, 2, 4, 5}};
    
    std::size_t ReferenceTriangularPrism::localNodesOnEdge_[9][2] = { {0, 1}, {0, 2}, {1, 2}, {3, 4}, {3, 5}, {4, 5}, {0, 3}, {1, 4}, {2, 5}};
    
    ReferenceTriangularPrism::ReferenceTriangularPrism()
            : ReferenceGeometry(6, 3, TRIANGULARPRISM), referenceGeometryCodim1TrianglePtr_(&ReferenceTriangle::Instance()), referenceGeometryCodim1SquarePtr_(&ReferenceSquare::Instance()), referenceGeometryCodim2Ptr_(&ReferenceLine::Instance())
    {
        PointReference p1(3), p2(3), p3(3), p4(3), p5(3), p6(3);
        
        p1[0] = +0.0;
        p1[1] = +0.0;
        p1[2] = -1.0;
        p2[0] = +1.0;
        p2[1] = +0.0;
        p2[2] = -1.0;
        p3[0] = +0.0;
        p3[1] = +1.0;
        p3[2] = -1.0;
        p4[0] = +0.0;
        p4[1] = +0.0;
        p4[2] = +1.0;
        p5[0] = +1.0;
        p5[1] = +0.0;
        p5[2] = +1.0;
        p6[0] = +0.0;
        p6[1] = +1.0;
        p6[2] = +1.0;
        
        points_[0] = p1;
        points_[1] = p2;
        points_[2] = p3;
        points_[3] = p4;
        points_[4] = p5;
        points_[5] = p6;
        
        /// Mappings between triangular prisms are not implemented
        mappingsTriangularPrismToTriangularPrism_[0] = 0;
        
        mappingsFaceToTriangularPrism_[0] = &MappingToRefFaceToTriangularPrism0::Instance();
        mappingsFaceToTriangularPrism_[1] = &MappingToRefFaceToTriangularPrism1::Instance();
        mappingsFaceToTriangularPrism_[2] = &MappingToRefFaceToTriangularPrism2::Instance();
        mappingsFaceToTriangularPrism_[3] = &MappingToRefFaceToTriangularPrism3::Instance();
        mappingsFaceToTriangularPrism_[4] = &MappingToRefFaceToTriangularPrism4::Instance();
    }
    
    ReferenceTriangularPrism::ReferenceTriangularPrism(const ReferenceTriangularPrism& copy)
            : ReferenceGeometry(copy), referenceGeometryCodim1TrianglePtr_(copy.referenceGeometryCodim1TrianglePtr_), referenceGeometryCodim1SquarePtr_(copy.referenceGeometryCodim1SquarePtr_), referenceGeometryCodim2Ptr_(copy.referenceGeometryCodim2Ptr_)
    {
    }
    
    bool ReferenceTriangularPrism::isInternalPoint(const PointReference& p) const
    {
        logger.assert(p.size()==3, "The dimension of the reference point is incorrect");
        return ((-1. <= p[2]) && (1. >= p[2]) && (p[0] >= 0.) && (p[0] <= 1.) && (p[1] >= 0.) && (p[1] <= 1. - p[0]));
    }
    
    PointReference ReferenceTriangularPrism::getCenter() const
    {
        PointReference p(3);
        p[0] = 1. / 3.;
        p[1] = 1. / 3.;
        p[2] = 0.;
        return p;
    }
    
    const PointReference& ReferenceTriangularPrism::getNode(const std::size_t& i) const
    {
        logger.assert(i<getNumberOfNodes(), "Asked for node %, but there are only % nodes", i, getNumberOfNodes());
        return points_[i];
    }
    
    std::ostream& operator<<(std::ostream& os, const ReferenceTriangularPrism& prism)
    {
        os << prism.getName() << " = ( ";
        ReferenceTriangularPrism::const_iterator it = prism.points_.begin();
        ReferenceTriangularPrism::const_iterator end = prism.points_.end();
        
        for (; it != end; ++it)
        {
            os << (*it) << '\t';
        }
        os << ')' << std::endl;
        
        return os;
    }
    
    // ================================== Codimension 0 ============================================
    
    std::size_t ReferenceTriangularPrism::getCodim0MappingIndex(const ListOfIndexesT& list1, const ListOfIndexesT& list2) const
    {
        /// \TODO: Implement tetrahedron to tetrahedron mappings.
        logger(FATAL, "ReferenceTriangularPrism::getCodim0MappingIndex: T.p to t.p mappings do not exist.\n");
        return 0;
    }
    
    const MappingReferenceToReference*
    ReferenceTriangularPrism::getCodim0MappingPtr(const std::size_t i) const
    {
        /// \TODO: Implement tetrahedron to tetrahedron mappings.
        logger(FATAL, "ReferenceTetrahedron::getCodim0MappingPtr: T.p to T.p mappings do not exist.\n");
        return 0;
    }
    
    // ================================== Codimension 1 ============================================
    
    std::vector<std::size_t> ReferenceTriangularPrism::getCodim1EntityLocalIndices(const std::size_t faceIndex) const
    {
        if (faceIndex < 2)
        {
            return std::vector<std::size_t>(localNodeIndexes_[faceIndex], localNodeIndexes_[faceIndex] + 3);
        }
        else if (faceIndex < 5)
        {
            return std::vector<std::size_t>(localNodeIndexes_[faceIndex], localNodeIndexes_[faceIndex] + 4);
        }
        else
        {
            logger(ERROR, "ReferenceTriangularPrism::getCodim1EntityLocalIndices: Index out of range. T.p has 5 faces.\n");
        }
        std::vector<std::size_t> dummy(1);
        return dummy;
    }
    
    const ReferenceGeometry*
    ReferenceTriangularPrism::getCodim1ReferenceGeometry(const std::size_t faceIndex) const
    {
        if (faceIndex < 2)
        {
            return referenceGeometryCodim1TrianglePtr_;
        }
        else if (faceIndex < 5)
        {
            return referenceGeometryCodim1SquarePtr_;
        }
        else
        {
            logger(ERROR, "ReferenceTriangularPrism::getCodim1ReferenceGeometry: Index out of range. T.p has 5 faces.\n");
        }
        return 0;
    }
    
    const MappingReferenceToReference*
    ReferenceTriangularPrism::getCodim1MappingPtr(const std::size_t faceIndex) const
    {
        logger.assert((faceIndex < 5), "Asked for a square point index larger than 3. There are only 4 nodes in a square!.\n");
        return mappingsFaceToTriangularPrism_[faceIndex];
    }
    
    // ================================== Codimension 2 ============================================
    
    std::vector<std::size_t> ReferenceTriangularPrism::getCodim2EntityLocalIndices(const std::size_t edgeIndex) const
    {
        logger.assert((edgeIndex < 9), "ReferenceTriangularPrism::getCodim2EntityLocalIndices Index out of range. T.p has only 9 edges.\n");
        return std::vector<std::size_t>(localNodesOnEdge_[edgeIndex], localNodesOnEdge_[edgeIndex] + 2);
    }
    
    const ReferenceGeometry*
    ReferenceTriangularPrism::getCodim2ReferenceGeometry(const std::size_t edgeIndex) const
    {
        logger.assert((edgeIndex < 9), "ReferenceTriangularPrism::getCodim2ReferenceGeometry Index out of range. T.p has only 9 edges.\n");
        return referenceGeometryCodim2Ptr_;
    }
    
    const MappingReferenceToReference*
    ReferenceTriangularPrism::getCodim2MappingPtr(const std::size_t faceIndex) const
    {
        /// \TODO: Implement line to t.p. mappings.
        logger(FATAL, "ReferenceTriangularPrism::getCodim2MappingPtr: Line to TP mappings do not exist.\n");
        return 0;
    }
    
    // ================================== Codimension 3 ============================================
    
    std::vector<std::size_t> ReferenceTriangularPrism::getCodim3EntityLocalIndices(const std::size_t nodeIndex) const
    {
        logger.assert((nodeIndex < 6), "ReferenceTriangularPrism::Index out of range. TP has only 6 nodes.\n");
        return std::vector<std::size_t>(1, nodeIndex);
    }
    
    // =============================== Refinement mappings =====================================
    
    void ReferenceTriangularPrism::refinementTransform(int refineType, std::size_t subElementIdx, const PointReference& p, PointReference& pMap) const
    {
        switch (refineType)
        {
            case 0:
            case 1:
            case 2:
            case 3:
                break;
                
            case 4:
                switch (subElementIdx)
                {
                    case 0:
                        pMap[0] = (p[0]) / 2.;
                        pMap[1] = (p[1]) / 2.;
                        pMap[2] = p[2];
                        break;
                    case 1:
                        pMap[0] = (p[0] + 1) / 2.;
                        pMap[1] = (p[1]) / 2.;
                        pMap[2] = p[2];
                        break;
                    case 2:
                        pMap[0] = (p[0]) / 2.;
                        pMap[1] = (p[1] + 1) / 2.;
                        pMap[2] = p[2];
                        break;
                    case 3:
                        pMap[0] = (1. - p[1]) / 2.;
                        pMap[1] = (1. - p[0]) / 2.;
                        pMap[2] = p[2];
                        break;
                }
                break;
                
            case 5: //????????????????????????
            case 6: //????????????????????????
                switch (subElementIdx)
                {
                    case 0:
                        pMap[0] = (1. + p[0]) / 4.;
                        pMap[1] = (3. - p[0] + 3. * p[1] - p[0] * p[1]) / 8.;
                        pMap[2] = p[2];
                        break;
                    case 1:
                        pMap[0] = (p[0] + 1) / 2.;
                        pMap[1] = (p[1]) / 2.;
                        pMap[2] = p[2];
                        break;
                }
                break;
                
            case 7: //????????????????????????
                switch (subElementIdx)
                {
                    case 0:
                        pMap[0] = (3. - p[1] + 3. * p[0] - p[0] * p[1]) / 8.;
                        pMap[1] = (1. + p[1]) / 4.;
                        pMap[2] = p[2];
                        break;
                    case 1:
                        pMap[0] = (p[0]) / 2.;
                        pMap[1] = (p[1] + 1) / 2.;
                        pMap[2] = p[2];
                        break;
                }
                break;
                
            case 20:
                // We use this for data transfer between a coarse element, which is
                // a hexahedron, and its subelements, which are three triangular-prisms.
                // The hexahedron itself is a subelement of a coarser triangular-prism
                // which element is refined in x_0 direction.
                switch (subElementIdx)
                {
                    case 0:
                        pMap[0] = -1. + p[0] * (1. - (-1.)) + p[1] * (-1. - (-1.));
                        pMap[1] = -1. + p[0] * (-1. - (-1.)) + p[1] * (0. - (-1.));
                        pMap[2] = p[2];
                        break;
                    case 1:
                        pMap[0] = -1. + p[0] * (1. - (-1.)) + p[1] * (-1. - (-1.));
                        pMap[1] = 0. + p[0] * (1. - (0.)) + p[1] * (1. - (0.));
                        pMap[2] = p[2];
                        break;
                    case 2:
                        pMap[0] = 1. + p[0] * (-1. - (1.)) + p[1] * (1. - (1.));
                        pMap[1] = 1. + p[0] * (0. - (1.)) + p[1] * (-1. - (1.));
                        pMap[2] = p[2];
                        break;
                }
                break;
                
            case 21:
                // We use this for data transfer between a coarse element, which is
                // a hexahedron, and its subelements, which are three triangular-prisms.
                // The hexahedron itself is a subelement of a coarser triangular-prism
                // element which is refined in x_1 direction.
                switch (subElementIdx)
                {
                    case 0:
                        pMap[0] = -1. + p[0] * (0. - (-1.)) + p[1] * (-1. - (-1.));
                        pMap[1] = -1. + p[0] * (-1. - (-1.)) + p[1] * (1. - (-1.));
                        pMap[2] = p[2];
                        break;
                    case 1:
                        pMap[0] = 0. + p[0] * (1. - (0.)) + p[1] * (1. - (0.));
                        pMap[1] = -1. + p[0] * (-1. - (-1.)) + p[1] * (1. - (-1.));
                        pMap[2] = p[2];
                        break;
                    case 2:
                        pMap[0] = 1. + p[0] * (-1. - (1.)) + p[1] * (0. - (1.));
                        pMap[1] = 1. + p[0] * (1. - (1.)) + p[1] * (-1. - (1.));
                        pMap[2] = p[2];
                        break;
                }
                break;
                
            case 22:
                // We use this for data transfer between a coarse element, which is
                // a hexahedron, and its triangular subelement. The coarse element 
                // is split into two hexahedrons and one triangular-prism.
                // The parent hexahedron itself is a subsubelement of a coarser 
                // triangular-prism element which is refined in x_0 direction.
                switch (subElementIdx)
                {
                    case 2:
                        pMap[0] = 1. + p[0] * (-1. - (1.)) + p[1] * (1. - (1.));
                        pMap[1] = 1. / 3. + p[0] * (0. - (1. / 3.)) + p[1] * (-1. / 3. - (1. / 3.));
                        pMap[2] = p[2];
                        break;
                        
                    case 0:
                    case 1:
                        logger(ERROR, "Subelement-0 and -1 are hexahedrons. \n");
                        break;
                }
                break;
                
            case 23:
                // We use this for data transfer between a coarse element, which is
                // a hexahedron, and its triangular subelement. The coarse element 
                // is split into two hexahedrons and one triangular-prism.
                // The parent hexahedron itself is a subsubelement of a coarser 
                // triangular-prism element which is refined in x_1 direction.
                switch (subElementIdx)
                {
                    case 2:
                        pMap[0] = 1. / 3. + p[0] * (-1. / 3. - (1. / 3.)) + p[1] * (0. - (1. / 3.));
                        pMap[1] = 1. + p[0] * (1. - (1.)) + p[1] * (-1. - (1.));
                        pMap[2] = p[2];
                        break;
                        
                    case 0:
                    case 1:
                        logger(ERROR, "subelement-0 and -1 are hexahedrons. \n");
                        break;
                }
                break;
                
            case 24:
                // We use this for data transfer between a coarse element, which is
                // a hexahedron, and its triangular subelement. The coarse element 
                // is split into one hexahedron and two triangular-prisms.
                // The parent hexahedron itself is a subsubelement of a coarser 
                // triangular-prism element which is refined in x_0 direction.
                switch (subElementIdx)
                {
                    case 0:
                        pMap[0] = -1. + p[0] * (1. - (-1.)) + p[1] * (-1. - (-1.));
                        pMap[1] = -1. + p[0] * (-1. - (-1.)) + p[1] * (-1. / 3. - (-1.));
                        pMap[2] = p[2];
                        break;
                        
                    case 1:
                        pMap[0] = -1. + p[0] * (1. - (-1.)) + p[1] * (-1. - (-1.));
                        pMap[1] = 1. / 3. + p[0] * (1. - (1. / 3.)) + p[1] * (1. - (1. / 3.));
                        pMap[2] = p[2];
                        break;
                        
                    case 2:
                        logger(ERROR, "subelement-2 is a hexahedron.\n");
                        break;
                }
                break;
                
            case 25:
                // We use this for data transfer between a coarse element, which is
                // a hexahedron, and its triangular subelement. The coarse element 
                // is split into one hexahedron and two triangular-prisms.
                // The parent hexahedron itself is a subsubelement of a coarser 
                // triangular-prism element which is refined in x_1 direction.
                switch (subElementIdx)
                {
                    case 0:
                        pMap[0] = -1. + p[0] * (-1. / 3. - (-1.)) + p[1] * (-1. - (-1.));
                        pMap[1] = -1. + p[0] * (-1. - (-1.)) + p[1] * (1. - (-1.));
                        pMap[2] = p[2];
                        break;
                        
                    case 1:
                        pMap[0] = 1. / 3. + p[0] * (1. - (1. / 3.)) + p[1] * (1. - (1. / 3.));
                        pMap[1] = -1. + p[0] * (-1. - (-1.)) + p[1] * (1. - (-1.));
                        pMap[2] = p[2];
                        break;
                        
                    case 2:
                        logger(ERROR, "Subelement-2 is a hexahedron.\n");
                        break;
                }
                break;
                
            default:
                pMap = p;
                break;
        }
    } // end of refinementTransform
    
    void ReferenceTriangularPrism::getRefinementMappingMatrixL(int refineType, std::size_t subElementIdx, LinearAlgebra::Matrix& Q) const
    {
        Q.resize(4, 4);
        Q = 0.;
        Q(3, 3) = 1.;
        switch (refineType)
        {
            case 0:
                switch (subElementIdx)
                {
                    case 0:
                        Q(0, 0) = .5;
                        Q(0, 3) = -.5;
                        Q(1, 1) = 1.;
                        Q(2, 2) = 1.;
                        break;
                    case 1:
                        Q(0, 0) = .5;
                        Q(0, 3) = .5;
                        Q(1, 1) = 1.;
                        Q(2, 2) = 1.;
                        break;
                }
                break;
                
            case 1:
                switch (subElementIdx)
                {
                    case 0:
                        Q(0, 0) = 1.;
                        Q(1, 1) = .5;
                        Q(1, 3) = -.5;
                        Q(2, 2) = 1.;
                        break;
                    case 1:
                        Q(0, 0) = 1.;
                        Q(1, 1) = .5;
                        Q(1, 3) = .5;
                        Q(2, 2) = 1.;
                        break;
                }
                break;
                
            case 2:
                switch (subElementIdx)
                {
                    case 0:
                        Q(0, 0) = 1.;
                        Q(1, 1) = 1.;
                        Q(2, 2) = .5;
                        Q(2, 3) = -.5;
                        break;
                    case 1:
                        Q(0, 0) = 1.;
                        Q(1, 1) = 1.;
                        Q(2, 2) = .5;
                        Q(2, 3) = .5;
                        break;
                }
                break;
                
            case 3:
                switch (subElementIdx)
                {
                    case 0:
                        Q(0, 0) = .5;
                        Q(0, 3) = -.5;
                        Q(1, 1) = .5;
                        Q(1, 3) = -.5;
                        Q(2, 2) = 1.;
                        break;
                    case 1:
                        Q(0, 0) = .5;
                        Q(0, 3) = .5;
                        Q(1, 1) = .5;
                        Q(1, 3) = -.5;
                        Q(2, 2) = 1.;
                        break;
                    case 2:
                        Q(0, 0) = .5;
                        Q(0, 3) = -.5;
                        Q(1, 1) = .5;
                        Q(1, 3) = .5;
                        Q(2, 2) = 1.;
                        break;
                    case 3:
                        Q(0, 0) = .5;
                        Q(0, 3) = .5;
                        Q(1, 1) = .5;
                        Q(1, 3) = .5;
                        Q(2, 2) = 1.;
                        break;
                }
                break;
                
            case 4:
                switch (subElementIdx)
                {
                    case 0:
                        Q(0, 0) = .5;
                        Q(0, 3) = -.5;
                        Q(1, 1) = 1.;
                        Q(2, 2) = .5;
                        Q(2, 3) = -.5;
                        break;
                    case 1:
                        Q(0, 0) = .5;
                        Q(0, 3) = .5;
                        Q(1, 1) = 1.;
                        Q(2, 2) = .5;
                        Q(2, 3) = -.5;
                        break;
                    case 2:
                        Q(0, 0) = .5;
                        Q(0, 3) = -.5;
                        Q(1, 1) = 1.;
                        Q(2, 2) = .5;
                        Q(2, 3) = .5;
                        break;
                    case 3:
                        Q(0, 0) = .5;
                        Q(0, 3) = .5;
                        Q(1, 1) = 1.;
                        Q(2, 2) = .5;
                        Q(2, 3) = .5;
                        break;
                }
                break;
                
            case 5:
                switch (subElementIdx)
                {
                    case 0:
                        Q(0, 0) = 1.;
                        Q(1, 1) = .5;
                        Q(1, 3) = -.5;
                        Q(2, 2) = .5;
                        Q(2, 3) = -.5;
                        break;
                    case 1:
                        Q(0, 0) = 1.;
                        Q(1, 1) = .5;
                        Q(1, 3) = .5;
                        Q(2, 2) = .5;
                        Q(2, 3) = -.5;
                        break;
                    case 2:
                        Q(0, 0) = 1.;
                        Q(1, 1) = .5;
                        Q(1, 3) = -.5;
                        Q(2, 2) = .5;
                        Q(2, 3) = .5;
                        break;
                    case 3:
                        Q(0, 0) = 1.;
                        Q(1, 1) = .5;
                        Q(1, 3) = .5;
                        Q(2, 2) = .5;
                        Q(2, 3) = .5;
                        break;
                }
                break;
                
            case 6:
                switch (subElementIdx)
                {
                    case 0:
                        Q(0, 0) = .5;
                        Q(0, 3) = -.5;
                        Q(1, 1) = .5;
                        Q(1, 3) = -.5;
                        Q(2, 2) = .5;
                        Q(2, 3) = -.5;
                        break;
                    case 1:
                        Q(0, 0) = .5;
                        Q(0, 3) = .5;
                        Q(1, 1) = .5;
                        Q(1, 3) = -.5;
                        Q(2, 2) = .5;
                        Q(2, 3) = -.5;
                        break;
                    case 2:
                        Q(0, 0) = .5;
                        Q(0, 3) = -.5;
                        Q(1, 1) = .5;
                        Q(1, 3) = .5;
                        Q(2, 2) = .5;
                        Q(2, 3) = -.5;
                        break;
                    case 3:
                        Q(0, 0) = .5;
                        Q(0, 3) = .5;
                        Q(1, 1) = .5;
                        Q(1, 3) = .5;
                        Q(2, 2) = .5;
                        Q(2, 3) = -.5;
                        break;
                    case 4:
                        Q(0, 0) = .5;
                        Q(0, 3) = -.5;
                        Q(1, 1) = .5;
                        Q(1, 3) = -.5;
                        Q(2, 2) = .5;
                        Q(2, 3) = .5;
                        break;
                    case 5:
                        Q(0, 0) = .5;
                        Q(0, 3) = .5;
                        Q(1, 1) = .5;
                        Q(1, 3) = -.5;
                        Q(2, 2) = .5;
                        Q(2, 3) = .5;
                        break;
                    case 6:
                        Q(0, 0) = .5;
                        Q(0, 3) = -.5;
                        Q(1, 1) = .5;
                        Q(1, 3) = .5;
                        Q(2, 2) = .5;
                        Q(2, 3) = .5;
                        break;
                    case 7:
                        Q(0, 0) = .5;
                        Q(0, 3) = .5;
                        Q(1, 1) = .5;
                        Q(1, 3) = .5;
                        Q(2, 2) = .5;
                        Q(2, 3) = .5;
                        break;
                }
                break;
                
            default:
                Q(0, 0) = 1.;
                Q(1, 1) = 1.;
                Q(2, 2) = 1.;
                break;
        }
    } // end of getRefinementMappingMatrixL
    
    void ReferenceTriangularPrism::getRefinementMappingMatrixR(int refineType, std::size_t subElementIdx, LinearAlgebra::Matrix& Q) const
    {
        Q.resize(4, 4);
        Q = 0.;
        Q(3, 3) = 1.;
        switch (refineType)
        {
            case 0:
                switch (subElementIdx)
                {
                    case 0:
                        Q(0, 0) = 2.;
                        Q(0, 3) = 1.;
                        Q(1, 1) = 1.;
                        Q(2, 2) = 1.;
                        break;
                    case 1:
                        Q(0, 0) = 2.;
                        Q(0, 3) = -1.;
                        Q(1, 1) = 1.;
                        Q(2, 2) = 1.;
                        break;
                }
                break;
                
            case 1:
                switch (subElementIdx)
                {
                    case 0:
                        Q(0, 0) = 1.;
                        Q(1, 1) = 2.;
                        Q(1, 3) = 1.;
                        Q(2, 2) = 1.;
                        break;
                    case 1:
                        Q(0, 0) = 1.;
                        Q(1, 1) = .5;
                        Q(1, 3) = -1.;
                        Q(2, 2) = 1.;
                        break;
                }
                break;
                
            case 2:
                switch (subElementIdx)
                {
                    case 0:
                        Q(0, 0) = 1.;
                        Q(1, 1) = 1.;
                        Q(2, 2) = 2.;
                        Q(2, 3) = 1.;
                        break;
                    case 1:
                        Q(0, 0) = 1.;
                        Q(1, 1) = 1.;
                        Q(2, 2) = 2.;
                        Q(2, 3) = -1.;
                        break;
                }
                break;
                
            case 3:
                switch (subElementIdx)
                {
                    case 0:
                        Q(0, 0) = 2.;
                        Q(0, 3) = 1.;
                        Q(1, 1) = 2.;
                        Q(1, 3) = 1.;
                        Q(2, 2) = 1.;
                        break;
                    case 1:
                        Q(0, 0) = 2.;
                        Q(0, 3) = -1.;
                        Q(1, 1) = 2.;
                        Q(1, 3) = 1.;
                        Q(2, 2) = 1.;
                        break;
                    case 2:
                        Q(0, 0) = 2.;
                        Q(0, 3) = 1.;
                        Q(1, 1) = 2.;
                        Q(1, 3) = -1.;
                        Q(2, 2) = 1.;
                        break;
                    case 3:
                        Q(0, 0) = 2.;
                        Q(0, 3) = -1.;
                        Q(1, 1) = 2.;
                        Q(1, 3) = -1.;
                        Q(2, 2) = 1.;
                        break;
                }
                break;
                
            case 4:
                switch (subElementIdx)
                {
                    case 0:
                        Q(0, 0) = 2.;
                        Q(0, 3) = 1.;
                        Q(1, 1) = 1.;
                        Q(2, 2) = 2.;
                        Q(2, 3) = 1.;
                        break;
                    case 1:
                        Q(0, 0) = 2.;
                        Q(0, 3) = -1.;
                        Q(1, 1) = 1.;
                        Q(2, 2) = 2.;
                        Q(2, 3) = 1.;
                        break;
                    case 2:
                        Q(0, 0) = 2.;
                        Q(0, 3) = 1.;
                        Q(1, 1) = 1.;
                        Q(2, 2) = 2.;
                        Q(2, 3) = -1.;
                        break;
                    case 3:
                        Q(0, 0) = 2.;
                        Q(0, 3) = -1.;
                        Q(1, 1) = 1.;
                        Q(2, 2) = 2.;
                        Q(2, 3) = -1.;
                        break;
                }
                break;
                
            case 5:
                switch (subElementIdx)
                {
                    case 0:
                        Q(0, 0) = 1.;
                        Q(1, 1) = 2.;
                        Q(1, 3) = 1.;
                        Q(2, 2) = 2.;
                        Q(2, 3) = 1.;
                        break;
                    case 1:
                        Q(0, 0) = 1.;
                        Q(1, 1) = 2.;
                        Q(1, 3) = -1.;
                        Q(2, 2) = 2.;
                        Q(2, 3) = 1.;
                        break;
                    case 2:
                        Q(0, 0) = 1.;
                        Q(1, 1) = 2.;
                        Q(1, 3) = 1.;
                        Q(2, 2) = 2.;
                        Q(2, 3) = -1.;
                        break;
                    case 3:
                        Q(0, 0) = 1.;
                        Q(1, 1) = 2.;
                        Q(1, 3) = -1.;
                        Q(2, 2) = 2.;
                        Q(2, 3) = -1.;
                        break;
                }
                break;
                
            case 6:
                switch (subElementIdx)
                {
                    case 0:
                        Q(0, 0) = 2.;
                        Q(0, 3) = 1.;
                        Q(1, 1) = 2.;
                        Q(1, 3) = 1.;
                        Q(2, 2) = 2.;
                        Q(2, 3) = 1.;
                        break;
                    case 1:
                        Q(0, 0) = 2.;
                        Q(0, 3) = -1.;
                        Q(1, 1) = 2.;
                        Q(1, 3) = 1.;
                        Q(2, 2) = 2.;
                        Q(2, 3) = 1.;
                        break;
                    case 2:
                        Q(0, 0) = 2.;
                        Q(0, 3) = 1.;
                        Q(1, 1) = 2.;
                        Q(1, 3) = -1.;
                        Q(2, 2) = 2.;
                        Q(2, 3) = 1.;
                        break;
                    case 3:
                        Q(0, 0) = 2.;
                        Q(0, 3) = -1.;
                        Q(1, 1) = 2.;
                        Q(1, 3) = -1.;
                        Q(2, 2) = 2.;
                        Q(2, 3) = 1.;
                        break;
                    case 4:
                        Q(0, 0) = 2.;
                        Q(0, 3) = 1.;
                        Q(1, 1) = 2.;
                        Q(1, 3) = 1.;
                        Q(2, 2) = 2.;
                        Q(2, 3) = -1.;
                        break;
                    case 5:
                        Q(0, 0) = 2.;
                        Q(0, 3) = -1.;
                        Q(1, 1) = 2.;
                        Q(1, 3) = 1.;
                        Q(2, 2) = 2.;
                        Q(2, 3) = -1.;
                        break;
                    case 6:
                        Q(0, 0) = 2.;
                        Q(0, 3) = 1.;
                        Q(1, 1) = 2.;
                        Q(1, 3) = -1.;
                        Q(2, 2) = 2.;
                        Q(2, 3) = -1.;
                        break;
                    case 7:
                        Q(0, 0) = 2.;
                        Q(0, 3) = -1.;
                        Q(1, 1) = 2.;
                        Q(1, 3) = -1.;
                        Q(2, 2) = 2.;
                        Q(2, 3) = -1.;
                        break;
                }
                break;
                
            default:
                Q(0, 0) = 1.;
                Q(1, 1) = 1.;
                Q(2, 2) = 1.;
                break;
        }
    } // end of getRefinementMappingMatrixR
    
    void ReferenceTriangularPrism::getCodim1RefinementMappingMatrixL(int refineType, std::size_t subElementIdx, std::size_t faLocalIndex, LinearAlgebra::Matrix& Q) const
    {
        int faRefinementType(-1);
        std::size_t subFaceIndex(0);
        
        Q.resize(3, 3);
        Q = 0.;
        Q(2, 2) = 1.;
        
        switch (refineType)
        {
            case 0:
                switch (faLocalIndex)
                {
                    case 2:
                        faRefinementType = 0;
                        subFaceIndex = ((subElementIdx + 1) % 2);
                        break;
                        
                    default:
                        faRefinementType = -1;
                        subFaceIndex = 100;
                        break;
                }
                getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixL(faRefinementType, subFaceIndex, Q);
                break;
                
            case 1:
                switch (faLocalIndex)
                {
                    case 3:
                        faRefinementType = 0;
                        subFaceIndex = subElementIdx % 2;
                        break;
                        
                    default:
                        faRefinementType = -1;
                        subFaceIndex = 100;
                        break;
                }
                getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixL(faRefinementType, subFaceIndex, Q);
                break;
                
            case 2:
                switch (faLocalIndex)
                {
                    case 4:
                        faRefinementType = 0;
                        subFaceIndex = subElementIdx % 2;
                        break;
                        
                    default:
                        faRefinementType = -1;
                        subFaceIndex = 100;
                        break;
                }
                getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixL(faRefinementType, subFaceIndex, Q);
                break;
                
            case 3:
                switch (faLocalIndex)
                {
                    case 2:
                    case 3:
                    case 4:
                        faRefinementType = 1;
                        subFaceIndex = subElementIdx;
                        break;
                        
                    default:
                        faRefinementType = -1;
                        subFaceIndex = 100;
                        break;
                }
                getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixL(faRefinementType, subFaceIndex, Q);
                break;
                
            case 4:
            case 5:
                switch (faLocalIndex)
                {
                    case 2:
                        faRefinementType = 0;
                        subFaceIndex = ((subElementIdx / 2) + 1) % 2;
                        break;
                        
                    case 3:
                        faRefinementType = 0;
                        subFaceIndex = subElementIdx % 2;
                        break;
                        
                    case 4:
                        faRefinementType = 0;
                        subFaceIndex = (subElementIdx + 1) % 2;
                        break;
                        
                    default:
                        faRefinementType = -1;
                        subFaceIndex = 100;
                        break;
                }
                getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixL(faRefinementType, subFaceIndex, Q);
                break;
                
            case 6:
                switch (faLocalIndex)
                {
                    case 2:
                        Q(0, 0) = -1.;
                        Q(1, 1) = 1.;
                        break;
                        
                    case 3:
                        faRefinementType = 0;
                        subFaceIndex = subElementIdx;
                        getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixL(faRefinementType, subFaceIndex, Q);
                        break;
                        
                    case 4:
                        faRefinementType = 0;
                        subFaceIndex = (subElementIdx + 1) % 2;
                        getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixL(faRefinementType, subFaceIndex, Q);
                        if (subElementIdx == 0)
                        {
                            Q(0, 0) = -Q(0, 0);
                        }
                        break;
                        
                    default:
                        Q(0, 0) = 1.;
                        Q(1, 1) = 1.;
                        break;
                }
                break;
                
            case 7:
                switch (faLocalIndex)
                {
                    case 2:
                        faRefinementType = 0;
                        subFaceIndex = (subElementIdx + 1) % 2;
                        getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixL(faRefinementType, subFaceIndex, Q);
                        if (subElementIdx == 0)
                        {
                            Q(0, 0) = -Q(0, 0);
                        }
                        break;
                        
                    case 4:
                        faRefinementType = 0;
                        subFaceIndex = subElementIdx;
                        getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixL(faRefinementType, subFaceIndex, Q);
                        break;
                        
                    default:
                        Q(0, 0) = 1.;
                        Q(1, 1) = 1.;
                        break;
                }
                break;
                
            default:
                Q(0, 0) = 1.;
                Q(1, 1) = 1.;
                break;
        }
    } // end of getCodim1RefinementMappingMatrixL
    
    void ReferenceTriangularPrism::getCodim1RefinementMappingMatrixR(int refineType, std::size_t subElementIdx, std::size_t faLocalIndex, LinearAlgebra::Matrix& Q) const
    {
        int faRefinementType(-1);
        std::size_t subFaceIndex(0);
        
        Q.resize(3, 3);
        Q = 0.;
        Q(2, 2) = 1.;
        
        switch (refineType)
        {
            case 0:
                switch (faLocalIndex)
                {
                    case 2:
                        faRefinementType = 0;
                        subFaceIndex = ((subElementIdx + 1) % 2);
                        break;
                        
                    default:
                        faRefinementType = -1;
                        subFaceIndex = 100;
                        break;
                }
                getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixR(faRefinementType, subFaceIndex, Q);
                break;
                
            case 1:
                switch (faLocalIndex)
                {
                    case 3:
                        faRefinementType = 0;
                        subFaceIndex = subElementIdx % 2;
                        break;
                        
                    default:
                        faRefinementType = -1;
                        subFaceIndex = 100;
                        break;
                }
                getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixR(faRefinementType, subFaceIndex, Q);
                break;
                
            case 2:
                switch (faLocalIndex)
                {
                    case 4:
                        faRefinementType = 0;
                        subFaceIndex = subElementIdx % 2;
                        break;
                        
                    default:
                        faRefinementType = -1;
                        subFaceIndex = 100;
                        break;
                }
                getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixR(faRefinementType, subFaceIndex, Q);
                break;
                
            case 3:
                switch (faLocalIndex)
                {
                    case 2:
                    case 3:
                    case 4:
                        faRefinementType = 1;
                        subFaceIndex = subElementIdx;
                        break;
                        
                    default:
                        faRefinementType = -1;
                        subFaceIndex = 100;
                        break;
                }
                getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixR(faRefinementType, subFaceIndex, Q);
                break;
                
            case 4:
            case 5:
                switch (faLocalIndex)
                {
                    case 2:
                        faRefinementType = 0;
                        subFaceIndex = ((subElementIdx / 2) + 1) % 2;
                        break;
                        
                    case 3:
                        faRefinementType = 0;
                        subFaceIndex = subElementIdx % 2;
                        break;
                        
                    case 4:
                        faRefinementType = 0;
                        subFaceIndex = (subElementIdx + 1) % 2;
                        break;
                        
                    default:
                        faRefinementType = -1;
                        subFaceIndex = 100;
                        break;
                }
                getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixR(faRefinementType, subFaceIndex, Q);
                break;
                
            case 6:
                switch (faLocalIndex)
                {
                    case 2:
                        Q(0, 0) = -1.;
                        Q(1, 1) = 1.;
                        break;
                        
                    case 3:
                        faRefinementType = 0;
                        subFaceIndex = subElementIdx;
                        getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixR(faRefinementType, subFaceIndex, Q);
                        break;
                        
                    case 4:
                        faRefinementType = 0;
                        subFaceIndex = (subElementIdx + 1) % 2;
                        getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixR(faRefinementType, subFaceIndex, Q);
                        if (subElementIdx == 0)
                        {
                            Q(0, 0) = -Q(0, 0);
                            Q(0, 1) = -Q(0, 1);
                            Q(0, 2) = -Q(0, 2);
                        }
                        break;
                        
                    default:
                        Q(0, 0) = 1.;
                        Q(1, 1) = 1.;
                        break;
                }
                break;
                
            case 7:
                switch (faLocalIndex)
                {
                    case 2:
                        faRefinementType = 0;
                        subFaceIndex = (subElementIdx + 1) % 2;
                        getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixR(faRefinementType, subFaceIndex, Q);
                        if (subElementIdx == 0)
                        {
                            Q(0, 0) = -Q(0, 0);
                            Q(0, 1) = -Q(0, 1);
                            Q(0, 2) = -Q(0, 2);
                        }
                        break;
                        
                    case 4:
                        faRefinementType = 0;
                        subFaceIndex = subElementIdx;
                        getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixR(faRefinementType, subFaceIndex, Q);
                        break;
                        
                    default:
                        Q(0, 0) = 1.;
                        Q(1, 1) = 1.;
                        break;
                }
                break;
                
            default:
                Q(0, 0) = 1.;
                Q(1, 1) = 1.;
                break;
        }
    } // end of getCodim1RefinementMappingMatrixR
    
}
;

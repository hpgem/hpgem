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

#include "ReferenceLine.hpp"
#include "Mappings/MappingToRefPointToLine.hpp"
#include "Mappings/MappingToRefLineToLine.hpp"
#include "Geometry/ReferencePoint.hpp"
#include "Geometry/PointReference.hpp"

namespace Geometry
{
    /* Behold the reference line:
     *
     * (-1) 0---1---1 (+1)
     *
     */
    std::size_t ReferenceLine::localNodeIndexes_[2][1] =
    {
        { 0 },
        { 1 }
    };

    ReferenceLine::ReferenceLine():
        ReferenceGeometry(2,1, LINE),/// Line has two points 1+1
        referenceGeometryCodim1Ptr_(&ReferencePoint::Instance())
    {
        PointReferenceT p1(1), p2(1);
        p1[0] = -1.0;
        p2[0] = 1.0;
        points_[0] = p1;
        points_[1] = p2;

        mappingsLineToLine_[0] = &MappingToRefLineToLine0::Instance();
        mappingsLineToLine_[1] = &MappingToRefLineToLine1::Instance();

        mappingsPointToLine_[0] = &MappingToRefPointToLine0::Instance();
        mappingsPointToLine_[1] = &MappingToRefPointToLine1::Instance();

    }

    ReferenceLine::ReferenceLine(const ReferenceLine& copy):
        ReferenceGeometry(copy),
        referenceGeometryCodim1Ptr_(&ReferencePoint::Instance())
    {
    }

    bool ReferenceLine::isInternalPoint(const PointReferenceT& p) const
    {
        return ((p[0] >= -1.) && (p[0] <= 1.));
    }

    void ReferenceLine::getCenter(PointReferenceT& p) const
    {
        p[0] = 0.;
    }

    void ReferenceLine::getNode(const IndexT& i, PointReferenceT& point) const
    {
        point = points_[i];
    }

    std::ostream& operator<<(std::ostream& os, const ReferenceLine& line)
    {
        os << line.getName() << " ={ ";
        ReferenceLine::const_iterator it = line.points_.begin();
        ReferenceLine::const_iterator end = line.points_.end();

        for ( ; it != end; ++it)
        {
            os << (*it);
        }
        os << '}' << std::endl;

        return os;
    }
    // ================================== Codimension 0 ============================================
    std::size_t ReferenceLine::
    getCodim0MappingIndex(const ListOfIndexesT& list1, const ListOfIndexesT& list2) const
    {
        if (list1.size() == 2 && list2.size() == 2)
        {
            if (list1[0] == list2[0])
                return 0;
            else
                return 1;
        }
        else
        {
            throw "ERROR: number of nodes of reference square was larger than 4.";
        }
    }

    const MappingReferenceToReference*
    ReferenceLine::getCodim0MappingPtr(const IndexT i) const
    {
        if (i < 2)
        {
            return mappingsLineToLine_[i];
        }
        else
        {
            throw "ERROR: Asked for a mappingSquareToSquare larger than 7. There are only 8!";
        }
    }

    // ================================== Codimension 1 ============================================

    void ReferenceLine::
    getCodim1EntityLocalIndices(const IndexT faceIndex, ListOfIndexesT& faceNodesLocal) const
    {
        if (faceIndex < 2)
        {
            faceNodesLocal.resize(1); // 2 nodes per face
            faceNodesLocal[0] = (IndexT) localNodeIndexes_[faceIndex][0];
        }
    }
    const ReferenceGeometry* ReferenceLine::getCodim1ReferenceGeometry(const IndexT faceIndex) const
    {
        if (faceIndex < 2)
        {
            return referenceGeometryCodim1Ptr_;
        }
        else
        {
            throw "ERROR: Asked for a line face index larger than 1. There are only 2 'faces' in a line!";
        }
    }
    const MappingReferenceToReference*
    ReferenceLine::getCodim1MappingPtr(const IndexT faceIndex) const
    {
        if (faceIndex < 2)
        {
            return mappingsPointToLine_[faceIndex];
        }
        else
        {
            throw "ERROR: Asked for a square point index larger than 3. There are only 4 nodes in a square!";
        }
    }

};

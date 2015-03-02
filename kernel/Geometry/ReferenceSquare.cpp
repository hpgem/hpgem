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

#include "ReferenceSquare.hpp"
#include "ReferenceLine.hpp"
#include "Mappings/MappingToRefLineToSquare.hpp"
#include "Mappings/MappingToRefSquareToSquare.hpp"
#include "Geometry/PointReference.hpp"
#include "Geometry/ReferencePoint.hpp"
#include "LinearAlgebra/Matrix.hpp"

namespace Geometry
{
    /* Behold the reference square:
     *
     * (-1,+1) 2---3---3 (+1,+1)
     *         |       |
     *         1       2
     *         |       |
     * (-1,-1) 0---0---1 (+1,-1)
     *
     */
    std::size_t ReferenceSquare::localNodeIndexes_[4][2] =
    {
        { 0, 1 },
        { 0, 2 },
        { 1, 3 },
        { 2, 3 },
    };

    ReferenceSquare::ReferenceSquare():
        ReferenceGeometry(4,2, SQUARE),
        referenceGeometryCodim1Ptr_(&ReferenceLine::Instance())
    {
        // See MappingLineToSquare.hpp for further info.                 Ref.Line     Ref.Sqr.Side
        mappingsLineToSquare_[0] = &MappingToRefLineToSquare0::Instance(); // x         -> 0:( x,-1)
        mappingsLineToSquare_[1] = &MappingToRefLineToSquare1::Instance(); // x         -> 1:(-1, x)
        mappingsLineToSquare_[2] = &MappingToRefLineToSquare2::Instance(); // x         -> 2:( 1, x)
        mappingsLineToSquare_[3] = &MappingToRefLineToSquare3::Instance(); // x         -> 3:( x, 1)

        // See Mapping SquareToSquare for further info.                      Ref. Sq.    Ref. Sq.
        mappingsSquareToSquare_[0] = &MappingToRefSquareToSquare0::Instance(); // (x,y)    -> (x,y)
        mappingsSquareToSquare_[1] = &MappingToRefSquareToSquare1::Instance(); // (x,y)    -> (-y,x)
        mappingsSquareToSquare_[2] = &MappingToRefSquareToSquare2::Instance(); // (x,y)    -> (-x,-y)
        mappingsSquareToSquare_[3] = &MappingToRefSquareToSquare3::Instance(); // (x,y)    -> (y,-x)
        mappingsSquareToSquare_[4] = &MappingToRefSquareToSquare4::Instance(); // (x,y)    -> (x,-y)
        mappingsSquareToSquare_[5] = &MappingToRefSquareToSquare5::Instance(); // (x,y)    -> (-x,y)
        mappingsSquareToSquare_[6] = &MappingToRefSquareToSquare6::Instance(); // (x,y)    -> (-y,-x)
        mappingsSquareToSquare_[7] = &MappingToRefSquareToSquare7::Instance(); // (x,y)    -> (y,x)

        // We set the actual coordinates (see top comment for drawing).
        PointReferenceT p1(2), p2(2), p3(2), p4(2);

        p1[0] = -1.0; p1[1] = -1.0;
        p2[0] = +1.0; p2[1] = -1.0;
        p3[0] = -1.0; p3[1] = +1.0;
        p4[0] = +1.0; p4[1] = +1.0;

        points_[0] = p1;
        points_[1] = p2;
        points_[2] = p3;
        points_[3] = p4;
    }
    
    ReferenceSquare::ReferenceSquare(const ReferenceSquare& copy):
        ReferenceGeometry(copy),
        referenceGeometryCodim1Ptr_(copy.referenceGeometryCodim1Ptr_)
    {
    }
    
    bool ReferenceSquare::isInternalPoint(const PointReferenceT& p) const
    {
        return ((p[0] >= -1.) && (p[0] <= 1.) && (p[1] >= -1.) && (p[1] <= 1.));
    }
    
    void ReferenceSquare::getCenter(PointReferenceT& p) const
    {
        p[0] = p[1] = 0.;
    }
    
    void ReferenceSquare::getNode(const IndexT& i, PointReferenceT& point) const
    {
        point = points_[i];
    }

    std::ostream& operator<<(std::ostream& os, const ReferenceSquare& square)
    {
        os <<square.getName()<<"={";
        ReferenceSquare::const_iterator it = square.points_.begin();
        ReferenceSquare::const_iterator end = square.points_.end();

        for ( ; it != end; ++it)
        {
            os << (*it) << '\t';
        }
        os <<'}'<<std::endl;

        return os;
    }

    // ================================== Codimension 0 ============================================
    std::size_t ReferenceSquare::
    getCodim0MappingIndex(const ListOfIndexesT& list1, const ListOfIndexesT& list2) const
    {
        if (list1.size() == 4 && list2.size() == 4)
        {
            if (list1[0] == list2[0])
            {
                if (list1[1] == list2[1])
                    return 0;
                else if (list1[3] == list2[3])
                    return 7;
                else
                    throw "ERROR: reference square indexes where given in impossible order.";
            }
            else if (list1[0] == list2[1])
            {
                if (list1[1] == list2[0])
                    return 5;
                else if (list1[1] == list2[3])
                    return 3;
                else
                    throw "ERROR: reference square indexes where given in impossible order.";
            }
            else if (list1[0] == list2[2])
            {
                if (list1[2] == list2[0])
                    return 4;
                else if (list1[2] == list2[3])
                    return 1;
                else
                    throw "ERROR: reference square indexes where given in impossible order.";
            }
            else
            {
                if (list1[1] == list2[1]) // (list1(3)==list2(0)) // Holds for both 2 and 6!
                    return 6;
                else if (list1[1] == list2[2])
                    return 2;
                else
                    throw "ERROR: reference square indexes where given in impossible order.";
            }
        }
        else
        {
            throw "ERROR: number of nodes of reference square was larger than 4.";
        }
    }

    const MappingReferenceToReference*
    ReferenceSquare::getCodim0MappingPtr(const IndexT i) const
    {
        if (i < 8)
        {
            return mappingsSquareToSquare_[i];
        }
        else
        {
            throw "ERROR: Asked for a mappingSquareToSquare larger than 7. There are only 8!";
        }
    }
    // ================================== Codimension 1 ============================================
    void ReferenceSquare::
    getCodim1EntityLocalIndices(const IndexT faceIndex, ListOfIndexesT& faceNodesLocal) const
    {
        if (faceIndex < 4)
        {
            faceNodesLocal.resize(2); // 2 nodes per face
            faceNodesLocal[0] = (IndexT) localNodeIndexes_[faceIndex][0];
            faceNodesLocal[1] = (IndexT) localNodeIndexes_[faceIndex][1];
        }
    }
    const ReferenceGeometry* ReferenceSquare::getCodim1ReferenceGeometry(const IndexT faceIndex) const
    {
        if (faceIndex < 4)
        {
            return referenceGeometryCodim1Ptr_;
        }
        else
        {
            throw "ERROR: Asked for a square face index larger than 3. There are only 4 faces in a square!";
        }
    }
    const MappingReferenceToReference*
    ReferenceSquare::getCodim1MappingPtr(const IndexT faceIndex) const
    {
        if (faceIndex < 4)
        {
            return mappingsLineToSquare_[faceIndex];
        }
        else
        {
            throw "ERROR: Asked for a square point index larger than 3. There are only 4 nodes in a square!";
        }
    }

	const ReferenceGeometry* ReferenceSquare::getCodim2ReferenceGeometry(const std::size_t) const {
		return &Geometry::ReferencePoint::Instance();
	}

    // =============================== Refinement mappings =====================================
    
    //! Transform a reference point using refinement mapping
    void ReferenceSquare::refinementTransform(std::size_t refineType, std::size_t subElementIdx, 
                  const PointReferenceT& p, PointReferenceT& pMap) const 
    {
        switch (refineType)
        {
          case 0:
            switch (subElementIdx)
            {
              case 0:
                pMap[0]= (p[0]-1)/2.;
                pMap[1]=p[1];
                break;
              case 1:
                pMap[0]= (p[0]+1)/2.;
                pMap[1]=p[1];
                break;
            }
            break;
            
          case 1:
            switch (subElementIdx)
            {
              case 0:
                pMap[0]=p[0];
                pMap[1]= (p[1]-1)/2.;
                break;
              case 1:
                pMap[0]=p[0];
                pMap[1]= (p[1]+1)/2.;
                break;
            }
            break;

          case 2:
            switch (subElementIdx)
            {
              case 0:
                pMap[0]= (p[0]-1)/2.;
                pMap[1]= (p[1]-1)/2.;
                break;
              case 1:
                pMap[0]= (p[0]+1)/2.;
                pMap[1]= (p[1]-1)/2.;
                break;
              case 2:
                pMap[0]= (p[0]-1)/2.;
                pMap[1]= (p[1]+1)/2.;
                break;
              case 3:
                pMap[0]= (p[0]+1)/2.;
                pMap[1]= (p[1]+1)/2.;
                break;
            }
            break;
            
          default:
            pMap = p;
            break;
        }
    }

    //! Transformation matrix of this refinement when located on the LEFT side
    void ReferenceSquare::getRefinementMappingMatrixL(std::size_t refineType, std::size_t subElementIdx, 
                LinearAlgebra::Matrix& Q) const 
    {
        Q.resize(3,3);
        Q = 0.;
        Q(2,2) = 1.;
        switch (refineType)
        {
          case 0:
            switch (subElementIdx)
            {
              case 0:
                Q(0,0) = .5;  Q(0,2) = -.5;
                Q(1,1) = 1.;
                break;
              case 1:
                Q(0,0) = .5;  Q(0,2) =  .5;
                Q(1,1) = 1.;
                break;
            }
            break;
            
          case 1:
            switch (subElementIdx)
            {
              case 0:
                Q(0,0) = 1.;
                Q(1,1) = .5;  Q(1,2) = -.5;
                break;
              case 1:
                Q(0,0) = 1.;
                Q(1,1) = .5;  Q(1,2) =  .5;
                break;
            }
            break;

          case 2:
            switch (subElementIdx)
            {
              case 0:
                Q(0,0) = .5;  Q(0,2) = -.5;
                Q(1,1) = .5;  Q(1,2) = -.5;
                break;
              case 1:
                Q(0,0) = .5;  Q(0,2) =  .5;
                Q(1,1) = .5;  Q(1,2) = -.5;
                break;
              case 2:
                Q(0,0) = .5;  Q(0,2) = -.5;
                Q(1,1) = .5;  Q(1,2) =  .5;
                break;
              case 3:
                Q(0,0) = .5;  Q(0,2) =  .5;
                Q(1,1) = .5;  Q(1,2) =  .5;
                break;
            }
            break;
            
          default:
            Q(0,0) = 1.;
            Q(1,1) = 1.;
            break;
            
        }
    }

    //! Transformation matrix of this refinement when located on the RIGHT side
    void ReferenceSquare::getRefinementMappingMatrixR(std::size_t refineType, std::size_t subElementIdx, 
                LinearAlgebra::Matrix& Q) const 
    {
        Q.resize(3,3);
        Q = 0.;
        Q(2,2) = 1.;
        switch (refineType)
        {
          case 0:
            switch (subElementIdx)
            {
              case 0:
                Q(0,0) = 2.;  Q(0,2) =  1.;
                Q(1,1) = 1.;
                break;
              case 1:
                Q(0,0) = 2.;  Q(0,2) = -1.;
                Q(1,1) = 1.;
                break;
            }
            break;
            
          case 1:
            switch (subElementIdx)
            {
              case 0:
                Q(0,0) = 1.;
                Q(1,1) = 2.;  Q(1,2) =  1.;
                break;
              case 1:
                Q(0,0) = 1.;
                Q(1,1) = 2.;  Q(1,2) = -1.;
                break;
            }
            break;

          case 2:
            switch (subElementIdx)
            {
              case 0:
                Q(0,0) = 2.;  Q(0,2) =  1.;
                Q(1,1) = 2.;  Q(1,2) =  1.;
                break;
              case 1:
                Q(0,0) = 2.;  Q(0,2) = -1.;
                Q(1,1) = 2.;  Q(1,2) =  1.;
                break;
              case 2:
                Q(0,0) = 2.;  Q(0,2) =  1.;
                Q(1,1) = 2.;  Q(1,2) = -1.;
                break;
              case 3:
                Q(0,0) = 2.;  Q(0,2) = -1.;
                Q(1,1) = 2.;  Q(1,2) = -1.;
                break;
            }
            break;
            
          default:
            Q(0,0) = 1.;
            Q(1,1) = 1.;
            break;
        }
    }

    //! Refinement mapping on codim1 for a given refinement on codim0
    //! Note: this should also applied on other dimensions
    void ReferenceSquare::getCodim1RefinementMappingMatrixL(std::size_t refineType, DimT subElementIdx, 
                            DimT faLocalIndex, LinearAlgebra::Matrix& Q) const 
    {}

    //! Refinement mapping on codim1 for a given refinement on codim0
    //! Note: this should also applied on other dimensions
    void ReferenceSquare::getCodim1RefinementMappingMatrixR(std::size_t refineType, DimT subElementIdx, 
                            DimT faLocalIndex, LinearAlgebra::Matrix& Q) const 
    {}

};

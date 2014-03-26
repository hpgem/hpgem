//  ReferenceHypercube.cpp
//
//  Created by Shavarsh Nurijanyan on 2/7/13.

#include "ReferenceHypercube.hpp"


namespace Geometry
{
    int ReferenceHypercube::localNodeIndexes_[8][8] =
    {
        { 0, 1, 2, 3, 4, 5, 6, 7 },
        { 0, 1, 2, 3, 8, 9, 10, 11},
        { 0, 1, 4, 5, 8, 9, 12, 13},
        { 0, 2, 4, 6, 8, 10, 12, 14},
        { 1, 3, 5, 7, 9, 11, 13, 15},
        { 2, 3, 6, 7, 10, 11, 14, 15},
        { 4, 5, 6, 7, 12, 13, 14, 15},
        { 8, 9, 10, 11, 12, 13, 14, 15},
    };

    ReferenceHypercube::ReferenceHypercube():
        ReferenceGeometry(FourD+12,4, HYPERCUBE),
        referenceGeometryCodim1Ptr_(&ReferenceCube::Instance()),
        referenceGeometryCodim2Ptr_(&ReferenceSquare::Instance()),
        referenceGeometryCodim3Ptr_(&ReferenceLine::Instance())

    {
        PointReferenceT p0(4), p1(4), p2(4), p3(4), p4(4), p5(4), p6(4), p7(4), p8(4);
        PointReferenceT p9(4), p10(4), p11(4), p12(4), p13(4), p14(4), p15(4);

        p0[0] = -1.0;   p0[1] = -1.0;   p0[2] = -1.0;   p0[3] = -1.0;
        p1[0] = +1.0;   p1[1] = -1.0;   p1[2] = -1.0;   p1[3] = -1.0;

        p2[0] = -1.0;   p2[1] = +1.0;   p2[2] = -1.0;   p2[3] = -1.0;
        p3[0] = +1.0;   p3[1] = +1.0;   p3[2] = -1.0;   p3[3] = -1.0;
        

        p4[0] = -1.0;   p4[1] = -1.0;   p4[2] = +1.0;   p4[3] = -1.0;
        p5[0] = +1.0;   p5[1] = -1.0;   p5[2] = +1.0;   p5[3] = -1.0;

        p6[0] = -1.0;   p6[1] = +1.0;   p6[2] = +1.0;   p6[3] = -1.0;
        p7[0] = +1.0;   p7[1] = +1.0;   p7[2] = +1.0;   p7[3] = -1.0;
        
        

        p8[0] = -1.0;   p8[1] = -1.0;   p8[2] = -1.0;   p8[3] = +1.0;
        p9[0] = +1.0;   p9[1] = -1.0;   p9[2] = -1.0;   p9[3] = +1.0;

        p10[0] = -1.0;  p10[1] = +1.0;  p10[2] = -1.0;   p10[3] = +1.0;
        p11[0] = +1.0;  p11[1] = +1.0;  p11[2] = -1.0;   p11[3] = +1.0;
        

        p12[0] = -1.0;  p12[1] = -1.0;  p12[2] = +1.0;   p12[3] = +1.0;
        p13[0] = +1.0;  p13[1] = -1.0;  p13[2] = +1.0;   p13[3] = +1.0;

        p14[0] = -1.0;  p14[1] = +1.0;  p14[2] = +1.0;   p14[3] = +1.0;
        p15[0] = +1.0;  p15[1] = +1.0;  p15[2] = +1.0;   p15[3] = +1.0;
        
        ///reference element
        
        points_[0]  = p0;  points_[1]  = p1;  points_[2]  = p2;  points_[3]  = p3;
        points_[4]  = p4;  points_[5]  = p5;  points_[6]  = p6;  points_[7]  = p7;
        points_[8]  = p8;  points_[9]  = p9;  points_[10] = p10; points_[11] = p11;
        points_[12] = p12; points_[13] = p13; points_[14] = p14; points_[15] = p15;

        mappingsCubeToHypercube_[0] = &MappingToRefCubeToHypercube0::Instance();
        mappingsCubeToHypercube_[1] = &MappingToRefCubeToHypercube1::Instance();
        mappingsCubeToHypercube_[2] = &MappingToRefCubeToHypercube2::Instance();
        mappingsCubeToHypercube_[3] = &MappingToRefCubeToHypercube3::Instance();
        mappingsCubeToHypercube_[4] = &MappingToRefCubeToHypercube4::Instance();
        mappingsCubeToHypercube_[5] = &MappingToRefCubeToHypercube5::Instance();
        mappingsCubeToHypercube_[6] = &MappingToRefCubeToHypercube6::Instance();
        mappingsCubeToHypercube_[7] = &MappingToRefCubeToHypercube7::Instance();
    }

    bool
    ReferenceHypercube::isInternalPoint(const PointReferenceT& p) const
    {
        return ((p[0] >= -1.) && (p[0] <= 1.) &&
                (p[1] >= -1.) && (p[1] <= 1.) &&
                (p[2] >= -1.) && (p[2] <= 1.) &&
                (p[3] >= -1.) && (p[3] <= 1.));
    }
    
    void
    ReferenceHypercube::getCenter(PointReferenceT& p) const
    {
        p[3] = p[2] = p[1] = p[0] = 0.;
    }
    
    void
    ReferenceHypercube::getNode(const IndexT& i, PointReferenceT& point) const
    {
        point = points_[i];
    }

    std::ostream& operator<<(std::ostream& os, const ReferenceHypercube& hypercube)
    {
        os <<hypercube.getName()<<"=( ";
        ReferenceHypercube::const_iterator cit = hypercube.points_.begin();
        ReferenceHypercube::const_iterator cend = hypercube.points_.end();

        for ( ; cit != cend; ++cit)
        {
            os << (*cit) << ' ';
        }
        os <<')'<<std::endl;

        return os;
    }

    // ================================== Codimension 0 ============================================

    int ReferenceHypercube::getCodim0MappingIndex(const ListOfIndexesT& list1, const ListOfIndexesT& list2) const
    {
        throw "ReferenceCube::getCodim0MappingIndex not implemented";
    }

    const MappingReferenceToReference*
    ReferenceHypercube::getCodim0MappingPtr(const IndexT i) const
    {
        throw "ReferenceCube::getCodim0MappingPtr not implemented";
    }

    // ================================== Codimension 1 ============================================

    const MappingReferenceToReference*
    ReferenceHypercube::getCodim1MappingPtr(const IndexT faceIndex) const
    {
        if (faceIndex < 8)
        {
            return mappingsCubeToHypercube_[faceIndex];
        }
        else
        {
            throw "ERROR: ReferenceHypercube::getCodim1MappingPtr requested face index does not exist";
        }
    }

    const ReferenceGeometry*
    ReferenceHypercube::getCodim1ReferenceGeometry(const IndexT faceIndex) const
    {
        if (faceIndex < 8)
        {
            return referenceGeometryCodim1Ptr_;
        }
        else
        {
            throw "ERROR: ReferenceHypercube::getCodim1ReferenceGeometry requested face index does not exist";
        }
    }

    void ReferenceHypercube::getCodim1EntityLocalIndices(const IndexT i, ListOfIndexesT& faceNodesLocal) const
    {
        if (i < 8)
        {
            faceNodesLocal.resize(8); // 2 nodes per edge
            faceNodesLocal[0] = localNodeIndexes_[i][0];
            faceNodesLocal[1] = localNodeIndexes_[i][1];
            faceNodesLocal[2] = localNodeIndexes_[i][2];
            faceNodesLocal[3] = localNodeIndexes_[i][3];
            faceNodesLocal[4] = localNodeIndexes_[i][4];
            faceNodesLocal[5] = localNodeIndexes_[i][5];
            faceNodesLocal[6] = localNodeIndexes_[i][6];
            faceNodesLocal[7] = localNodeIndexes_[i][7];
        }
        else
        {
            throw "ERROR: ReferenceHypercube::getCodim1EntityLocalIndices requested face index does not exist";
        }

    }

    // ================================== Codimension 2 ============================================

    const MappingReferenceToReference*
    ReferenceHypercube::getCodim2MappingPtr(const IndexT lineIndex) const
    {
        /// TODO: Implement face to hypercube mappings.
        throw "ERROR: ReferenceHypercube::getCodim2MappingPtr: face to hypercube mappings not implemented";
    }

    const ReferenceGeometry*
    ReferenceHypercube::getCodim2ReferenceGeometry(const IndexT e) const
    {
        if (e < 12)
        {
            return referenceGeometryCodim2Ptr_;
        }
        else
        {
            throw "ERROR: ReferenceHypercube::getCodim2ReferenceGeometry requested side index does not exist";
        }
    }

    void ReferenceHypercube::getCodim2EntityLocalIndices(const IndexT i, ListOfIndexesT& edgeNodesLocal) const
    {
        if (i < 24)
        {
            throw "ReferenceHypercube::getCodim2EntityLocalIndices: not implemented";
        }
        else
        {
            throw "ReferenceHypercube::getCodim2EntityLocalIndices requested side index does not exist";
        }
    }

    // ================================== Codimension 3 ============================================

    void ReferenceHypercube::getCodim3EntityLocalIndices(const IndexT i, ListOfIndexesT& edgeNodesLocal) const
    {
        if (i < 16)
        {
            throw "ReferenceHypercube::getCodim3EntityLocalIndices not implemented";
        }
        else
        {
            throw "ReferenceHypercube::getCodim2EntityLocalIndices requested side index does not exist";
        }

    }
    
    
    // ================================== Quadrature rules =====================================

    /// Add a quadrature rule into the list of valid quadrature rules for this geometry.
    void ReferenceHypercube::addGaussQuadratureRule(QuadratureRules::GaussQuadratureRule* const qr)
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
    QuadratureRules::GaussQuadratureRule* const ReferenceHypercube::getGaussQuadratureRule(int order) const
    {
        for (std::list<QuadratureRules::GaussQuadratureRule*>::const_iterator it = lstGaussQuadratureRules_.begin();
              it != lstGaussQuadratureRules_.end(); ++it)
          if ((*it)->order() >= order) return *it;

        return NULL;
    }

};

//
//  ReferenceLine.hpp
//  
//
//  Created by Shavarsh Nurijanyan on 2/4/13.
//
//

#include "ReferenceLine.hpp"
#include "Mappings/MappingToRefPointToLine.hpp"
#include "Mappings/MappingToRefLineToLine.hpp"

namespace Geometry
{
    /* Behold the reference line:
     *
     * (-1) 0---1---1 (+1)
     *
     */
    int ReferenceLine::localNodeIndexes_[2][1] =
    {
        { 0 },
        { 1 }
    };

    ReferenceLine::ReferenceLine():
        ReferenceGeometry(OneD+1,1, LINE),/// Line has two points 1+1
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
    int ReferenceLine::
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


    // ================================== Quadrature rules =====================================

    /// Add a quadrature rule into the list of valid quadrature rules for this geometry.
    void ReferenceLine::addGaussQuadratureRule(QuadratureRules::GaussQuadratureRule* const qr)
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
    QuadratureRules::GaussQuadratureRule* const ReferenceLine::getGaussQuadratureRule(int order) const
    {
        for (std::list<QuadratureRules::GaussQuadratureRule*>::const_iterator it = lstGaussQuadratureRules_.begin();
              it != lstGaussQuadratureRules_.end(); ++it)
          if ((*it)->order() >= order) return *it;

        return NULL;
    }

};

//
//  ReferenceTriangle.cpp
//  
//
//  Created by Shavarsh Nurijanyan on 2/5/13.
//
//

#include "ReferenceTriangle.hpp"

namespace Geometry
{
    /* The ordering of the vertex and faces in a triangle:
     *
     *   (0,1) 2
     *         | \
     *         1   2
     *         |     \
     *   (0,0) 0---0---1 (1,0)
     *
     */
    int ReferenceTriangle::localNodeIndexes_[3][2] =
    {
        { 0, 1 },
        { 0, 2 },
        { 1, 2 }
    };

    ReferenceTriangle::ReferenceTriangle():
        ReferenceGeometry<TwoD>(TwoD+1,TRIANGLE),
        referenceGeometryCodim1Ptr_(&ReferenceLine::Instance())
    {

        // See MappingLineToTriangle.hpp for further info.                   Ref.Line->Ref.Tr.Side
        mappingsLineToTriangle_[0] = &MappingToRefLineToTriangle0::Instance(); // x -> 0:((1+x)/2,0)
        mappingsLineToTriangle_[1] = &MappingToRefLineToTriangle1::Instance(); // x -> 1:(0,(1+x)/2)
        mappingsLineToTriangle_[2] = &MappingToRefLineToTriangle2::Instance(); // x -> 2:((1-x)/2,(1+x)/2)

        // See Mapping TriangleToTriangle for further info.                          Ref.Tr. Ref. Tr.
        mappingsTriangleToTriangle_[0] = &MappingToRefTriangleToTriangle0::Instance(); // (x,y) -> (x,y)
        mappingsTriangleToTriangle_[1] = &MappingToRefTriangleToTriangle1::Instance(); // (x,y) -> (-y,x)
        mappingsTriangleToTriangle_[2] = &MappingToRefTriangleToTriangle2::Instance(); // (x,y  -> (-x,-y)
        mappingsTriangleToTriangle_[3] = &MappingToRefTriangleToTriangle3::Instance(); // (x,y) -> (y,x)
        mappingsTriangleToTriangle_[4] = &MappingToRefTriangleToTriangle4::Instance(); // (x,y) -> (x,-y)
        mappingsTriangleToTriangle_[5] = &MappingToRefTriangleToTriangle5::Instance(); // (x,y) -> (-x,y)

        PointReferenceT p1, p2, p3;
        p1[0] = 0.0; p1[1] = 0.0;
        p2[0] = 1.0; p2[1] = 0.0;
        p3[0] = 0.0; p3[1] = 1.0;

        points_[0] = p1;
        points_[1] = p2;
        points_[2] = p3;
    }
    
    ReferenceTriangle::ReferenceTriangle(const ReferenceTriangle& copy):
        ReferenceGeometry<TwoD>(copy),
        referenceGeometryCodim1Ptr_(&ReferenceLine::Instance())
    {}
    
    bool ReferenceTriangle::isInternalPoint(const PointReferenceT& p) const
    {
        return ((p[0] >= 0.) && (p[0] <= 1.) && (p[1] >= 0.) && (p[1] <= 1. - p[0]));
    }
    
    void ReferenceTriangle::getCenter(PointReferenceT& p) const
    {
        p[0] = p[1] = 1. / 3.;
    }

    void ReferenceTriangle::getNode(const IndexT& i, PointReferenceT& point) const
    {
        point = points_[i];
    }

    std::ostream& operator<<(std::ostream& os, const ReferenceTriangle& triangle)
    {
        os <<triangle.getName()<<" =( ";
        ReferenceTriangle::const_iterator it = triangle.points_.begin();
        ReferenceTriangle::const_iterator end = triangle.points_.end();

        for ( ; it != end; ++it)
        {
            os << (*it) << ' ';
        }
        os <<')'<<std::endl;

        return os;
    }
    // ================================== Codimension 0 ============================================
    int ReferenceTriangle::
    getCodim0MappingIndex(const ListOfIndexesT& list1, const ListOfIndexesT& list2) const
    {
        if (list1.size() == 3 && list2.size() == 3)
        {
            if (list1[0] == list2[0])
            {
                if (list1[1] == list2[1])
                    return 0; // 0.1.2.
                else
                    return 1; // 0.2.1.
            }
            else
            {
                if (list1[0] == list2[1])
                {
                    if (list1[1] == list2[2])
                        return 5; // 1.2.0.
                    else
                        return 3; // 1.0.2.
                }
                else
                {
                    if (list1[1] == list2[1])
                        return 4; // 2.1.0.
                    else
                        return 2; // 2.0.1.
                }
            }
        }
        else
        {
            throw "ERROR: number of nodes of reference triangle was larger than 3.";
        }
    }

    const MappingReferenceToReference<2, 2>*
    ReferenceTriangle::getCodim0MappingPtr(const IndexT i) const
    {
        if (i < 6)
        {
            return mappingsTriangleToTriangle_[i];
        }
        else
        {
            throw "ERROR: Asked for a mappingTriangleToTriangle larger than 5. There are only 6!";
        }
    }
    // ================================== Codimension 1 ============================================
    void ReferenceTriangle::
    getCodim1EntityLocalIndices(const IndexT faceIndex, ListOfIndexesT& faceNodesLocal) const
    {
        if (faceIndex < 3)
        {
            faceNodesLocal.resize(2); // 2 nodes per face
            faceNodesLocal[0] = (IndexT) localNodeIndexes_[faceIndex][0];
            faceNodesLocal[1] = (IndexT) localNodeIndexes_[faceIndex][1];
        }
    }
    const ReferenceGeometry<1>* ReferenceTriangle::getCodim1ReferenceGeometry(const IndexT faceIndex) const
    {
        if (faceIndex < 3)
        {
            return referenceGeometryCodim1Ptr_;
        }
        else
        {
            throw "ERROR: Asked for a triangle face index larger than 2. There are only 3 faces in a triangle!";
        }
    }
    const MappingReferenceToReference<1, 2>*
    ReferenceTriangle::getCodim1MappingPtr(const IndexT faceIndex) const
    {
        if (faceIndex < 3)
        {
            return mappingsLineToTriangle_[faceIndex];
        }
        else
        {
            throw "ERROR: Asked for a triangle point index larger than 3. There are only 3 nodes in a triangle!";
        }
    }


    // ================================== Quadrature rules =====================================

    /// Add a quadrature rule into the list of valid quadrature rules for this geometry.
    void ReferenceTriangle::addGaussQuadratureRule(QuadratureRules::GaussQuadratureRule<2>* const qr)
    {
        std::list<QuadratureRules::GaussQuadratureRule<2>*>::iterator it = lstGaussQuadratureRules_.begin();
        while (it != lstGaussQuadratureRules_.end())
        {
          if ((*it)->order() < qr->order()) ++it;
          else break;
        }
        lstGaussQuadratureRules_.insert(it,qr);
    }

    /// Get a valid quadrature for this geometry.
    QuadratureRules::GaussQuadratureRule<2>* const ReferenceTriangle::getGaussQuadratureRule(int order) const
    {
        for (std::list<QuadratureRules::GaussQuadratureRule<2>*>::const_iterator it = lstGaussQuadratureRules_.begin();
              it != lstGaussQuadratureRules_.end(); ++it)
          if ((*it)->order() >= order) return *it;

        return NULL;
    }
            
};

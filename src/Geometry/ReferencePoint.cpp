//
//  ReferencePoint.cpp
//  
//
//  Created by Shavarsh Nurijanyan on 2/4/13.
//
//

#include "ReferencePoint.hpp"
#include "Mappings/MappingToRefPointToPoint.hpp"

namespace Geometry
{
    /* Behold the reference point, ruler of them all:
     *
     *                  0.
     *
     */
    ReferencePoint::ReferencePoint():
         ReferenceGeometry(POINT)
    { mappingsPointToPoint_=&MappingToRefPointToPoint::Instance();}

    ReferencePoint::ReferencePoint(const ReferencePoint& copy):
        ReferenceGeometry(copy),mappingsPointToPoint_(copy.mappingsPointToPoint_)
    { }

    bool ReferencePoint::isInternalPoint(const PointReference& p) const
    {
        return true;
    }

    void ReferencePoint::getCenter(PointReference& p) const { }

    void ReferencePoint::getNode(const IndexT& i, PointReference& point) const { }

    int ReferencePoint::getCodim0MappingIndex(const ListOfIndexesT&, const ListOfIndexesT&) const
    {
        return 0;
    }

    const MappingReferenceToReference* ReferencePoint::getCodim0MappingPtr(const IndexT a) const
    {
        return mappingsPointToPoint_;
    }

    // ================================== Quadrature rules =====================================

    /// Add a quadrature rule into the list of valid quadrature rules for this geometry.
    void ReferencePoint::addGaussQuadratureRule(QuadratureRules::GaussQuadratureRule* const qr)
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
    QuadratureRules::GaussQuadratureRule* const ReferencePoint::getGaussQuadratureRule(int order) const
    {
        for (std::list<QuadratureRules::GaussQuadratureRule*>::const_iterator it = lstGaussQuadratureRules_.begin();
              it != lstGaussQuadratureRules_.end(); ++it)
          if ((*it)->order() >= order) return *it;

        return NULL;
    }

};

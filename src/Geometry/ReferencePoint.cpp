//
//  ReferencePoint.cpp
//  
//
//  Created by Shavarsh Nurijanyan on 2/4/13.
//
//

#include "ReferencePoint.hpp"

namespace Geometry
{
    /* Behold the reference point, ruler of them all:
     *
     *                  0.
     *
     */
    ReferencePoint::ReferencePoint():
         ReferenceGeometry<0>(POINT)
    { }

    ReferencePoint::ReferencePoint(const ReferencePoint& copy):
        ReferenceGeometry<0>(copy)
    { }

    bool ReferencePoint::isInternalPoint(const PointReference<0>& p) const
    {
        return true;
    }

    void ReferencePoint::getCenter(PointReference<0>& p) const { }

    void ReferencePoint::getNode(const IndexT& i, PointReference<0>& point) const { }

    int ReferencePoint::getCodim0MappingIndex(const ListOfIndexesT&, const ListOfIndexesT&) const
    {
        return 0;
    }

    const MappingReferenceToReference<0,0>* ReferencePoint::getCodim0MappingPtr(const IndexT a) const
    {
        return 0;
    }

    // ================================== Quadrature rules =====================================

    /// Add a quadrature rule into the list of valid quadrature rules for this geometry.
    void ReferencePoint::addGaussQuadratureRule(QuadratureRules::GaussQuadratureRule<0>* const qr) 
    {
        std::list<QuadratureRules::GaussQuadratureRule<0>*>::iterator it = lstGaussQuadratureRules_.begin();
        while (it != lstGaussQuadratureRules_.end())
        {
          if ((*it)->order() < qr->order()) ++it;
          else break;
        }
        lstGaussQuadratureRules_.insert(it,qr);
    }

    /// Get a valid quadrature for this geometry.
    QuadratureRules::GaussQuadratureRule<0>* const ReferencePoint::getGaussQuadratureRule(int order) const
    {
        for (std::list<QuadratureRules::GaussQuadratureRule<0>*>::const_iterator it = lstGaussQuadratureRules_.begin();
              it != lstGaussQuadratureRules_.end(); ++it)
          if ((*it)->order() >= order) return *it;

        return NULL;
    }

};

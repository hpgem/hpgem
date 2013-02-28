//------------------------------------------------------------------------------
// ~OC~
// File: outwardNormalVectorSign.hh
// Hack to avoid rearranging the face generation wrt
// ReferenceGeometries and Mappings. This will vanish with the geometry remake.
// Lars Pesch, Wed Aug 17 15:05:17 CEST 2005
//------------------------------------------------------------------------------
#ifndef OUTWARDNORMALVECTORSIGN_HH
#define OUTWARDNORMALVECTORSIGN_HH
#include "MappingReferenceToReference.hpp"

namespace Geometry
{
    /*!
     *  Welcome!, to the most ugly piece of code of hpGEM. Warning!, mental sanity in danger!
     *  ~OC~
     *  This function is a temporary fix: some of Henk's old face-to-element mappings have inward
     *  normal vectors and thus need a sign to conform to the convention that they all must be
     *  outward. This sign is applied in the function getNormalVector of Face (the only client of
     *  outwardNormalVectorSign). For all newly added reference geometries all face-to-element
     *  mappings were designed from the outset so that the normal vector is outward.
     *
     * TODO: Switch old mappings to be outward.
     */
    template <unsigned int DIM>
    double OutwardNormalVectorSign(const MappingReferenceToReference<DIM-1, DIM>* const map);
}
#endif

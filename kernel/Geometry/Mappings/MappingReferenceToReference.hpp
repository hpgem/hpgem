/*
 * MappingReferenceToReference.hpp
 *
 *  Created on: Feb 11, 2013
 *      Author: nicorivas
 */

#ifndef REFERENCETOREFERENCEM_H_
#define REFERENCETOREFERENCEM_H_


#include "../PointReference.hpp"
#include "MappingInterface.hpp"

namespace Geometry
{

/*! \brief Intermediate ABC for reference space to reference space mappings.

    Mappings from and to reference space are used in two contexts:

        o from a face of a ReferenceGeometry onto the Geometry, hence (dim-1)->dim; this is needed
          when evaluating the trace of a function (defined on the elements).

        o from a ReferenceGeometry onto itself, needed for the face-to-face mappings (because of
          the possible permutations of the global node indices of the face nodes when seen from the
          two elements on either side of the face).

    Note that the implementations have constant state, as they concern only reference coordinates.
    There is no need for reinit function, as the ReferenceToPhysicalMappings. Also transform can be
    declared here because, in contrast to Ref2PhysSpaceMapping, it concerns RefSpacePoint objects
    as both arguments. */
    
        //typedef MappingInterface<unsigned int dimFrom, unsigned int dimTo> shit;
    
    class MappingReferenceToReference: public MappingInterface
    {
        public:
            virtual void transform(const Geometry::PointReference&,
                                         Geometry::PointReference&) const = 0;
    };
};
#endif /* REFERENCETOPHYSICALM_H_ */

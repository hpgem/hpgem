/*
 * ReferenceToPhysicalM.hpp
 *
 *  Created on: Feb 5, 2013
 *      Author: nicorivas
 */

#ifndef REFERENCETOPHYSICALM_H_
#define REFERENCETOPHYSICALM_H_

#include "MappingInterface.hpp"
#include "../PointReference.hpp"
#include "../PointPhysical.hpp"
#include "../PhysicalGeometry.hpp"

namespace Geometry
{
/*! ~OC~
    Second layer abstract base class (derived from Mapping) for mappings
    that go from a reference element to physical space (hence different point
    types in the two spaces, cf. member function transform).

    In the current design of mappings, these do not store the point
    coordinates, but this can neither be enforced, nor am I sure that this is
    essential. This should be reconsidered at a later stage. The access to
    physical space information is via a reference to the NodeContainer;
    mappings can index into that one with the global node numbers of the
    elements.

    The rein(*it)-function is meant to alert an object of a change of the layout
    of the Element in physical space. Since the reference geometry of the
    Element does not change, but rather only (some of) the vertex positions,
    the Mapping can be adjusted to the new layout.*/

    template <unsigned int dimFrom, unsigned int dimTo = dimFrom>
    class MappingReferenceToPhysical: public MappingInterface<dimFrom,dimTo>
    {
        public:
            /// \bug This is a work around for g++ bug 14258 which is fixed in modern compliers so at some point change back
            typedef typename MappingInterface<dimFrom,dimTo>::PointToT           PointToT;
            typedef typename MappingInterface<dimFrom,dimTo>::PointPhysicalT     PointPhysicalT;
            typedef typename MappingInterface<dimFrom,dimTo>::PointReferenceT    PointReferenceT;
            typedef typename MappingInterface<dimFrom,dimTo>::JacobianT          JacobianT;

            typedef PhysicalGeometry<dimTo>  PhysicalGeometryT;

            typedef std::vector<PointPhysicalT>*    VectorOfPointsT;
        

        public:
            MappingReferenceToPhysical(): MappingInterface<dimFrom,dimTo>(){ }

            // Sets.
            void setNodesPtr(VectorOfPointsT* nodes) {nodes_ = nodes;}

            // Methods.
            //! ~OC~ Transform a point from reference space to physical space.
            virtual void transform(const PointReferenceT&, PointPhysicalT&) const = 0;
            //! ~OC~ Recompute mapping after physical nodes have moved.
            virtual void reinit(const PhysicalGeometryT* const) = 0;
            void getNodeCoordinates(const int index, PointPhysicalT& coords) const
                {coords = (*nodes_)[index].getCoordinates();}

        private:
            static std::vector<PointPhysical<dimFrom> >* nodes_; /// Pointer to the global node container.
    };

    template <unsigned int dimFrom, unsigned int dimTo>
    std::vector<Geometry::PointPhysical<dimFrom> >* MappingReferenceToPhysical<dimFrom, dimTo>::nodes_;
};

#endif /* REFERENCETOPHYSICALM_H_ */

/*
 * Map.hpp
 *
 *  Created on: Feb 3, 2013
 *      Author: nicorivas
 */

#ifndef MAPPINGINTERFACE_H_
#define MAPPINGINTERFACE_H_

#include "../Point.hpp"
#include "../PointReference.hpp"
#include "../ReferenceGeometry.hpp"
#include "../Jacobian.hpp"

namespace Geometry
{
    /*!
    \brief (OC): Base class that defines the functionality of a mapping between two coordinate systems.

    Mappings transform coordinates of a point with respect to one system to
    coordinates in another system. This can be to obtain physical coordinates
    for points in a reference geometry, or to map points from a face to one of
    the neighbouring elements' coordinate systems, etc. <BR>

    Using mappings, calculations can be carried out in reference space;
    e.g. quadrature rules are defined on the reference geometry, while some of
    the functions that are to be integrated live in physical space.<BR>

    Apart from mapping points to (physical) space, a Mapping gives information
    about the relative layout of the connected spaces, namely in its Jacobian.
    This information also has to be provided by implementation classes.<BR>

    End of interface description.

    Some remarks:
    <ul>
        <li> At the moment I expect the assignment of a certain Mapping
             implementation to an Element to be done by a function
             makeMapping() in ElementFactory. This will have to be
             user-implemented as the choice of the order, use of curved
             elements near the boundary etc. are problem-dependent.
        <li> In general one might expect a Mapping to be fully
             characterized by the dimension of the space(s) and the order.
             In that respect one may wonder that the Mapping implemenations
             are adapted to the reference geometries. This however allows
             them to be tailored to their actual purpose. Connected to this
             is the choice NOT to implement the transform and calcJacobian
             functions using shape functions. While they could make the
             representation (especially for higher order cases) simpler
             it would probably lead to slower computation: The
             \f$ \sum_{i=1}^N b_i(\xi) \vec{x}_i\f$ would need \f$ N\f$
             calls to basis function objects. (Note that such an
             implementation would require a Mapping to store the
             NodeList-indices of the vertices (and additional points) it uses
             for its mapping, so it would have to be rebuild after a
             compaction step, but would look up the \f$\vec{x}_i\f$
             coordinates for each evaluation, so it would not matter if these
             coordinates changed; with the current implementation
             it is just the other way round: it does not care about the
             compaction step, but it has to be reinitialized after changing
             vertex coordinates. I would expect most problems to work with
             fixed elements, so I deem this the more efficient approach.
             Additionally I expect the performance to be better since the
             necessary vectors are stored in Mapping, hence no need to
             interfere with a list that is somewhere else in memory; another
             minor pro is that the transformation formulae can be written
             (partially) in Horner form, so that some multiplications are
             saved. This effect is small though, as these are multiplications
             of scalars, while a MUCH bigger proportion will go into the
             operations on Point objects. I would expect an
             improvement of runtime performance using expression templates
             (especially if higher than first order is used).
        <li> Also about the previous item: possibly the call to a ValidPoint
             function of the reference geometry could hence be avoided, but I
             am not sure at the current stage (depends on whether some
             mappings can be used for different kinds of geometries).
    </ul> */
    template <unsigned int dimFrom, unsigned int dimTo>
    class MappingInterface
    {
        public:
            typedef Point<dimTo>                        PointToT;
            typedef PointPhysical<dimTo>                PointPhysicalT;
            typedef PointReference<dimFrom>             PointReferenceT;
            typedef Jacobian<dimFrom, dimTo>            JacobianT;

        public:
            /*! (OC): Jacobian has a gradient in each line, hence as many lines as target
                space (DIM2) and as many columns as original space (DIM1),
                \frac{\partial x_i}{\partial \xi_j}. */
            virtual void calcJacobian(const PointReferenceT&, JacobianT&) const = 0;
            MappingInterface() { }
            virtual ~MappingInterface() { };
    };
};
#endif /* MAP_H_ */

/*
 * MappingToPhysSimplexLinear.hpp
 *
 *  Created on: Feb 5, 2013
 *      Author: nicorivas
 */

#ifndef MAPPINGSIMPLEXLINEAR_H_
#define MAPPINGSIMPLEXLINEAR_H_

#include "MappingReferenceToPhysical.hpp"
#include "../GlobalNamespaceGeometry.hpp"
#include "../Jacobian.hpp"

namespace Geometry
{
    /*!
     * "In geometry, a simplex (plural simplexes or simplices) is a generalization of the notion of
     *  a triangle or tetrahedron to arbitrary dimension." -Wikipedia.
     *
     * This class defines the linear mappings between simplexes. See the comments in the Physical<Simplex>.cpp files to know the order of the
     * vertex of each simplex, an order which is kept by the mappings.
     * No specialization is needed because the mapping is general for geometries of vertex number
     * (vn) one greater than dimension (d), that is, vn = d+1;
     */

    template <unsigned int DIM>
    class MappingToPhysSimplexLinear: public MappingReferenceToPhysical<DIM,DIM>
    {
        private:
            typedef Geometry::PhysicalGeometry<DIM> PhysicalGeometryT;
            typedef Geometry::PointReference<DIM> PointReferenceT;
            typedef Geometry::PointPhysical<DIM> PointPhysicalT;
            typedef Geometry::Jacobian<DIM,DIM> JacobianT;

        public:
            MappingToPhysSimplexLinear(const PhysicalGeometryT*const& pG) { reinit(pG); };
            virtual void transform(const PointReferenceT&, PointPhysicalT&) const;
            virtual void calcJacobian(const PointReferenceT&, JacobianT&) const;
            virtual void reinit(const PhysicalGeometryT*const);

        private:
            //bool isValidPoint(const PointReferenceT&) const; //TODO: Implement this function.
            //! ~OC~
            //! In this case it is worth using an array for the mapping factors,
            //! since they are just difference vectors (see loop in reinit)
            PointPhysicalT a[DIM+1];
    };
    #include "MappingToPhysSimplexLinear_Impl.hpp"
};
#endif

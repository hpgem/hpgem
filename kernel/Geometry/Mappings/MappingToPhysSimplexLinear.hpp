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
     * \BUG this mapping does not work for DIM=1 (use MappingToPhysHypercubeLinear<1> instead)
     */

    template <unsigned int DIM>
    class MappingToPhysSimplexLinear: public MappingReferenceToPhysical
    {
        private:
            typedef Geometry::PhysicalGeometry PhysicalGeometryT;
            typedef Geometry::PointReference PointReferenceT;
            typedef Geometry::PointPhysical PointPhysicalT;
            typedef Geometry::Jacobian JacobianT;

        public:
            MappingToPhysSimplexLinear(const PhysicalGeometryT*const& pG):a(DIM+1,DIM){
                MappingReferenceToPhysical::setNodesPtr(&pG->getNodes()); reinit(pG); };
            virtual void transform(const PointReferenceT&, PointPhysicalT&) const;
            virtual void calcJacobian(const PointReferenceT&, JacobianT&) const;
            virtual void reinit(const PhysicalGeometryT*const);
            virtual int getTargetDimension() const {return DIM;}

        private:
            //bool isValidPoint(const PointReferenceT&) const; //TODO: Implement this function.
            //! ~OC~
            //! In this case it is worth using an array for the mapping factors,
            //! since they are just difference vectors (see loop in reinit)
            std::vector<PointPhysicalT> a;
    };
    #include "MappingToPhysSimplexLinear_Impl.hpp"
};
#endif

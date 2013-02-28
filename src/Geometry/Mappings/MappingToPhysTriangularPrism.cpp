/*
 * MappingToPhysTriangularPrism.hpp
 *
 *  Created on: Feb 5, 2013
 *      Author: nicorivas
 */

#ifndef MAPPINGTRIANGULARPRISM_CPP_
#define MAPPINGTRIANGULARPRISM_CPP_

#include "MappingToPhysTriangularPrism.hpp"

namespace Geometry
{
    MappingToPhysTriangularPrism::MappingToPhysTriangularPrism(const PhysicalGeometryT* const physicalGeometry)
    {
        reinit(physicalGeometry);
    }

    void MappingToPhysTriangularPrism::
    transform(const PointReferenceT& pR, PointPhysicalT& pP) const
    {
        if (isValidPoint(pR))
        {
#if SAVECOEFFS
            pointPhysical = a0
                          + a1 * xi[0]
                          + a2 * xi[1]
                          + a3 * xi[2]
                          + a4 * xi[0] * xi[2]
                          + a5 * xi[1] * xi[2];
#else
            const double t1 = pR[0] * pR[2];
            const double t2 = pR[1] * pR[2];
            LinearAlgebra::NumericalVector f2(6);

            f2[0] = 0.5 * (1.0 - pR[0] - pR[1] - pR[2] + t1 + t2);
            f2[1] = 0.5 * (pR[0] - t1);
            f2[2] = 0.5 * (pR[1] - t2);
            f2[3] = 0.5 * (-t2 + 1.0 - pR[0] - t1 - pR[1] + pR[2]);
            f2[4] = 0.5 * (t1 + pR[0]);
            f2[5] = 0.5 * (t2 + pR[1]);

            PointPhysicalT p;

            pP[0] = pP[1] = pP[2] = 0.;

            for (int i = 0; i < 6; ++i)
            {
                    // physicalGeometry_->getNodeCoordinates(i,p);
                    //pP += f2[i] * p;
            }
#endif
        }
        else
        {
            throw "ERROR: MappingToPhysTriangularPrism::transform, mapping point outside geometry.";
        }
    }

    void MappingToPhysTriangularPrism::
    calcJacobian(const PointReferenceT& pR, JacobianT& jacobian) const
    {
        if (isValidPoint(pR))
        {
#ifdef SAVECOEFFS
            Geometry::PointPhysical<3> d_dxi0(a1 + xi[2] * a4);
            Geometry::PointPhysical<3> d_dxi1(a2 + xi[2] * a5);
            Geometry::PointPhysical<3> d_dxi2(a3 + xi[0] * a4 + xi[1] * a5);

            for (int i = 0; i < 3; ++i)
            {
                jacobian(i,0) = d_dxi0[i];
                jacobian(i,1) = d_dxi1[i];
                jacobian(i,2) = d_dxi2[i];
            }
#else
            LinearAlgebra::NumericalVector df_dxi0(6), df_dxi1(6), df_dxi2(6);

            df_dxi0[0] = +0.5 * (-1. + pR[2]);
            df_dxi0[1] = +0.5 * (+1. - pR[2]);
            df_dxi0[2] = +0.0;
            df_dxi0[3] = +0.5 * (-1. - pR[2]);
            df_dxi0[4] = +0.5 * (+1. + pR[2]);
            df_dxi0[5] = +0.0;

            df_dxi1[0] = +0.5 * (-1. + pR[2]);
            df_dxi1[1] = +0.0;
            df_dxi1[2] = +0.5 * (+1. - pR[2]);
            df_dxi1[3] = +0.5 * (-1. - pR[2]);
            df_dxi1[4] = +0.0;
            df_dxi1[5] = +0.5 * (+1. + pR[2]);

            df_dxi2[0] = +0.5 * (-1. + pR[0] + pR[1]);
            df_dxi2[1] = -0.5 * pR[0];
            df_dxi2[2] = -0.5 * pR[1];
            df_dxi2[3] = +0.5 * (+1. - pR[1] - pR[0]);
            df_dxi2[4] = +0.5 * pR[0];
            df_dxi2[5] = +0.5 * pR[1];

            Geometry::PointPhysical<3> d_dxi0;
            Geometry::PointPhysical<3> d_dxi1;
            Geometry::PointPhysical<3> d_dxi2;

            for (int i = 0; i < 3; ++i)
            {
                d_dxi0[i] = 0.;
                d_dxi1[i] = 0.;
                d_dxi2[i] = 0.;
            }

            Geometry::PointPhysical<3> p;

            for (int i = 0; i < 6; ++i)
            {
                    //physicalGeometry_->getNodeCoordinates(i,p);

                d_dxi0 += df_dxi0[i] * p;
                d_dxi1 += df_dxi1[i] * p;
                d_dxi2 += df_dxi2[i] * p;
            }

            for (int i = 0; i < 3; ++i)
            {
                jacobian(i,0) = d_dxi0[i];
                jacobian(i,1) = d_dxi1[i];
                jacobian(i,2) = d_dxi2[i];
            }
        }
        else
        {
            // ERROR
        }
#endif
    }

    void
    MappingToPhysTriangularPrism::reinit(const PhysicalGeometryT*const physicalGeometry)
    {
#if SAVECOEFFS
        FixedVector<Geometry::PointPhysical<3>, 6> p;

        for (unsigned int i = 0; i < 6; ++i)
        {
            physicalGeometry->getVertexPoint(i, p[i]);
        }

        // see Maple/triangularprismMapping1.mws
        a0 = 0.5 * (p[0] + p[3]);
        a1 = 0.5 * (-p[0] + p[4] - p[3] + p[1]);
        a2 = 0.5 * (-p[0] + p[5] - p[3] + p[2]);
        a3 = 0.5 * (-p[0] + p[3]);
        a4 = 0.5 * (p[0] + p[4] - p[3] - p[1]);
        a5 = 0.5 * (p[0] + p[5] - p[3] - p[2]);
#else
            //  physicalGeometry_ = physicalGeometry;
#endif
    }

    bool MappingToPhysTriangularPrism::isValidPoint(const PointReferenceT& pointReference) const
    {
        if ((0. <= pointReference[0]) && (pointReference[0] <= 1.) &&
            (0. <= pointReference[1]) && (pointReference[1] <= 1.) &&
            (pointReference[0] + pointReference[1] - 1. <= 1.e-16) &&
            (pointReference[2] >= -1.) && (pointReference[2] <=  1.))
            return true;
        else
            return false;
    }
};
#endif

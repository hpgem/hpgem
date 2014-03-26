/*
 * MappingToPhysPyramid.cpp
 *
 *  Created on: Feb 5, 2013
 *      Author: nicorivas
 */

#ifndef MAPPINGPYRAMID_CPP_
#define MAPPINGPYRAMID_CPP_

#include "MappingToPhysPyramid.hpp"

namespace Geometry
{
    MappingToPhysPyramid::MappingToPhysPyramid(const PhysicalGeometryT*const physicalGeometry)
    {
        reinit(physicalGeometry);
    }

    void MappingToPhysPyramid::
    transform(const PointReferenceT& pR, PointPhysicalT& pP) const
    {
        if (isValidPoint(pR))
        {
            const double t1 = pR[0] * pR[1];
            const double t2 = pR[0] * pR[1]  * pR[2] / (1 - pR[2]);

            std::vector<double> f8;
            f8.resize(5);
            f8[0] = pR[2];
            f8[1] = 0.25 * ( 1. - pR[0] - pR[1] + t1 - pR[2] + t2 );
            f8[2] = 0.25 * ( 1. + pR[0] - pR[1] - t1 - pR[2] - t2 );
            f8[3] = 0.25 * ( 1. - pR[0] + pR[1] - t1 - pR[2] - t2 );
            f8[4] = 0.25 * ( 1. + pR[0] + pR[1] + t1 - pR[2] + t2 );

            PointPhysicalT p(3);

            pP[0] = pP[1] = pP[2] = 0.0;

            for (int i=0; i<5; ++i)
            {
                getNodeCoordinates(globalNodeIndices_[i], p);
                pP += f8[i] * p;
            }
        }
        else
        {
            throw "ERROR: MappingToPhysPyramid::transform, mapping point outside geometry.";
        }
    }
    void MappingToPhysPyramid::calcJacobian(const PointReferenceT& pR, JacobianT& jacobian) const
    {
        if (isValidPoint(pR))
        {
            std::vector<double> df_dxi0, df_dxi1, df_dxi2;

            const double dt6dx0 = pR[1] * pR[2] / (1. - pR[2]);
            const double dt6dx1 = pR[0] * pR[2] / (1. - pR[2]);
            const double dt6dx2 = pR[0] * pR[1] / ((1. - pR[2]) * (1. - pR[2]));

            df_dxi0[0] = 0.0;
            df_dxi0[1] = 0.25 * (-1. + pR[1] + dt6dx0);
            df_dxi0[2] = 0.25 * ( 1. - pR[1] - dt6dx0);
            df_dxi0[3] = 0.25 * (-1. - pR[1] - dt6dx0);
            df_dxi0[4] = 0.25 * ( 1. + pR[1] + dt6dx0);

            df_dxi1[0] = 0.0;
            df_dxi1[1] = 0.25 * (-1. + pR[0] + dt6dx1);
            df_dxi1[2] = 0.25 * (-1. - pR[0] - dt6dx1);
            df_dxi1[3] = 0.25 * ( 1. - pR[0] - dt6dx1);
            df_dxi1[4] = 0.25 * ( 1. + pR[0] + dt6dx1);

            df_dxi2[0] = 1.0;
            df_dxi2[1] = 0.25 * (-1. + dt6dx2);
            df_dxi2[2] = 0.25 * (-1. - dt6dx2);
            df_dxi2[3] = 0.25 * (-1. - dt6dx2);
            df_dxi2[4] = 0.25 * (-1. + dt6dx2);

            PointPhysicalT d_dxi0(3);
            PointPhysicalT d_dxi1(3);
            PointPhysicalT d_dxi2(3);

            for (unsigned int i = 0; i < 3; ++i)
            {
                d_dxi0[i] = 0.;
                d_dxi1[i] = 0.;
                d_dxi2[i] = 0.;
            }

            PointPhysicalT p(3);

            for (int i = 0; i < 5; ++i)
            {
                getNodeCoordinates(globalNodeIndices_[i], p);

                d_dxi0 += df_dxi0[i] * p;
                d_dxi1 += df_dxi1[i] * p;
                d_dxi2 += df_dxi2[i] * p;
            }

            for (unsigned int i = 0; i < 3; ++i)
            {
                jacobian(i,0) = d_dxi0[i];
                jacobian(i,1) = d_dxi1[i];
                jacobian(i,2) = d_dxi2[i];
            }
        }
        else
        {
            throw "ERROR: MappingToPhysPyramid::calcJacobian, mapping point outside geometry.";
        }
    }

    void MappingToPhysPyramid::reinit(const PhysicalGeometryT*const physicalGeometry)
    {
        for (int i = 0; i < 5; ++i)
        {
            globalNodeIndices_[i] = physicalGeometry->getNodeIndex(i);
        }
    }

    bool MappingToPhysPyramid::isValidPoint(const PointReferenceT& pointReference) const
    {
        static const double eps = 1.e-14;
        const double z = pointReference[2];
        if ((std::abs(pointReference[0]) <= 1. - z + eps)
                && (std::abs(pointReference[1]) <= 1. - z + eps)
                && (z >= 0. - eps)
                && (z <= 1. + eps))
        {
            return true;
        }
        else
        {
            return false;
        }
    }

};
#endif

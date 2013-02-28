/*
 * MappingSimpleCubeNLinear.hpp
 *
 *  Created on: Feb 5, 2013
 *      Author: nicorivas
 */

#ifndef MAPPINGHYPERCUBELINEAR_CPP_
#define MAPPINGHYPERCUBELINEAR_CPP_

#include "MappingToPhysHypercubeLinear.hpp"

namespace Geometry
{
    // =============================================================================================
    // ~~~ Dimension 1 ~~~==========================================================================
    // =============================================================================================
    MappingToPhysHypercubeLinear<1>::
    MappingToPhysHypercubeLinear(const PhysicalGeometryT*const& physicalGeometry)
    {
        mid = slope = 0.0;
        reinit(physicalGeometry);
    }

    void MappingToPhysHypercubeLinear<1>::
    transform(const PointReferenceT& pointReference, PointPhysicalT& pointPhysical) const
    {
        if (isValidPoint(pointReference))
        {
            pointPhysical[0] = mid + pointReference[0] * slope;
        }
        else
        {
            // ERROR
        }
    }

    void MappingToPhysHypercubeLinear<1>::reinit(const PhysicalGeometryT*const physicalGeometry)
    {
        PointPhysicalT p0;
        PointPhysicalT p1;
        physicalGeometry->getNodeCoordinates(0, p0);
        physicalGeometry->getNodeCoordinates(1, p1);
        mid   = 0.5 * (p1[0] + p0[0]);
        slope = 0.5 * (p1[0] - p0[0]);
    }

    bool MappingToPhysHypercubeLinear<1>::isValidPoint(const PointReferenceT& pointReference) const
    {
        if ((pointReference[0] < -1.) || (pointReference[0] > 1.))
            return false;
        else
            return true;
    }

    void MappingToPhysHypercubeLinear<1>::calcJacobian(const PointReferenceT&, JacobianT&) const
    {
        // No Jacobian in 1D.
    }

    // =============================================================================================
    // ~~~ Dimension 2 ~~~==========================================================================
    // =============================================================================================

    MappingToPhysHypercubeLinear<2>::
    MappingToPhysHypercubeLinear(const PhysicalGeometryT*const& physicalGeometry):
    MappingReferenceToPhysical<2,2>()
    {
        reinit(physicalGeometry);
    }
    void MappingToPhysHypercubeLinear<2>::
    transform(const PointReferenceT& pointReference, PointPhysicalT& pointPhysical) const
    {
        if (isValidPoint(pointReference))
        {
            pointPhysical = a0 + pointReference[0] * a1
                               + pointReference[1] * (a2 + pointReference[0] * a12);
        }
        else
        {
            // ERROR
        }
    }

    void MappingToPhysHypercubeLinear<2>::
    calcJacobian(const PointReferenceT& pointReference, JacobianT& jacobian) const
    {
        if (isValidPoint(pointReference))
        {
            jacobian(0,0) = a1[0] + pointReference[1] * a12[0];
            jacobian(0,1) = a2[0] + pointReference[0] * a12[0];
            jacobian(1,0) = a1[1] + pointReference[1] * a12[1];
            jacobian(1,1) = a2[1] + pointReference[0] * a12[1];
        }
        else
        {
            // ERROR
        }
    }

    void MappingToPhysHypercubeLinear<2>::reinit(const PhysicalGeometryT*const physicalGeometry)
    {
        PointPhysicalT p0,p1,p2,p3;
        physicalGeometry->getNodeCoordinates(0, p0);
        physicalGeometry->getNodeCoordinates(1, p1);
        physicalGeometry->getNodeCoordinates(2, p2);
        physicalGeometry->getNodeCoordinates(3, p3);

        a0  = 0.25 * (p0 + p1 + p2 + p3);
        a1  = 0.25 * (p1 - p0 + p3 - p2);
        a2  = 0.25 * (p2 - p0 + p3 - p1);
        a12 = 0.25 * (p3 - p1 + p0 - p2);

    }

    bool MappingToPhysHypercubeLinear<2>::isValidPoint(const PointReferenceT& pointReference) const
    {
        if ((pointReference[0] < -1.) || (pointReference[0] > 1.) ||
            (pointReference[1] < -1.) || (pointReference[1] > 1.))
            return false;
        else
            return true;
    }

    // =============================================================================================
    // ~~~ Dimension 3 ~~~==========================================================================
    // =============================================================================================

    MappingToPhysHypercubeLinear<3>
    ::MappingToPhysHypercubeLinear(const PhysicalGeometryT*const& physicalGeometry)
    {
        reinit(physicalGeometry);
    }

    void MappingToPhysHypercubeLinear<3>::
    transform(const PointReferenceT& pR, PointPhysicalT& pointPhysical) const
    {
        if (isValidPoint(pR))
        {
            pointPhysical = a0 + pR[0] * (a1 + pR[1] * (a12 +  pR[2] * a123) + pR[2] * a13)
                               + pR[1] * (a2 + pR[2] * a23)
                               + pR[2] * a3;
        }
        else
        {
            // ERROR
        }
    }

    void MappingToPhysHypercubeLinear<3>::
    calcJacobian(const PointReferenceT& pR, JacobianT& jacobian) const
    {
        if (isValidPoint(pR))
        {
            for(int i = 0; i < 3; ++i)
            {
                double pR01 = pR[0] * pR[1];
                double pR02 = pR[0] * pR[2];
                double pR12 = pR[1] * pR[2];

                jacobian(0,0) = a1[0] + pR[1] * a12[0] + pR[2] * a13[0] + pR12 * a123[0];
                jacobian(0,1) = a2[0] + pR[0] * a12[0] + pR[2] * a23[0] + pR02 * a123[0];
                jacobian(0,2) = a3[0] + pR[0] * a13[0] + pR[1] * a23[0] + pR01 * a123[0];

                jacobian(1,0) = a1[1] + pR[1] * a12[1] + pR[2] * a13[1] + pR12 * a123[1];
                jacobian(1,1) = a2[1] + pR[0] * a12[1] + pR[2] * a23[1] + pR02 * a123[1];
                jacobian(1,2) = a3[1] + pR[0] * a13[1] + pR[1] * a23[1] + pR01 * a123[1];

                jacobian(2,0) = a1[2] + pR[1] * a12[2] + pR[2] * a13[2] + pR12 * a123[2];
                jacobian(2,1) = a2[2] + pR[0] * a12[2] + pR[2] * a23[2] + pR02 * a123[2];
                jacobian(2,2) = a3[2] + pR[0] * a13[2] + pR[1] * a23[2] + pR01 * a123[2];
            }
        }
        else
        {
            // ERROR
        }
    }

    void MappingToPhysHypercubeLinear<3>::
    reinit(const PhysicalGeometryT*const physicalGeometry)
    {
        PointPhysicalT p0,p1,p2,p3,p4,p5,p6,p7;

        physicalGeometry->getNodeCoordinates(0, p0);
        physicalGeometry->getNodeCoordinates(1, p1);
        physicalGeometry->getNodeCoordinates(2, p2);
        physicalGeometry->getNodeCoordinates(3, p3);
        physicalGeometry->getNodeCoordinates(4, p4);
        physicalGeometry->getNodeCoordinates(5, p5);
        physicalGeometry->getNodeCoordinates(6, p6);
        physicalGeometry->getNodeCoordinates(7, p7);

        a0   = 0.125 * (p0 + p1 + p2 + p3 + p4 + p5 + p6 + p7);
        a1   = 0.125 * (p1 - p0 + p3 - p2 + p5 - p4 + p7 - p6);
        a2   = 0.125 * (p2 - p0 + p3 - p1 + p6 - p4 + p7 - p5);
        a3   = 0.125 * (p4 - p0 + p5 - p1 + p6 - p2 + p7 - p3);
        a12  = 0.125 * (p0 - p1 + p3 - p2 + p4 - p5 + p7 - p6);
        a13  = 0.125 * (p0 - p1 + p2 - p3 + p5 - p4 + p7 - p6);
        a23  = 0.125 * (p0 - p2 + p1 - p3 + p6 - p4 + p7 - p5);
        a123 = 0.125 * (p1 - p0 + p2 - p3 + p4 - p5 + p7 - p6);
    }

    bool MappingToPhysHypercubeLinear<3>::
    isValidPoint(const PointReferenceT& pointReference) const
    {
        if ((pointReference[0] < -1.) || (pointReference[0] > 1.) ||
            (pointReference[1] < -1.) || (pointReference[1] > 1.) ||
            (pointReference[2] < -1.) || (pointReference[2] > 1.))
                return false;
            else
                return true;
    }

    // =============================================================================================
    // ~~~ Dimension 4 ~~~==========================================================================
    // =============================================================================================

    MappingToPhysHypercubeLinear<4>::
    MappingToPhysHypercubeLinear(const PhysicalGeometryT* const& physicalGeometry)
    {
        reinit(physicalGeometry);
    }

    void MappingToPhysHypercubeLinear<4>::
    transform(const PointReferenceT& pR, PointPhysicalT& pointPhysical) const
    {
        if (isValidPoint(pR))
        {
            pointPhysical = abar + pR[0] * a0 + pR[1] * a1 + pR[2] * a2 + pR[3] * a3;
        }
        else
        {
            // ERROR
        }
    }

    void MappingToPhysHypercubeLinear<4>::
    calcJacobian(const PointReferenceT& pR, JacobianT& jacobian) const
    {
        if (isValidPoint(pR))
        {
            for(int i = 0; i < 4; ++i)
            {
                jacobian(i,0) = a0[i]   + a01[i] * pR[1] + a02[i]   * pR[2] + a03[i] * pR[3]
                              + a012[i] * pR[1]  * pR[2] + a013[i]  * pR[1] * pR[3]
                              + a230[i] * pR[2]  * pR[3] + a0123[i] * pR[1] * pR[2]  * pR[3];

                jacobian(i,1) = a1[i]   + a01[i] * pR[0] + a12[i]   * pR[2] + a13[i] * pR[3]
                              + a012[i] * pR[0]  * pR[2] + a013[i]  * pR[0] * pR[3]
                              + a123[i] * pR[2]  * pR[3] + a0123[i] * pR[0] * pR[2]  * pR[3];

                jacobian(i,2) = a2[i]   + a02[i] * pR[0] + a12[i]   * pR[1] + a23[i] * pR[3]
                              + a012[i] * pR[0]  * pR[1] + a230[i]  * pR[0] * pR[3]
                              + a123[i] * pR[1]  * pR[3] + a0123[i] * pR[0] * pR[1]  * pR[3];

                jacobian(i,3) = a3[i]   + a03[i] * pR[0] + a13[i]   * pR[1] + a23[i] * pR[2]
                              + a013[i] * pR[0]  * pR[1] + a230[i]  * pR[0] * pR[2]
                              + a123[i] * pR[1]  * pR[2] + a0123[i] * pR[0] * pR[1]  * pR[2];
            }
        }
        else
        {
            // ERROR
        }
    }

    void MappingToPhysHypercubeLinear<4>::reinit(const PhysicalGeometryT* const physicalGeometry)
    {
        PointPhysicalT P[16];
        for (int i = 0; i < 16; ++i) physicalGeometry->getNodeCoordinates(i, P[i]);

        abar = 0.0625 *
            ( P[0] + P[1] + P[2]  + P[3]  + P[4]  + P[5]  + P[6]  + P[7]
             +P[8] + P[9] + P[10] + P[11] + P[12] + P[13] + P[14] + P[15]);

        a0 = 0.0625 *
            (-P[0] + P[1] - P[2]  + P[3]  - P[4]  + P[5]  - P[6]  + P[7]
             -P[8] + P[9] - P[10] + P[11] - P[12] + P[13] - P[14] + P[15]);

        a1 = 0.0625 *
            (-P[0] - P[1] + P[2]  + P[3]  - P[4]  - P[5]  + P[6]  + P[7]
             -P[8] - P[9] + P[10] + P[11] - P[12] - P[13] + P[14] + P[15]);

        a2 = 0.0625 *
            (-P[0] - P[1] - P[2]  - P[3]  + P[4]  + P[5]  + P[6]  + P[7]
             -P[8] - P[9] - P[10] - P[11] + P[12] + P[13] + P[14] + P[15]);

        a3 = 0.0625 *
            (-P[0] - P[1] - P[2]  - P[3]  - P[4]  - P[5]  - P[6]  - P[7]
             +P[8] + P[9] + P[10] + P[11] + P[12] + P[13] + P[14] + P[15]);
        // Only supporting 4-cubes for the moment
    }

    bool MappingToPhysHypercubeLinear<4>::isValidPoint(const PointReferenceT& pointReference) const
    {
        if ((pointReference[0] < -1.) || (pointReference[0] > 1.) ||
            (pointReference[1] < -1.) || (pointReference[1] > 1.) ||
            (pointReference[2] < -1.) || (pointReference[2] > 1.) ||
            (pointReference[3] < -1.) || (pointReference[3] > 1.))
                return false;
            else
                return true;
    }
};
#endif /* MAPPINGSIMPLECUBENLINEAR_H_ */

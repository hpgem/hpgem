/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found below.


 Copyright (c) 2014, University of Twente
 All rights reserved.

 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */


#ifndef MAPPINGHYPERCUBELINEAR_CPP_
#define MAPPINGHYPERCUBELINEAR_CPP_

#include "MappingToPhysHypercubeLinear.hpp"
#include <vector>
#include "Geometry/PhysicalGeometry.hpp"
#include "Geometry/PointReference.hpp"
#include "Geometry/Jacobian.hpp"

namespace Geometry
{
    // =============================================================================================
    // ~~~ Dimension 1 ~~~==========================================================================
    // =============================================================================================
    MappingToPhysHypercubeLinear<1>::
    MappingToPhysHypercubeLinear(const PhysicalGeometryT*const& physicalGeometry)
    {
        mid = slope = 0.0;
        MappingReferenceToPhysical::setNodesPtr(&physicalGeometry->getNodes());
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
            // ERROR ///\TODO emit a warning at a non-deadly priority level
            pointPhysical[0] = mid + pointReference[0] * slope;
        }
    }

    void MappingToPhysHypercubeLinear<1>::reinit(const PhysicalGeometryT*const physicalGeometry)
    {
        PointPhysicalT p0(1);
        PointPhysicalT p1(1);
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

    void MappingToPhysHypercubeLinear<1>::calcJacobian(const PointReferenceT&, JacobianT& jac) const
    {
        jac[0] = slope;
    }

    // =============================================================================================
    // ~~~ Dimension 2 ~~~==========================================================================
    // =============================================================================================

    MappingToPhysHypercubeLinear<2>::
    MappingToPhysHypercubeLinear(const PhysicalGeometryT*const& physicalGeometry):
    MappingReferenceToPhysical(),a0(2),a1(2),a12(2),a2(2)
    {
        MappingReferenceToPhysical::setNodesPtr(&physicalGeometry->getNodes());
        reinit(physicalGeometry);
    }
    void MappingToPhysHypercubeLinear<2>::
    transform(const PointReferenceT& pointReference, PointPhysicalT& pointPhysical) const
    {
    	//assert(L2norm(a1)>1e-14&&L2norm(a2 + pointReference[0] * a12)>1e-14);
        if (isValidPoint(pointReference))
        {
            pointPhysical = a0 + pointReference[0] * a1
                               + pointReference[1] * (a2 + pointReference[0] * a12);
        }
        else
        {
            // ERROR///\TODO emit a warning at a non-deadly priority level
            pointPhysical = a0 + pointReference[0] * a1
                               + pointReference[1] * (a2 + pointReference[0] * a12);
        }
    }

    void MappingToPhysHypercubeLinear<2>::
    calcJacobian(const PointReferenceT& pointReference, JacobianT& jacobian) const
    {
        //if (isValidPoint(pointReference))
        //{
            jacobian(0,0) = a1[0] + pointReference[1] * a12[0];
            jacobian(0,1) = a2[0] + pointReference[0] * a12[0];
            jacobian(1,0) = a1[1] + pointReference[1] * a12[1];
            jacobian(1,1) = a2[1] + pointReference[0] * a12[1];
        //}
        //else
        //{
        //    // ERROR///\TODO emit a warning at a non-deadly priority level
        //}
    }

    void MappingToPhysHypercubeLinear<2>::reinit(const PhysicalGeometryT*const physicalGeometry)
    {
        PointPhysicalT p0(2),p1(2),p2(2),p3(2);
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
    ::MappingToPhysHypercubeLinear(const PhysicalGeometryT*const& physicalGeometry):
     a1(3),a12(3),a13(3),a123(3),a0(3),a2(3),a23(3),a3(3)
    {
        MappingReferenceToPhysical::setNodesPtr(&physicalGeometry->getNodes());
        reinit(physicalGeometry);
    }

    void MappingToPhysHypercubeLinear<3>::
    transform(const PointReferenceT& pR, PointPhysicalT& pointPhysical) const
    {
    	//assert(L2norm(a3)>1e-14&&L2norm(a2 + pR[2] * a23)>1e-14&&L2norm(a1 + pR[1]*(a12+pR[2]*a123) pR[2] * a23)>1e-14);
        if (isValidPoint(pR))
        {
            pointPhysical = a0 + pR[0] * (a1 + pR[1] * (a12 +  pR[2] * a123) + pR[2] * a13)
                               + pR[1] * (a2 + pR[2] * a23)
                               + pR[2] * a3;
        }
        else
        {
            // ERROR
            pointPhysical = a0 + pR[0] * (a1 + pR[1] * (a12 +  pR[2] * a123) + pR[2] * a13)
                               + pR[1] * (a2 + pR[2] * a23)
                               + pR[2] * a3;
        }
    }

    void MappingToPhysHypercubeLinear<3>::
    calcJacobian(const PointReferenceT& pR, JacobianT& jacobian) const
    {
        //if (isValidPoint(pR))
        //{
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
        //}
        //else
        //{
        //    // ERROR
        //}
    }

    void MappingToPhysHypercubeLinear<3>::
    reinit(const PhysicalGeometryT*const physicalGeometry)
    {
        PointPhysicalT p0(3),p1(3),p2(3),p3(3),p4(3),p5(3),p6(3),p7(3);

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
    MappingToPhysHypercubeLinear(const PhysicalGeometryT* const& physicalGeometry):
    abar(4),a0(4),a01(4),a02(4),a03(4),a012(4),a013(4),a0123(4),a1(4),a12(4),a13(4),a123(4),a2(4),a23(4),a230(4),a3(4)
    {
        MappingReferenceToPhysical::setNodesPtr(&physicalGeometry->getNodes());
        reinit(physicalGeometry);
    }

    void MappingToPhysHypercubeLinear<4>::
    transform(const PointReferenceT& pR, PointPhysicalT& pointPhysical) const
    {
        //if (isValidPoint(pR))
        //{
            pointPhysical = abar + pR[0] * a0 + pR[1] * a1 + pR[2] * a2 + pR[3] * a3;
        //}
        //else
        //{
            // ERROR
        //}
    }

    void MappingToPhysHypercubeLinear<4>::
    calcJacobian(const PointReferenceT& pR, JacobianT& jacobian) const
    {
    	//assert(...)
        //if (isValidPoint(pR))
        //{
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
        //}
        //else
        //{
        //    // ERROR
        //}
    }

    void MappingToPhysHypercubeLinear<4>::reinit(const PhysicalGeometryT* const physicalGeometry)
    {
        std::vector<PointPhysicalT> P(16,4);
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

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

#include "MappingToPhysHypercubeLinear.h"
#include <vector>
#include "Geometry/PhysicalGeometry.h"
#include "Geometry/PointReference.h"
#include "Geometry/Jacobian.h"
#include "Base/L2Norm.h"

namespace Geometry
{
    // =============================================================================================
    // ~~~ Dimension 1 ~~~==========================================================================
    // =============================================================================================
    MappingToPhysHypercubeLinear<1>::MappingToPhysHypercubeLinear(const PhysicalGeometry<1>* const & physicalGeometry)
            : MappingReferenceToPhysical(physicalGeometry)
    {
        logger.assert(physicalGeometry!=nullptr, "Invalid physical geometry passed");
        mid = slope = 0.0;
        reinit();
    }
    
    PointPhysical<1> MappingToPhysHypercubeLinear<1>::transform(const PointReference<1>& pointReference) const
    {
        PointPhysical<1> pointPhysical;
        if (isValidPoint(pointReference))
        {
            pointPhysical[0] = mid + pointReference[0] * slope;
        }
        else
        {
            logger(WARN, "In MappingToPhysHypercubeLinear<1>::transform, the given PointReference is not between -1 and 1.");
            pointPhysical[0] = mid + pointReference[0] * slope;
        }
        return pointPhysical;
    }
    
    PointReference<1> MappingToPhysHypercubeLinear<1>::inverseTransform(const PointPhysical<1>& pointPhysical) const
    {
        return {(pointPhysical[0] - mid) / slope};
    }

    void MappingToPhysHypercubeLinear<1>::reinit()
    {
        PointPhysical<1> p0 = this->geometry->getLocalNodeCoordinates(0);
        PointPhysical<1> p1 = this->geometry->getLocalNodeCoordinates(1);
        mid = 0.5 * (p1[0] + p0[0]);
        slope = 0.5 * (p1[0] - p0[0]);
    }
    
    bool MappingToPhysHypercubeLinear<1>::isValidPoint(const PointReference<1>& pointReference) const
    {
        if ((pointReference[0] < -1.) || (pointReference[0] > 1.))
            return false;
        else
            return true;
    }
    
    Jacobian<1, 1> MappingToPhysHypercubeLinear<1>::calcJacobian(const PointReference<1>& pointReference) const
    {
        Jacobian<1, 1> jac;
        jac[0] = slope;
        return jac;
    }
    
    // =============================================================================================
    // ~~~ Dimension 2 ~~~==========================================================================
    // =============================================================================================
    
    MappingToPhysHypercubeLinear<2>::MappingToPhysHypercubeLinear(const PhysicalGeometry<2>* const & physicalGeometry)
            : MappingReferenceToPhysical(physicalGeometry)
    {
        logger.assert(physicalGeometry!=nullptr, "Invalid physical geometry passed");
        reinit();
    }

    PointPhysical<2> MappingToPhysHypercubeLinear<2>::transform(const PointReference<2>& pointReference) const
    {
        PointPhysical<2> pointPhysical;
        if (isValidPoint(pointReference))
        {
            //In the part below, we compute: pointPhysical = a0 + pointReference[0] * a1
            //                   + pointReference[1] * (a2 + pointReference[0] * a12);
            pointPhysical = a2;
            pointPhysical.axpy(pointReference[0], a12);
            pointPhysical *= pointReference[1];
            pointPhysical += a0;
            pointPhysical.axpy(pointReference[0], a1);
        }
        else
        {
            //logger(WARN, "In MappingToPhysHypercubeLinear<2>::transform, the "
            //        "given pointReference is not in the square [-1,1]x[-1,1]");
            //In the part below, we compute: pointPhysical = a0 + pointReference[0] * a1
            //                   + pointReference[1] * (a2 + pointReference[0] * a12);
            pointPhysical = a2;
            pointPhysical.axpy(pointReference[0], a12);
            pointPhysical *= pointReference[1];
            pointPhysical += a0;
            pointPhysical.axpy(pointReference[0], a1);
        }
        return pointPhysical;
    }

    PointReference<2> MappingToPhysHypercubeLinear<2>::inverseTransform(const PointPhysical<2>& pointPhysical) const
    {
        Geometry::PointReference<2> result;
        Geometry::PointPhysical<2> comparison = transform(result);
        LinearAlgebra::SmallVector<2> correction;
        double error = Base::L2Norm(pointPhysical - comparison);
        while(error > 1e-14)
        {
            correction = (pointPhysical - comparison).getCoordinates();
            calcJacobian(result).solve(correction);
            result = PointReference<2>(result + correction);
            comparison = transform(result);
            error = Base::L2Norm(pointPhysical - comparison);
        }
        return result;
    }
    
    Jacobian<2, 2> MappingToPhysHypercubeLinear<2>::calcJacobian(const PointReference<2>& pointReference) const
    {
        logger.assert(pointReference.size()==2, "Reference point has the wrong dimension");
        Jacobian<2, 2> jacobian;
        jacobian(0, 0) = a1[0] + pointReference[1] * a12[0];
        jacobian(0, 1) = a2[0] + pointReference[0] * a12[0];
        jacobian(1, 0) = a1[1] + pointReference[1] * a12[1];
        jacobian(1, 1) = a2[1] + pointReference[0] * a12[1];
        return jacobian;
    }
    
    void MappingToPhysHypercubeLinear<2>::reinit()
    {
        
        //this routine computes the following quantities
        // a0  = 0.25 * (p0 + p1 + p2 + p3);
        // a1  = 0.25 * (p1 - p0 + p3 - p2);
        // a2  = 0.25 * (p2 - p0 + p3 - p1);
        // a12 = 0.25 * (p3 - p1 + p0 - p2);
        

        a0 = geometry->getLocalNodeCoordinates(0);
        a1 = geometry->getLocalNodeCoordinates(1);
        a2 = geometry->getLocalNodeCoordinates(2);
        PointPhysical<2> temp = geometry->getLocalNodeCoordinates(3);
        a0 += temp;
        a12 = a1;
        a12 += a2;
        a12 -= a0;
        a12 *= -0.25;
        a0 += a1;
        a0 += a2;
        a1 += temp;
        a2 += temp;
        a1.axpy(-0.5, a0);
        a2.axpy(-0.5, a0);
        a0 *= 0.25;
        a1 *= 0.5;
        a2 *= 0.5;
        
    }
    
    bool MappingToPhysHypercubeLinear<2>::isValidPoint(const PointReference<2>& pointReference) const
    {
        if ((pointReference[0] < -1.) || (pointReference[0] > 1.) || (pointReference[1] < -1.) || (pointReference[1] > 1.))
            return false;
        else
            return true;
    }
    
    // =============================================================================================
    // ~~~ Dimension 3 ~~~==========================================================================
    // =============================================================================================
    
    MappingToPhysHypercubeLinear<3>::MappingToPhysHypercubeLinear(const PhysicalGeometry<3>* const & physicalGeometry)
            : MappingReferenceToPhysical(physicalGeometry)
    {
        logger.assert(physicalGeometry!=nullptr, "Invalid physical geometry passed");
        reinit();
    }
    
    PointPhysical<3> MappingToPhysHypercubeLinear<3>::transform(const PointReference<3>& pR) const
    {
        if (isValidPoint(pR))
        {
            return a0 + pR[0] * (a1 + pR[1] * (a12 + pR[2] * a123) + pR[2] * a13) + pR[1] * (a2 + pR[2] * a23) + pR[2] * a3;
        }
        else
        {
            // ERROR
            return a0 + pR[0] * (a1 + pR[1] * (a12 + pR[2] * a123) + pR[2] * a13) + pR[1] * (a2 + pR[2] * a23) + pR[2] * a3;
        }
    }
    
    PointReference<3> MappingToPhysHypercubeLinear<3>::inverseTransform(const PointPhysical<3>& pointPhysical) const
    {
        Geometry::PointReference<3> result;
        Geometry::PointPhysical<3> comparison = transform(result);
        LinearAlgebra::SmallVector<3> correction;
        double error = Base::L2Norm(pointPhysical - comparison);
        while(error > 1e-14)
        {
            correction = (pointPhysical - comparison).getCoordinates();
            calcJacobian(result).solve(correction);
            result = PointReference<3>(result + correction);
            comparison = transform(result);
            error = Base::L2Norm(pointPhysical - comparison);
        }
        return result;
    }

    Jacobian<3, 3> MappingToPhysHypercubeLinear<3>::calcJacobian(const PointReference<3>& pR) const
    {
        Jacobian<3, 3> jacobian;
        for (std::size_t i = 0; i < 3; ++i)
        {
            double pR01 = pR[0] * pR[1];
            double pR02 = pR[0] * pR[2];
            double pR12 = pR[1] * pR[2];
            
            jacobian(0, 0) = a1[0] + pR[1] * a12[0] + pR[2] * a13[0] + pR12 * a123[0];
            jacobian(0, 1) = a2[0] + pR[0] * a12[0] + pR[2] * a23[0] + pR02 * a123[0];
            jacobian(0, 2) = a3[0] + pR[0] * a13[0] + pR[1] * a23[0] + pR01 * a123[0];
            
            jacobian(1, 0) = a1[1] + pR[1] * a12[1] + pR[2] * a13[1] + pR12 * a123[1];
            jacobian(1, 1) = a2[1] + pR[0] * a12[1] + pR[2] * a23[1] + pR02 * a123[1];
            jacobian(1, 2) = a3[1] + pR[0] * a13[1] + pR[1] * a23[1] + pR01 * a123[1];
            
            jacobian(2, 0) = a1[2] + pR[1] * a12[2] + pR[2] * a13[2] + pR12 * a123[2];
            jacobian(2, 1) = a2[2] + pR[0] * a12[2] + pR[2] * a23[2] + pR02 * a123[2];
            jacobian(2, 2) = a3[2] + pR[0] * a13[2] + pR[1] * a23[2] + pR01 * a123[2];
        }
        return jacobian;
    }
    
    void MappingToPhysHypercubeLinear<3>::reinit()
    {
        PointPhysical<3> p0 = geometry->getLocalNodeCoordinates(0);
        PointPhysical<3> p1 = geometry->getLocalNodeCoordinates(1);
        PointPhysical<3> p2 = geometry->getLocalNodeCoordinates(2);
        PointPhysical<3> p3 = geometry->getLocalNodeCoordinates(3);
        PointPhysical<3> p4 = geometry->getLocalNodeCoordinates(4);
        PointPhysical<3> p5 = geometry->getLocalNodeCoordinates(5);
        PointPhysical<3> p6 = geometry->getLocalNodeCoordinates(6);
        PointPhysical<3> p7 = geometry->getLocalNodeCoordinates(7);
        
        a0 = 0.125 * (p0 + p1 + p2 + p3 + p4 + p5 + p6 + p7);
        a1 = 0.125 * (p1 - p0 + p3 - p2 + p5 - p4 + p7 - p6);
        a2 = 0.125 * (p2 - p0 + p3 - p1 + p6 - p4 + p7 - p5);
        a3 = 0.125 * (p4 - p0 + p5 - p1 + p6 - p2 + p7 - p3);
        a12 = 0.125 * (p0 - p1 + p3 - p2 + p4 - p5 + p7 - p6);
        a13 = 0.125 * (p0 - p1 + p2 - p3 + p5 - p4 + p7 - p6);
        a23 = 0.125 * (p0 - p2 + p1 - p3 + p6 - p4 + p7 - p5);
        a123 = 0.125 * (p1 - p0 + p2 - p3 + p4 - p5 + p7 - p6);
    }
    
    bool MappingToPhysHypercubeLinear<3>::isValidPoint(const PointReference<3>& pointReference) const
    {
        logger.assert(pointReference.size()==3, "Reference point has the wrong dimension");
        if ((pointReference[0] < -1.) || (pointReference[0] > 1.) || (pointReference[1] < -1.) || (pointReference[1] > 1.) || (pointReference[2] < -1.) || (pointReference[2] > 1.))
            return false;
        else
            return true;
    }
    
    // =============================================================================================
    // ~~~ Dimension 4 ~~~==========================================================================
    // =============================================================================================
    
    MappingToPhysHypercubeLinear<4>::MappingToPhysHypercubeLinear(const PhysicalGeometry<4>* const & physicalGeometry)
            : MappingReferenceToPhysical(physicalGeometry)
    {
        logger.assert(physicalGeometry!=nullptr, "Invalid physical geometry passed");
        reinit();
    }
    
    PointPhysical<4> MappingToPhysHypercubeLinear<4>::transform(const PointReference<4>& pR) const
    {
        logger.assert(pR.size()==4, "Reference point has the wrong dimension");
        return abar + pR[0] * a0 + pR[1] * a1 + pR[2] * a2 + pR[3] * a3;
    }
    
    PointReference<4> MappingToPhysHypercubeLinear<4>::inverseTransform(const PointPhysical<4>& pointPhysical) const
    {
        Geometry::PointReference<4> result;
        Geometry::PointPhysical<4> comparison = transform(result);
        LinearAlgebra::SmallVector<4> correction;
        double error = Base::L2Norm(pointPhysical - comparison);
        while(error > 1e-14)
        {
            correction = (pointPhysical - comparison).getCoordinates();
            calcJacobian(result).solve(correction);
            result = PointReference<4>(result + correction);
            comparison = transform(result);
            error = Base::L2Norm(pointPhysical - comparison);
        }
        return result;
    }

    Jacobian<4, 4> MappingToPhysHypercubeLinear<4>::calcJacobian(const PointReference<4>& pR) const
    {
        logger.assert(pR.size()==4, "Reference point has the wrong dimension");
        Jacobian<4, 4> jacobian;
        for (std::size_t i = 0; i < 4; ++i)
        {
            jacobian(i, 0) = a0[i];// + a01[i] * pR[1] + a02[i] * pR[2] + a03[i] * pR[3] + a012[i] * pR[1] * pR[2] + a013[i] * pR[1] * pR[3] + a230[i] * pR[2] * pR[3] + a0123[i] * pR[1] * pR[2] * pR[3];
            
            jacobian(i, 1) = a1[i];// + a01[i] * pR[0] + a12[i] * pR[2] + a13[i] * pR[3] + a012[i] * pR[0] * pR[2] + a013[i] * pR[0] * pR[3] + a123[i] * pR[2] * pR[3] + a0123[i] * pR[0] * pR[2] * pR[3];
            
            jacobian(i, 2) = a2[i];// + a02[i] * pR[0] + a12[i] * pR[1] + a23[i] * pR[3] + a012[i] * pR[0] * pR[1] + a230[i] * pR[0] * pR[3] + a123[i] * pR[1] * pR[3] + a0123[i] * pR[0] * pR[1] * pR[3];
            
            jacobian(i, 3) = a3[i];// + a03[i] * pR[0] + a13[i] * pR[1] + a23[i] * pR[2] + a013[i] * pR[0] * pR[1] + a230[i] * pR[0] * pR[2] + a123[i] * pR[1] * pR[2] + a0123[i] * pR[0] * pR[1] * pR[2];
        }
        return jacobian;
    }
    
    void MappingToPhysHypercubeLinear<4>::reinit()
    {
        std::array<PointPhysical<4>, 16> P;
        for (std::size_t i = 0; i < 16; ++i)
        {
            P[i] = geometry->getLocalNodeCoordinates(i);
        }
        
        abar = 0.0625 * (P[0] + P[1] + P[2] + P[3] + P[4] + P[5] + P[6] + P[7] + P[8] + P[9] + P[10] + P[11] + P[12] + P[13] + P[14] + P[15]);
        
        a0 = 0.0625 * (-P[0] + P[1] - P[2] + P[3] - P[4] + P[5] - P[6] + P[7] - P[8] + P[9] - P[10] + P[11] - P[12] + P[13] - P[14] + P[15]);
        
        a1 = 0.0625 * (-P[0] - P[1] + P[2] + P[3] - P[4] - P[5] + P[6] + P[7] - P[8] - P[9] + P[10] + P[11] - P[12] - P[13] + P[14] + P[15]);
        
        a2 = 0.0625 * (-P[0] - P[1] - P[2] - P[3] + P[4] + P[5] + P[6] + P[7] - P[8] - P[9] - P[10] - P[11] + P[12] + P[13] + P[14] + P[15]);
        
        a3 = 0.0625 * (-P[0] - P[1] - P[2] - P[3] - P[4] - P[5] - P[6] - P[7] + P[8] + P[9] + P[10] + P[11] + P[12] + P[13] + P[14] + P[15]);
        // Only supporting 4-cubes for the moment
    }
    
    bool MappingToPhysHypercubeLinear<4>::isValidPoint(const PointReference<4>& pointReference) const
    {
        logger.assert(pointReference.size()==4, "Reference point has the wrong dimension");
        if ((pointReference[0] < -1.) || (pointReference[0] > 1.) || (pointReference[1] < -1.) || (pointReference[1] > 1.) || (pointReference[2] < -1.) || (pointReference[2] > 1.) || (pointReference[3] < -1.) || (pointReference[3] > 1.))
            return false;
        else
            return true;
    }
}

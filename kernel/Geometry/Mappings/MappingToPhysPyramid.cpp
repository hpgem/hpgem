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


#ifndef MAPPINGPYRAMID_CPP_
#define MAPPINGPYRAMID_CPP_

#include "MappingToPhysPyramid.hpp"

#include "Geometry/PhysicalGeometry.hpp"
#include "Geometry/PointReference.hpp"
#include "Geometry/Jacobian.hpp"
#include "Geometry/PointPhysical.hpp"
#include <cmath>

namespace Geometry
{
    MappingToPhysPyramid::MappingToPhysPyramid(const PhysicalGeometryT*const physicalGeometry)
    {
        MappingReferenceToPhysical::setNodesPtr(&physicalGeometry->getNodes());
        reinit(physicalGeometry);
    }

    void MappingToPhysPyramid::
    transform(const PointReferenceT& pR, PointPhysicalT& pP) const
    {
        //if (isValidPoint(pR))
        //{
            const double t1 = pR[0] * pR[1];
            const double t2 = pR[0] * pR[1]  * pR[2] / (1 - pR[2]+1e-50);//prevents trouble at the tip of the pyramid

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
        //}
        //else  //while I agree this is a bad situation I dont want a mapping deciding for me that a crappy quadrature rule is illegal
        //{
        //    throw "ERROR: MappingToPhysPyramid::transform, mapping point outside geometry.";
        //}
    }
    void MappingToPhysPyramid::calcJacobian(const PointReferenceT& pR, JacobianT& jacobian) const
    {
        //if (isValidPoint(pR))
        //{
            std::vector<double> df_dxi0(5), df_dxi1(5), df_dxi2(5);

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
        //}
        //else
        //{
        //    throw "ERROR: MappingToPhysPyramid::calcJacobian, mapping point outside geometry.";
        //}
    }

    void MappingToPhysPyramid::reinit(const PhysicalGeometryT*const physicalGeometry)
    {
    	globalNodeIndices_.resize(5);
        for (int i = 0; i < 5; ++i)
        {
            globalNodeIndices_[i] = physicalGeometry->getNodeIndex(i);
        }
    }

    //* \bug Changed this function to use abs not std::abs is this gave the error error:
    /*  call to 'abs' is ambiguous, but only if we are using STL vector for the LinearAlgebrea
    /*  The problem is prob the use of math.h instead of cmath some where in the code but I cannot find it at the moment [Ant]
     *  (resolved) you forgot to #include <cmath>, but something #include <cstdlib> (where the integer type std::abs is defined) -FB
     */
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

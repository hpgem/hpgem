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

#include "MappingToPhysPyramid.h"

#include "Geometry/PhysicalGeometry.h"
#include "Geometry/PointReference.h"
#include "Geometry/Jacobian.h"
#include "Geometry/PointPhysical.h"
#include <cmath>

namespace Geometry
{
    MappingToPhysPyramid::MappingToPhysPyramid(const PhysicalGeometry* const physicalGeometry)
    {
        logger.assert(physicalGeometry!=nullptr, "Invalid physical geometry passed");
        MappingReferenceToPhysical::setNodesPtr(&physicalGeometry->getNodes());
        reinit(physicalGeometry);
    }
    
    PointPhysical MappingToPhysPyramid::transform(const PointReference& pR) const
    {
        logger.assert(pR.size()==3, "Reference point has the wrong dimension");
        PointPhysical pP(3);
        const double t1 = pR[0] * pR[1];
        const double t2 = pR[0] * pR[1] * pR[2] / (1 - pR[2] + 1e-50); //prevents trouble at the tip of the pyramid
                
        std::vector<double> f8;
        f8.resize(5);
        f8[0] = pR[2];
        f8[1] = 0.25 * (1. - pR[0] - pR[1] + t1 - pR[2] + t2);
        f8[2] = 0.25 * (1. + pR[0] - pR[1] - t1 - pR[2] - t2);
        f8[3] = 0.25 * (1. - pR[0] + pR[1] - t1 - pR[2] - t2);
        f8[4] = 0.25 * (1. + pR[0] + pR[1] + t1 - pR[2] + t2);
        
        PointPhysical p(3);
        
        pP[0] = pP[1] = pP[2] = 0.0;
        
        for (std::size_t i = 0; i < 5; ++i)
        {
            p = getNodeCoordinates(globalNodeIndices_[i]);
            pP += f8[i] * p;
        }
        return pP;
    }
    Jacobian MappingToPhysPyramid::calcJacobian(const PointReference& pR) const
    {
        logger.assert(pR.size()==3, "Reference point has the wrong dimension");
        Jacobian jacobian(3, 3);
        std::vector<double> df_dxi0(5), df_dxi1(5), df_dxi2(5);
        
        const double dt6dx0 = pR[1] * pR[2] / (1. - pR[2]);
        const double dt6dx1 = pR[0] * pR[2] / (1. - pR[2]);
        const double dt6dx2 = pR[0] * pR[1] / ((1. - pR[2]) * (1. - pR[2]));
        
        df_dxi0[0] = 0.0;
        df_dxi0[1] = 0.25 * (-1. + pR[1] + dt6dx0);
        df_dxi0[2] = 0.25 * (1. - pR[1] - dt6dx0);
        df_dxi0[3] = 0.25 * (-1. - pR[1] - dt6dx0);
        df_dxi0[4] = 0.25 * (1. + pR[1] + dt6dx0);
        
        df_dxi1[0] = 0.0;
        df_dxi1[1] = 0.25 * (-1. + pR[0] + dt6dx1);
        df_dxi1[2] = 0.25 * (-1. - pR[0] - dt6dx1);
        df_dxi1[3] = 0.25 * (1. - pR[0] - dt6dx1);
        df_dxi1[4] = 0.25 * (1. + pR[0] + dt6dx1);
        
        df_dxi2[0] = 1.0;
        df_dxi2[1] = 0.25 * (-1. + dt6dx2);
        df_dxi2[2] = 0.25 * (-1. - dt6dx2);
        df_dxi2[3] = 0.25 * (-1. - dt6dx2);
        df_dxi2[4] = 0.25 * (-1. + dt6dx2);
        
        PointPhysical d_dxi0(3);
        PointPhysical d_dxi1(3);
        PointPhysical d_dxi2(3);
        
        for (std::size_t i = 0; i < 3; ++i)
        {
            d_dxi0[i] = 0.;
            d_dxi1[i] = 0.;
            d_dxi2[i] = 0.;
        }
        
        PointPhysical p(3);
        
        for (std::size_t i = 0; i < 5; ++i)
        {
            p = getNodeCoordinates(globalNodeIndices_[i]);
            
            d_dxi0 += df_dxi0[i] * p;
            d_dxi1 += df_dxi1[i] * p;
            d_dxi2 += df_dxi2[i] * p;
        }
        
        for (std::size_t i = 0; i < 3; ++i)
        {
            jacobian(i, 0) = d_dxi0[i];
            jacobian(i, 1) = d_dxi1[i];
            jacobian(i, 2) = d_dxi2[i];
        }
        return jacobian;
    }
    
    void MappingToPhysPyramid::reinit(const PhysicalGeometry* const physicalGeometry)
    {
        logger.assert(physicalGeometry!=nullptr, "Invalid physical geometry passed");
        globalNodeIndices_.resize(5);
        for (std::size_t i = 0; i < 5; ++i)
        {
            globalNodeIndices_[i] = physicalGeometry->getNodeIndex(i);
        }
    }
    
    bool MappingToPhysPyramid::isValidPoint(const PointReference& pointReference) const
    {
        logger.assert(pointReference.size()==3, "Reference point has the wrong dimension");
        static const double eps = 1.e-14;
        const double z = pointReference[2];
        if ((std::abs(pointReference[0]) <= 1. - z + eps) && (std::abs(pointReference[1]) <= 1. - z + eps) && (z >= 0. - eps) && (z <= 1. + eps))
        {
            return true;
        }
        else
        {
            return false;
        }
    }

}
;
#endif

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

#ifndef MESHMOVERCONTRACTION_HPP_
#define MESHMOVERCONTRACTION_HPP_

#include "Base/MeshMoverBase.h"

/// MeshMoverContraction moves the grid such that a contraction in the domain arises
/// If you want only the first half of the contraction, set xMiddle to the end of the
/// domain.
class MeshMoverContraction : public Base::MeshMoverBase<DIM>
{
    
public:
    
    using PointPhysicalT = Geometry::PointPhysical<DIM>;
    
    MeshMoverContraction()
    {
    }
    
    virtual ~MeshMoverContraction()
    {
    }

    void movePoint(PointPhysicalT& point) const override final
    {
        logger.assert(2 == DIM, "Called mesh mover for contraction while DIM != 2");
        const double xBegin = 1;
        const double xMiddle = 5;
        const double xEnd = 100;
        const double contractionWidth = 0.8;        
        
        const double distFirst = xMiddle - xBegin;
        const double distSecond = xEnd - xMiddle;
        const double indentationWidth = (1. - contractionWidth) / 2;
        if (point[0] > xBegin && point[0] < xMiddle)
        {
            point[1] = (indentationWidth + point[1] * contractionWidth)  - (point[1] - (indentationWidth + point[1] * contractionWidth))* ((point[0] - xMiddle)/distFirst );
        }
        if (point[0] >= xMiddle && point[0] <= xEnd)
        {
            point[1] = (indentationWidth + point[1] * contractionWidth)  + (point[1] - (indentationWidth + point[1] * contractionWidth))* ((point[0] - xMiddle)/distSecond );
        }
    }
};

#endif /* MESHMOVER_HPP_ */

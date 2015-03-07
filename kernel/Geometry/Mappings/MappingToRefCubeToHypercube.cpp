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

#include "MappingToRefCubeToHypercube.h"
#include "Geometry/PointReference.h"
#include "Geometry/Jacobian.h"

namespace Geometry
{
    // ~~~ index 0 ~~~==============================================================================
    
    const MappingToRefCubeToHypercube0& MappingToRefCubeToHypercube0::Instance()
    {
        static const MappingToRefCubeToHypercube0 theInstance;
        return theInstance;
    }
    
    PointReference MappingToRefCubeToHypercube0::transform(const Geometry::PointReference& p1) const
    {
        PointReference p2(4);
        p2[0] = p1[0];
        p2[1] = p1[1];
        p2[2] = p1[2];
        p2[3] = -1.0;
        return p2;
    }
    
    Jacobian MappingToRefCubeToHypercube0::calcJacobian(const Geometry::PointReference& p1) const
    {
        Jacobian jacobian(4, 3);
        jacobian(0, 0) = 1.0;
        jacobian(1, 0) = 0.0;
        jacobian(2, 0) = 0.0;
        jacobian(3, 0) = 0.0;
        
        jacobian(0, 1) = 0.0;
        jacobian(1, 1) = 1.0;
        jacobian(2, 1) = 0.0;
        jacobian(3, 1) = 0.0;
        
        jacobian(0, 2) = 0.0;
        jacobian(1, 2) = 0.0;
        jacobian(2, 2) = 1.0;
        jacobian(3, 2) = 0.0;
        return jacobian;
    }
    
    MappingToRefCubeToHypercube0::MappingToRefCubeToHypercube0()
    {
    }
    MappingToRefCubeToHypercube0::~MappingToRefCubeToHypercube0()
    {
    }
    
    // ~~~ index 1 ~~~==============================================================================
    
    const MappingToRefCubeToHypercube1& MappingToRefCubeToHypercube1::Instance()
    {
        static const MappingToRefCubeToHypercube1 theInstance;
        return theInstance;
    }
    
    PointReference MappingToRefCubeToHypercube1::transform(const Geometry::PointReference& p1) const
    {
        PointReference p2(4);
        p2[0] = p1[0];
        p2[1] = p1[1];
        p2[2] = -1.0;
        p2[3] = p1[2];
        return p2;
    }
    
    Jacobian MappingToRefCubeToHypercube1::calcJacobian(const Geometry::PointReference& p1) const
    {
        Jacobian jacobian(4, 3);
        jacobian(0, 0) = 1.0;
        jacobian(1, 0) = 0.0;
        jacobian(2, 0) = 0.0;
        jacobian(3, 0) = 0.0;
        
        jacobian(0, 1) = 0.0;
        jacobian(1, 1) = 1.0;
        jacobian(2, 1) = 0.0;
        jacobian(3, 1) = 0.0;
        
        jacobian(0, 2) = 0.0;
        jacobian(1, 2) = 0.0;
        jacobian(2, 2) = 0.0;
        jacobian(3, 2) = 1.0;
        return jacobian;
    }
    
    MappingToRefCubeToHypercube1::MappingToRefCubeToHypercube1()
    {
    }
    MappingToRefCubeToHypercube1::~MappingToRefCubeToHypercube1()
    {
    }
    
    // ~~~ index 2 ~~~==============================================================================
    
    const MappingToRefCubeToHypercube2& MappingToRefCubeToHypercube2::Instance()
    {
        static const MappingToRefCubeToHypercube2 theInstance;
        return theInstance;
    }
    
    PointReference MappingToRefCubeToHypercube2::transform(const Geometry::PointReference& p1) const
    {
        PointReference p2(4);
        p2[0] = p1[0];
        p2[1] = -1.0;
        p2[2] = p1[1];
        p2[3] = p1[2];
        return p2;
    }
    
    Jacobian MappingToRefCubeToHypercube2::calcJacobian(const Geometry::PointReference& p1) const
    {
        Jacobian jacobian(4, 3);
        jacobian(0, 0) = 1.0;
        jacobian(1, 0) = 0.0;
        jacobian(2, 0) = 0.0;
        jacobian(3, 0) = 0.0;
        
        jacobian(0, 1) = 0.0;
        jacobian(1, 1) = 0.0;
        jacobian(2, 1) = 1.0;
        jacobian(3, 1) = 0.0;
        
        jacobian(0, 2) = 0.0;
        jacobian(1, 2) = 0.0;
        jacobian(2, 2) = 0.0;
        jacobian(3, 2) = 1.0;
        return jacobian;
    }
    
    MappingToRefCubeToHypercube2::MappingToRefCubeToHypercube2()
    {
    }
    MappingToRefCubeToHypercube2::~MappingToRefCubeToHypercube2()
    {
    }
    
    // ~~~ index 3 ~~~==============================================================================
    
    const MappingToRefCubeToHypercube3& MappingToRefCubeToHypercube3::Instance()
    {
        static const MappingToRefCubeToHypercube3 theInstance;
        return theInstance;
    }
    
    PointReference MappingToRefCubeToHypercube3::transform(const Geometry::PointReference& p1) const
    {
        PointReference p2(4);
        p2[0] = -1.0;
        p2[1] = p1[0];
        p2[2] = p1[1];
        p2[3] = p1[2];
        return p2;
    }
    
    Jacobian MappingToRefCubeToHypercube3::calcJacobian(const Geometry::PointReference& p1) const
    {
        Jacobian jacobian(4, 3);
        jacobian(0, 0) = 0.0;
        jacobian(1, 0) = 1.0;
        jacobian(2, 0) = 0.0;
        jacobian(3, 0) = 0.0;
        
        jacobian(0, 1) = 0.0;
        jacobian(1, 1) = 0.0;
        jacobian(2, 1) = 1.0;
        jacobian(3, 1) = 0.0;
        
        jacobian(0, 2) = 0.0;
        jacobian(1, 2) = 0.0;
        jacobian(2, 2) = 0.0;
        jacobian(3, 2) = 1.0;
        return jacobian;
    }
    
    MappingToRefCubeToHypercube3::MappingToRefCubeToHypercube3()
    {
    }
    MappingToRefCubeToHypercube3::~MappingToRefCubeToHypercube3()
    {
    }
    
    // ~~~ index 4 ~~~==============================================================================
    
    const MappingToRefCubeToHypercube4& MappingToRefCubeToHypercube4::Instance()
    {
        static const MappingToRefCubeToHypercube4 theInstance;
        return theInstance;
    }
    
    PointReference MappingToRefCubeToHypercube4::transform(const Geometry::PointReference& p1) const
    {
        PointReference p2(4);
        p2[0] = +1.0;
        p2[1] = p1[0];
        p2[2] = p1[1];
        p2[3] = p1[2];
        return p2;
    }
    
    Jacobian MappingToRefCubeToHypercube4::calcJacobian(const Geometry::PointReference& p1) const
    {
        Jacobian jacobian(4, 3);
        jacobian(0, 0) = 0.0;
        jacobian(1, 0) = 1.0;
        jacobian(2, 0) = 0.0;
        jacobian(3, 0) = 0.0;
        
        jacobian(0, 1) = 0.0;
        jacobian(1, 1) = 0.0;
        jacobian(2, 1) = 1.0;
        jacobian(3, 1) = 0.0;
        
        jacobian(0, 2) = 0.0;
        jacobian(1, 2) = 0.0;
        jacobian(2, 2) = 0.0;
        jacobian(3, 2) = 1.0;
        return jacobian;
    }
    
    MappingToRefCubeToHypercube4::MappingToRefCubeToHypercube4()
    {
    }
    MappingToRefCubeToHypercube4::~MappingToRefCubeToHypercube4()
    {
    }
    
    // ~~~ index 5 ~~~==============================================================================
    
    const MappingToRefCubeToHypercube5& MappingToRefCubeToHypercube5::Instance()
    {
        static const MappingToRefCubeToHypercube5 theInstance;
        return theInstance;
    }
    
    PointReference MappingToRefCubeToHypercube5::transform(const Geometry::PointReference& p1) const
    {
        PointReference p2(4);
        p2[0] = p1[0];
        p2[1] = +1.0;
        p2[2] = p1[1];
        p2[3] = p1[2];
        return p2;
    }
    
    Jacobian MappingToRefCubeToHypercube5::calcJacobian(const Geometry::PointReference& p1) const
    {
        Jacobian jacobian(4, 3);
        jacobian(0, 0) = 1.0;
        jacobian(1, 0) = 0.0;
        jacobian(2, 0) = 0.0;
        jacobian(3, 0) = 0.0;
        
        jacobian(0, 1) = 0.0;
        jacobian(1, 1) = 0.0;
        jacobian(2, 1) = 1.0;
        jacobian(3, 1) = 0.0;
        
        jacobian(0, 2) = 0.0;
        jacobian(1, 2) = 0.0;
        jacobian(2, 2) = 0.0;
        jacobian(3, 2) = 1.0;
        return jacobian;
    }
    
    MappingToRefCubeToHypercube5::MappingToRefCubeToHypercube5()
    {
    }
    MappingToRefCubeToHypercube5::~MappingToRefCubeToHypercube5()
    {
    }
    
    // ~~~ index 6 ~~~==============================================================================
    
    const MappingToRefCubeToHypercube6& MappingToRefCubeToHypercube6::Instance()
    {
        static const MappingToRefCubeToHypercube6 theInstance;
        return theInstance;
    }
    
    PointReference MappingToRefCubeToHypercube6::transform(const Geometry::PointReference& p1) const
    {
        PointReference p2(4);
        p2[0] = p1[0];
        p2[1] = p1[1];
        p2[2] = +1.0;
        p2[3] = p1[2];
        return p2;
    }
    
    Jacobian MappingToRefCubeToHypercube6::calcJacobian(const Geometry::PointReference& p1) const
    {
        Jacobian jacobian(4, 3);
        jacobian(0, 0) = 1.0;
        jacobian(1, 0) = 0.0;
        jacobian(2, 0) = 0.0;
        jacobian(3, 0) = 0.0;
        
        jacobian(0, 1) = 0.0;
        jacobian(1, 1) = 1.0;
        jacobian(2, 1) = 0.0;
        jacobian(3, 1) = 0.0;
        
        jacobian(0, 2) = 0.0;
        jacobian(1, 2) = 0.0;
        jacobian(2, 2) = 0.0;
        jacobian(3, 2) = 1.0;
        return jacobian;
    }
    
    MappingToRefCubeToHypercube6::MappingToRefCubeToHypercube6()
    {
    }
    MappingToRefCubeToHypercube6::~MappingToRefCubeToHypercube6()
    {
    }
    
    // ~~~ index 7 ~~~==============================================================================
    
    const MappingToRefCubeToHypercube7& MappingToRefCubeToHypercube7::Instance()
    {
        static const MappingToRefCubeToHypercube7 theInstance;
        return theInstance;
    }
    
    PointReference MappingToRefCubeToHypercube7::transform(const Geometry::PointReference& p1) const
    {
        PointReference p2(4);
        p2[0] = p1[0];
        p2[1] = p1[1];
        p2[2] = p1[2];
        p2[3] = +1.0;
        return p2;
        
    }
    
    Jacobian MappingToRefCubeToHypercube7::calcJacobian(const Geometry::PointReference& p1) const
    {
        Jacobian jacobian(4, 3);
        jacobian(0, 0) = 1.0;
        jacobian(1, 0) = 0.0;
        jacobian(2, 0) = 0.0;
        jacobian(3, 0) = 0.0;
        
        jacobian(0, 1) = 0.0;
        jacobian(1, 1) = 1.0;
        jacobian(2, 1) = 0.0;
        jacobian(3, 1) = 0.0;
        
        jacobian(0, 2) = 0.0;
        jacobian(1, 2) = 0.0;
        jacobian(2, 2) = 1.0;
        jacobian(3, 2) = 0.0;
        return jacobian;
    }
    
    MappingToRefCubeToHypercube7::MappingToRefCubeToHypercube7()
    {
    }
    MappingToRefCubeToHypercube7::~MappingToRefCubeToHypercube7()
    {
    }
}

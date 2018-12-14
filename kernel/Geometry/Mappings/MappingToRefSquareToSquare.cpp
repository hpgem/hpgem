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


#include "MappingToRefSquareToSquare.h"
#include "Geometry/Jacobian.h"
#include "Geometry/PointReference.h"

namespace Geometry
{
    // ~~~~~~~~~~~~~~~==============================================================================
    // ~~~ index 0 ~~~==============================================================================
    // ~~~~~~~~~~~~~~~==============================================================================
    
    const MappingToRefSquareToSquare0& MappingToRefSquareToSquare0::Instance()
    {
        static const MappingToRefSquareToSquare0 theInstance;
        return theInstance;
    }
    
    PointReference<2> MappingToRefSquareToSquare0::transform(const Geometry::PointReference<2>& p1) const
    {
        logger.assert_debug(p1.size() == 2, "Reference point has the wrong dimension");
        return {p1[0], p1[1]};
    }
    
    Jacobian<2, 2> MappingToRefSquareToSquare0::calcJacobian(const Geometry::PointReference<2>& p1) const
    {
        logger.assert_debug(p1.size() == 2, "Reference point has the wrong dimension");
        Jacobian<2, 2> jacobian;
        jacobian(0, 0) = 1.0;
        jacobian(0, 1) = 0.0;
        jacobian(1, 0) = 0.0;
        jacobian(1, 1) = 1.0;
        return jacobian;
    }
    
    MappingToRefSquareToSquare0::MappingToRefSquareToSquare0()
    {
    }
    
    // ~~~~~~~~~~~~~~~==============================================================================
    // ~~~ index 1 ~~~==============================================================================
    // ~~~~~~~~~~~~~~~==============================================================================
    
    const MappingToRefSquareToSquare1& MappingToRefSquareToSquare1::Instance()
    {
        static const MappingToRefSquareToSquare1 theInstance;
        return theInstance;
    }
    
    PointReference<2> MappingToRefSquareToSquare1::transform(const Geometry::PointReference<2>& p1) const
    {
        logger.assert_debug(p1.size() == 2, "Reference point has the wrong dimension");
        return {-p1[1], p1[0]};
    }
    
    Jacobian<2, 2> MappingToRefSquareToSquare1::calcJacobian(const Geometry::PointReference<2>& p1) const
    {
        logger.assert_debug(p1.size() == 2, "Reference point has the wrong dimension");
        Jacobian<2, 2> jacobian;
        jacobian(0, 0) = 0.0;
        jacobian(0, 1) = -1.0;
        jacobian(1, 0) = 1.0;
        jacobian(1, 1) = 0.0;
        return jacobian;
    }
    
    MappingToRefSquareToSquare1::MappingToRefSquareToSquare1()
    {
    }
    
    // ~~~~~~~~~~~~~~~==============================================================================
    // ~~~ index 2 ~~~==============================================================================
    // ~~~~~~~~~~~~~~~==============================================================================
    
    const MappingToRefSquareToSquare2& MappingToRefSquareToSquare2::Instance()
    {
        static const MappingToRefSquareToSquare2 theInstance;
        return theInstance;
    }
    
    PointReference<2> MappingToRefSquareToSquare2::transform(const Geometry::PointReference<2>& p1) const
    {
        logger.assert_debug(p1.size() == 2, "Reference point has the wrong dimension");
        return {-p1[0], -p1[1]};
    }
    
    Jacobian<2, 2> MappingToRefSquareToSquare2::calcJacobian(const Geometry::PointReference<2>& p1) const
    {
        logger.assert_debug(p1.size() == 2, "Reference point has the wrong dimension");
        Jacobian<2, 2> jacobian;
        jacobian(0, 0) = -1.0;
        jacobian(0, 1) = 0.0;
        jacobian(1, 0) = 0.0;
        jacobian(1, 1) = -1.0;
        return jacobian;
    }
    
    MappingToRefSquareToSquare2::MappingToRefSquareToSquare2()
    {
    }
    
    // ~~~~~~~~~~~~~~~==============================================================================
    // ~~~ index 3 ~~~==============================================================================
    // ~~~~~~~~~~~~~~~==============================================================================
    
    const MappingToRefSquareToSquare3& MappingToRefSquareToSquare3::Instance()
    {
        static const MappingToRefSquareToSquare3 theInstance;
        return theInstance;
    }
    
    PointReference<2> MappingToRefSquareToSquare3::transform(const Geometry::PointReference<2>& p1) const
    {
        logger.assert_debug(p1.size() == 2, "Reference point has the wrong dimension");
        return {p1[1], -p1[0]};
    }
    
    Jacobian<2, 2> MappingToRefSquareToSquare3::calcJacobian(const Geometry::PointReference<2>& p1) const
    {
        logger.assert_debug(p1.size() == 2, "Reference point has the wrong dimension");
        Jacobian<2, 2> jacobian;
        jacobian(0, 0) = 0.0;
        jacobian(0, 1) = 1.0;
        jacobian(1, 0) = -1.0;
        jacobian(1, 1) = 0.0;
        return jacobian;
    }
    
    MappingToRefSquareToSquare3::MappingToRefSquareToSquare3()
    {
    }
    
    // ~~~~~~~~~~~~~~~==============================================================================
    // ~~~ index 4 ~~~==============================================================================
    // ~~~~~~~~~~~~~~~==============================================================================
    
    const MappingToRefSquareToSquare4& MappingToRefSquareToSquare4::Instance()
    {
        static const MappingToRefSquareToSquare4 theInstance;
        return theInstance;
    }
    
    PointReference<2> MappingToRefSquareToSquare4::transform(const Geometry::PointReference<2>& p1) const
    {
        logger.assert_debug(p1.size() == 2, "Reference point has the wrong dimension");
        return {p1[0], -p1[1]};
    }
    
    Jacobian<2, 2> MappingToRefSquareToSquare4::calcJacobian(const Geometry::PointReference<2>& p1) const
    {
        logger.assert_debug(p1.size() == 2, "Reference point has the wrong dimension");
        Jacobian<2, 2> jacobian;
        jacobian(0, 0) = 1.0;
        jacobian(0, 1) = 0.0;
        jacobian(1, 0) = 0.0;
        jacobian(1, 1) = -1.0;
        return jacobian;
    }
    
    MappingToRefSquareToSquare4::MappingToRefSquareToSquare4()
    {
    }
    
    // ~~~~~~~~~~~~~~~==============================================================================
    // ~~~ index 5 ~~~==============================================================================
    // ~~~~~~~~~~~~~~~==============================================================================
    
    const MappingToRefSquareToSquare5& MappingToRefSquareToSquare5::Instance()
    {
        static const MappingToRefSquareToSquare5 theInstance;
        return theInstance;
    }
    
    PointReference<2> MappingToRefSquareToSquare5::transform(const Geometry::PointReference<2>& p1) const
    {
        logger.assert_debug(p1.size() == 2, "Reference point has the wrong dimension");
        return {-p1[0], p1[1]};
    }
    
    Jacobian<2, 2> MappingToRefSquareToSquare5::calcJacobian(const Geometry::PointReference<2>& p1) const
    {
        logger.assert_debug(p1.size() == 2, "Reference point has the wrong dimension");
        Jacobian<2, 2> jacobian;
        jacobian(0, 0) = -1.0;
        jacobian(0, 1) = 0.0;
        jacobian(1, 0) = 0.0;
        jacobian(1, 1) = 1.0;
        return jacobian;
    }
    
    MappingToRefSquareToSquare5::MappingToRefSquareToSquare5()
    {
    }
    
    // ~~~~~~~~~~~~~~~==============================================================================
    // ~~~ index 6 ~~~==============================================================================
    // ~~~~~~~~~~~~~~~==============================================================================
    
    const MappingToRefSquareToSquare6& MappingToRefSquareToSquare6::Instance()
    {
        static const MappingToRefSquareToSquare6 theInstance;
        return theInstance;
    }
    
    PointReference<2> MappingToRefSquareToSquare6::transform(const Geometry::PointReference<2>& p1) const
    {
        logger.assert_debug(p1.size() == 2, "Reference point has the wrong dimension");
        return {-p1[1], -p1[0]};
    }
    
    Jacobian<2, 2> MappingToRefSquareToSquare6::calcJacobian(const Geometry::PointReference<2>& p1) const
    {
        logger.assert_debug(p1.size() == 2, "Reference point has the wrong dimension");
        Jacobian<2, 2> jacobian;
        jacobian(0, 0) = 0.0;
        jacobian(0, 1) = -1.0;
        jacobian(1, 0) = -1.0;
        jacobian(1, 1) = 0.0;
        return jacobian;
    }
    
    MappingToRefSquareToSquare6::MappingToRefSquareToSquare6()
    {
    }
    
    // ~~~~~~~~~~~~~~~==============================================================================
    // ~~~ index 7 ~~~==============================================================================
    // ~~~~~~~~~~~~~~~==============================================================================
    
    const MappingToRefSquareToSquare7& MappingToRefSquareToSquare7::Instance()
    {
        static const MappingToRefSquareToSquare7 theInstance;
        return theInstance;
    }
    
    PointReference<2> MappingToRefSquareToSquare7::transform(const Geometry::PointReference<2>& p1) const
    {
        logger.assert_debug(p1.size() == 2, "Reference point has the wrong dimension");
        return {p1[1], p1[0]};
    }
    
    Jacobian<2, 2> MappingToRefSquareToSquare7::calcJacobian(const Geometry::PointReference<2>& p1) const
    {
        logger.assert_debug(p1.size() == 2, "Reference point has the wrong dimension");
        Jacobian<2, 2> jacobian;
        jacobian(0, 0) = 0.0;
        jacobian(0, 1) = 1.0;
        jacobian(1, 0) = 1.0;
        jacobian(1, 1) = 0.0;
        return jacobian;
    }
    
    MappingToRefSquareToSquare7::MappingToRefSquareToSquare7()
    {
    }
}

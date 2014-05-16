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

#include "MappingToRefSquareToCube.hpp"
#include "Geometry/PointReference.hpp"
#include "Geometry/Jacobian.hpp"

namespace Geometry
{
    // ~~~ index 0 ~~~==============================================================================

    const MappingToRefSquareToCube0& MappingToRefSquareToCube0::Instance()
    {
        static const MappingToRefSquareToCube0 theInstance;
        return theInstance;
    }

    void MappingToRefSquareToCube0::transform(const Geometry::PointReference& p1,
                                               Geometry::PointReference& p2) const
    {
        p2[0] = p1[0];
        p2[1] = p1[1];
        p2[2] = -1.0;
    }

    void MappingToRefSquareToCube0::calcJacobian(const Geometry::PointReference& p1,
                                                  Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) = 1.0;
        jacobian(1,0) = 0.0;
        jacobian(2,0) = 0.0;

        jacobian(0,1) = 0.0;
        jacobian(1,1) = 1.0;
        jacobian(2,1) = 0.0;
    }

    MappingToRefSquareToCube0::MappingToRefSquareToCube0() { }
    MappingToRefSquareToCube0::~MappingToRefSquareToCube0() { }

    // ~~~ index 1 ~~~==============================================================================

    const MappingToRefSquareToCube1& MappingToRefSquareToCube1::Instance()
    {
        static const MappingToRefSquareToCube1 theInstance;
        return theInstance;
    }

    void MappingToRefSquareToCube1::transform(const Geometry::PointReference& p1,
                                               Geometry::PointReference& p2) const
    {
        p2[0] = p1[0];
        p2[1] = -1.0;
        p2[2] = p1[1];
    }

    void MappingToRefSquareToCube1::calcJacobian(const Geometry::PointReference& p1,
                                                  Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) = 1.0;
        jacobian(1,0) = 0.0;
        jacobian(2,0) = 0.0;

        jacobian(0,1) = 0.0;
        jacobian(1,1) = 0.0;
        jacobian(2,1) = 1.0;
    }

    MappingToRefSquareToCube1::MappingToRefSquareToCube1() { }
    MappingToRefSquareToCube1::~MappingToRefSquareToCube1() { }

    // ~~~ index 2 ~~~==============================================================================

    const MappingToRefSquareToCube2& MappingToRefSquareToCube2::Instance()
    {
        static const MappingToRefSquareToCube2 theInstance;
        return theInstance;
    }

    void MappingToRefSquareToCube2::transform(const Geometry::PointReference& p1,
                                               Geometry::PointReference& p2) const
    {
        p2[0] = -1.0;
        p2[1] = p1[0];
        p2[2] = p1[1];
    }

    void MappingToRefSquareToCube2::calcJacobian(const Geometry::PointReference& p1,
                                                  Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) = 0.0;
        jacobian(1,0) = 1.0;
        jacobian(2,0) = 0.0;

        jacobian(0,1) = 0.0;
        jacobian(1,1) = 0.0;
        jacobian(2,1) = 1.0;
    }

    MappingToRefSquareToCube2::MappingToRefSquareToCube2() { }
    MappingToRefSquareToCube2::~MappingToRefSquareToCube2() { }

    // ~~~ index 3 ~~~==============================================================================

    const MappingToRefSquareToCube3& MappingToRefSquareToCube3::Instance()
    {
        static const MappingToRefSquareToCube3 theInstance;
        return theInstance;
    }

    void MappingToRefSquareToCube3::transform(const Geometry::PointReference& p1,
                                               Geometry::PointReference& p2) const
    {
        p2[0] = 1.0;
        p2[1] = p1[0];
        p2[2] = p1[1];
    }

    void MappingToRefSquareToCube3::calcJacobian(const Geometry::PointReference& p1,
                                                  Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) = 0.0;
        jacobian(1,0) = 1.0;
        jacobian(2,0) = 0.0;

        jacobian(0,1) = 0.0;
        jacobian(1,1) = 0.0;
        jacobian(2,1) = 1.0;
    }

    MappingToRefSquareToCube3::MappingToRefSquareToCube3() { }
    MappingToRefSquareToCube3::~MappingToRefSquareToCube3() { }

    // ~~~ index 4 ~~~==============================================================================

    const MappingToRefSquareToCube4& MappingToRefSquareToCube4::Instance()
    {
        static const MappingToRefSquareToCube4 theInstance;
        return theInstance;
    }

    void MappingToRefSquareToCube4::transform(const Geometry::PointReference& p1,
                                               Geometry::PointReference& p2) const
    {
        p2[0] = p1[0];
        p2[1] = 1.0;
        p2[2] = p1[1];
    }

    void MappingToRefSquareToCube4::calcJacobian(const Geometry::PointReference& p1,
                                                  Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) = 1.0;
        jacobian(1,0) = 0.0;
        jacobian(2,0) = 0.0;

        jacobian(0,1) = 0.0;
        jacobian(1,1) = 0.0;
        jacobian(2,1) = 1.0;
    }

    MappingToRefSquareToCube4::MappingToRefSquareToCube4() { }
    MappingToRefSquareToCube4::~MappingToRefSquareToCube4() { }

    // ~~~ index 5 ~~~==============================================================================

    const MappingToRefSquareToCube5& MappingToRefSquareToCube5::Instance()
    {
        static const MappingToRefSquareToCube5 theInstance;
        return theInstance;
    }

    void MappingToRefSquareToCube5::transform(const Geometry::PointReference& p1,
                                               Geometry::PointReference& p2) const
    {
        p2[0] = p1[0];
        p2[1] = p1[1];
        p2[2] = 1.0;
    }

    void MappingToRefSquareToCube5::calcJacobian(const Geometry::PointReference& p1,
                                                  Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) = 1.0;
        jacobian(1,0) = 0.0;
        jacobian(2,0) = 0.0;

        jacobian(0,1) = 0.0;
        jacobian(1,1) = 1.0;
        jacobian(2,1) = 0.0;
    }

    MappingToRefSquareToCube5::MappingToRefSquareToCube5() { }
    MappingToRefSquareToCube5::~MappingToRefSquareToCube5() { }


}

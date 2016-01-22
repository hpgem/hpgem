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

#ifndef MAPPINGSQUARETOSQUARE_H_
#define MAPPINGSQUARETOSQUARE_H_

#include "MappingReferenceToReference.h"

namespace Geometry
{
    /*
     * The ordering of the vertex and faces in a square:
     *
     * (-1,+1) 2---3---3 (+1,+1)
     *         |       |
     *         1       2
     *         |       |
     * (-1,-1) 0---0---1 (1,-1)
     *
     * This implements the linear mappings of a Square [-1,1]^2 onto itself.
     * There are 8 possible mappings:
     *
     *      index 1: (x,y)->(x,y)    (Identity)
     *      index 2: (x,y)->(-y,x)   (Right rotation)
     *      index 3: (x,y)->(-x,-y)  (2x Right rotation)
     *      index 4: (x,y)->(y,-x)   (3x Right rotation (or Left rotation))
     *      index 5: (x,y)->(x,-y)   (Mirror in the x-axis)
     *      index 6: (x,y)->(-x,y)   (Mirror in the y-axis)
     *      index 7: (x,y)->(-y,-x)  (Mirror in the line y=-x)
     *      index 8: (x,y)->(y,x)    (Mirror in the line y=x)
     *
     * These correspond to the node mutations : (1)(2)(3)(4), (1243), (14)(23),
     * (1342), (13)(24), (12)(34), (14)(2)(3) and (1)(4)(23). The first 4
     * are rotations, while the last 4 are mirrorings.
     *
     * The node numbering can be found in the ReferenceSquare class definition.
     */

    // ~~~ index 0 ~~~==============================================================================
    class MappingToRefSquareToSquare0 : public MappingReferenceToReference<0>
    {
    public:
        static const MappingToRefSquareToSquare0& Instance();
        PointReference<2> transform(const Geometry::PointReference<2>& p1) const override final;
        Jacobian<2, 2> calcJacobian(const Geometry::PointReference<2>&) const override final;
        std::size_t getTargetDimension() const override final
        {
            return 2;
        }
        MappingToRefSquareToSquare0(const MappingToRefSquareToSquare0&) = delete;
        MappingToRefSquareToSquare0& operator=(const MappingToRefSquareToSquare0&) = delete;
    private:
        MappingToRefSquareToSquare0();
        std::map<const PointReference<2>*, const PointReference<2>*> transformedCoordinates;
    };
    // ~~~ index 1 ~~~==============================================================================
    class MappingToRefSquareToSquare1 : public MappingReferenceToReference<0>
    {
    public:
        static const MappingToRefSquareToSquare1& Instance();
        PointReference<2> transform(const Geometry::PointReference<2>& p1) const override final;
        Jacobian<2, 2> calcJacobian(const Geometry::PointReference<2>&) const override final;
        std::size_t getTargetDimension() const override final
        {
            return 2;
        }
        MappingToRefSquareToSquare1(const MappingToRefSquareToSquare1&) = delete;
        MappingToRefSquareToSquare1& operator=(const MappingToRefSquareToSquare1&) = delete;
    private:
        MappingToRefSquareToSquare1();
        std::map<const PointReference<2>*, const PointReference<2>*> transformedCoordinates;
    };
    // ~~~ index 2 ~~~==============================================================================
    class MappingToRefSquareToSquare2 : public MappingReferenceToReference<0>
    {
    public:
        static const MappingToRefSquareToSquare2& Instance();
        PointReference<2> transform(const Geometry::PointReference<2>& p1) const override final;
        Jacobian<2, 2> calcJacobian(const Geometry::PointReference<2>&) const override final;
        std::size_t getTargetDimension() const override final
        {
            return 2;
        }
        MappingToRefSquareToSquare2(const MappingToRefSquareToSquare2&) = delete;
        MappingToRefSquareToSquare2& operator=(const MappingToRefSquareToSquare2&) = delete;
    private:
        MappingToRefSquareToSquare2();
        std::map<const PointReference<2>*, const PointReference<2>*> transformedCoordinates;
    };
    // ~~~ index 3 ~~~==============================================================================
    class MappingToRefSquareToSquare3 : public MappingReferenceToReference<0>
    {
    public:
        static const MappingToRefSquareToSquare3& Instance();
        PointReference<2> transform(const Geometry::PointReference<2>& p1) const override final;
        Jacobian<2, 2> calcJacobian(const Geometry::PointReference<2>&) const override final;
        std::size_t getTargetDimension() const override final
        {
            return 2;
        }
        MappingToRefSquareToSquare3(const MappingToRefSquareToSquare3&) = delete;
        MappingToRefSquareToSquare3& operator=(const MappingToRefSquareToSquare3&) = delete;
    private:
        MappingToRefSquareToSquare3();
        std::map<const PointReference<2>*, const PointReference<2>*> transformedCoordinates;
    };
    // ~~~ index 4 ~~~==============================================================================
    class MappingToRefSquareToSquare4 : public MappingReferenceToReference<0>
    {
    public:
        static const MappingToRefSquareToSquare4& Instance();
        PointReference<2> transform(const Geometry::PointReference<2>& p1) const override final;
        Jacobian<2, 2> calcJacobian(const Geometry::PointReference<2>&) const override final;
        std::size_t getTargetDimension() const override final
        {
            return 2;
        }
        MappingToRefSquareToSquare4(const MappingToRefSquareToSquare4&) = delete;
        MappingToRefSquareToSquare4& operator=(const MappingToRefSquareToSquare4&) = delete;
    private:
        MappingToRefSquareToSquare4();
        std::map<const PointReference<2>*, const PointReference<2>*> transformedCoordinates;
    };
    // ~~~ index 5 ~~~==============================================================================
    class MappingToRefSquareToSquare5 : public MappingReferenceToReference<0>
    {
    public:
        static const MappingToRefSquareToSquare5& Instance();
        PointReference<2> transform(const Geometry::PointReference<2>& p1) const override final;
        Jacobian<2, 2> calcJacobian(const Geometry::PointReference<2>&) const override final;
        std::size_t getTargetDimension() const override final
        {
            return 2;
        }
        MappingToRefSquareToSquare5(const MappingToRefSquareToSquare5&) = delete;
        MappingToRefSquareToSquare5& operator=(const MappingToRefSquareToSquare5&) = delete;
    private:
        MappingToRefSquareToSquare5();
        std::map<const PointReference<2>*, const PointReference<2>*> transformedCoordinates;
    };
    // ~~~ index 6 ~~~==============================================================================
    class MappingToRefSquareToSquare6 : public MappingReferenceToReference<0>
    {
    public:
        static const MappingToRefSquareToSquare6& Instance();
        PointReference<2> transform(const Geometry::PointReference<2>& p1) const override final;
        Jacobian<2, 2> calcJacobian(const Geometry::PointReference<2>&) const override final;
        std::size_t getTargetDimension() const override final
        {
            return 2;
        }
        MappingToRefSquareToSquare6(const MappingToRefSquareToSquare6&) = delete;
        MappingToRefSquareToSquare6& operator=(const MappingToRefSquareToSquare6&) = delete;
    private:
        MappingToRefSquareToSquare6();
        std::map<const PointReference<2>*, const PointReference<2>*> transformedCoordinates;
    };
    // ~~~ index 7 ~~~==============================================================================
    class MappingToRefSquareToSquare7 : public MappingReferenceToReference<0>
    {
    public:
        static const MappingToRefSquareToSquare7& Instance();
        PointReference<2> transform(const Geometry::PointReference<2>& p1) const override final;
        Jacobian<2, 2> calcJacobian(const Geometry::PointReference<2>&) const override final;
        std::size_t getTargetDimension() const override final
        {
            return 2;
        }
        MappingToRefSquareToSquare7(const MappingToRefSquareToSquare7&) = delete;
        MappingToRefSquareToSquare7& operator=(const MappingToRefSquareToSquare7&) = delete;
    private:
        MappingToRefSquareToSquare7();
        std::map<const PointReference<2>*, const PointReference<2>*> transformedCoordinates;
    };
}
#endif

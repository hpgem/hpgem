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

#ifndef MAPPINGTOREFCUBETOCUBE_H_
#define MAPPINGTOREFCUBETOCUBE_H_

#include "MappingReferenceToReference.h"

namespace Geometry
{
    /*!
     *
     * This implements the linear mappings of a Cube [-1,1]^3 onto itself.
     * The ordering of the vertex and faces in a cube:
     *
     *     6o---------o7
     *     /|        /|
     *    / |       / |
     *  4o---------o5 |
     *   | 2o------|--o3
     *   | /       | /
     *   |/        |/
     *  0o---------o1
     *
     * This is only used as a face-to-face map in space-time discretisation.
     *
     * There are 8 implemented mappings:
     *
     *      index 1: (t,x,y)->(t,x,y)    (Identity)
     *      index 2: (t,x,y)->(t,-y,x)   (Right rotation)
     *      index 3: (t,x,y)->(t,-x,-y)  (2x Right rotation)
     *      index 4: (t,x,y)->(t,y,-x)   (3x Right rotation (or Left rotation))
     *      index 5: (t,x,y)->(t,x,-y)   (Mirror in the x-axis)
     *      index 6: (t,x,y)->(t,-x,y)   (Mirror in the y-axis)
     *      index 7: (t,x,y)->(t,-y,-x)  (Mirror in the line y=-x)
     *      index 8: (t,x,y)->(t,y,x)    (Mirror in the line y=x)
     *
     * These correspond to the node mutations : (1)(2)(3)(4)(5)(6)(7)(0), (0264)(1375),
     * (06)(17)(24)(35), (0462)(1573), (04)(15)(26)(37), (02)(13)(46)(57),
     * (06)(2)(3)(17)(4)(5)   and (0)(6)(24)(1)(7)(35). The first 4 are rotations,
     * while the last 4 are mirrorings.
     */

    // ~~~ index 0 ~~~==============================================================================
    class MappingToRefCubeToCube0 : public MappingReferenceToReference
    {
    public:
        static const MappingToRefCubeToCube0& Instance();
        PointReference transform(const Geometry::PointReference& p1) const override final;
        Jacobian calcJacobian(const Geometry::PointReference&) const override final;
        std::size_t getTargetDimension() const override final
        {
            return 3;
        }        
        MappingToRefCubeToCube0(const MappingToRefCubeToCube0&) = delete;
        MappingToRefCubeToCube0& operator=(const MappingToRefCubeToCube0&) = delete;
    private:
        MappingToRefCubeToCube0();
    };
    // ~~~ index 1 ~~~==============================================================================
    class MappingToRefCubeToCube1 : public MappingReferenceToReference
    {
    public:
        static const MappingToRefCubeToCube1& Instance();
        PointReference transform(const Geometry::PointReference& p1) const override final;
        Jacobian calcJacobian(const Geometry::PointReference&) const override final;
        std::size_t getTargetDimension() const override final
        {
            return 3;
        }
        MappingToRefCubeToCube1(const MappingToRefCubeToCube1&) = delete;
        MappingToRefCubeToCube1& operator=(const MappingToRefCubeToCube1&) = delete;
    private:
        MappingToRefCubeToCube1();
    };
    // ~~~ index 2 ~~~==============================================================================
    class MappingToRefCubeToCube2 : public MappingReferenceToReference
    {
    public:
        static const MappingToRefCubeToCube2& Instance();
        PointReference transform(const Geometry::PointReference& p1) const override final;
        Jacobian calcJacobian(const Geometry::PointReference&) const override final;
        std::size_t getTargetDimension() const override final
        {
            return 3;
        }
        MappingToRefCubeToCube2(const MappingToRefCubeToCube2&) = delete;
        MappingToRefCubeToCube2& operator=(const MappingToRefCubeToCube2&) = delete;
    private:
        MappingToRefCubeToCube2();
    };
    // ~~~ index 3 ~~~==============================================================================
    class MappingToRefCubeToCube3 : public MappingReferenceToReference
    {
    public:
        static const MappingToRefCubeToCube3& Instance();
        PointReference transform(const Geometry::PointReference& p1) const override final;
        Jacobian calcJacobian(const Geometry::PointReference&) const override final;
        std::size_t getTargetDimension() const override final
        {
            return 3;
        }
        MappingToRefCubeToCube3(const MappingToRefCubeToCube3&) = delete;
        MappingToRefCubeToCube3& operator=(const MappingToRefCubeToCube3&) = delete;
    private:
        MappingToRefCubeToCube3();
    };
    // ~~~ index 4 ~~~==============================================================================
    class MappingToRefCubeToCube4 : public MappingReferenceToReference
    {
    public:
        static const MappingToRefCubeToCube4& Instance();
        PointReference transform(const Geometry::PointReference& p1) const override final;
        Jacobian calcJacobian(const Geometry::PointReference&) const override final;
        std::size_t getTargetDimension() const override final
        {
            return 3;
        }
        
        MappingToRefCubeToCube4(const MappingToRefCubeToCube4&) = delete;
        MappingToRefCubeToCube4& operator=(const MappingToRefCubeToCube4&) = delete;
    private:
        MappingToRefCubeToCube4();
    };
    // ~~~ index 5 ~~~==============================================================================
    class MappingToRefCubeToCube5 : public MappingReferenceToReference
    {
    public:
        static const MappingToRefCubeToCube5& Instance();
        PointReference transform(const Geometry::PointReference& p1) const override final;
        Jacobian calcJacobian(const Geometry::PointReference&) const override final;
        std::size_t getTargetDimension() const override final
        {
            return 3;
        }
        MappingToRefCubeToCube5(const MappingToRefCubeToCube5&) = delete;
        MappingToRefCubeToCube5& operator=(const MappingToRefCubeToCube5&) = delete;
    private:
        MappingToRefCubeToCube5();
    };
    // ~~~ index 6 ~~~==============================================================================
    class MappingToRefCubeToCube6 : public MappingReferenceToReference
    {
    public:
        static const MappingToRefCubeToCube6& Instance();
        PointReference transform(const Geometry::PointReference& p1) const override final;
        Jacobian calcJacobian(const Geometry::PointReference&) const override final;
        std::size_t getTargetDimension() const override final
        {
            return 3;
        }
        MappingToRefCubeToCube6(const MappingToRefCubeToCube6&) = delete;
        MappingToRefCubeToCube6& operator=(const MappingToRefCubeToCube6&) = delete;
    private:
        MappingToRefCubeToCube6();
    };
    // ~~~ index 7 ~~~==============================================================================
    class MappingToRefCubeToCube7 : public MappingReferenceToReference
    {
    public:
        static const MappingToRefCubeToCube7& Instance();
        PointReference transform(const Geometry::PointReference& p1) const override final;
        Jacobian calcJacobian(const Geometry::PointReference&) const override final;
        std::size_t getTargetDimension() const override final
        {
            return 3;
        }
        MappingToRefCubeToCube7(const MappingToRefCubeToCube7&) = delete;
        MappingToRefCubeToCube7& operator=(const MappingToRefCubeToCube7&) = delete;
    private:
        MappingToRefCubeToCube7();
    };
}
#endif

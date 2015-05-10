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
     * This implements the linear mappings of a Cube [-1,1]^2 onto itself.
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
     * ~OC~
     * \todo I don't understand what is 't', and why do we use the same mappings as in a square
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
        virtual PointReference transform(const Geometry::PointReference& p1) const;
        virtual Jacobian calcJacobian(const Geometry::PointReference&) const;
        virtual std::size_t getTargetDimension() const
        {
            return 3;
        }
    private:
        MappingToRefCubeToCube0();
        MappingToRefCubeToCube0(const MappingToRefCubeToCube0&);
        MappingToRefCubeToCube0& operator=(const MappingToRefCubeToCube0&);
        virtual ~MappingToRefCubeToCube0();
    };
    // ~~~ index 1 ~~~==============================================================================
    class MappingToRefCubeToCube1 : public MappingReferenceToReference
    {
    public:
        static const MappingToRefCubeToCube1& Instance();
        virtual PointReference transform(const Geometry::PointReference& p1) const;
        virtual Jacobian calcJacobian(const Geometry::PointReference&) const;
        virtual std::size_t getTargetDimension() const
        {
            return 3;
        }
    private:
        MappingToRefCubeToCube1();
        MappingToRefCubeToCube1(const MappingToRefCubeToCube1&);
        MappingToRefCubeToCube1& operator=(const MappingToRefCubeToCube1&);
        virtual ~MappingToRefCubeToCube1();
    };
    // ~~~ index 2 ~~~==============================================================================
    class MappingToRefCubeToCube2 : public MappingReferenceToReference
    {
    public:
        static const MappingToRefCubeToCube2& Instance();
        virtual PointReference transform(const Geometry::PointReference& p1) const;
        virtual Jacobian calcJacobian(const Geometry::PointReference&) const;
        virtual std::size_t getTargetDimension() const
        {
            return 3;
        }
    private:
        MappingToRefCubeToCube2();
        MappingToRefCubeToCube2(const MappingToRefCubeToCube2&);
        MappingToRefCubeToCube2& operator=(const MappingToRefCubeToCube2&);
        virtual ~MappingToRefCubeToCube2();
    };
    // ~~~ index 3 ~~~==============================================================================
    class MappingToRefCubeToCube3 : public MappingReferenceToReference
    {
    public:
        static const MappingToRefCubeToCube3& Instance();
        virtual PointReference transform(const Geometry::PointReference& p1) const;
        virtual Jacobian calcJacobian(const Geometry::PointReference&) const;
        virtual std::size_t getTargetDimension() const
        {
            return 3;
        }
    private:
        MappingToRefCubeToCube3();
        MappingToRefCubeToCube3(const MappingToRefCubeToCube3&);
        MappingToRefCubeToCube3& operator=(const MappingToRefCubeToCube3&);
        virtual ~MappingToRefCubeToCube3();
    };
    // ~~~ index 4 ~~~==============================================================================
    class MappingToRefCubeToCube4 : public MappingReferenceToReference
    {
    public:
        static const MappingToRefCubeToCube4& Instance();
        virtual PointReference transform(const Geometry::PointReference& p1) const;
        virtual Jacobian calcJacobian(const Geometry::PointReference&) const;
        virtual std::size_t getTargetDimension() const
        {
            return 3;
        }
    private:
        MappingToRefCubeToCube4();
        MappingToRefCubeToCube4(const MappingToRefCubeToCube4&);
        MappingToRefCubeToCube4& operator=(const MappingToRefCubeToCube4&);
        virtual ~MappingToRefCubeToCube4();
    };
    // ~~~ index 5 ~~~==============================================================================
    class MappingToRefCubeToCube5 : public MappingReferenceToReference
    {
    public:
        static const MappingToRefCubeToCube5& Instance();
        virtual PointReference transform(const Geometry::PointReference& p1) const;
        virtual Jacobian calcJacobian(const Geometry::PointReference&) const;
        virtual std::size_t getTargetDimension() const
        {
            return 3;
        }
    private:
        MappingToRefCubeToCube5();
        MappingToRefCubeToCube5(const MappingToRefCubeToCube5&);
        MappingToRefCubeToCube5& operator=(const MappingToRefCubeToCube5&);
        virtual ~MappingToRefCubeToCube5();
    };
    // ~~~ index 6 ~~~==============================================================================
    class MappingToRefCubeToCube6 : public MappingReferenceToReference
    {
    public:
        static const MappingToRefCubeToCube6& Instance();
        virtual PointReference transform(const Geometry::PointReference& p1) const;
        virtual Jacobian calcJacobian(const Geometry::PointReference&) const;
        virtual std::size_t getTargetDimension() const
        {
            return 3;
        }
    private:
        MappingToRefCubeToCube6();
        MappingToRefCubeToCube6(const MappingToRefCubeToCube6&);
        MappingToRefCubeToCube6& operator=(const MappingToRefCubeToCube6&);
        virtual ~MappingToRefCubeToCube6();
    };
    // ~~~ index 7 ~~~==============================================================================
    class MappingToRefCubeToCube7 : public MappingReferenceToReference
    {
    public:
        static const MappingToRefCubeToCube7& Instance();
        virtual PointReference transform(const Geometry::PointReference& p1) const;
        virtual Jacobian calcJacobian(const Geometry::PointReference&) const;
        virtual std::size_t getTargetDimension() const
        {
            return 3;
        }
    private:
        MappingToRefCubeToCube7();
        MappingToRefCubeToCube7(const MappingToRefCubeToCube7&);
        MappingToRefCubeToCube7& operator=(const MappingToRefCubeToCube7&);
        virtual ~MappingToRefCubeToCube7();
    };
}
;
#endif

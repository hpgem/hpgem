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

#include "MappingReferenceToReference.hpp"
#include "../Jacobian.hpp"

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
    class MappingToRefSquareToSquare0: public MappingReferenceToReference
    {
        public:
            static const MappingToRefSquareToSquare0& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 2;}
        private:
            MappingToRefSquareToSquare0();
            MappingToRefSquareToSquare0(const MappingToRefSquareToSquare0&);
            MappingToRefSquareToSquare0& operator=(const MappingToRefSquareToSquare0&);
            virtual ~MappingToRefSquareToSquare0();
    };
    // ~~~ index 1 ~~~==============================================================================
    class MappingToRefSquareToSquare1: public MappingReferenceToReference
    {
        public:
            static const MappingToRefSquareToSquare1& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 2;}
        private:
            MappingToRefSquareToSquare1();
            MappingToRefSquareToSquare1(const MappingToRefSquareToSquare1&);
            MappingToRefSquareToSquare1& operator=(const MappingToRefSquareToSquare1&);
            virtual ~MappingToRefSquareToSquare1();
    };
    // ~~~ index 2 ~~~==============================================================================
    class MappingToRefSquareToSquare2: public MappingReferenceToReference
    {
        public:
            static const MappingToRefSquareToSquare2& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 2;}
        private:
            MappingToRefSquareToSquare2();
            MappingToRefSquareToSquare2(const MappingToRefSquareToSquare2&);
            MappingToRefSquareToSquare2& operator=(const MappingToRefSquareToSquare2&);
            virtual ~MappingToRefSquareToSquare2();
    };
    // ~~~ index 3 ~~~==============================================================================
    class MappingToRefSquareToSquare3: public MappingReferenceToReference
    {
        public:
            static const MappingToRefSquareToSquare3& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 2;}
        private:
            MappingToRefSquareToSquare3();
            MappingToRefSquareToSquare3(const MappingToRefSquareToSquare3&);
            MappingToRefSquareToSquare3& operator=(const MappingToRefSquareToSquare3&);
            virtual ~MappingToRefSquareToSquare3();
    };
    // ~~~ index 4 ~~~==============================================================================
    class MappingToRefSquareToSquare4: public MappingReferenceToReference
    {
        public:
            static const MappingToRefSquareToSquare4& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 2;}
        private:
            MappingToRefSquareToSquare4();
            MappingToRefSquareToSquare4(const MappingToRefSquareToSquare4&);
            MappingToRefSquareToSquare4& operator=(const MappingToRefSquareToSquare4&);
            virtual ~MappingToRefSquareToSquare4();
    };
    // ~~~ index 5 ~~~==============================================================================
    class MappingToRefSquareToSquare5: public MappingReferenceToReference
    {
        public:
            static const MappingToRefSquareToSquare5& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 2;}
        private:
            MappingToRefSquareToSquare5();
            MappingToRefSquareToSquare5(const MappingToRefSquareToSquare5&);
            MappingToRefSquareToSquare5& operator=(const MappingToRefSquareToSquare5&);
            virtual ~MappingToRefSquareToSquare5();
    };
    // ~~~ index 6 ~~~==============================================================================
    class MappingToRefSquareToSquare6: public MappingReferenceToReference
    {
        public:
            static const MappingToRefSquareToSquare6& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 2;}
        private:
            MappingToRefSquareToSquare6();
            MappingToRefSquareToSquare6(const MappingToRefSquareToSquare6&);
            MappingToRefSquareToSquare6& operator=(const MappingToRefSquareToSquare6&);
            virtual ~MappingToRefSquareToSquare6();
    };
    // ~~~ index 7 ~~~==============================================================================
    class MappingToRefSquareToSquare7: public MappingReferenceToReference
    {
        public:
            static const MappingToRefSquareToSquare7& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 2;}
        private:
            MappingToRefSquareToSquare7();
            MappingToRefSquareToSquare7(const MappingToRefSquareToSquare7&);
            MappingToRefSquareToSquare7& operator=(const MappingToRefSquareToSquare7&);
            virtual ~MappingToRefSquareToSquare7();
    };
};
#endif

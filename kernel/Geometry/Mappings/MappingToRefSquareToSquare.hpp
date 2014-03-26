/*
 * MappingToRefSquareToSquare.hpp
 *
 *  Created on: Feb 11, 2013
 *      Author: nicorivas
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
     * These correspond to the node mutations : (1)(2)(3)(4), (1234), (14)(23),
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

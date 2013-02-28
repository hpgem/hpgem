/*
 * MappingToRefCubeToCube.hpp
 *
 *  Created on: Feb 15, 2013
 *      Author: nicorivas
 */

#ifndef MAPPINGTOREFCUBETOCUBE_H_
#define MAPPINGTOREFCUBETOCUBE_H_

#include "MappingReferenceToReference.hpp"
#include "../Jacobian.hpp"

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
     * There are 8 possible mappings:
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
     * These correspond to the node mutations : (1)(2)(3)(4)(5)(6)(7)(8), (1234)(5678),
     * (14)(23)(58)(67), (1342) (5786), (13)(24)(57)(68), (12)(34)(56)(78),
     * (14)(2)(3)(58)(6)(7)   and (1)(4)(23) (5)(8)(67). The first 4 are rotations,
     * while the last 4 are mirrorings.
     */

    // ~~~ index 0 ~~~==============================================================================
    class MappingToRefCubeToCube0: public MappingReferenceToReference<3,3>
    {
        public:
            static const MappingToRefCubeToCube0& Instance();
            virtual void transform(const Geometry::PointReference<3>& p1,
                                         Geometry::PointReference<3>& p2) const;
            virtual void calcJacobian(const Geometry::PointReference<3>&,
                                            Geometry::Jacobian<3,3>&) const;
        private:
            MappingToRefCubeToCube0();
            MappingToRefCubeToCube0(const MappingToRefCubeToCube0&);
            MappingToRefCubeToCube0& operator=(const MappingToRefCubeToCube0&);
            virtual ~MappingToRefCubeToCube0();
    };
    // ~~~ index 1 ~~~==============================================================================
    class MappingToRefCubeToCube1: public MappingReferenceToReference<3,3>
    {
        public:
            static const MappingToRefCubeToCube1& Instance();
            virtual void transform(const Geometry::PointReference<3>& p1,
                                         Geometry::PointReference<3>& p2) const;
            virtual void calcJacobian(const Geometry::PointReference<3>&,
                                            Geometry::Jacobian<3,3>&) const;
        private:
            MappingToRefCubeToCube1();
            MappingToRefCubeToCube1(const MappingToRefCubeToCube1&);
            MappingToRefCubeToCube1& operator=(const MappingToRefCubeToCube1&);
            virtual ~MappingToRefCubeToCube1();
    };
    // ~~~ index 2 ~~~==============================================================================
    class MappingToRefCubeToCube2: public MappingReferenceToReference<3,3>
    {
        public:
            static const MappingToRefCubeToCube2& Instance();
            virtual void transform(const Geometry::PointReference<3>& p1,
                                         Geometry::PointReference<3>& p2) const;
            virtual void calcJacobian(const Geometry::PointReference<3>&,
                                            Geometry::Jacobian<3,3>&) const;
        private:
            MappingToRefCubeToCube2();
            MappingToRefCubeToCube2(const MappingToRefCubeToCube2&);
            MappingToRefCubeToCube2& operator=(const MappingToRefCubeToCube2&);
            virtual ~MappingToRefCubeToCube2();
    };
    // ~~~ index 3 ~~~==============================================================================
    class MappingToRefCubeToCube3: public MappingReferenceToReference<3,3>
    {
        public:
            static const MappingToRefCubeToCube3& Instance();
            virtual void transform(const Geometry::PointReference<3>& p1,
                                         Geometry::PointReference<3>& p2) const;
            virtual void calcJacobian(const Geometry::PointReference<3>&,
                                            Geometry::Jacobian<3,3>&) const;
        private:
            MappingToRefCubeToCube3();
            MappingToRefCubeToCube3(const MappingToRefCubeToCube3&);
            MappingToRefCubeToCube3& operator=(const MappingToRefCubeToCube3&);
            virtual ~MappingToRefCubeToCube3();
    };
    // ~~~ index 4 ~~~==============================================================================
    class MappingToRefCubeToCube4: public MappingReferenceToReference<3,3>
    {
        public:
            static const MappingToRefCubeToCube4& Instance();
            virtual void transform(const Geometry::PointReference<3>& p1,
                                         Geometry::PointReference<3>& p2) const;
            virtual void calcJacobian(const Geometry::PointReference<3>&,
                                            Geometry::Jacobian<3,3>&) const;
        private:
            MappingToRefCubeToCube4();
            MappingToRefCubeToCube4(const MappingToRefCubeToCube4&);
            MappingToRefCubeToCube4& operator=(const MappingToRefCubeToCube4&);
            virtual ~MappingToRefCubeToCube4();
    };
    // ~~~ index 5 ~~~==============================================================================
    class MappingToRefCubeToCube5: public MappingReferenceToReference<3,3>
    {
        public:
            static const MappingToRefCubeToCube5& Instance();
            virtual void transform(const Geometry::PointReference<3>& p1,
                                         Geometry::PointReference<3>& p2) const;
            virtual void calcJacobian(const Geometry::PointReference<3>&,
                                            Geometry::Jacobian<3,3>&) const;
        private:
            MappingToRefCubeToCube5();
            MappingToRefCubeToCube5(const MappingToRefCubeToCube5&);
            MappingToRefCubeToCube5& operator=(const MappingToRefCubeToCube5&);
            virtual ~MappingToRefCubeToCube5();
    };
    // ~~~ index 6 ~~~==============================================================================
    class MappingToRefCubeToCube6: public MappingReferenceToReference<3,3>
    {
        public:
            static const MappingToRefCubeToCube6& Instance();
            virtual void transform(const Geometry::PointReference<3>& p1,
                                         Geometry::PointReference<3>& p2) const;
            virtual void calcJacobian(const Geometry::PointReference<3>&,
                                            Geometry::Jacobian<3,3>&) const;
        private:
            MappingToRefCubeToCube6();
            MappingToRefCubeToCube6(const MappingToRefCubeToCube6&);
            MappingToRefCubeToCube6& operator=(const MappingToRefCubeToCube6&);
            virtual ~MappingToRefCubeToCube6();
    };
    // ~~~ index 7 ~~~==============================================================================
    class MappingToRefCubeToCube7: public MappingReferenceToReference<3,3>
    {
        public:
            static const MappingToRefCubeToCube7& Instance();
            virtual void transform(const Geometry::PointReference<3>& p1,
                                         Geometry::PointReference<3>& p2) const;
            virtual void calcJacobian(const Geometry::PointReference<3>&,
                                            Geometry::Jacobian<3,3>&) const;
        private:
            MappingToRefCubeToCube7();
            MappingToRefCubeToCube7(const MappingToRefCubeToCube7&);
            MappingToRefCubeToCube7& operator=(const MappingToRefCubeToCube7&);
            virtual ~MappingToRefCubeToCube7();
    };
};
#endif

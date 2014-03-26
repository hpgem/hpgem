/*
 * MappingToRefSquareToCube.hpp
 *
 *  Created on: Feb 15, 2013
 *      Author: nicorivas
 */

#ifndef MAPPINGSQUARETOCUBE_H_
#define MAPPINGSQUARETOCUBE_H_

#include "MappingReferenceToReference.hpp"
#include "../Jacobian.hpp"

namespace Geometry
{
    /* The ordering of the vertex in a cube:
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
     *  faces indexes:
     *              0: (0,1,2,3)
     *              1: (0,1,4,5)
     *              2: (0,2,4,6)
     *              3: (1,3,5,7)
     *              4: (2,3,6,7)
     *              5: (4,5,6,7)
     *
     * ~OC~
     *
     * This implements the mappings of a Square [-1,1]^2 onto a cube [-1,1]^3 depending on the
     * faceindex. The mappings are defined as:
     *
     *      faceindex 0: (x,y)->(x,y,-1)
     *      faceindex 1: (x,y)->(x,-1,y)
     *      faceindex 2: (x,y)->(-1,x,y)
     *      faceindex 3: (x,y)->(1,x,y)
     *      faceindex 4: (x,y)->(x,1,y)
     *      faceindex 5: (x,y)->(x,y,1)
     *
     * The mappings are chosen to preserve the ordering of the vertices. (Ordering by coordinate).
     */

    // ~~~ index 0 ~~~==============================================================================

    class MappingToRefSquareToCube0: public MappingReferenceToReference
    {
        public:
            static const MappingToRefSquareToCube0& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 3;}
        private:
            MappingToRefSquareToCube0();
            MappingToRefSquareToCube0(const MappingToRefSquareToCube0&);
            MappingToRefSquareToCube0& operator=(const MappingToRefSquareToCube0&);
            virtual ~MappingToRefSquareToCube0();
    };

    // ~~~ index 1 ~~~==============================================================================

    class MappingToRefSquareToCube1: public MappingReferenceToReference
    {
        public:
            static const MappingToRefSquareToCube1& Instance();
            virtual void transform(const Geometry::PointReference&,
                                         Geometry::PointReference&) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 3;}
        private:
            MappingToRefSquareToCube1();
            MappingToRefSquareToCube1(const MappingToRefSquareToCube1&);
            MappingToRefSquareToCube1& operator=(const MappingToRefSquareToCube1&);
            virtual ~MappingToRefSquareToCube1();
    };

    // ~~~ index 2 ~~~==============================================================================

    class MappingToRefSquareToCube2: public MappingReferenceToReference
    {
        public:
            static const MappingToRefSquareToCube2& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 3;}
        private:
            MappingToRefSquareToCube2();
            MappingToRefSquareToCube2(const MappingToRefSquareToCube2&);
            MappingToRefSquareToCube1& operator=(const MappingToRefSquareToCube2&);
            virtual ~MappingToRefSquareToCube2();
    };

    // ~~~ index 3 ~~~==============================================================================

    class MappingToRefSquareToCube3: public MappingReferenceToReference
    {
        public:
            static const MappingToRefSquareToCube3& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 3;}
        private:
            MappingToRefSquareToCube3();
            MappingToRefSquareToCube3(const MappingToRefSquareToCube3&);
            MappingToRefSquareToCube3& operator=(const MappingToRefSquareToCube3&);
            virtual ~MappingToRefSquareToCube3();
    };

    // ~~~ index 4 ~~~==============================================================================

    class MappingToRefSquareToCube4: public MappingReferenceToReference
    {
        public:
            static const MappingToRefSquareToCube4& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 3;}
        private:
            MappingToRefSquareToCube4();
            MappingToRefSquareToCube4(const MappingToRefSquareToCube4&);
            MappingToRefSquareToCube4& operator=(const MappingToRefSquareToCube4&);
            virtual ~MappingToRefSquareToCube4();
    };

    // ~~~ index 5 ~~~==============================================================================

    class MappingToRefSquareToCube5: public MappingReferenceToReference
    {
        public:
            static const MappingToRefSquareToCube5& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 3;}
        private:
            MappingToRefSquareToCube5();
            MappingToRefSquareToCube5(const MappingToRefSquareToCube5&);
            MappingToRefSquareToCube5& operator=(const MappingToRefSquareToCube5&);
            virtual ~MappingToRefSquareToCube5();
    };

};
#endif /* MAPPINGSIMPLECUBENLINEAR_H_ */

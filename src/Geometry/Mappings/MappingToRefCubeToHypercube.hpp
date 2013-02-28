/*
 * MappingToRefCubeToHypercube.hpp
 *
 *  Created on: Feb 19, 2013
 *      Author: nicorivas
 */

#ifndef MappingToRefCubeToHypercube_H_
#define MappingToRefCubeToHypercube_H_

#include "MappingReferenceToReference.hpp"
#include "../Jacobian.hpp"

namespace Geometry
{
    /*
     *
     */

    // ~~~ index 0 ~~~==============================================================================

    class MappingToRefCubeToHypercube0: public MappingReferenceToReference<3,4>
    {
        public:
            static const MappingToRefCubeToHypercube0& Instance();
            virtual void transform(const Geometry::PointReference<3>& p1,
                                         Geometry::PointReference<4>& p2) const;
            virtual void calcJacobian(const Geometry::PointReference<3>&,
                                            Geometry::Jacobian<3,4>&) const;
        private:
            MappingToRefCubeToHypercube0();
            MappingToRefCubeToHypercube0(const MappingToRefCubeToHypercube0&);
            MappingToRefCubeToHypercube0& operator=(const MappingToRefCubeToHypercube0&);
            virtual ~MappingToRefCubeToHypercube0();
    };

    // ~~~ index 1 ~~~==============================================================================

    class MappingToRefCubeToHypercube1: public MappingReferenceToReference<3,4>
    {
        public:
            static const MappingToRefCubeToHypercube1& Instance();
            virtual void transform(const Geometry::PointReference<3>&,
                                         Geometry::PointReference<4>&) const;
            virtual void calcJacobian(const Geometry::PointReference<3>&,
                                            Geometry::Jacobian<3,4>&) const;
        private:
            MappingToRefCubeToHypercube1();
            MappingToRefCubeToHypercube1(const MappingToRefCubeToHypercube1&);
            MappingToRefCubeToHypercube1& operator=(const MappingToRefCubeToHypercube1&);
            virtual ~MappingToRefCubeToHypercube1();
    };

    // ~~~ index 2 ~~~==============================================================================

    class MappingToRefCubeToHypercube2: public MappingReferenceToReference<3,4>
    {
        public:
            static const MappingToRefCubeToHypercube2& Instance();
            virtual void transform(const Geometry::PointReference<3>& p1,
                                         Geometry::PointReference<4>& p2) const;
            virtual void calcJacobian(const Geometry::PointReference<3>&,
                                            Geometry::Jacobian<3,4>&) const;
        private:
            MappingToRefCubeToHypercube2();
            MappingToRefCubeToHypercube2(const MappingToRefCubeToHypercube2&);
            MappingToRefCubeToHypercube1& operator=(const MappingToRefCubeToHypercube2&);
            virtual ~MappingToRefCubeToHypercube2();
    };

    // ~~~ index 3 ~~~==============================================================================

    class MappingToRefCubeToHypercube3: public MappingReferenceToReference<3,4>
    {
        public:
            static const MappingToRefCubeToHypercube3& Instance();
            virtual void transform(const Geometry::PointReference<3>& p1,
                                         Geometry::PointReference<4>& p2) const;
            virtual void calcJacobian(const Geometry::PointReference<3>&,
                                            Geometry::Jacobian<3,4>&) const;
        private:
            MappingToRefCubeToHypercube3();
            MappingToRefCubeToHypercube3(const MappingToRefCubeToHypercube3&);
            MappingToRefCubeToHypercube3& operator=(const MappingToRefCubeToHypercube3&);
            virtual ~MappingToRefCubeToHypercube3();
    };

    // ~~~ index 4 ~~~==============================================================================

    class MappingToRefCubeToHypercube4: public MappingReferenceToReference<3,4>
    {
        public:
            static const MappingToRefCubeToHypercube4& Instance();
            virtual void transform(const Geometry::PointReference<3>& p1,
                                         Geometry::PointReference<4>& p2) const;
            virtual void calcJacobian(const Geometry::PointReference<3>&,
                                            Geometry::Jacobian<3,4>&) const;
        private:
            MappingToRefCubeToHypercube4();
            MappingToRefCubeToHypercube4(const MappingToRefCubeToHypercube4&);
            MappingToRefCubeToHypercube4& operator=(const MappingToRefCubeToHypercube4&);
            virtual ~MappingToRefCubeToHypercube4();
    };

    // ~~~ index 5 ~~~==============================================================================

    class MappingToRefCubeToHypercube5: public MappingReferenceToReference<3,4>
    {
        public:
            static const MappingToRefCubeToHypercube5& Instance();
            virtual void transform(const Geometry::PointReference<3>& p1,
                                         Geometry::PointReference<4>& p2) const;
            virtual void calcJacobian(const Geometry::PointReference<3>&,
                                            Geometry::Jacobian<3,4>&) const;
        private:
            MappingToRefCubeToHypercube5();
            MappingToRefCubeToHypercube5(const MappingToRefCubeToHypercube5&);
            MappingToRefCubeToHypercube5& operator=(const MappingToRefCubeToHypercube5&);
            virtual ~MappingToRefCubeToHypercube5();
    };

    // ~~~ index 6 ~~~==============================================================================

    class MappingToRefCubeToHypercube6: public MappingReferenceToReference<3,4>
    {
        public:
            static const MappingToRefCubeToHypercube6& Instance();
            virtual void transform(const Geometry::PointReference<3>& p1,
                                         Geometry::PointReference<4>& p2) const;
            virtual void calcJacobian(const Geometry::PointReference<3>&,
                                            Geometry::Jacobian<3,4>&) const;
        private:
            MappingToRefCubeToHypercube6();
            MappingToRefCubeToHypercube6(const MappingToRefCubeToHypercube6&);
            MappingToRefCubeToHypercube6& operator=(const MappingToRefCubeToHypercube6&);
            virtual ~MappingToRefCubeToHypercube6();
    };

    // ~~~ index 7 ~~~==============================================================================

    class MappingToRefCubeToHypercube7: public MappingReferenceToReference<3,4>
    {
        public:
            static const MappingToRefCubeToHypercube7& Instance();
            virtual void transform(const Geometry::PointReference<3>& p1,
                                         Geometry::PointReference<4>& p2) const;
            virtual void calcJacobian(const Geometry::PointReference<3>&,
                                            Geometry::Jacobian<3,4>&) const;
        private:
            MappingToRefCubeToHypercube7();
            MappingToRefCubeToHypercube7(const MappingToRefCubeToHypercube7&);
            MappingToRefCubeToHypercube7& operator=(const MappingToRefCubeToHypercube7&);
            virtual ~MappingToRefCubeToHypercube7();
    };
};
#endif

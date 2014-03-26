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

    class MappingToRefCubeToHypercube0: public MappingReferenceToReference
    {
        public:
            static const MappingToRefCubeToHypercube0& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 4;}
        private:
            MappingToRefCubeToHypercube0();
            MappingToRefCubeToHypercube0(const MappingToRefCubeToHypercube0&);
            MappingToRefCubeToHypercube0& operator=(const MappingToRefCubeToHypercube0&);
            virtual ~MappingToRefCubeToHypercube0();
    };

    // ~~~ index 1 ~~~==============================================================================

    class MappingToRefCubeToHypercube1: public MappingReferenceToReference
    {
        public:
            static const MappingToRefCubeToHypercube1& Instance();
            virtual void transform(const Geometry::PointReference&,
                                         Geometry::PointReference&) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 4;}
        private:
            MappingToRefCubeToHypercube1();
            MappingToRefCubeToHypercube1(const MappingToRefCubeToHypercube1&);
            MappingToRefCubeToHypercube1& operator=(const MappingToRefCubeToHypercube1&);
            virtual ~MappingToRefCubeToHypercube1();
    };

    // ~~~ index 2 ~~~==============================================================================

    class MappingToRefCubeToHypercube2: public MappingReferenceToReference
    {
        public:
            static const MappingToRefCubeToHypercube2& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 4;}
        private:
            MappingToRefCubeToHypercube2();
            MappingToRefCubeToHypercube2(const MappingToRefCubeToHypercube2&);
            MappingToRefCubeToHypercube1& operator=(const MappingToRefCubeToHypercube2&);
            virtual ~MappingToRefCubeToHypercube2();
    };

    // ~~~ index 3 ~~~==============================================================================

    class MappingToRefCubeToHypercube3: public MappingReferenceToReference
    {
        public:
            static const MappingToRefCubeToHypercube3& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 4;}
        private:
            MappingToRefCubeToHypercube3();
            MappingToRefCubeToHypercube3(const MappingToRefCubeToHypercube3&);
            MappingToRefCubeToHypercube3& operator=(const MappingToRefCubeToHypercube3&);
            virtual ~MappingToRefCubeToHypercube3();
    };

    // ~~~ index 4 ~~~==============================================================================

    class MappingToRefCubeToHypercube4: public MappingReferenceToReference
    {
        public:
            static const MappingToRefCubeToHypercube4& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 4;}
        private:
            MappingToRefCubeToHypercube4();
            MappingToRefCubeToHypercube4(const MappingToRefCubeToHypercube4&);
            MappingToRefCubeToHypercube4& operator=(const MappingToRefCubeToHypercube4&);
            virtual ~MappingToRefCubeToHypercube4();
    };

    // ~~~ index 5 ~~~==============================================================================

    class MappingToRefCubeToHypercube5: public MappingReferenceToReference
    {
        public:
            static const MappingToRefCubeToHypercube5& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 4;}
        private:
            MappingToRefCubeToHypercube5();
            MappingToRefCubeToHypercube5(const MappingToRefCubeToHypercube5&);
            MappingToRefCubeToHypercube5& operator=(const MappingToRefCubeToHypercube5&);
            virtual ~MappingToRefCubeToHypercube5();
    };

    // ~~~ index 6 ~~~==============================================================================

    class MappingToRefCubeToHypercube6: public MappingReferenceToReference
    {
        public:
            static const MappingToRefCubeToHypercube6& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 4;}
        private:
            MappingToRefCubeToHypercube6();
            MappingToRefCubeToHypercube6(const MappingToRefCubeToHypercube6&);
            MappingToRefCubeToHypercube6& operator=(const MappingToRefCubeToHypercube6&);
            virtual ~MappingToRefCubeToHypercube6();
    };

    // ~~~ index 7 ~~~==============================================================================

    class MappingToRefCubeToHypercube7: public MappingReferenceToReference
    {
        public:
            static const MappingToRefCubeToHypercube7& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 4;}
        private:
            MappingToRefCubeToHypercube7();
            MappingToRefCubeToHypercube7(const MappingToRefCubeToHypercube7&);
            MappingToRefCubeToHypercube7& operator=(const MappingToRefCubeToHypercube7&);
            virtual ~MappingToRefCubeToHypercube7();
    };
};
#endif

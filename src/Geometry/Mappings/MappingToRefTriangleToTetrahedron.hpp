/*
 * MappingToRefTriangleToTetrahedron.hpp
 *
 *  Created on: Feb 17, 2013
 *      Author: nicorivas
 */

#ifndef MappingToRefTriangleToTetrahedron_H_
#define MappingToRefTriangleToTetrahedron_H_

#include "MappingReferenceToReference.hpp"
#include "../Jacobian.hpp"

namespace Geometry
{
    /* The ordering of the vertex in a cube:
     *
     *  (0,0,1) 3
     *          |\
     *          |  \2 (0,1,0)
     *          |  / \
     *          |/     \
     *  (0,0,0) 0--------1 (1,0,0)
     *
     *  faces indexes:
     *              0: (0,3,2)
     *              1: (0,1,3)
     *              2: (0,2,1)
     *              3: (1,2,3)
     */

    // ~~~ index 0 ~~~==============================================================================

    class MappingToRefTriangleToTetrahedron0: public MappingReferenceToReference<2,3>
    {
        public:
            static const MappingToRefTriangleToTetrahedron0& Instance();
            virtual void transform(const Geometry::PointReference<2>& p1,
                                         Geometry::PointReference<3>& p2) const;
            virtual void calcJacobian(const Geometry::PointReference<2>&,
                                            Geometry::Jacobian<2,3>&) const;
        private:
            MappingToRefTriangleToTetrahedron0();
            MappingToRefTriangleToTetrahedron0(const MappingToRefTriangleToTetrahedron0&);
            MappingToRefTriangleToTetrahedron0& operator=(const MappingToRefTriangleToTetrahedron0&);
            virtual ~MappingToRefTriangleToTetrahedron0();
    };

    // ~~~ index 1 ~~~==============================================================================

    class MappingToRefTriangleToTetrahedron1: public MappingReferenceToReference<2,3>
    {
        public:
            static const MappingToRefTriangleToTetrahedron1& Instance();
            virtual void transform(const Geometry::PointReference<2>&,
                                         Geometry::PointReference<3>&) const;
            virtual void calcJacobian(const Geometry::PointReference<2>&,
                                            Geometry::Jacobian<2,3>&) const;
        private:
            MappingToRefTriangleToTetrahedron1();
            MappingToRefTriangleToTetrahedron1(const MappingToRefTriangleToTetrahedron1&);
            MappingToRefTriangleToTetrahedron1& operator=(const MappingToRefTriangleToTetrahedron1&);
            virtual ~MappingToRefTriangleToTetrahedron1();
    };

    // ~~~ index 2 ~~~==============================================================================

    class MappingToRefTriangleToTetrahedron2: public MappingReferenceToReference<2,3>
    {
        public:
            static const MappingToRefTriangleToTetrahedron2& Instance();
            virtual void transform(const Geometry::PointReference<2>& p1,
                                         Geometry::PointReference<3>& p2) const;
            virtual void calcJacobian(const Geometry::PointReference<2>&,
                                            Geometry::Jacobian<2,3>&) const;
        private:
            MappingToRefTriangleToTetrahedron2();
            MappingToRefTriangleToTetrahedron2(const MappingToRefTriangleToTetrahedron2&);
            MappingToRefTriangleToTetrahedron1& operator=(const MappingToRefTriangleToTetrahedron2&);
            virtual ~MappingToRefTriangleToTetrahedron2();
    };

    // ~~~ index 3 ~~~==============================================================================

    class MappingToRefTriangleToTetrahedron3: public MappingReferenceToReference<2,3>
    {
        public:
            static const MappingToRefTriangleToTetrahedron3& Instance();
            virtual void transform(const Geometry::PointReference<2>& p1,
                                         Geometry::PointReference<3>& p2) const;
            virtual void calcJacobian(const Geometry::PointReference<2>&,
                                            Geometry::Jacobian<2,3>&) const;
        private:
            MappingToRefTriangleToTetrahedron3();
            MappingToRefTriangleToTetrahedron3(const MappingToRefTriangleToTetrahedron3&);
            MappingToRefTriangleToTetrahedron3& operator=(const MappingToRefTriangleToTetrahedron3&);
            virtual ~MappingToRefTriangleToTetrahedron3();
    };
};
#endif

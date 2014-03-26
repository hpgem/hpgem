/*
 * MappingToRefTriangleToTetrahedron.cpp
 *
 *  Created on: Feb 17, 2013
 *      Author: nicorivas
 */
#include "MappingToRefTriangleToTetrahedron.hpp"

namespace Geometry
{
    // ~~~ index 0 ~~~==============================================================================

    const MappingToRefTriangleToTetrahedron0& MappingToRefTriangleToTetrahedron0::Instance()
    {
        static const MappingToRefTriangleToTetrahedron0 theInstance;
        return theInstance;
    }

    void MappingToRefTriangleToTetrahedron0::transform(const Geometry::PointReference& p1,
                                                        Geometry::PointReference& p2) const
    {
        p2[0] = 0.0;
        p2[1] = p1[1];
        p2[2] = p1[0];
    }

    void MappingToRefTriangleToTetrahedron0::calcJacobian(const Geometry::PointReference& p1,
                                                           Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) = 0.0;
        jacobian(1,0) = 0.0;
        jacobian(2,0) = 1.0;

        jacobian(0,1) = 0.0;
        jacobian(1,1) = 1.0;
        jacobian(2,1) = 0.0;
    }

    MappingToRefTriangleToTetrahedron0::MappingToRefTriangleToTetrahedron0() { }
    MappingToRefTriangleToTetrahedron0::~MappingToRefTriangleToTetrahedron0() { }

    // ~~~ index 1 ~~~==============================================================================

    const MappingToRefTriangleToTetrahedron1& MappingToRefTriangleToTetrahedron1::Instance()
    {
        static const MappingToRefTriangleToTetrahedron1 theInstance;
        return theInstance;
    }

    void MappingToRefTriangleToTetrahedron1::transform(const Geometry::PointReference& p1,
                                                        Geometry::PointReference& p2) const
    {
        p2[0] = p1[0];
        p2[1] = 0.0;
        p2[2] = p1[1];
    }

    void MappingToRefTriangleToTetrahedron1::calcJacobian(const Geometry::PointReference& p1,
                                                           Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) = 1.0;
        jacobian(1,0) = 0.0;
        jacobian(2,0) = 0.0;

        jacobian(0,1) = 0.0;
        jacobian(1,1) = 0.0;
        jacobian(2,1) = 1.0;
    }

    MappingToRefTriangleToTetrahedron1::MappingToRefTriangleToTetrahedron1() { }
    MappingToRefTriangleToTetrahedron1::~MappingToRefTriangleToTetrahedron1() { }

    // ~~~ index 2 ~~~==============================================================================

    const MappingToRefTriangleToTetrahedron2& MappingToRefTriangleToTetrahedron2::Instance()
    {
        static const MappingToRefTriangleToTetrahedron2 theInstance;
        return theInstance;
    }

    void MappingToRefTriangleToTetrahedron2::transform(const Geometry::PointReference& p1,
                                                        Geometry::PointReference& p2) const
    {
        p2[0] = p1[1];
        p2[1] = p1[0];
        p2[2] = 0.0;
    }

    void MappingToRefTriangleToTetrahedron2::calcJacobian(const Geometry::PointReference& p1,
                                                           Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) = 0.0;
        jacobian(1,0) = 1.0;
        jacobian(2,0) = 0.0;

        jacobian(0,1) = 1.0;
        jacobian(1,1) = 0.0;
        jacobian(2,1) = 0.0;
    }

    MappingToRefTriangleToTetrahedron2::MappingToRefTriangleToTetrahedron2() { }
    MappingToRefTriangleToTetrahedron2::~MappingToRefTriangleToTetrahedron2() { }

    // ~~~ index 3 ~~~==============================================================================

    const MappingToRefTriangleToTetrahedron3& MappingToRefTriangleToTetrahedron3::Instance()
    {
        static const MappingToRefTriangleToTetrahedron3 theInstance;
        return theInstance;
    }

    void MappingToRefTriangleToTetrahedron3::transform(const Geometry::PointReference& p1,
                                                        Geometry::PointReference& p2) const
    {
        p2[0] = 1.0 - p1[0] - p1[1];
        p2[1] = p1[0];
        p2[2] = p1[1];
    }

    void MappingToRefTriangleToTetrahedron3::calcJacobian(const Geometry::PointReference& p1,
                                                           Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) = -1.0;
        jacobian(1,0) =  1.0;
        jacobian(2,0) =  0.0;

        jacobian(0,1) = -1.0;
        jacobian(1,1) =  0.0;
        jacobian(2,1) =  1.0;
    }

    MappingToRefTriangleToTetrahedron3::MappingToRefTriangleToTetrahedron3() { }
    MappingToRefTriangleToTetrahedron3::~MappingToRefTriangleToTetrahedron3() { }
}

/*
 * MappingToRefTriangleToTriangle.cpp
 *
 *  Created on: Feb 17, 2013
 *      Author: nicorivas
 */

#ifndef MappingToRefTriangleToTriangle_C_
#define MappingToRefTriangleToTriangle_C_

#include "MappingToRefTriangleToTriangle.hpp"
#include "../Jacobian.hpp"

namespace Geometry
{
    // ~~~~~~~~~~~~~~~==============================================================================
    // ~~~ index 0 ~~~==============================================================================
    // ~~~~~~~~~~~~~~~==============================================================================

    const MappingToRefTriangleToTriangle0& MappingToRefTriangleToTriangle0::Instance()
    {
        static const MappingToRefTriangleToTriangle0 theInstance;
        return theInstance;
    }

    void MappingToRefTriangleToTriangle0::transform(const Geometry::PointReference& p1,
                                                 Geometry::PointReference& p2) const
    {
        p2[0] = p1[0];
        p2[1] = p1[1];
    }

    void MappingToRefTriangleToTriangle0::calcJacobian(const Geometry::PointReference&,
                                                    Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) = 1.0; jacobian(0,1) = 0.0;
        jacobian(1,0) = 0.0; jacobian(1,1) = 1.0;
    }

    MappingToRefTriangleToTriangle0::MappingToRefTriangleToTriangle0() { }
    MappingToRefTriangleToTriangle0::~MappingToRefTriangleToTriangle0() { }

    // ~~~~~~~~~~~~~~~==============================================================================
    // ~~~ index 1 ~~~==============================================================================
    // ~~~~~~~~~~~~~~~==============================================================================

    const MappingToRefTriangleToTriangle1& MappingToRefTriangleToTriangle1::Instance()
    {
        static const MappingToRefTriangleToTriangle1 theInstance;
        return theInstance;
    }

    void MappingToRefTriangleToTriangle1::transform(const Geometry::PointReference& p1,
                                                 Geometry::PointReference& p2) const
    {
        p2[0] = p1[1];
        p2[1] = p1[0];
    }

    void MappingToRefTriangleToTriangle1::calcJacobian(const Geometry::PointReference&,
                                                    Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) =  0.0; jacobian(0,1) = 1.0;
        jacobian(1,0) =  1.0; jacobian(1,1) = 0.0;
    }

    MappingToRefTriangleToTriangle1::MappingToRefTriangleToTriangle1() { }
    MappingToRefTriangleToTriangle1::~MappingToRefTriangleToTriangle1() { }

    // ~~~~~~~~~~~~~~~==============================================================================
    // ~~~ index 2 ~~~==============================================================================
    // ~~~~~~~~~~~~~~~==============================================================================

    const MappingToRefTriangleToTriangle2& MappingToRefTriangleToTriangle2::Instance()
    {
        static const MappingToRefTriangleToTriangle2 theInstance;
        return theInstance;
    }

    void MappingToRefTriangleToTriangle2::transform(const Geometry::PointReference& p1,
                                                 Geometry::PointReference& p2) const
    {
        p2[0] = 1.0 - p1[0] - p1[1];
        p2[1] = p1[0];
    }

    void MappingToRefTriangleToTriangle2::calcJacobian(const Geometry::PointReference&,
                                                    Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) = -1.0; jacobian(0,1) = -1.0;
        jacobian(1,0) =  1.0; jacobian(1,1) =  0.0;
    }

    MappingToRefTriangleToTriangle2::MappingToRefTriangleToTriangle2() { }
    MappingToRefTriangleToTriangle2::~MappingToRefTriangleToTriangle2() { }

    // ~~~~~~~~~~~~~~~==============================================================================
    // ~~~ index 3 ~~~==============================================================================
    // ~~~~~~~~~~~~~~~==============================================================================

    const MappingToRefTriangleToTriangle3& MappingToRefTriangleToTriangle3::Instance()
    {
        static const MappingToRefTriangleToTriangle3 theInstance;
        return theInstance;
    }

    void MappingToRefTriangleToTriangle3::transform(const Geometry::PointReference& p1,
                                                 Geometry::PointReference& p2) const
    {
        p2[0] = 1.0 - p1[0] - p1[1];
        p2[1] = p1[1];
    }

    void MappingToRefTriangleToTriangle3::calcJacobian(const Geometry::PointReference&,
                                                    Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) = -1.0; jacobian(0,1) = -1.0;
        jacobian(1,0) =  0.0; jacobian(1,1) =  1.0;
    }

    MappingToRefTriangleToTriangle3::MappingToRefTriangleToTriangle3() { }
    MappingToRefTriangleToTriangle3::~MappingToRefTriangleToTriangle3() { }

    // ~~~~~~~~~~~~~~~==============================================================================
    // ~~~ index 4 ~~~==============================================================================
    // ~~~~~~~~~~~~~~~==============================================================================

    const MappingToRefTriangleToTriangle4& MappingToRefTriangleToTriangle4::Instance()
    {
        static const MappingToRefTriangleToTriangle4 theInstance;
        return theInstance;
    }

    void MappingToRefTriangleToTriangle4::transform(const Geometry::PointReference& p1,
                                                 Geometry::PointReference& p2) const
    {
        p2[0] = p1[0];
        p2[1] = 1.0 - p1[0] - p1[1];
    }

    void MappingToRefTriangleToTriangle4::calcJacobian(const Geometry::PointReference&,
                                                    Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) =  1.0; jacobian(0,1) =  0.0;
        jacobian(1,0) = -1.0; jacobian(1,1) = -1.0;
    }

    MappingToRefTriangleToTriangle4::MappingToRefTriangleToTriangle4() { }
    MappingToRefTriangleToTriangle4::~MappingToRefTriangleToTriangle4() { }

    // ~~~~~~~~~~~~~~~==============================================================================
    // ~~~ index 5 ~~~==============================================================================
    // ~~~~~~~~~~~~~~~==============================================================================

    const MappingToRefTriangleToTriangle5& MappingToRefTriangleToTriangle5::Instance()
    {
        static const MappingToRefTriangleToTriangle5 theInstance;
        return theInstance;
    }

    void MappingToRefTriangleToTriangle5::transform(const Geometry::PointReference& p1,
                                                 Geometry::PointReference& p2) const
    {
        p2[0] = p1[1];
        p2[1] = 1.0 - p1[0] - p1[1];
    }

    void MappingToRefTriangleToTriangle5::calcJacobian(const Geometry::PointReference&,
                                                    Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) =  0.0; jacobian(0,1) =  1.0;
        jacobian(1,0) = -1.0; jacobian(1,1) = -1.0;
    }

    MappingToRefTriangleToTriangle5::MappingToRefTriangleToTriangle5() { }
    MappingToRefTriangleToTriangle5::~MappingToRefTriangleToTriangle5() { }
};
#endif /* MAPPINGSIMPLECUBENLINEAR_C_ */

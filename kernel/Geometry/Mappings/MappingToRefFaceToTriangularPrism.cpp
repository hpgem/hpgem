/*
 * MappingToRefFaceToTriangularPrism.cpp
 *
 *  Created on: Feb 17, 2013
 *      Author: nicorivas
 */
#include "MappingToRefFaceToTriangularPrism.hpp"

namespace Geometry
{
    // ~~~ index 0 ~~~==============================================================================

    const MappingToRefFaceToTriangularPrism0& MappingToRefFaceToTriangularPrism0::Instance()
    {
        static const MappingToRefFaceToTriangularPrism0 theInstance;
        return theInstance;
    }

    void MappingToRefFaceToTriangularPrism0::transform(const Geometry::PointReference& p1,
                                                        Geometry::PointReference& p2) const
    {
        p2[0] =  p1[1];
        p2[1] =  p1[0];
        p2[2] = -1.0;
    }

    void MappingToRefFaceToTriangularPrism0::calcJacobian(const Geometry::PointReference& p1,
                                                           Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) = 0.0;
        jacobian(1,0) = 1.0;
        jacobian(2,0) = 0.0;

        jacobian(0,1) = 1.0;
        jacobian(1,1) = 0.0;
        jacobian(2,1) = 0.0;
    }

    MappingToRefFaceToTriangularPrism0::MappingToRefFaceToTriangularPrism0() { }
    MappingToRefFaceToTriangularPrism0::~MappingToRefFaceToTriangularPrism0() { }

    // ~~~ index 1 ~~~==============================================================================

    const MappingToRefFaceToTriangularPrism1& MappingToRefFaceToTriangularPrism1::Instance()
    {
        static const MappingToRefFaceToTriangularPrism1 theInstance;
        return theInstance;
    }

    void MappingToRefFaceToTriangularPrism1::transform(const Geometry::PointReference& p1,
                                                        Geometry::PointReference& p2) const
    {
        p2[0] = p1[0];
        p2[1] = p1[1];
        p2[2] = 1.0;
    }

    void MappingToRefFaceToTriangularPrism1::calcJacobian(const Geometry::PointReference& p1,
                                                           Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) = 1.0;
        jacobian(1,0) = 0.0;
        jacobian(2,0) = 0.0;

        jacobian(0,1) = 0.0;
        jacobian(1,1) = 1.0;
        jacobian(2,1) = 0.0;
    }

    MappingToRefFaceToTriangularPrism1::MappingToRefFaceToTriangularPrism1() { }
    MappingToRefFaceToTriangularPrism1::~MappingToRefFaceToTriangularPrism1() { }

    // ~~~ index 2 ~~~==============================================================================

    const MappingToRefFaceToTriangularPrism2& MappingToRefFaceToTriangularPrism2::Instance()
    {
        static const MappingToRefFaceToTriangularPrism2 theInstance;
        return theInstance;
    }

    void MappingToRefFaceToTriangularPrism2::transform(const Geometry::PointReference& p1,
                                                        Geometry::PointReference& p2) const
    {
        p2[0] = 0.0;
        p2[1] = 0.5 * (1.0 - p1[0]);
        p2[2] = p1[1];
    }

    void MappingToRefFaceToTriangularPrism2::calcJacobian(const Geometry::PointReference& p1,
                                                           Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) =  0.0;
        jacobian(1,0) = -0.5;
        jacobian(2,0) =  0.0;

        jacobian(0,1) =  0.0;
        jacobian(1,1) =  0.0;
        jacobian(2,1) =  1.0;
    }

    MappingToRefFaceToTriangularPrism2::MappingToRefFaceToTriangularPrism2() { }
    MappingToRefFaceToTriangularPrism2::~MappingToRefFaceToTriangularPrism2() { }

    // ~~~ index 3 ~~~==============================================================================

    const MappingToRefFaceToTriangularPrism3& MappingToRefFaceToTriangularPrism3::Instance()
    {
        static const MappingToRefFaceToTriangularPrism3 theInstance;
        return theInstance;
    }

    void MappingToRefFaceToTriangularPrism3::transform(const Geometry::PointReference& p1,
                                                        Geometry::PointReference& p2) const
    {
        p2[0] = 0.5 * (1.0 + p1[0]);
        p2[1] = 0.0;
        p2[2] = p1[1];
    }

    void MappingToRefFaceToTriangularPrism3::calcJacobian(const Geometry::PointReference& p1,
                                                           Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) =  0.5;
        jacobian(1,0) =  0.0;
        jacobian(2,0) =  0.0;

        jacobian(0,1) =  0.0;
        jacobian(1,1) =  0.0;
        jacobian(2,1) =  1.0;
    }

    MappingToRefFaceToTriangularPrism3::MappingToRefFaceToTriangularPrism3() { }
    MappingToRefFaceToTriangularPrism3::~MappingToRefFaceToTriangularPrism3() { }

    // ~~~ index 4 ~~~==============================================================================

    const MappingToRefFaceToTriangularPrism4& MappingToRefFaceToTriangularPrism4::Instance()
    {
        static const MappingToRefFaceToTriangularPrism4 theInstance;
        return theInstance;
    }

    void MappingToRefFaceToTriangularPrism4::transform(const Geometry::PointReference& p1,
                                                Geometry::PointReference& p2) const
    {
        p2[0] = 0.5 * (1.0 - p1[0]);
        p2[1] = 0.5 * (1.0 + p1[0]);
        p2[2] = p1[1];
    }

    void MappingToRefFaceToTriangularPrism4::calcJacobian(const Geometry::PointReference& p1,
                                                   Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) = -0.5;
        jacobian(1,0) =  0.5;
        jacobian(2,0) =  0.0;

        jacobian(0,1) =  0.0;
        jacobian(1,1) =  0.0;
        jacobian(2,1) =  1.0;
    }

    MappingToRefFaceToTriangularPrism4::MappingToRefFaceToTriangularPrism4() { }
    MappingToRefFaceToTriangularPrism4::~MappingToRefFaceToTriangularPrism4() { }

}

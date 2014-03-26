/*
 * MappingToRefCubeToHypercube.cpp
 *
 *  Created on: Feb 17, 2013
 *      Author: nicorivas
 */
#include "MappingToRefCubeToHypercube.hpp"

namespace Geometry
{
    // ~~~ index 0 ~~~==============================================================================

    const MappingToRefCubeToHypercube0& MappingToRefCubeToHypercube0::Instance()
    {
        static const MappingToRefCubeToHypercube0 theInstance;
        return theInstance;
    }

    void MappingToRefCubeToHypercube0::transform(const Geometry::PointReference& p1,
                                                  Geometry::PointReference& p2) const
    {
        p2[0] =  p1[0];
        p2[1] =  p1[1];
        p2[2] =  p1[2];
        p2[3] = -1.0;
    }

    void MappingToRefCubeToHypercube0::calcJacobian(const Geometry::PointReference& p1,
                                                     Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) = 1.0;
        jacobian(1,0) = 0.0;
        jacobian(2,0) = 0.0;
        jacobian(3,0) = 0.0;

        jacobian(0,1) = 0.0;
        jacobian(1,1) = 1.0;
        jacobian(2,1) = 0.0;
        jacobian(3,1) = 0.0;

        jacobian(0,1) = 0.0;
        jacobian(1,1) = 0.0;
        jacobian(2,1) = 1.0;
        jacobian(3,1) = 0.0;
    }

    MappingToRefCubeToHypercube0::MappingToRefCubeToHypercube0() { }
    MappingToRefCubeToHypercube0::~MappingToRefCubeToHypercube0() { }

    // ~~~ index 1 ~~~==============================================================================

    const MappingToRefCubeToHypercube1& MappingToRefCubeToHypercube1::Instance()
    {
        static const MappingToRefCubeToHypercube1 theInstance;
        return theInstance;
    }

    void MappingToRefCubeToHypercube1::transform(const Geometry::PointReference& p1,
                                                  Geometry::PointReference& p2) const
    {
        p2[0] =  p1[0];
        p2[1] =  p1[1];
        p2[2] = -1.0;
        p2[3] =  p1[2];
    }

    void MappingToRefCubeToHypercube1::calcJacobian(const Geometry::PointReference& p1,
                                                     Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) = 1.0;
        jacobian(1,0) = 0.0;
        jacobian(2,0) = 0.0;
        jacobian(3,0) = 0.0;

        jacobian(0,1) = 0.0;
        jacobian(1,1) = 1.0;
        jacobian(2,1) = 0.0;
        jacobian(3,1) = 0.0;

        jacobian(0,1) = 0.0;
        jacobian(1,1) = 0.0;
        jacobian(2,1) = 0.0;
        jacobian(3,1) = 1.0;
    }

    MappingToRefCubeToHypercube1::MappingToRefCubeToHypercube1() { }
    MappingToRefCubeToHypercube1::~MappingToRefCubeToHypercube1() { }

    // ~~~ index 2 ~~~==============================================================================

    const MappingToRefCubeToHypercube2& MappingToRefCubeToHypercube2::Instance()
    {
        static const MappingToRefCubeToHypercube2 theInstance;
        return theInstance;
    }

    void MappingToRefCubeToHypercube2::transform(const Geometry::PointReference& p1,
                                                  Geometry::PointReference& p2) const
    {
        p2[0] =  p1[0];
        p2[1] =  -1.0;
        p2[2] =  p1[1];
        p2[3] =  p1[2];
    }

    void MappingToRefCubeToHypercube2::calcJacobian(const Geometry::PointReference& p1,
                                                     Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) = 1.0;
        jacobian(1,0) = 0.0;
        jacobian(2,0) = 0.0;
        jacobian(3,0) = 0.0;

        jacobian(0,1) = 0.0;
        jacobian(1,1) = 0.0;
        jacobian(2,1) = 1.0;
        jacobian(3,1) = 0.0;

        jacobian(0,1) = 0.0;
        jacobian(1,1) = 0.0;
        jacobian(2,1) = 0.0;
        jacobian(3,1) = 1.0;
    }

    MappingToRefCubeToHypercube2::MappingToRefCubeToHypercube2() { }
    MappingToRefCubeToHypercube2::~MappingToRefCubeToHypercube2() { }

    // ~~~ index 3 ~~~==============================================================================

    const MappingToRefCubeToHypercube3& MappingToRefCubeToHypercube3::Instance()
    {
        static const MappingToRefCubeToHypercube3 theInstance;
        return theInstance;
    }

    void MappingToRefCubeToHypercube3::transform(const Geometry::PointReference& p1,
                                                  Geometry::PointReference& p2) const
    {
        p2[0] = -1.0;
        p2[1] =  p1[0];
        p2[2] =  p1[1];
        p2[3] =  p1[2];
    }

    void MappingToRefCubeToHypercube3::calcJacobian(const Geometry::PointReference& p1,
                                                     Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) = 0.0;
        jacobian(1,0) = 1.0;
        jacobian(2,0) = 0.0;
        jacobian(3,0) = 0.0;

        jacobian(0,1) = 0.0;
        jacobian(1,1) = 0.0;
        jacobian(2,1) = 1.0;
        jacobian(3,1) = 0.0;

        jacobian(0,1) = 0.0;
        jacobian(1,1) = 0.0;
        jacobian(2,1) = 0.0;
        jacobian(3,1) = 1.0;
    }

    MappingToRefCubeToHypercube3::MappingToRefCubeToHypercube3() { }
    MappingToRefCubeToHypercube3::~MappingToRefCubeToHypercube3() { }

    // ~~~ index 4 ~~~==============================================================================

    const MappingToRefCubeToHypercube4& MappingToRefCubeToHypercube4::Instance()
    {
        static const MappingToRefCubeToHypercube4 theInstance;
        return theInstance;
    }

    void MappingToRefCubeToHypercube4::transform(const Geometry::PointReference& p1,
                                                  Geometry::PointReference& p2) const
    {
        p2[0] = +1.0;
        p2[1] =  p1[0];
        p2[2] =  p1[1];
        p2[3] =  p1[2];
    }

    void MappingToRefCubeToHypercube4::calcJacobian(const Geometry::PointReference& p1,
                                                     Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) = 0.0;
        jacobian(1,0) = 1.0;
        jacobian(2,0) = 0.0;
        jacobian(3,0) = 0.0;

        jacobian(0,1) = 0.0;
        jacobian(1,1) = 0.0;
        jacobian(2,1) = 1.0;
        jacobian(3,1) = 0.0;

        jacobian(0,1) = 0.0;
        jacobian(1,1) = 0.0;
        jacobian(2,1) = 0.0;
        jacobian(3,1) = 1.0;
    }

    MappingToRefCubeToHypercube4::MappingToRefCubeToHypercube4() { }
    MappingToRefCubeToHypercube4::~MappingToRefCubeToHypercube4() { }

    // ~~~ index 5 ~~~==============================================================================

    const MappingToRefCubeToHypercube5& MappingToRefCubeToHypercube5::Instance()
    {
        static const MappingToRefCubeToHypercube5 theInstance;
        return theInstance;
    }

    void MappingToRefCubeToHypercube5::transform(const Geometry::PointReference& p1,
                                                  Geometry::PointReference& p2) const
    {
        p2[0] =  p1[0];
        p2[1] = +1.0;
        p2[2] =  p1[1];
        p2[3] =  p1[2];
    }

    void MappingToRefCubeToHypercube5::calcJacobian(const Geometry::PointReference& p1,
                                                     Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) = 1.0;
        jacobian(1,0) = 0.0;
        jacobian(2,0) = 0.0;
        jacobian(3,0) = 0.0;

        jacobian(0,1) = 0.0;
        jacobian(1,1) = 0.0;
        jacobian(2,1) = 1.0;
        jacobian(3,1) = 0.0;

        jacobian(0,1) = 0.0;
        jacobian(1,1) = 0.0;
        jacobian(2,1) = 0.0;
        jacobian(3,1) = 1.0;
    }

    MappingToRefCubeToHypercube5::MappingToRefCubeToHypercube5() { }
    MappingToRefCubeToHypercube5::~MappingToRefCubeToHypercube5() { }

    // ~~~ index 6 ~~~==============================================================================

    const MappingToRefCubeToHypercube6& MappingToRefCubeToHypercube6::Instance()
    {
        static const MappingToRefCubeToHypercube6 theInstance;
        return theInstance;
    }

    void MappingToRefCubeToHypercube6::transform(const Geometry::PointReference& p1,
                                                  Geometry::PointReference& p2) const
    {
        p2[0] =  p1[0];
        p2[1] =  p1[1];
        p2[2] = +1.0;
        p2[3] =  p1[2];
    }

    void MappingToRefCubeToHypercube6::calcJacobian(const Geometry::PointReference& p1,
                                                     Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) = 1.0;
        jacobian(1,0) = 0.0;
        jacobian(2,0) = 0.0;
        jacobian(3,0) = 0.0;

        jacobian(0,1) = 0.0;
        jacobian(1,1) = 1.0;
        jacobian(2,1) = 0.0;
        jacobian(3,1) = 0.0;

        jacobian(0,1) = 0.0;
        jacobian(1,1) = 0.0;
        jacobian(2,1) = 0.0;
        jacobian(3,1) = 1.0;
    }

    MappingToRefCubeToHypercube6::MappingToRefCubeToHypercube6() { }
    MappingToRefCubeToHypercube6::~MappingToRefCubeToHypercube6() { }

    // ~~~ index 7 ~~~==============================================================================

    const MappingToRefCubeToHypercube7& MappingToRefCubeToHypercube7::Instance()
    {
        static const MappingToRefCubeToHypercube7 theInstance;
        return theInstance;
    }

    void MappingToRefCubeToHypercube7::transform(const Geometry::PointReference& p1,
                                                  Geometry::PointReference& p2) const
    {
        p2[0] =  p1[0];
        p2[1] =  p1[1];
        p2[2] =  p1[2];
        p2[3] = +1.0;

    }

    void MappingToRefCubeToHypercube7::calcJacobian(const Geometry::PointReference& p1,
                                                     Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) = 1.0;
        jacobian(1,0) = 0.0;
        jacobian(2,0) = 0.0;
        jacobian(3,0) = 0.0;

        jacobian(0,1) = 0.0;
        jacobian(1,1) = 1.0;
        jacobian(2,1) = 0.0;
        jacobian(3,1) = 0.0;

        jacobian(0,1) = 0.0;
        jacobian(1,1) = 0.0;
        jacobian(2,1) = 1.0;
        jacobian(3,1) = 0.0;
    }

    MappingToRefCubeToHypercube7::MappingToRefCubeToHypercube7() { }
    MappingToRefCubeToHypercube7::~MappingToRefCubeToHypercube7() { }
}

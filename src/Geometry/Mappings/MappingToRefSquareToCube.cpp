/*
 * MappingToRefSquareToCube.cpp
 *
 *  Created on: Feb 15, 2013
 *      Author: nicorivas
 */
#include "MappingToRefSquareToCube.hpp"

namespace Geometry
{
    // ~~~ index 0 ~~~==============================================================================

    const MappingToRefSquareToCube0& MappingToRefSquareToCube0::Instance()
    {
        static const MappingToRefSquareToCube0 theInstance;
        return theInstance;
    }

    void MappingToRefSquareToCube0::transform(const Geometry::PointReference<2>& p1,
                                               Geometry::PointReference<3>& p2) const
    {
        p2[0] = p1[0];
        p2[1] = p1[1];
        p2[2] = -1.0;
    }

    void MappingToRefSquareToCube0::calcJacobian(const Geometry::PointReference<2>& p1,
                                                  Geometry::Jacobian<2,3>& jacobian) const
    {
        jacobian(0,0) = 1.0;
        jacobian(1,0) = 0.0;
        jacobian(2,0) = 0.0;

        jacobian(0,1) = 0.0;
        jacobian(1,1) = 1.0;
        jacobian(2,1) = 0.0;
    }

    MappingToRefSquareToCube0::MappingToRefSquareToCube0() { }
    MappingToRefSquareToCube0::~MappingToRefSquareToCube0() { }

    // ~~~ index 1 ~~~==============================================================================

    const MappingToRefSquareToCube1& MappingToRefSquareToCube1::Instance()
    {
        static const MappingToRefSquareToCube1 theInstance;
        return theInstance;
    }

    void MappingToRefSquareToCube1::transform(const Geometry::PointReference<2>& p1,
                                               Geometry::PointReference<3>& p2) const
    {
        p2[0] = p1[0];
        p2[1] = -1.0;
        p2[2] = p1[1];
    }

    void MappingToRefSquareToCube1::calcJacobian(const Geometry::PointReference<2>& p1,
                                                  Geometry::Jacobian<2,3>& jacobian) const
    {
        jacobian(0,0) = 1.0;
        jacobian(1,0) = 0.0;
        jacobian(2,0) = 0.0;

        jacobian(0,1) = 0.0;
        jacobian(1,1) = 0.0;
        jacobian(2,1) = 1.0;
    }

    MappingToRefSquareToCube1::MappingToRefSquareToCube1() { }
    MappingToRefSquareToCube1::~MappingToRefSquareToCube1() { }

    // ~~~ index 2 ~~~==============================================================================

    const MappingToRefSquareToCube2& MappingToRefSquareToCube2::Instance()
    {
        static const MappingToRefSquareToCube2 theInstance;
        return theInstance;
    }

    void MappingToRefSquareToCube2::transform(const Geometry::PointReference<2>& p1,
                                               Geometry::PointReference<3>& p2) const
    {
        p2[0] = -1.0;
        p2[1] = p1[0];
        p2[2] = p1[1];
    }

    void MappingToRefSquareToCube2::calcJacobian(const Geometry::PointReference<2>& p1,
                                                  Geometry::Jacobian<2,3>& jacobian) const
    {
        jacobian(0,0) = 0.0;
        jacobian(1,0) = 1.0;
        jacobian(2,0) = 0.0;

        jacobian(0,1) = 0.0;
        jacobian(1,1) = 0.0;
        jacobian(2,1) = 1.0;
    }

    MappingToRefSquareToCube2::MappingToRefSquareToCube2() { }
    MappingToRefSquareToCube2::~MappingToRefSquareToCube2() { }

    // ~~~ index 3 ~~~==============================================================================

    const MappingToRefSquareToCube3& MappingToRefSquareToCube3::Instance()
    {
        static const MappingToRefSquareToCube3 theInstance;
        return theInstance;
    }

    void MappingToRefSquareToCube3::transform(const Geometry::PointReference<2>& p1,
                                               Geometry::PointReference<3>& p2) const
    {
        p2[0] = 1.0;
        p2[1] = p1[0];
        p2[2] = p1[1];
    }

    void MappingToRefSquareToCube3::calcJacobian(const Geometry::PointReference<2>& p1,
                                                  Geometry::Jacobian<2,3>& jacobian) const
    {
        jacobian(0,0) = 0.0;
        jacobian(1,0) = 1.0;
        jacobian(2,0) = 0.0;

        jacobian(0,1) = 0.0;
        jacobian(1,1) = 0.0;
        jacobian(2,1) = 1.0;
    }

    MappingToRefSquareToCube3::MappingToRefSquareToCube3() { }
    MappingToRefSquareToCube3::~MappingToRefSquareToCube3() { }

    // ~~~ index 4 ~~~==============================================================================

    const MappingToRefSquareToCube4& MappingToRefSquareToCube4::Instance()
    {
        static const MappingToRefSquareToCube4 theInstance;
        return theInstance;
    }

    void MappingToRefSquareToCube4::transform(const Geometry::PointReference<2>& p1,
                                               Geometry::PointReference<3>& p2) const
    {
        p2[0] = p1[0];
        p2[1] = 1.0;
        p2[2] = p1[1];
    }

    void MappingToRefSquareToCube4::calcJacobian(const Geometry::PointReference<2>& p1,
                                                  Geometry::Jacobian<2,3>& jacobian) const
    {
        jacobian(0,0) = 1.0;
        jacobian(1,0) = 0.0;
        jacobian(2,0) = 0.0;

        jacobian(0,1) = 0.0;
        jacobian(1,1) = 0.0;
        jacobian(2,1) = 1.0;
    }

    MappingToRefSquareToCube4::MappingToRefSquareToCube4() { }
    MappingToRefSquareToCube4::~MappingToRefSquareToCube4() { }

    // ~~~ index 5 ~~~==============================================================================

    const MappingToRefSquareToCube5& MappingToRefSquareToCube5::Instance()
    {
        static const MappingToRefSquareToCube5 theInstance;
        return theInstance;
    }

    void MappingToRefSquareToCube5::transform(const Geometry::PointReference<2>& p1,
                                               Geometry::PointReference<3>& p2) const
    {
        p2[0] = p1[0];
        p2[1] = p1[1];
        p2[2] = 1.0;
    }

    void MappingToRefSquareToCube5::calcJacobian(const Geometry::PointReference<2>& p1,
                                                  Geometry::Jacobian<2,3>& jacobian) const
    {
        jacobian(0,0) = 1.0;
        jacobian(1,0) = 0.0;
        jacobian(2,0) = 0.0;

        jacobian(0,1) = 0.0;
        jacobian(1,1) = 1.0;
        jacobian(2,1) = 0.0;
    }

    MappingToRefSquareToCube5::MappingToRefSquareToCube5() { }
    MappingToRefSquareToCube5::~MappingToRefSquareToCube5() { }


}

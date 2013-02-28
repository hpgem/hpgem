/*
 * MappingToRefLineToSquare.cpp
 *
 *  Created on: Feb 11, 2013
 *      Author: nicorivas
 */
#include "MappingToRefLineToSquare.hpp"

namespace Geometry
{
    // ~~~~~~~~~~~~~~~==============================================================================
    // ~~~ index 0 ~~~==============================================================================
    // ~~~~~~~~~~~~~~~==============================================================================

    const MappingToRefLineToSquare0& MappingToRefLineToSquare0::Instance()
    {
        static const MappingToRefLineToSquare0 theInstance;
        return theInstance;
    }

    void MappingToRefLineToSquare0::transform(const Geometry::PointReference<1>& p1,
                                               Geometry::PointReference<2>& p2) const
    {
        p2[0] = p1[0];
        p2[1] = -1.0;
    }

    void MappingToRefLineToSquare0::calcJacobian(const Geometry::PointReference<1>& p1,
                                                  Geometry::Jacobian<1,2>& jacobian) const
    {
        jacobian(0,0) = 1.0;
        jacobian(1,0) = 0.0;
    }

    MappingToRefLineToSquare0::MappingToRefLineToSquare0() { }
    MappingToRefLineToSquare0::~MappingToRefLineToSquare0() { }

    // ~~~~~~~~~~~~~~~==============================================================================
    // ~~~ index 1 ~~~==============================================================================
    // ~~~~~~~~~~~~~~~==============================================================================

    const MappingToRefLineToSquare1& MappingToRefLineToSquare1::Instance()
    {
        static const MappingToRefLineToSquare1 theInstance;
        return theInstance;
    }

    void MappingToRefLineToSquare1::transform(const Geometry::PointReference<1>& p1,
                                               Geometry::PointReference<2>& p2) const
    {
        p2[0] = -1.0;
        p2[1] = p1[0];
    }

    void MappingToRefLineToSquare1::calcJacobian(const Geometry::PointReference<1>& p1,
                                                  Geometry::Jacobian<1,2>& jacobian) const
    {
        jacobian(0,0) = 0.0;
        jacobian(1,0) = 1.0;
    }

    MappingToRefLineToSquare1::MappingToRefLineToSquare1() { }
    MappingToRefLineToSquare1::~MappingToRefLineToSquare1() { }

    // ~~~~~~~~~~~~~~~==============================================================================
    // ~~~ index 2 ~~~==============================================================================
    // ~~~~~~~~~~~~~~~==============================================================================

    const MappingToRefLineToSquare2& MappingToRefLineToSquare2::Instance()
    {
        static const MappingToRefLineToSquare2 theInstance;
        return theInstance;
    }

    void MappingToRefLineToSquare2::transform(const Geometry::PointReference<1>& p1,
                                               Geometry::PointReference<2>& p2) const
    {
        p2[0] = 1.0;
        p2[1] = p1[0];
    }

    void MappingToRefLineToSquare2::calcJacobian(const Geometry::PointReference<1>& p1,
                                                  Geometry::Jacobian<1,2>& jacobian) const
    {
        jacobian(0,0) = 0.0;
        jacobian(1,0) = 1.0;
    }

    MappingToRefLineToSquare2::MappingToRefLineToSquare2() { }
    MappingToRefLineToSquare2::~MappingToRefLineToSquare2() { }

    // ~~~~~~~~~~~~~~~==============================================================================
    // ~~~ index 3 ~~~==============================================================================
    // ~~~~~~~~~~~~~~~==============================================================================

    const MappingToRefLineToSquare3& MappingToRefLineToSquare3::Instance()
    {
        static const MappingToRefLineToSquare3 theInstance;
        return theInstance;
    }

    void MappingToRefLineToSquare3::transform(const Geometry::PointReference<1>& p1,
                                               Geometry::PointReference<2>& p2) const
    {
        p2[0] = p1[0];
        p2[1] = 1.0;
    }

    void MappingToRefLineToSquare3::calcJacobian(const Geometry::PointReference<1>& p1,
                                                  Geometry::Jacobian<1,2>& jacobian) const
    {
        jacobian(0,0) = 1.0;
        jacobian(1,0) = 0.0;
    }

    MappingToRefLineToSquare3::MappingToRefLineToSquare3() { }
    MappingToRefLineToSquare3::~MappingToRefLineToSquare3() { }
}

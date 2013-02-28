/*
 * MappingToRefSquareToSquare.cpp
 *
 *  Created on: Feb 11, 2013
 *      Author: nicorivas
 */

#ifndef MAPPINGSQUARETOSQUARE_C_
#define MAPPINGSQUARETOSQUARE_C_

#include "MappingToRefSquareToSquare.hpp"
#include "../Jacobian.hpp"

namespace Geometry
{
    // ~~~~~~~~~~~~~~~==============================================================================
    // ~~~ index 0 ~~~==============================================================================
    // ~~~~~~~~~~~~~~~==============================================================================

    const MappingToRefSquareToSquare0& MappingToRefSquareToSquare0::Instance()
    {
        static const MappingToRefSquareToSquare0 theInstance;
        return theInstance;
    }

    void MappingToRefSquareToSquare0::transform(const Geometry::PointReference<2>& p1,
                                                 Geometry::PointReference<2>& p2) const
    {
        p2[0] = p1[0]; p2[1] = p1[1];
    }

    void MappingToRefSquareToSquare0::calcJacobian(const Geometry::PointReference<2>&,
                                                    Geometry::Jacobian<2,2>& jacobian) const
    {
        jacobian(0,0) = 1.0; jacobian(0,1) = 0.0;
        jacobian(1,0) = 0.0; jacobian(1,1) = 1.0;
    }

    MappingToRefSquareToSquare0::MappingToRefSquareToSquare0() { }
    MappingToRefSquareToSquare0::~MappingToRefSquareToSquare0() { }

    // ~~~~~~~~~~~~~~~==============================================================================
    // ~~~ index 1 ~~~==============================================================================
    // ~~~~~~~~~~~~~~~==============================================================================

    const MappingToRefSquareToSquare1& MappingToRefSquareToSquare1::Instance()
    {
        static const MappingToRefSquareToSquare1 theInstance;
        return theInstance;
    }

    void MappingToRefSquareToSquare1::transform(const Geometry::PointReference<2>& p1,
                                                 Geometry::PointReference<2>& p2) const
    {
        p2[0] = -p1[0]; p2[1] = p1[1];
    }

    void MappingToRefSquareToSquare1::calcJacobian(const Geometry::PointReference<2>&,
                                                    Geometry::Jacobian<2,2>& jacobian) const
    {
        jacobian(0,0) =  0.0; jacobian(0,1) = 1.0;
        jacobian(1,0) = -1.0; jacobian(1,1) = 0.0;
    }

    MappingToRefSquareToSquare1::MappingToRefSquareToSquare1() { }
    MappingToRefSquareToSquare1::~MappingToRefSquareToSquare1() { }

    // ~~~~~~~~~~~~~~~==============================================================================
    // ~~~ index 2 ~~~==============================================================================
    // ~~~~~~~~~~~~~~~==============================================================================

    const MappingToRefSquareToSquare2& MappingToRefSquareToSquare2::Instance()
    {
        static const MappingToRefSquareToSquare2 theInstance;
        return theInstance;
    }

    void MappingToRefSquareToSquare2::transform(const Geometry::PointReference<2>& p1,
                                                 Geometry::PointReference<2>& p2) const
    {
        p2[0] = -p1[0]; p2[1] = -p1[1];
    }

    void MappingToRefSquareToSquare2::calcJacobian(const Geometry::PointReference<2>&,
                                                    Geometry::Jacobian<2,2>& jacobian) const
    {
        jacobian(0,0) = -1.0; jacobian(0,1) =  0.0;
        jacobian(1,0) =  0.0; jacobian(1,1) = -1.0;
    }

    MappingToRefSquareToSquare2::MappingToRefSquareToSquare2() { }
    MappingToRefSquareToSquare2::~MappingToRefSquareToSquare2() { }

    // ~~~~~~~~~~~~~~~==============================================================================
    // ~~~ index 3 ~~~==============================================================================
    // ~~~~~~~~~~~~~~~==============================================================================

    const MappingToRefSquareToSquare3& MappingToRefSquareToSquare3::Instance()
    {
        static const MappingToRefSquareToSquare3 theInstance;
        return theInstance;
    }

    void MappingToRefSquareToSquare3::transform(const Geometry::PointReference<2>& p1,
                                                 Geometry::PointReference<2>& p2) const
    {
        p2[0] = p1[0]; p2[1] = -p1[1];
    }

    void MappingToRefSquareToSquare3::calcJacobian(const Geometry::PointReference<2>&,
                                                    Geometry::Jacobian<2,2>& jacobian) const
    {
        jacobian(0,0) = 0.0; jacobian(0,1) = -1.0;
        jacobian(1,0) = 1.0; jacobian(1,1) =  0.0;
    }

    MappingToRefSquareToSquare3::MappingToRefSquareToSquare3() { }
    MappingToRefSquareToSquare3::~MappingToRefSquareToSquare3() { }

    // ~~~~~~~~~~~~~~~==============================================================================
    // ~~~ index 4 ~~~==============================================================================
    // ~~~~~~~~~~~~~~~==============================================================================

    const MappingToRefSquareToSquare4& MappingToRefSquareToSquare4::Instance()
    {
        static const MappingToRefSquareToSquare4 theInstance;
        return theInstance;
    }

    void MappingToRefSquareToSquare4::transform(const Geometry::PointReference<2>& p1,
                                                 Geometry::PointReference<2>& p2) const
    {
        p2[0] = p1[0]; p2[1] = -p1[1];
    }

    void MappingToRefSquareToSquare4::calcJacobian(const Geometry::PointReference<2>&,
                                                    Geometry::Jacobian<2,2>& jacobian) const
    {
        jacobian(0,0) = 1.0; jacobian(0,1) =  0.0;
        jacobian(1,0) = 0.0; jacobian(1,1) = -1.0;
    }

    MappingToRefSquareToSquare4::MappingToRefSquareToSquare4() { }
    MappingToRefSquareToSquare4::~MappingToRefSquareToSquare4() { }

    // ~~~~~~~~~~~~~~~==============================================================================
    // ~~~ index 5 ~~~==============================================================================
    // ~~~~~~~~~~~~~~~==============================================================================

    const MappingToRefSquareToSquare5& MappingToRefSquareToSquare5::Instance()
    {
        static const MappingToRefSquareToSquare5 theInstance;
        return theInstance;
    }

    void MappingToRefSquareToSquare5::transform(const Geometry::PointReference<2>& p1,
                                                 Geometry::PointReference<2>& p2) const
    {
        p2[0] = -p1[0]; p2[1] = p1[1];
    }

    void MappingToRefSquareToSquare5::calcJacobian(const Geometry::PointReference<2>&,
                                                    Geometry::Jacobian<2,2>& jacobian) const
    {
        jacobian(0,0) = -1.0; jacobian(0,1) = 0.0;
        jacobian(1,0) =  0.0; jacobian(1,1) = 1.0;
    }

    MappingToRefSquareToSquare5::MappingToRefSquareToSquare5() { }
    MappingToRefSquareToSquare5::~MappingToRefSquareToSquare5() { }

    // ~~~~~~~~~~~~~~~==============================================================================
    // ~~~ index 6 ~~~==============================================================================
    // ~~~~~~~~~~~~~~~==============================================================================

    const MappingToRefSquareToSquare6& MappingToRefSquareToSquare6::Instance()
    {
        static const MappingToRefSquareToSquare6 theInstance;
        return theInstance;
    }

    void MappingToRefSquareToSquare6::transform(const Geometry::PointReference<2>& p1,
                                                 Geometry::PointReference<2>& p2) const
    {
        p2[0] = -p1[0]; p2[1] = -p1[1];
    }

    void MappingToRefSquareToSquare6::calcJacobian(const Geometry::PointReference<2>&,
                                                    Geometry::Jacobian<2,2>& jacobian) const
    {
        jacobian(0,0) =  1.0; jacobian(0,1) = -1.0;
        jacobian(1,0) = -1.0; jacobian(1,1) =  0.0;
    }

    MappingToRefSquareToSquare6::MappingToRefSquareToSquare6() { }
    MappingToRefSquareToSquare6::~MappingToRefSquareToSquare6() { }

    // ~~~~~~~~~~~~~~~==============================================================================
    // ~~~ index 7 ~~~==============================================================================
    // ~~~~~~~~~~~~~~~==============================================================================

    const MappingToRefSquareToSquare7& MappingToRefSquareToSquare7::Instance()
    {
        static const MappingToRefSquareToSquare7 theInstance;
        return theInstance;
    }

    void MappingToRefSquareToSquare7::transform(const Geometry::PointReference<2>& p1,
                                                 Geometry::PointReference<2>& p2) const
    {
        p2[0] = p1[0]; p2[1] = p1[1];
    }

    void MappingToRefSquareToSquare7::calcJacobian(const Geometry::PointReference<2>&,
                                                    Geometry::Jacobian<2,2>& jacobian) const
    {
        jacobian(0,0) = 0.0; jacobian(0,1) = 1.0;
        jacobian(1,0) = 1.0; jacobian(1,1) = 0.0;
    }

    MappingToRefSquareToSquare7::MappingToRefSquareToSquare7() { }
    MappingToRefSquareToSquare7::~MappingToRefSquareToSquare7() { }
};
#endif

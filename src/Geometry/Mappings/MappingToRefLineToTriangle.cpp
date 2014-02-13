/*
 * MappingToRefLineToTriangle.cpp
 *
 *  Created on: Feb 11, 2013
 *      Author: nicorivas
 */
#include "MappingToRefLineToTriangle.hpp"

namespace Geometry
{
    // ~~~~~~~~~~~~~~~==============================================================================
    // ~~~ index 0 ~~~==============================================================================
    // ~~~~~~~~~~~~~~~==============================================================================

    const MappingToRefLineToTriangle0& MappingToRefLineToTriangle0::Instance()
    {
        static const MappingToRefLineToTriangle0 theInstance;
        return theInstance;
    }

    void MappingToRefLineToTriangle0::transform(const Geometry::PointReference& p1,
                                                 Geometry::PointReference& p2) const
    {
        p2[0] = 0.5 * (p1[0] + 1.0);
        p2[1] = 0.0;
    }

    void MappingToRefLineToTriangle0::calcJacobian(const Geometry::PointReference& p1,
                                                    Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) = 0.5;
        jacobian(1,0) = 0.0;
    }

    MappingToRefLineToTriangle0::MappingToRefLineToTriangle0() { }
    MappingToRefLineToTriangle0::~MappingToRefLineToTriangle0() { }

    // ~~~~~~~~~~~~~~~==============================================================================
    // ~~~ index 1 ~~~==============================================================================
    // ~~~~~~~~~~~~~~~==============================================================================

    const MappingToRefLineToTriangle1& MappingToRefLineToTriangle1::Instance()
    {
        static const MappingToRefLineToTriangle1 theInstance;
        return theInstance;
    }

    void MappingToRefLineToTriangle1::transform(const Geometry::PointReference& p1,
                                               Geometry::PointReference& p2) const
    {
        p2[0] = 0.0;
        p2[1] = 0.5 * (p1[0] + 1.0);
    }

    void MappingToRefLineToTriangle1::calcJacobian(const Geometry::PointReference& p1,
                                                  Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) = 0.0;
        jacobian(1,0) = 0.5;
    }

    MappingToRefLineToTriangle1::MappingToRefLineToTriangle1() { }
    MappingToRefLineToTriangle1::~MappingToRefLineToTriangle1() { }

    // ~~~~~~~~~~~~~~~==============================================================================
    // ~~~ index 2 ~~~==============================================================================
    // ~~~~~~~~~~~~~~~==============================================================================

    const MappingToRefLineToTriangle2& MappingToRefLineToTriangle2::Instance()
    {
        static const MappingToRefLineToTriangle2 theInstance;
        return theInstance;
    }

    void MappingToRefLineToTriangle2::transform(const Geometry::PointReference& p1,
                                               Geometry::PointReference& p2) const
    {
        p2[0] = 0.5 * (-p1[0] + 1.0);
        p2[1] = 0.5 * ( p1[0] + 1.0);
    }

    void MappingToRefLineToTriangle2::calcJacobian(const Geometry::PointReference& p1,
                                                  Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) = -0.5;
        jacobian(1,0) =  0.5;
    }

    MappingToRefLineToTriangle2::MappingToRefLineToTriangle2() { }
    MappingToRefLineToTriangle2::~MappingToRefLineToTriangle2() { }
}

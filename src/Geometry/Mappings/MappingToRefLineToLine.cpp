/*
 * MappingToRefLineToLine.cpp
 *
 *  Created on: Feb 14, 2013
 *      Author: nicorivas
 */
#include "MappingToRefLineToLine.hpp"

namespace Geometry
{
    /*
     * The reference line:
     *
     * (-1) 0-------1 (+1)
     *
     * Linear maps of a line into itself. There are only two possible mappings:
     *
     *      index 0: x -> x
     *      index 1: x -> -x
     *
     */

    // ~~~~~~~~~~~~~~~==============================================================================
    // ~~~ index 0 ~~~==============================================================================
    // ~~~~~~~~~~~~~~~==============================================================================

    const MappingToRefLineToLine0& MappingToRefLineToLine0::Instance()
    {
        static const MappingToRefLineToLine0 theInstance;
        return theInstance;
    }

    void MappingToRefLineToLine0::transform(const Geometry::PointReference& p1,
                                             Geometry::PointReference& p2) const
    { p2[0] = p1[0]; }

    void MappingToRefLineToLine0::calcJacobian(const Geometry::PointReference& p1,
                                                Geometry::Jacobian& jacobian) const
    { jacobian(0,0) = 1.0; }

    MappingToRefLineToLine0::MappingToRefLineToLine0() { }
    MappingToRefLineToLine0::~MappingToRefLineToLine0() { }

    // ~~~~~~~~~~~~~~~==============================================================================
    // ~~~ index 1 ~~~==============================================================================
    // ~~~~~~~~~~~~~~~==============================================================================

    const MappingToRefLineToLine1& MappingToRefLineToLine1::Instance()
    {
        static const MappingToRefLineToLine1 theInstance;
        return theInstance;
    }

    void MappingToRefLineToLine1::transform(const Geometry::PointReference& p1,
                                             Geometry::PointReference& p2) const
    { p2[0] = -p1[0]; }

    void MappingToRefLineToLine1::calcJacobian(const Geometry::PointReference& p1,
                                                Geometry::Jacobian& jacobian) const
    { jacobian(0,0) = -1.0; }

    MappingToRefLineToLine1::MappingToRefLineToLine1() { }
    MappingToRefLineToLine1::~MappingToRefLineToLine1() { }

}

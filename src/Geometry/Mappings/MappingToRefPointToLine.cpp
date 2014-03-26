/*
 * MappingToRefPointToLine.cpp
 *
 *  Created on: Feb 14, 2013
 *      Author: nicorivas
 */
#include "MappingToRefPointToLine.hpp"

namespace Geometry
{
    /*
     * The reference line:
     *
     * (-1) 0-------1 (+1)
     *
     * Linear map a point into a line. There are only two possible mappings:
     *
     *      index 0: () -> -1.0
     *      index 1: () -> 1.0
     *
     * \todo I don't quite get this.
     *
     */

    // ~~~~~~~~~~~~~~~==============================================================================
    // ~~~ index 0 ~~~==============================================================================
    // ~~~~~~~~~~~~~~~==============================================================================

    const MappingToRefPointToLine0& MappingToRefPointToLine0::Instance()
    {
        static const MappingToRefPointToLine0 theInstance;
        return theInstance;
    }

    void MappingToRefPointToLine0::transform(const Geometry::PointReference&,
                                             Geometry::PointReference& p2) const
    { p2[0] = -1.0; }

    void MappingToRefPointToLine0::calcJacobian(const Geometry::PointReference&,
                                                Geometry::Jacobian&) const
    { }

    MappingToRefPointToLine0::MappingToRefPointToLine0() { }
    MappingToRefPointToLine0::~MappingToRefPointToLine0() { }

    // ~~~~~~~~~~~~~~~==============================================================================
    // ~~~ index 1 ~~~==============================================================================
    // ~~~~~~~~~~~~~~~==============================================================================

    const MappingToRefPointToLine1& MappingToRefPointToLine1::Instance()
    {
        static const MappingToRefPointToLine1 theInstance;
        return theInstance;
    }

    void MappingToRefPointToLine1::transform(const Geometry::PointReference&,
                                              Geometry::PointReference& p2) const
    { p2[0] = 1.0; }

    void MappingToRefPointToLine1::calcJacobian(const Geometry::PointReference&,
                                                 Geometry::Jacobian&) const
    { }

    MappingToRefPointToLine1::MappingToRefPointToLine1() { }
    MappingToRefPointToLine1::~MappingToRefPointToLine1() { }

}

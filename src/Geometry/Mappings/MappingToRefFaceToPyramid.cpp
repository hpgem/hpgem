/*
 * MappingToRefFaceToPyramid.cpp
 *
 *  Created on: Feb 17, 2013
 *      Author: nicorivas
 */
#include "MappingToRefFaceToPyramid.hpp"

namespace Geometry
{
    // ~~~ index 0 ~~~==============================================================================

    const MappingToRefFaceToPyramid0& MappingToRefFaceToPyramid0::Instance()
    {
        static const MappingToRefFaceToPyramid0 theInstance;
        return theInstance;
    }

    void MappingToRefFaceToPyramid0::transform(const Geometry::PointReference& p1,
                                                        Geometry::PointReference& p2) const
    {
        p2[0] =  p1[0];
        p2[1] = -p1[1];
        p2[2] =  0.0;
    }

    void MappingToRefFaceToPyramid0::calcJacobian(const Geometry::PointReference& p1,
                                                           Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) = 0.0;
        jacobian(1,0) = 1.0;
        jacobian(2,0) = 0.0;

        jacobian(0,1) = 1.0;
        jacobian(1,1) = 0.0;
        jacobian(2,1) = 0.0;
    }

    MappingToRefFaceToPyramid0::MappingToRefFaceToPyramid0() { }
    MappingToRefFaceToPyramid0::~MappingToRefFaceToPyramid0() { }

    // ~~~ index 1 ~~~==============================================================================

    const MappingToRefFaceToPyramid1& MappingToRefFaceToPyramid1::Instance()
    {
        static const MappingToRefFaceToPyramid1 theInstance;
        return theInstance;
    }

    void MappingToRefFaceToPyramid1::transform(const Geometry::PointReference& p1,
                                                        Geometry::PointReference& p2) const
    {
        p2[0] = -1.0 + p1[1];
        p2[1] = +1.0 - 2.0 * p1[0] - p1[1];
        p2[2] = p1[1];
    }

    void MappingToRefFaceToPyramid1::calcJacobian(const Geometry::PointReference& p1,
                                                           Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) =  0.0;
        jacobian(1,0) = -2.0;
        jacobian(2,0) =  0.0;

        jacobian(0,1) =  1.0;
        jacobian(1,1) = -1.0;
        jacobian(2,1) =  1.0;
    }

    MappingToRefFaceToPyramid1::MappingToRefFaceToPyramid1() { }
    MappingToRefFaceToPyramid1::~MappingToRefFaceToPyramid1() { }

    // ~~~ index 2 ~~~==============================================================================

    const MappingToRefFaceToPyramid2& MappingToRefFaceToPyramid2::Instance()
    {
        static const MappingToRefFaceToPyramid2 theInstance;
        return theInstance;
    }

    void MappingToRefFaceToPyramid2::transform(const Geometry::PointReference& p1,
                                                        Geometry::PointReference& p2) const
    {
        p2[0] =  1.0 - p1[1];
        p2[1] = -1.0 + 2.0 * p1[0] + p1[1];
        p2[2] =  p1[1];
    }

    void MappingToRefFaceToPyramid2::calcJacobian(const Geometry::PointReference& p1,
                                                           Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) = 0.0;
        jacobian(1,0) = 2.0;
        jacobian(2,0) = 0.0;

        jacobian(0,1) = -1.0;
        jacobian(1,1) =  1.0;
        jacobian(2,1) =  1.0;
    }

    MappingToRefFaceToPyramid2::MappingToRefFaceToPyramid2() { }
    MappingToRefFaceToPyramid2::~MappingToRefFaceToPyramid2() { }

    // ~~~ index 3 ~~~==============================================================================

    const MappingToRefFaceToPyramid3& MappingToRefFaceToPyramid3::Instance()
    {
        static const MappingToRefFaceToPyramid3 theInstance;
        return theInstance;
    }

    void MappingToRefFaceToPyramid3::transform(const Geometry::PointReference& p1,
                                                        Geometry::PointReference& p2) const
    {
        p2[0] = -1.0 + 2.0 * p1[0] + p1[1];
        p2[1] = -1.0 + p1[1];
        p2[2] = p1[1];
    }

    void MappingToRefFaceToPyramid3::calcJacobian(const Geometry::PointReference& p1,
                                                           Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) =  2.0;
        jacobian(1,0) =  0.0;
        jacobian(2,0) =  0.0;

        jacobian(0,1) =  1.0;
        jacobian(1,1) =  1.0;
        jacobian(2,1) =  1.0;
    }

    MappingToRefFaceToPyramid3::MappingToRefFaceToPyramid3() { }
    MappingToRefFaceToPyramid3::~MappingToRefFaceToPyramid3() { }

    // ~~~ index 4 ~~~==============================================================================

    const MappingToRefFaceToPyramid4& MappingToRefFaceToPyramid4::Instance()
    {
        static const MappingToRefFaceToPyramid4 theInstance;
        return theInstance;
    }

    void MappingToRefFaceToPyramid4::transform(const Geometry::PointReference& p1,
                                                Geometry::PointReference& p2) const
    {
        p2[0] = 1.0 - 2.0 * p1[0] - p1[1];
        p2[1] = 1.0 - p1[1];
        p2[2] = p1[1];
    }

    void MappingToRefFaceToPyramid4::calcJacobian(const Geometry::PointReference& p1,
                                                   Geometry::Jacobian& jacobian) const
    {
        jacobian(0,0) = -2.0;
        jacobian(1,0) =  0.0;
        jacobian(2,0) =  0.0;

        jacobian(0,1) = -1.0;
        jacobian(1,1) = -1.0;
        jacobian(2,1) =  1.0;
    }

    MappingToRefFaceToPyramid4::MappingToRefFaceToPyramid4() { }
    MappingToRefFaceToPyramid4::~MappingToRefFaceToPyramid4() { }

}

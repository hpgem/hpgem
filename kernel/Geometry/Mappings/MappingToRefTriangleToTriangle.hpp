/*
 * MappingToRefTriangleToTriangle.hpp
 *
 *  Created on: Feb 17, 2013
 *      Author: nicorivas
 */

#ifndef MappingToRefTriangleToTriangle_H_
#define MappingToRefTriangleToTriangle_H_

#include "MappingReferenceToReference.hpp"
#include "../Jacobian.hpp"

namespace Geometry
{
    /*
     * The ordering of the vertex and faces in a triangle:
     *
     *   (0,1) 2
     *         | \
     *         1   2
     *         |     \
     *   (0,0) 0---0---1 (1,0)
     *
     *
     * This implements the linear mappings of a ReferenceTriangle [0,1]^2 onto itself.
     * There are 6 possible mappings.
     *
     * TODO: Proper mapping documentation missing here
     *
     *      index 0: (x,y)->(x,y)     (Identity)
     *      index 1: (x,y)->(y,x)     ((021) Rotation with respect to x = y)
     *      index 2: (x,y)->(1-x-y,x) ((120) Counter-clockwise rotation)
     *      index 3: (x,y)->(1-x-y,y) ((102) Something weird)
     *      index 4: (x,y)->(x,1-x-y) ((210) Something weird)
     *      index 5: (x,y)->(y,1-x-y) ((201) Clockwise rotation)
     */

    // ~~~ index 0 ~~~==============================================================================
    class MappingToRefTriangleToTriangle0: public MappingReferenceToReference
    {
        public:
            static const MappingToRefTriangleToTriangle0& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 2;}
        private:
            MappingToRefTriangleToTriangle0();
            MappingToRefTriangleToTriangle0(const MappingToRefTriangleToTriangle0&);
            MappingToRefTriangleToTriangle0& operator=(const MappingToRefTriangleToTriangle0&);
            virtual ~MappingToRefTriangleToTriangle0();
    };
    // ~~~ index 1 ~~~==============================================================================
    class MappingToRefTriangleToTriangle1: public MappingReferenceToReference
    {
        public:
            static const MappingToRefTriangleToTriangle1& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 2;}
        private:
            MappingToRefTriangleToTriangle1();
            MappingToRefTriangleToTriangle1(const MappingToRefTriangleToTriangle1&);
            MappingToRefTriangleToTriangle1& operator=(const MappingToRefTriangleToTriangle1&);
            virtual ~MappingToRefTriangleToTriangle1();
    };
    // ~~~ index 2 ~~~==============================================================================
    class MappingToRefTriangleToTriangle2: public MappingReferenceToReference
    {
        public:
            static const MappingToRefTriangleToTriangle2& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 2;}
        private:
            MappingToRefTriangleToTriangle2();
            MappingToRefTriangleToTriangle2(const MappingToRefTriangleToTriangle2&);
            MappingToRefTriangleToTriangle2& operator=(const MappingToRefTriangleToTriangle2&);
            virtual ~MappingToRefTriangleToTriangle2();
    };
    // ~~~ index 3 ~~~==============================================================================
    class MappingToRefTriangleToTriangle3: public MappingReferenceToReference
    {
        public:
            static const MappingToRefTriangleToTriangle3& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 2;}
        private:
            MappingToRefTriangleToTriangle3();
            MappingToRefTriangleToTriangle3(const MappingToRefTriangleToTriangle3&);
            MappingToRefTriangleToTriangle3& operator=(const MappingToRefTriangleToTriangle3&);
            virtual ~MappingToRefTriangleToTriangle3();
    };
    // ~~~ index 4 ~~~==============================================================================
    class MappingToRefTriangleToTriangle4: public MappingReferenceToReference
    {
        public:
            static const MappingToRefTriangleToTriangle4& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 2;}
        private:
            MappingToRefTriangleToTriangle4();
            MappingToRefTriangleToTriangle4(const MappingToRefTriangleToTriangle4&);
            MappingToRefTriangleToTriangle4& operator=(const MappingToRefTriangleToTriangle4&);
            virtual ~MappingToRefTriangleToTriangle4();
    };
    // ~~~ index 5 ~~~==============================================================================
    class MappingToRefTriangleToTriangle5: public MappingReferenceToReference
    {
        public:
            static const MappingToRefTriangleToTriangle5& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 2;}
        private:
            MappingToRefTriangleToTriangle5();
            MappingToRefTriangleToTriangle5(const MappingToRefTriangleToTriangle5&);
            MappingToRefTriangleToTriangle5& operator=(const MappingToRefTriangleToTriangle5&);
            virtual ~MappingToRefTriangleToTriangle5();
    };
};
#endif

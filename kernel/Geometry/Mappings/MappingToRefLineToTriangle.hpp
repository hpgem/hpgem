/*
 * MappingToRefLineToTriangle.hpp
 *
 *  Created on: Feb 11, 2013
 *      Author: nicorivas
 */

#ifndef MappingToRefLineToTriangle_H_
#define MappingToRefLineToTriangle_H_

#include "MappingReferenceToReference.hpp"
#include "../Jacobian.hpp"

namespace Geometry
{
    /* The ordering of the vertex and faces in a triangle:
     *
     *   (0,1) 2
     *         | \
     *         1   2
     *         |     \
     *   (0,0) 0---0---1 (1,0)
     *
     *
     * This maps the reference line [-1,1] to the triangle shown above. The mappings are defined as
     * follows:
     *
     *      faceindex 0: x -> ((1+x)/2,0)
     *      faceindex 1: x -> (0,(1+x)/2)
     *      faceindex 2: x -> ((1-x)/2,(1+x)/2)
     *
     */

    // ~~~ index 0 ~~~==============================================================================
    class MappingToRefLineToTriangle0: public MappingReferenceToReference
    {
        public:
            static const MappingToRefLineToTriangle0& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 2;}
        private:
            MappingToRefLineToTriangle0();
            MappingToRefLineToTriangle0(const MappingToRefLineToTriangle0&);
            MappingToRefLineToTriangle0& operator=(const MappingToRefLineToTriangle0&);
            virtual ~MappingToRefLineToTriangle0();
    };
    // ~~~ index 1 ~~~==============================================================================
    class MappingToRefLineToTriangle1: public MappingReferenceToReference
    {
        public:
            static const MappingToRefLineToTriangle1& Instance();
            virtual void transform(const Geometry::PointReference&,
                                         Geometry::PointReference&) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 2;}
        private:
            MappingToRefLineToTriangle1();
            MappingToRefLineToTriangle1(const MappingToRefLineToTriangle1&);
            MappingToRefLineToTriangle1& operator=(const MappingToRefLineToTriangle1&);
            virtual ~MappingToRefLineToTriangle1();
    };
    // ~~~ index 2 ~~~==============================================================================
    class MappingToRefLineToTriangle2: public MappingReferenceToReference
    {
        public:
            static const MappingToRefLineToTriangle2& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 2;}
        private:
            MappingToRefLineToTriangle2();
            MappingToRefLineToTriangle2(const MappingToRefLineToTriangle2&);
            MappingToRefLineToTriangle1& operator=(const MappingToRefLineToTriangle2&);
            virtual ~MappingToRefLineToTriangle2();
    };
};
#endif /* MAPPINGSIMPLECUBENLINEAR_H_ */

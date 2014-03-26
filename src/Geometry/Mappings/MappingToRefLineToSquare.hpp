/*
 * MappingToRefLineToSquare.hpp
 *
 *  Created on: Feb 11, 2013
 *      Author: nicorivas
 */

#ifndef MAPPINGLINETOSQUARE_H_
#define MAPPINGLINETOSQUARE_H_

#include "MappingReferenceToReference.hpp"
#include "../Jacobian.hpp"

namespace Geometry
{
    /*
     * The ordering of the vertex and faces in a square:
     *
     * (-1,+1) 2---3---3 (+1,+1)
     *         |       |
     *         1       2
     *         |       |
     * (-1,-1) 0---0---1 (1,-1)
     *
     * This maps the reference line [-1,1] to the square shown above. The mappings are defined as
     * follows:
     *
     *      faceindex 0: x -> (x,-1)
     *      faceindex 1: x -> (-1,x)
     *      faceindex 2: x -> (1,x)
     *      faceindex 3: x -> (x,1)
     *
     * The mapping can in principle be defined freely, but to simplify finding the right
     * RefFace2RefFaceMapping they are chosen to preserve the ordering of the vertices (Ordering by
     * coordinate). This makes it possible to use the FaceNodeList from the Reference elements to
     * find the RF2RFMapping. This will also eliminate the need for mirror type of RF2RFMappings,
     * since these will not occur when the ordering is preserved (I think).
     */

    // ~~~ index 0 ~~~==============================================================================
    class MappingToRefLineToSquare0: public MappingReferenceToReference
    {
        public:
            static const MappingToRefLineToSquare0& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 2;}
        private:
            MappingToRefLineToSquare0();
            MappingToRefLineToSquare0(const MappingToRefLineToSquare0&);
            MappingToRefLineToSquare0& operator=(const MappingToRefLineToSquare0&);
            virtual ~MappingToRefLineToSquare0();
    };
    // ~~~ index 1 ~~~==============================================================================
    class MappingToRefLineToSquare1: public MappingReferenceToReference
    {
        public:
            static const MappingToRefLineToSquare1& Instance();
            virtual void transform(const Geometry::PointReference&,
                                         Geometry::PointReference&) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 2;}
        private:
            MappingToRefLineToSquare1();
            MappingToRefLineToSquare1(const MappingToRefLineToSquare1&);
            MappingToRefLineToSquare1& operator=(const MappingToRefLineToSquare1&);
            virtual ~MappingToRefLineToSquare1();
    };
    // ~~~ index 2 ~~~==============================================================================
    class MappingToRefLineToSquare2: public MappingReferenceToReference
    {
        public:
            static const MappingToRefLineToSquare2& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 2;}
        private:
            MappingToRefLineToSquare2();
            MappingToRefLineToSquare2(const MappingToRefLineToSquare2&);
            MappingToRefLineToSquare1& operator=(const MappingToRefLineToSquare2&);
            virtual ~MappingToRefLineToSquare2();
    };
    // ~~~ index 3 ~~~==============================================================================
    class MappingToRefLineToSquare3: public MappingReferenceToReference
    {
        public:
            static const MappingToRefLineToSquare3& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 2;}
        private:
            MappingToRefLineToSquare3();
            MappingToRefLineToSquare3(const MappingToRefLineToSquare3&);
            MappingToRefLineToSquare3& operator=(const MappingToRefLineToSquare3&);
            virtual ~MappingToRefLineToSquare3();
    };
};
#endif /* MAPPINGSIMPLECUBENLINEAR_H_ */

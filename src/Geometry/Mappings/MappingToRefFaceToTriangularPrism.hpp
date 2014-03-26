/*
 * MappingToRefFaceToTriangularPrism.hpp
 *
 *  Created on: Feb 19, 2013
 *      Author: nicorivas
 */

#ifndef MappingToRefFaceToTriangularPrism_H_
#define MappingToRefFaceToTriangularPrism_H_

#include "MappingReferenceToReference.hpp"
#include "../Jacobian.hpp"

namespace Geometry
{
    /*
     *
     */

    // ~~~ index 0 ~~~==============================================================================

    class MappingToRefFaceToTriangularPrism0: public MappingReferenceToReference
    {
        public:
            static const MappingToRefFaceToTriangularPrism0& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 3;}
        private:
            MappingToRefFaceToTriangularPrism0();
            MappingToRefFaceToTriangularPrism0(const MappingToRefFaceToTriangularPrism0&);
            MappingToRefFaceToTriangularPrism0& operator=(const MappingToRefFaceToTriangularPrism0&);
            virtual ~MappingToRefFaceToTriangularPrism0();
    };

    // ~~~ index 1 ~~~==============================================================================

    class MappingToRefFaceToTriangularPrism1: public MappingReferenceToReference
    {
        public:
            static const MappingToRefFaceToTriangularPrism1& Instance();
            virtual void transform(const Geometry::PointReference&,
                                         Geometry::PointReference&) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 3;}
        private:
            MappingToRefFaceToTriangularPrism1();
            MappingToRefFaceToTriangularPrism1(const MappingToRefFaceToTriangularPrism1&);
            MappingToRefFaceToTriangularPrism1& operator=(const MappingToRefFaceToTriangularPrism1&);
            virtual ~MappingToRefFaceToTriangularPrism1();
    };

    // ~~~ index 2 ~~~==============================================================================

    class MappingToRefFaceToTriangularPrism2: public MappingReferenceToReference
    {
        public:
            static const MappingToRefFaceToTriangularPrism2& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 3;}
        private:
            MappingToRefFaceToTriangularPrism2();
            MappingToRefFaceToTriangularPrism2(const MappingToRefFaceToTriangularPrism2&);
            MappingToRefFaceToTriangularPrism1& operator=(const MappingToRefFaceToTriangularPrism2&);
            virtual ~MappingToRefFaceToTriangularPrism2();
    };

    // ~~~ index 3 ~~~==============================================================================

    class MappingToRefFaceToTriangularPrism3: public MappingReferenceToReference
    {
        public:
            static const MappingToRefFaceToTriangularPrism3& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 3;}
        private:
            MappingToRefFaceToTriangularPrism3();
            MappingToRefFaceToTriangularPrism3(const MappingToRefFaceToTriangularPrism3&);
            MappingToRefFaceToTriangularPrism3& operator=(const MappingToRefFaceToTriangularPrism3&);
            virtual ~MappingToRefFaceToTriangularPrism3();
    };

    // ~~~ index 4 ~~~==============================================================================

    class MappingToRefFaceToTriangularPrism4: public MappingReferenceToReference
    {
        public:
            static const MappingToRefFaceToTriangularPrism4& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 3;}
        private:
            MappingToRefFaceToTriangularPrism4();
            MappingToRefFaceToTriangularPrism4(const MappingToRefFaceToTriangularPrism4&);
            MappingToRefFaceToTriangularPrism4& operator=(const MappingToRefFaceToTriangularPrism4&);
            virtual ~MappingToRefFaceToTriangularPrism4();
    };

};
#endif

/*
 * MappingToRefLineToLine.hpp
 *
 *  Created on: Feb 11, 2013
 *      Author: nicorivas
 */

#ifndef MAPPINGLINETOLINE_H_
#define MAPPINGLINETOLINE_H_

#include "MappingReferenceToReference.hpp"
#include "../Jacobian.hpp"

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

    // ~~~ index 0 ~~~=========================================================================== //
    class MappingToRefLineToLine0: public MappingReferenceToReference
    {
        public:
            static const MappingToRefLineToLine0& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 1;}
        private:
            MappingToRefLineToLine0();
            MappingToRefLineToLine0(const MappingToRefLineToLine0&);
            MappingToRefLineToLine0& operator=(const MappingToRefLineToLine0&);
            virtual ~MappingToRefLineToLine0();
    };

    // ~~~ index 1 ~~~=========================================================================== //
    class MappingToRefLineToLine1: public MappingReferenceToReference
    {
        public:
            static const MappingToRefLineToLine1& Instance();
            virtual void transform(const Geometry::PointReference& p1,
                                         Geometry::PointReference& p2) const;
            virtual void calcJacobian(const Geometry::PointReference&,
                                            Geometry::Jacobian&) const;
            virtual int getTargetDimension() const {return 1;}
        private:
            MappingToRefLineToLine1();
            MappingToRefLineToLine1(const MappingToRefLineToLine1&);
            MappingToRefLineToLine1& operator=(const MappingToRefLineToLine1&);
            virtual ~MappingToRefLineToLine1();
    };
};
#endif

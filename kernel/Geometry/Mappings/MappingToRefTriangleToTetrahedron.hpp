/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found below.


 Copyright (c) 2014, University of Twente
 All rights reserved.

 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef MappingToRefTriangleToTetrahedron_H_
#define MappingToRefTriangleToTetrahedron_H_

#include "MappingReferenceToReference.hpp"

namespace Geometry
{
    /* The ordering of the vertex in a cube:
     *
     *  (0,0,1) 3
     *          |\
     *          |  \2 (0,1,0)
     *          |  / \
     *          |/     \
     *  (0,0,0) 0--------1 (1,0,0)
     *
     *  faces indexes:
     *              0: (0,3,2)
     *              1: (0,1,3)
     *              2: (0,2,1)
     *              3: (1,2,3)
     */

    // ~~~ index 0 ~~~==============================================================================

    class MappingToRefTriangleToTetrahedron0: public MappingReferenceToReference
    {
        public:
            static const MappingToRefTriangleToTetrahedron0& Instance();
            virtual PointReference transform(const Geometry::PointReference& p1) const;
            virtual Jacobian calcJacobian(const Geometry::PointReference&) const;
            virtual std::size_t getTargetDimension() const {return 3;}
        private:
            MappingToRefTriangleToTetrahedron0();
            MappingToRefTriangleToTetrahedron0(const MappingToRefTriangleToTetrahedron0&);
            MappingToRefTriangleToTetrahedron0& operator=(const MappingToRefTriangleToTetrahedron0&);
            virtual ~MappingToRefTriangleToTetrahedron0();
    };

    // ~~~ index 1 ~~~==============================================================================

    class MappingToRefTriangleToTetrahedron1: public MappingReferenceToReference
    {
        public:
            static const MappingToRefTriangleToTetrahedron1& Instance();
            virtual PointReference transform(const Geometry::PointReference& p1) const;
            virtual Jacobian calcJacobian(const Geometry::PointReference&) const;
            virtual std::size_t getTargetDimension() const {return 3;}
        private:
            MappingToRefTriangleToTetrahedron1();
            MappingToRefTriangleToTetrahedron1(const MappingToRefTriangleToTetrahedron1&);
            MappingToRefTriangleToTetrahedron1& operator=(const MappingToRefTriangleToTetrahedron1&);
            virtual ~MappingToRefTriangleToTetrahedron1();
    };

    // ~~~ index 2 ~~~==============================================================================

    class MappingToRefTriangleToTetrahedron2: public MappingReferenceToReference
    {
        public:
            static const MappingToRefTriangleToTetrahedron2& Instance();
            virtual PointReference transform(const Geometry::PointReference& p1) const;
            virtual Jacobian calcJacobian(const Geometry::PointReference&) const;
            virtual std::size_t getTargetDimension() const {return 3;}
        private:
            MappingToRefTriangleToTetrahedron2();
            MappingToRefTriangleToTetrahedron2(const MappingToRefTriangleToTetrahedron2&);
            MappingToRefTriangleToTetrahedron1& operator=(const MappingToRefTriangleToTetrahedron2&);
            virtual ~MappingToRefTriangleToTetrahedron2();
    };

    // ~~~ index 3 ~~~==============================================================================

    class MappingToRefTriangleToTetrahedron3: public MappingReferenceToReference
    {
        public:
            static const MappingToRefTriangleToTetrahedron3& Instance();
            virtual PointReference transform(const Geometry::PointReference& p1) const;
            virtual Jacobian calcJacobian(const Geometry::PointReference&) const;
            virtual std::size_t getTargetDimension() const {return 3;}
        private:
            MappingToRefTriangleToTetrahedron3();
            MappingToRefTriangleToTetrahedron3(const MappingToRefTriangleToTetrahedron3&);
            MappingToRefTriangleToTetrahedron3& operator=(const MappingToRefTriangleToTetrahedron3&);
            virtual ~MappingToRefTriangleToTetrahedron3();
    };
};
#endif

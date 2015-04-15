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

#ifndef MappingToRefLineToTriangle_H_
#define MappingToRefLineToTriangle_H_

#include "MappingReferenceToReference.h"

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
    class MappingToRefLineToTriangle0 : public MappingReferenceToReference
    {
    public:
        static const MappingToRefLineToTriangle0& Instance();
        PointReference transform(const Geometry::PointReference& p1) const override final;
        Jacobian calcJacobian(const Geometry::PointReference&) const override final;
        std::size_t getTargetDimension() const override final
        {
            return 2;
        }
        MappingToRefLineToTriangle0(const MappingToRefLineToTriangle0&) = delete;
        MappingToRefLineToTriangle0& operator=(const MappingToRefLineToTriangle0&) = delete;
    private:
        MappingToRefLineToTriangle0();
    };
    // ~~~ index 1 ~~~==============================================================================
    class MappingToRefLineToTriangle1 : public MappingReferenceToReference
    {
    public:
        static const MappingToRefLineToTriangle1& Instance();
        PointReference transform(const Geometry::PointReference& p1) const override final;
        Jacobian calcJacobian(const Geometry::PointReference&) const override final;
        std::size_t getTargetDimension() const override final
        {
            return 2;
        }
        MappingToRefLineToTriangle1(const MappingToRefLineToTriangle1&) = delete;
        MappingToRefLineToTriangle1& operator=(const MappingToRefLineToTriangle1&) = delete;
    private:
        MappingToRefLineToTriangle1();
    };
    // ~~~ index 2 ~~~==============================================================================
    class MappingToRefLineToTriangle2 : public MappingReferenceToReference
    {
    public:
        static const MappingToRefLineToTriangle2& Instance();
        PointReference transform(const Geometry::PointReference& p1) const override final;
        Jacobian calcJacobian(const Geometry::PointReference&) const override final;
        std::size_t getTargetDimension() const override final
        {
            return 2;
        }
        MappingToRefLineToTriangle2(const MappingToRefLineToTriangle2&) = delete;
        MappingToRefLineToTriangle1& operator=(const MappingToRefLineToTriangle2&) = delete;
    private:
        MappingToRefLineToTriangle2();
    };
}
#endif /* MAPPINGSIMPLECUBENLINEAR_H_ */

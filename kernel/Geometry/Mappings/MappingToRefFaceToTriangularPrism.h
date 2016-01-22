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

#ifndef MappingToRefFaceToTriangularPrism_H_
#define MappingToRefFaceToTriangularPrism_H_

#include "MappingReferenceToReference.h"

namespace Geometry
{

    // ~~~ index 0 ~~~==============================================================================
    class MappingToRefFaceToTriangularPrism0 : public MappingReferenceToReference<1>
    {
    public:
        static const MappingToRefFaceToTriangularPrism0& Instance();
        PointReference<3> transform(const Geometry::PointReference<2>& p1) const override final;
        Jacobian<2, 3> calcJacobian(const Geometry::PointReference<2>&) const override final;
        std::size_t getTargetDimension() const override final
        {
            return 3;
        }
        MappingToRefFaceToTriangularPrism0(const MappingToRefFaceToTriangularPrism0&) = delete;
        MappingToRefFaceToTriangularPrism0& operator=(const MappingToRefFaceToTriangularPrism0&) = delete;
    private:
        MappingToRefFaceToTriangularPrism0();
        std::map<const PointReference<2>*, const PointReference<3>*> transformedCoordinates;
    };
    
    // ~~~ index 1 ~~~==============================================================================
    
    class MappingToRefFaceToTriangularPrism1 : public MappingReferenceToReference<1>
    {
    public:
        static const MappingToRefFaceToTriangularPrism1& Instance();
        PointReference<3> transform(const Geometry::PointReference<2>& p1) const override final;
        Jacobian<2, 3> calcJacobian(const Geometry::PointReference<2>&) const override final;
        std::size_t getTargetDimension() const override final
        {
            return 3;
        }
        MappingToRefFaceToTriangularPrism1(const MappingToRefFaceToTriangularPrism1&) = delete;
        MappingToRefFaceToTriangularPrism1& operator=(const MappingToRefFaceToTriangularPrism1&) = delete;
    private:
        MappingToRefFaceToTriangularPrism1();
        std::map<const PointReference<2>*, const PointReference<3>*> transformedCoordinates;
    };
    
    // ~~~ index 2 ~~~==============================================================================
    
    class MappingToRefFaceToTriangularPrism2 : public MappingReferenceToReference<1>
    {
    public:
        static const MappingToRefFaceToTriangularPrism2& Instance();
        PointReference<3> transform(const Geometry::PointReference<2>& p1) const override final;
        Jacobian<2, 3> calcJacobian(const Geometry::PointReference<2>&) const override final;
        std::size_t getTargetDimension() const override final
        {
            return 3;
        }
        MappingToRefFaceToTriangularPrism2(const MappingToRefFaceToTriangularPrism2&) = delete;
        MappingToRefFaceToTriangularPrism1& operator=(const MappingToRefFaceToTriangularPrism2&) = delete;
    private:
        MappingToRefFaceToTriangularPrism2();
        std::map<const PointReference<2>*, const PointReference<3>*> transformedCoordinates;
    };
    
    // ~~~ index 3 ~~~==============================================================================
    
    class MappingToRefFaceToTriangularPrism3 : public MappingReferenceToReference<1>
    {
    public:
        static const MappingToRefFaceToTriangularPrism3& Instance();
        PointReference<3> transform(const Geometry::PointReference<2>& p1) const override final;
        Jacobian<2, 3> calcJacobian(const Geometry::PointReference<2>&) const override final;
        std::size_t getTargetDimension() const override final
        {
            return 3;
        }
        MappingToRefFaceToTriangularPrism3(const MappingToRefFaceToTriangularPrism3&) = delete;
        MappingToRefFaceToTriangularPrism3& operator=(const MappingToRefFaceToTriangularPrism3&) = delete;
    private:
        MappingToRefFaceToTriangularPrism3();
        std::map<const PointReference<2>*, const PointReference<3>*> transformedCoordinates;
    };
    
    // ~~~ index 4 ~~~==============================================================================
    
    class MappingToRefFaceToTriangularPrism4 : public MappingReferenceToReference<1>
    {
    public:
        static const MappingToRefFaceToTriangularPrism4& Instance();
        PointReference<3> transform(const Geometry::PointReference<2>& p1) const override final;
        Jacobian<2, 3> calcJacobian(const Geometry::PointReference<2>&) const override final;
        std::size_t getTargetDimension() const override final
        {
            return 3;
        }
        MappingToRefFaceToTriangularPrism4(const MappingToRefFaceToTriangularPrism4&) = delete;
        MappingToRefFaceToTriangularPrism4& operator=(const MappingToRefFaceToTriangularPrism4&) = delete;
    private:
        MappingToRefFaceToTriangularPrism4();
        std::map<const PointReference<2>*, const PointReference<3>*> transformedCoordinates;
    };

}
#endif

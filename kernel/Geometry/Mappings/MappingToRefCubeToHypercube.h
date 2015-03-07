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

#ifndef MappingToRefCubeToHypercube_H_
#define MappingToRefCubeToHypercube_H_

#include "MappingReferenceToReference.h"

namespace Geometry
{
    /*
     *
     */

    // ~~~ index 0 ~~~==============================================================================
    class MappingToRefCubeToHypercube0 : public MappingReferenceToReference
    {
    public:
        static const MappingToRefCubeToHypercube0& Instance();
        virtual PointReference transform(const Geometry::PointReference& p1) const;
        virtual Jacobian calcJacobian(const Geometry::PointReference&) const;
        virtual std::size_t getTargetDimension() const
        {
            return 4;
        }
    private:
        MappingToRefCubeToHypercube0();
        MappingToRefCubeToHypercube0(const MappingToRefCubeToHypercube0&);
        MappingToRefCubeToHypercube0& operator=(const MappingToRefCubeToHypercube0&);
        virtual ~MappingToRefCubeToHypercube0();
    };
    
    // ~~~ index 1 ~~~==============================================================================
    
    class MappingToRefCubeToHypercube1 : public MappingReferenceToReference
    {
    public:
        static const MappingToRefCubeToHypercube1& Instance();
        virtual PointReference transform(const Geometry::PointReference& p1) const;
        virtual Jacobian calcJacobian(const Geometry::PointReference&) const;
        virtual std::size_t getTargetDimension() const
        {
            return 4;
        }
    private:
        MappingToRefCubeToHypercube1();
        MappingToRefCubeToHypercube1(const MappingToRefCubeToHypercube1&);
        MappingToRefCubeToHypercube1& operator=(const MappingToRefCubeToHypercube1&);
        virtual ~MappingToRefCubeToHypercube1();
    };
    
    // ~~~ index 2 ~~~==============================================================================
    
    class MappingToRefCubeToHypercube2 : public MappingReferenceToReference
    {
    public:
        static const MappingToRefCubeToHypercube2& Instance();
        virtual PointReference transform(const Geometry::PointReference& p1) const;
        virtual Jacobian calcJacobian(const Geometry::PointReference&) const;
        virtual std::size_t getTargetDimension() const
        {
            return 4;
        }
    private:
        MappingToRefCubeToHypercube2();
        MappingToRefCubeToHypercube2(const MappingToRefCubeToHypercube2&);
        MappingToRefCubeToHypercube1& operator=(const MappingToRefCubeToHypercube2&);
        virtual ~MappingToRefCubeToHypercube2();
    };
    
    // ~~~ index 3 ~~~==============================================================================
    
    class MappingToRefCubeToHypercube3 : public MappingReferenceToReference
    {
    public:
        static const MappingToRefCubeToHypercube3& Instance();
        virtual PointReference transform(const Geometry::PointReference& p1) const;
        virtual Jacobian calcJacobian(const Geometry::PointReference&) const;
        virtual std::size_t getTargetDimension() const
        {
            return 4;
        }
    private:
        MappingToRefCubeToHypercube3();
        MappingToRefCubeToHypercube3(const MappingToRefCubeToHypercube3&);
        MappingToRefCubeToHypercube3& operator=(const MappingToRefCubeToHypercube3&);
        virtual ~MappingToRefCubeToHypercube3();
    };
    
    // ~~~ index 4 ~~~==============================================================================
    
    class MappingToRefCubeToHypercube4 : public MappingReferenceToReference
    {
    public:
        static const MappingToRefCubeToHypercube4& Instance();
        virtual PointReference transform(const Geometry::PointReference& p1) const;
        virtual Jacobian calcJacobian(const Geometry::PointReference&) const;
        virtual std::size_t getTargetDimension() const
        {
            return 4;
        }
    private:
        MappingToRefCubeToHypercube4();
        MappingToRefCubeToHypercube4(const MappingToRefCubeToHypercube4&);
        MappingToRefCubeToHypercube4& operator=(const MappingToRefCubeToHypercube4&);
        virtual ~MappingToRefCubeToHypercube4();
    };
    
    // ~~~ index 5 ~~~==============================================================================
    
    class MappingToRefCubeToHypercube5 : public MappingReferenceToReference
    {
    public:
        static const MappingToRefCubeToHypercube5& Instance();
        virtual PointReference transform(const Geometry::PointReference& p1) const;
        virtual Jacobian calcJacobian(const Geometry::PointReference&) const;
        virtual std::size_t getTargetDimension() const
        {
            return 4;
        }
    private:
        MappingToRefCubeToHypercube5();
        MappingToRefCubeToHypercube5(const MappingToRefCubeToHypercube5&);
        MappingToRefCubeToHypercube5& operator=(const MappingToRefCubeToHypercube5&);
        virtual ~MappingToRefCubeToHypercube5();
    };
    
    // ~~~ index 6 ~~~==============================================================================
    
    class MappingToRefCubeToHypercube6 : public MappingReferenceToReference
    {
    public:
        static const MappingToRefCubeToHypercube6& Instance();
        virtual PointReference transform(const Geometry::PointReference& p1) const;
        virtual Jacobian calcJacobian(const Geometry::PointReference&) const;
        virtual std::size_t getTargetDimension() const
        {
            return 4;
        }
    private:
        MappingToRefCubeToHypercube6();
        MappingToRefCubeToHypercube6(const MappingToRefCubeToHypercube6&);
        MappingToRefCubeToHypercube6& operator=(const MappingToRefCubeToHypercube6&);
        virtual ~MappingToRefCubeToHypercube6();
    };
    
    // ~~~ index 7 ~~~==============================================================================
    
    class MappingToRefCubeToHypercube7 : public MappingReferenceToReference
    {
    public:
        static const MappingToRefCubeToHypercube7& Instance();
        virtual PointReference transform(const Geometry::PointReference& p1) const;
        virtual Jacobian calcJacobian(const Geometry::PointReference&) const;
        virtual std::size_t getTargetDimension() const
        {
            return 4;
        }
    private:
        MappingToRefCubeToHypercube7();
        MappingToRefCubeToHypercube7(const MappingToRefCubeToHypercube7&);
        MappingToRefCubeToHypercube7& operator=(const MappingToRefCubeToHypercube7&);
        virtual ~MappingToRefCubeToHypercube7();
    };
}
;
#endif

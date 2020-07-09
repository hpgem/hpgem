/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2014, University of Twente
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors
 may be used to endorse or promote products derived from this software without
 specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef HPGEM_KERNEL_POINTPHYSICALBASE_H
#define HPGEM_KERNEL_POINTPHYSICALBASE_H

#include <cstdlib>

namespace Geometry {
template <std::size_t DIM>
class PointPhysical;

/// \brief evades a compile-time conflict between the virtual and template
/// keywords \details for technical reasons c++ does not support abstract
/// templated functions. For some functions this means we get into trouble with
/// the return type. By defining a base class that can be converted back into
/// every possible templated PointReference we can return the base class
/// instead. Be warned that this reduces type safety. (In Debug mode, there are
/// run-time checks to make sure everything in done correctly)
// first define the base class, then define the subclass, then implement the
// type conversion operators
class PointPhysicalBase {
   public:
    operator PointPhysical<0> &();

    operator PointPhysical<1> &();

    operator PointPhysical<2> &();

    operator PointPhysical<3> &();

    operator PointPhysical<4> &();

    operator const PointPhysical<0> &() const;

    operator const PointPhysical<1> &() const;

    operator const PointPhysical<2> &() const;

    operator const PointPhysical<3> &() const;

    operator const PointPhysical<4> &() const;

    virtual std::size_t size() const = 0;

   protected:
    PointPhysicalBase() = default;
    virtual ~PointPhysicalBase() = default;
};
std::ostream& operator<<(std::ostream& out, const PointPhysicalBase& point);
}  // namespace Geometry

#endif // HPGEM_KERNEL_POINTPHYSICALBASE_H

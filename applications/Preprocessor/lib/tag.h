/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2021, University of Twente
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
#ifndef HPGEM_TAG_H
#define HPGEM_TAG_H

#include <memory>

namespace Preprocessor {

/**
 * Simple tag to statically pass an integer.
 *
 * For example used when one needs to specialize templates of the form:
 *
 * template<std::size_t dimension>
 * void function(non-template-arguments)
 *
 * Which can't be specialized, as the template argument is not present in the
 * argument list, nor return type. To allow specialization add the tag:
 *
 * template<std::size_t dimension>
 * void function(non-template-arguments, tag<dimension>)
 *
 * @tparam i The integer to pass.
 */
template <std::size_t i>
struct tag {};

/**
 * @brief Same as tag, but with integer parameter.
 *
 * Example use is similar to that of tag, but when using recursion on the
 * dimension this allows having itag<-1> as base case. With tag, the base case
 * should be tag<0>, which usually means that the logic of tag<d> needs to be
 * duplicated.
 */
template <int>
struct itag {};

}  // namespace Preprocessor

#endif  // HPGEM_TAG_H

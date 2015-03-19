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

#ifndef REFERENCETOREFERENCEM_H_
#define REFERENCETOREFERENCEM_H_

#include "MappingInterface.h"

namespace Geometry
{
    
    /*! \brief Intermediate ABC for reference space to reference space mappings.

     Mappings from and to reference space are used in two contexts:

     o from a face of a ReferenceGeometry onto the Geometry, hence (dim-1)->dim; this is needed
     when evaluating the trace of a function (defined on the elements).

     o from a ReferenceGeometry onto itself, needed for the face-to-face mappings (because of
     the possible permutations of the global node indices of the face nodes when seen from the
     two elements on either side of the face).

     Note that the implementations have constant state, as they concern only reference coordinates.
     There is no need for reinit function, as the ReferenceToPhysicalMappings. Also transform can be
     declared here because, in contrast to Ref2PhysSpaceMapping, it concerns RefSpacePoint objects
     as both arguments. */

    class MappingReferenceToReference : public MappingInterface
    {
    public:
        virtual PointReference transform(const Geometry::PointReference&) const = 0;
    };
}
;
#endif /* REFERENCETOPHYSICALM_H_ */

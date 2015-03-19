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

#ifndef ____ConcatenatedMapping__
#define ____ConcatenatedMapping__

#include "MappingReferenceToReference.h"
//------------------------------------------------------------------------------

namespace Geometry
{
    /*! ConcatenatedMapping allows to compose two mappings of type
     *  Ref2RefSpaceMapping and treat the result as a mapping itself. This
     *  functionality is needed for the face-to-reference element
     *  mappings. There, for the mapping to the right (R) element, the
     *  face-to-face mapping (from left to right side) and the right side's
     *  face-to-element mapping have to be combined. The resulting mapping can
     *  be handed out by the face, and that's where ConcatenatedMapping comes
     *  in. Note that it is necessary to allow different dimensions, since a
     *  face-to-element mapping is always (dim-1) -> dim.
     *
     *  Internally, ConcatenatedMapping keeps references only, so that it is
     *  clear that this class will not delete anything. Hence the lifetime of
     *  the ConcatenatedMapping must not exceed the one of the component
     *  mappings, but this is hardly possible for reasons in the mathematical
     *  usage. */
    class ConcatenatedMapping : public MappingReferenceToReference
    {
    public:
        //! Ctor gets two references to existing mappings.
        ConcatenatedMapping(const MappingReferenceToReference& m1, const MappingReferenceToReference& m2)
                : map1_(m1), map2_(m2)
        {
        }
        
        //! Transformation is simply via the intermediate space.
        virtual PointReference transform(const PointReference& pIn) const;

        //! To compute the Jacobian, the two component ones have to multiplied.
        virtual Jacobian calcJacobian(const PointReference& p) const;

        virtual std::size_t getTargetDimension() const;

    private:
        const MappingReferenceToReference& map1_;
        const MappingReferenceToReference& map2_;
    };

} // close namespace Geometry

#endif /* defined(____ConcatenatedMapping__) */

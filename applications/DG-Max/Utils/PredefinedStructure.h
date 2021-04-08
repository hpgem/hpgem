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
#ifndef HPGEM_PREDEFINEDSTRUCTURE_H
#define HPGEM_PREDEFINEDSTRUCTURE_H

#include "StructureDescription.h"

namespace DGMax {


// Predefined structures
enum class PredefinedStructure : std::size_t {
    VACUUM = 0,
    BRAGG_STACK = 1,           // [-inf,0.5] material 1, [0.5, inf] material 2
    CYLINDER = 2,              // Cylinder at x,y (0.5, 0.5) radius 0.2
    SQUARE_HOLE = 3,           // Square hole for x,y in [0.1, 0.9] x [0.1, 0.9]
    INVERSE_WOODPILE_OLD = 4,  // Old IW definition with a=1, c=1/sqrt(2)
    INVERSE_WOODPILE_NEW = 5,  // New IW definition with a=sqrt(2),c=1
};

/// 'Parse' an integer to a PredefinedStructure
PredefinedStructure structureFromInt(std::size_t value);

/// Some predefined structures.
///
/// Both for historical meshes without zone information and structured meshes
class PredefinedStructureDescription : public StructureDescription {
   public:
    PredefinedStructureDescription(PredefinedStructure structure,
                                  std::size_t dimension)
        : structure_(structure), dimension(dimension){};

    ElementInfos* createElementInfo(const Base::Element* element) final;

   private:
    template<std::size_t DIM>
    ElementInfos* createElementInfoDim(const Base::Element* element) const;

    PredefinedStructure structure_;
    std::size_t dimension;
};
}

#endif  // HPGEM_PREDEFINEDSTRUCTURE_H

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
#ifndef HPGEM_VTKELEMENT_H
#define HPGEM_VTKELEMENT_H

#include <cstdint>
#include <vector>

#include "Geometry/PointReference.h"

namespace hpgem {

namespace Output {

/**
 * \brief VTK description of a type of mesh elements
 *
 * Essential information about a single type of reference element that can be
 * used in a VTK file (in VTK unstructured grid these are known as cells).
 *
 * @tparam DIM The dimension of the geometry
 */
template <std::size_t DIM>
class VTKElement {
   public:
    virtual ~VTKElement() = default;
    /**
     * Magic number identifying the element/cell type.
     *
     * @see
     * https://github.com/Kitware/VTK/blob/0ce0d74e67927fd964a27c045d68e2f32b5f65f7/Common/DataModel/vtkCellType.h
     * @return
     */
    virtual std::uint8_t vtkId() const = 0;

    /**
     * The reference coordinates on the corresponding ReferenceGeometry. The
     * ordering of these is used to for the corresponding physical coordinates
     * in the output file.
     *
     * @return A vector with reference coordinates of the points
     */
    virtual const std::vector<Geometry::PointReference<DIM>>& getPoints()
        const = 0;
};

}  // namespace Output
}  // namespace hpgem

#endif  // HPGEM_VTKELEMENT_H

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

#ifndef KERNEL_GEOMETRY_MAPPINGS_REFINEMENTMAPSFORTETRAHEDRON_H_
#define KERNEL_GEOMETRY_MAPPINGS_REFINEMENTMAPSFORTETRAHEDRON_H_


#include "RefinementMapping.h"
#include "RefinementMapsForTriangle.h"
#include "Geometry/PointReference.h"

namespace Geometry
{
    /* Stores all the refinement mappings for a tetrahedron.
     *
     * 3 o
     *   |\
     *   |  \2
     *   |  / \
     *   |/     \
     * 0 o--------o 1
     *
     *     index 0
     *
     * The legacy implementation did not have any refinement maps for tetrahedra, so only the identity mapping is provided
     *
     */

    class RefinementMapForTetrahedron0 : public RefinementMapping
    {
    public:
        static const RefinementMapForTetrahedron0* instance()
        {
            static RefinementMapForTetrahedron0 theInstance;
            return &theInstance;
        }

        std::string getName() const override final
        {
            return "Identity map";
        }

        PointReference<3> refinementTransform(std::size_t subElementIndex, const PointReference<3>& p) const override final
        {
            logger.assert(subElementIndex == 0, "asked for subElement %, but the % has only % subElements", subElementIndex, getName(), getNumberOfSubElements());
            return p;
        }

        LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixL(std::size_t subElementIndex, const PointReference<3>& p) const override final
        {
            logger.assert(subElementIndex == 0, "asked for subElement %, but the % has only % subElements", subElementIndex, getName(), getNumberOfSubElements());
            //the bonus column is most likely {0, 0, 0, 1}
            return LinearAlgebra::SmallMatrix<3, 3>{{{{1., 0., 0.}}, {{0., 1., 0.}}, {{0., 0., 1.}}}};
        }

        LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixR(std::size_t subElementIndex, const PointReference<3>& p) const override final
        {
            logger.assert(subElementIndex == 0, "asked for subElement %, but the % has only % subElements", subElementIndex, getName(), getNumberOfSubElements());
            //the bonus column is most likely {0, 0, 0, 1}
            return LinearAlgebra::SmallMatrix<3, 3>{{{{1., 0., 0.}}, {{0., 1., 0.}}, {{0., 0., 1.}}}};
        }

        std::size_t getNumberOfNewNodes() const override final
        {
            return 0;
        }

        std::size_t getNumberOfSubElements() const override final
        {
            return 1;
        }

        std::vector<std::size_t> getSubElementLocalNodeIndices(std::size_t subElementIndex) const override final
        {
            logger.assert(subElementIndex == 0, "asked for subElement %, but the % has only % subElements", subElementIndex, getName(), getNumberOfSubElements());
            return std::vector<std::size_t>{0, 1, 2, 3};
        }

        std::vector<const RefinementMapping*> getCodim1RefinementMaps() const override final
        {
            return std::vector<const RefinementMapping*>{RefinementMapForTriangle0::instance(), RefinementMapForTriangle0::instance(), RefinementMapForTriangle0::instance(), RefinementMapForTriangle0::instance()};
        }
    private:
        RefinementMapForTetrahedron0() = default;
        RefinementMapForTetrahedron0(const RefinementMapForTetrahedron0&) = delete;
    };
}



#endif /* KERNEL_GEOMETRY_MAPPINGS_REFINEMENTMAPSFORTETRAHEDRON_H_ */

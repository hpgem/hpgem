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
#ifndef RefinementMapping_h
#define RefinementMapping_h

#include <iostream>
#include <string>
#include "LinearAlgebra/SmallMatrix.h"
#include "Geometry/PointReference.h"
#include "Geometry/ReferenceGeometry.h"

namespace Geometry {
template <std::size_t DIM>
class PointReference;

class RefinementMapping {
   public:
    virtual ~RefinementMapping() = default;

    virtual std::string getName() const = 0;

    //! Transform a reference point using refinement mapping
    virtual PointReference<0> refinementTransform(
        std::size_t subElementIdx, const PointReference<0>& p) const {
        logger(ERROR, "You passed a coordinate with the wrong dimension!");
        return p;
    }

    virtual PointReference<1> refinementTransform(
        std::size_t subElementIdx, const PointReference<1>& p) const {
        logger(ERROR, "You passed a coordinate with the wrong dimension!");
        return p;
    }

    virtual PointReference<2> refinementTransform(
        std::size_t subElementIdx, const PointReference<2>& p) const {
        logger(ERROR, "You passed a coordinate with the wrong dimension!");
        return p;
    }

    virtual PointReference<3> refinementTransform(
        std::size_t subElementIdx, const PointReference<3>& p) const {
        logger(ERROR, "You passed a coordinate with the wrong dimension!");
        return p;
    }

    virtual PointReference<4> refinementTransform(
        std::size_t subElementIdx, const PointReference<4>& p) const {
        logger(ERROR, "You passed a coordinate with the wrong dimension!");
        return p;
    }

    /// you have to pass one reference point so the refinement mapping can get
    /// the dimension right \todo think of a solution that doesn't need a dummy
    ///argument
    virtual std::vector<PointReference<0>> getNewNodeLocations(
        const PointReference<0>&) const {
        logger(ERROR, "You passed a coordinate with the wrong dimension!");
        return {};
    }

    /// you have to pass one reference point so the refinement mapping can get
    /// the dimension right \todo think of a solution that doesn't need a dummy
    ///argument
    virtual std::vector<PointReference<1>> getNewNodeLocations(
        const PointReference<1>&) const {
        logger(ERROR, "You passed a coordinate with the wrong dimension!");
        return {};
    }

    /// you have to pass one reference point so the refinement mapping can get
    /// the dimension right \todo think of a solution that doesn't need a dummy
    ///argument
    virtual std::vector<PointReference<2>> getNewNodeLocations(
        const PointReference<2>&) const {
        logger(ERROR, "You passed a coordinate with the wrong dimension!");
        return {};
    }

    /// you have to pass one reference point so the refinement mapping can get
    /// the dimension right \todo think of a solution that doesn't need a dummy
    ///argument
    virtual std::vector<PointReference<3>> getNewNodeLocations(
        const PointReference<3>&) const {
        logger(ERROR, "You passed a coordinate with the wrong dimension!");
        return {};
    }

    /// you have to pass one reference point so the refinement mapping can get
    /// the dimension right \todo think of a solution that doesn't need a dummy
    ///argument
    virtual std::vector<PointReference<4>> getNewNodeLocations(
        const PointReference<4>&) const {
        logger(ERROR, "You passed a coordinate with the wrong dimension!");
        return {};
    }

    //! Transformation matrix of this refinement when located on the LEFT side.
    //! This is also the jacobian of the refinementTransform
    virtual LinearAlgebra::SmallMatrix<0, 0> getRefinementMappingMatrixL(
        std::size_t subElementIdx, const PointReference<0>& p) const {
        logger(ERROR, "You passed a coordinate with the wrong dimension!");
        return LinearAlgebra::SmallMatrix<0, 0>();
    }

    virtual LinearAlgebra::SmallMatrix<1, 1> getRefinementMappingMatrixL(
        std::size_t subElementIdx, const PointReference<1>& p) const {
        logger(ERROR, "You passed a coordinate with the wrong dimension!");
        return LinearAlgebra::SmallMatrix<1, 1>();
    }

    virtual LinearAlgebra::SmallMatrix<2, 2> getRefinementMappingMatrixL(
        std::size_t subElementIdx, const PointReference<2>& p) const {
        logger(ERROR, "You passed a coordinate with the wrong dimension!");
        return LinearAlgebra::SmallMatrix<2, 2>();
    }

    virtual LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixL(
        std::size_t subElementIdx, const PointReference<3>& p) const {
        logger(ERROR, "You passed a coordinate with the wrong dimension!");
        return LinearAlgebra::SmallMatrix<3, 3>();
    }

    virtual LinearAlgebra::SmallMatrix<4, 4> getRefinementMappingMatrixL(
        std::size_t subElementIdx, const PointReference<4>& p) const {
        logger(ERROR, "You passed a coordinate with the wrong dimension!");
        return LinearAlgebra::SmallMatrix<4, 4>();
    }

    //! Transformation matrix of this refinement when located on the RIGHT side.
    //! This is the inverse of the LEFT transformation
    // it would be nicer if the return type was a smallMatrix, but it is illegal
    // to change only the return type and nothing else
    virtual LinearAlgebra::SmallMatrix<0, 0> getRefinementMappingMatrixR(
        std::size_t subElementIdx, const PointReference<0>& p) const {
        logger(ERROR, "You passed a coordinate with the wrong dimension!");
        return LinearAlgebra::SmallMatrix<0, 0>();
    }

    virtual LinearAlgebra::SmallMatrix<1, 1> getRefinementMappingMatrixR(
        std::size_t subElementIdx, const PointReference<1>& p) const {
        logger(ERROR, "You passed a coordinate with the wrong dimension!");
        return LinearAlgebra::SmallMatrix<1, 1>();
    }

    virtual LinearAlgebra::SmallMatrix<2, 2> getRefinementMappingMatrixR(
        std::size_t subElementIdx, const PointReference<2>& p) const {
        logger(ERROR, "You passed a coordinate with the wrong dimension!");
        return LinearAlgebra::SmallMatrix<2, 2>();
    }

    virtual LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixR(
        std::size_t subElementIdx, const PointReference<3>& p) const {
        logger(ERROR, "You passed a coordinate with the wrong dimension!");
        return LinearAlgebra::SmallMatrix<3, 3>();
    }

    virtual LinearAlgebra::SmallMatrix<4, 4> getRefinementMappingMatrixR(
        std::size_t subElementIdx, const PointReference<4>& p) const {
        logger(ERROR, "You passed a coordinate with the wrong dimension!");
        return LinearAlgebra::SmallMatrix<4, 4>();
    }

    //! \brief Number of new nodes due to the refinement that should be added to
    //! the vector of localNodeIndices
    virtual std::size_t getNumberOfNewNodes() const = 0;

    //! \brief Number of sub-elements due to the refinement
    virtual std::size_t getNumberOfSubElements() const = 0;

    //! \brief Assembly nodes for sub-element
    virtual std::vector<std::size_t> getSubElementLocalNodeIndices(
        std::size_t subElement) const = 0;

    virtual Geometry::ReferenceGeometry* getBigElementReferenceGeometry()
        const = 0;

    virtual Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const = 0;

    virtual std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const = 0;

    ///\brief returns the index of the sub-element and the local index of the
    ///specified subface in the sub-element
    virtual std::tuple<std::size_t, std::size_t> getSubElementAndLocalFaceIndex(
        std::size_t face, std::size_t subFaceIndex) const = 0;

    ///\brief return the nodes that belong on the refined face; in the ordering
    ///dictated by the refined face
    virtual std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const = 0;
};
}  // namespace Geometry
#endif

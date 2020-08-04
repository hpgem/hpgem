/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2020, University of Twente
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
#ifndef HPGEM_CGDGMATRIXKPHASESHIFTBUILDER_H
#define HPGEM_CGDGMATRIXKPHASESHIFTBUILDER_H

#include "KPhaseShift.h"

namespace DGMax {
/// \brief Builder for KPhaseShift for a matrix with conforming basis function
/// rows and discontinuous basis functions columns, both having the
/// k-phase-shifted periodic boundary conditions.
///
/// For conforming basis functions there will be basis functions that have
/// support on two or more elements that are separated by k-phase-shifted
/// boundary conditions. To keep these basis functions conforming with
/// k-phase-shifted boundary conditions they need a phase factor e^{ika_T} that
/// is dependent on the element T. Thus for element matrices formed between such
/// conforming basis functions and DG basis functions, they gain the phase
/// factor of the element T on which the DG basis function has support.
///
/// To fix a_T for a conforming basis function we look at the geometric element
/// (element, face, edge, node) to which it is associated. This geometrical
/// element has a list of elements to which it belongs, the first element from
/// this list is taken and on this element a_T is set to zero. The value of a_T
/// on the other elements T' on which the basis function has support can then be
/// derived. For this take a point x that is shared by T and T' and let x_T and
/// x_T' be the coordinate as observed from the respective elements, then a_T =
/// x_T - x_T'.
///
/// This builder builds the KPhaseShift for the element matrices of elements
/// that touch the k-phase-shifted periodic boundary. As on these elements the
/// discontinuous basis functions overlap with conforming basis functions that
/// have support on the periodic boundary. For compatibility with the extra
/// phase shift vector b used in FaceMatrixKPhaseShiftBuilder, there is the
/// possibility to subtract a vector b from a_T, where b depends on the element.
///
/// \tparam DIM The dimension of the domain.
template <std::size_t DIM>
class CGDGMatrixKPhaseShiftBuilder {
   public:
    using MatrixExtractor =
        std::function<LinearAlgebra::MiddleSizeMatrix(const Base::Element*)>;

    CGDGMatrixKPhaseShiftBuilder()
        : matrixExtractor_(nullptr),
          extraShift_(nullptr),
          cgIndexing_(nullptr),
          dgIndexing_(nullptr),
          cgUnknowns_(),
          dgUnknowns_(),
          hermitian_(false){};

    KPhaseShifts<DIM> build() const;

    /// Set the function that provides the element matrix for an element that
    /// touches the boundary.
    void setMatrixExtractor(MatrixExtractor extractor) {
        matrixExtractor_ = std::move(extractor);
    }

    /// Extra shift b, that is to be subtracted from a_T.
    void setExtraShift(
        std::function<LinearAlgebra::SmallVector<DIM>(const Base::Element*)>
            extraShift) {
        extraShift_ = extraShift;
    }

    void setIndices(const Utilities::GlobalIndexing* cgIndexing,
                    const Utilities::GlobalIndexing* dgIndexing) {
        logger.assert_debug(cgIndexing != nullptr, "Null cg index");
        logger.assert_debug(cgIndexing != nullptr, "Null dg index");
        cgIndexing_ = cgIndexing;
        dgIndexing_ = dgIndexing;
        cgUnknowns_ = cgIndexing_->getIncludedUnknowns();
        dgUnknowns_ = dgIndexing_->getIncludedUnknowns();
    }

    void setCGUnknowns(const std::vector<std::size_t>& cgUnknowns);
    void setDGUnknowns(const std::vector<std::size_t>& dgUnknowns);

    /// Whether this to phase shift just CG rows and DG columns, or also the DG
    /// rows and CG columns. If set to true it is expected that cgIndexing ==
    /// dgIndexing and that cgUnknowns and dgUnknowns are disjoint.
    void setHermitian(bool hermitian) { hermitian_ = hermitian; }

   private:
    MatrixExtractor matrixExtractor_;
    std::function<LinearAlgebra::SmallVector<DIM>(const Base::Element*)>
        extraShift_;

    const Utilities::GlobalIndexing* cgIndexing_;
    const Utilities::GlobalIndexing* dgIndexing_;

    std::vector<std::size_t> cgUnknowns_;
    std::vector<std::size_t> dgUnknowns_;

    bool hermitian_;

    void addFacePhaseShifts(
        const Base::Face* face,
        std::vector<DGMax::KPhaseShiftBlock<DIM>>& out) const;

    void addEdgePhaseShifts(
        const Base::Edge* edge,
        std::vector<DGMax::KPhaseShiftBlock<DIM>>& out) const;
    void addNodePhaseShifts(
        const Base::Node* node,
        std::vector<DGMax::KPhaseShiftBlock<DIM>>& out) const;

    template <typename GEOM>
    void addElementPhaseShift(
        const GEOM* geom, const Geometry::PointPhysical<DIM>& owningCoord,
        const Base::Element* element,
        std::vector<DGMax::KPhaseShiftBlock<DIM>>& out) const;
};
}  // namespace DGMax

class CGDGMatrixKPhaseShiftBuilder {};

#endif  // HPGEM_CGDGMATRIXKPHASESHIFTBUILDER_H

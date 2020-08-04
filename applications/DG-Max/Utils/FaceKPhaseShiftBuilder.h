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
#ifndef HPGEM_FACEKPHASESHIFTBUILDER_H
#define HPGEM_FACEKPHASESHIFTBUILDER_H

#include "KPhaseShift.h"

namespace DGMax {

/// \brief Builder for KPhaseShifts for a DG-matrix where the k-phase-shift
/// happens through the face coupling at periodic boundary faces.
///
/// For a DG the k-phase-shifted boundary conditions for a matrix S(k) result in
/// that the face matrices for the faces on the periodic boundary gain some
/// phase factors e^{ika}. The lattice vector a here is the difference x_l -
/// x_r, where x_l and x_r are the coordinates for the same point on the
/// periodic boundary face, as observed from the left or right element
/// respectively. On the face matrices this introduces a phase factor e^{ika}
/// for the rows corresponding to the right element, and e^{-ika} for the
/// columns of that element. The effect is thus that the two off diagonal blocks
/// gain a phase factor e^{+-ika}.
///
/// This builder creates a KPhaseShifts that represents these off-diagonal
/// blocks of face matrices for faces on the periodic boundary. In addition to
/// this standard phase shift it also supports the possibility for an extra
/// vector b that is added to a.
///
/// \tparam DIM The dimension of the problem
template <std::size_t DIM>
class FaceMatrixKPhaseShiftBuilder {
   public:
    using MatrixExtractor =
        std::function<std::pair<LinearAlgebra::MiddleSizeMatrix,
                                LinearAlgebra::MiddleSizeMatrix>(
            const Base::Face*)>;

    KPhaseShifts<DIM> build() const;

    /// Set the function to get the two off diagonal face matrix blocks for a
    /// periodic boundary face. The first matrix should correspond to the test
    /// functions of the left element with the trial functions of the right
    /// element, with the second matrix vice versa. If the GlobalIndex does not
    /// include all unknowns, then the matrices rows and columns should only be
    /// for DoFs that correspond to included unknowns.
    ///
    /// \param extractor The function used to extract the two matrix blocks
    void setMatrixExtractor(MatrixExtractor extractor) {
        matrixExtractor_ = std::move(extractor);
    }

    /// Set function to compute the extra phase shift vector b that is to be
    /// added to a.
    ///
    /// \param extraShift The function to compute the shift.
    void setExtraShift(
        std::function<LinearAlgebra::SmallVector<DIM>(const Base::Face*)>
            extraShift) {
        extraShift_ = extraShift;
    }

    void setIndexing(const Utilities::GlobalIndexing* indexing) {
        logger.assert_debug(indexing != nullptr, "Null index");
        indexing_ = indexing;
        unknowns_ = indexing->getIncludedUnknowns();
    }

    void setUnknowns(const std::vector<std::size_t>& unknowns);

   private:
    KPhaseShiftBlock<DIM> facePhaseShift(const Base::Face* face) const;

    const Utilities::GlobalIndexing* indexing_;
    std::vector<std::size_t> unknowns_;

    MatrixExtractor matrixExtractor_;
    std::function<LinearAlgebra::SmallVector<DIM>(const Base::Face*)>
        extraShift_;
};
}  // namespace DGMax

#endif  // HPGEM_FACEKPHASESHIFTBUILDER_H

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
#ifndef HPGEM_KPHASESHIFT_H
#define HPGEM_KPHASESHIFT_H

#include <memory>
#include <utility>

#include "LinearAlgebra/SmallVector.h"
#include "Utilities/GlobalIndexing.h"

#include "MatrixBlocks.h"

namespace hpgem {

namespace LinearAlgebra {
class MiddleSizeMatrix;
}

namespace Base {
class Element;
class Face;
class Edge;
class Node;
}  // namespace Base

namespace Geometry {
template <std::size_t DIM>
class PointPhysical;
}

namespace Utilities {
class GlobalIndexing;
}
}  // namespace hpgem

namespace DGMax {

using namespace hpgem;

/// A MatrixBlock which incurs a phase factor e^{ikx} (e^{-ikx} for the
/// symmetric counter part).
/// \tparam DIM The dimension of the problem
template <std::size_t DIM>
class KPhaseShiftBlock {
   public:
    KPhaseShiftBlock(MatrixBlocks blocks, LinearAlgebra::SmallVector<DIM> dx)
        : blocks_(std::move(blocks)), dx_(dx){};

    /// Perform the actual phase shift
    ///
    /// \param k The wave vector
    /// \param storage Temporary storage space to use
    /// \param mat The matrix in which to insert the block(s).
    void apply(LinearAlgebra::SmallVector<DIM> k,
               std::vector<PetscScalar>& storage, Mat mat) const;

   private:
    /// Check for a k-vector if a phase factor dx_ * k is present
    bool shiftNeeded(LinearAlgebra::SmallVector<DIM> k) const {
        return std::abs(dx_ * k) > 1e-12;
    }

    /// The blocks that need to be shifted
    MatrixBlocks blocks_;
    /// The distance x for the phase factor e^{ikx}.
    LinearAlgebra::SmallVector<DIM> dx_;
};

/// \brief Phase Shifts in a GlobalMatrix which has a k-shifted periodic
/// boundary.
///
/// Background:
/// With standard periodic boundary conditions a solution u satisfies the
/// original PDE when traversing over the periodic boundary. Thus if x and x+a
/// are the coordinates on two sides of a periodic boundary then x and x+a are
/// the same point. Periodic boundary conditions then identify u(x) and u(x+a),
/// including their derivatives. For k-shifted boundary conditions we do not
/// identify u(x) with u(x+a), but u(x) e^{ika} with u(x+a), where k is the wave
/// vector. As result a global FEM matrix A becomes a parameterized matrix A(k).
/// In practice only very few entries, those corresponding to elements on the
/// periodic boundary, depend on k. Where this k-dependence is of the form
/// e^{ika} for some a.
///
/// This class represents a way to efficiently implement these phase shifts for
/// a k-shifted periodic boundary. This is done under a twofold assumption:
///   - For each entry the phase shift is of the form e^{ika}, where the vector
///     `a` may depend on the entry but not on `k`.
///   - The entries for k=0 are known and fixed, the k-dependent entries for
///     S(k) can therefore be replaced by the known entries multiplied with
///     e^{ika}.
///
/// Implementation note. The alternative to replacing the block is taking it
/// out and multiplying it with exp(i dk*x), with dk the change with the
/// previous k-vector. However, after replacing a block a call to assemble is
/// needed, resulting in possibly significant synchronization overhead. The
/// current approach only needs a single assembly call.
/// \tparam DIM The dimension
template <std::size_t DIM>
class KPhaseShifts {
   public:
    KPhaseShifts() = default;

    explicit KPhaseShifts(std::vector<KPhaseShiftBlock<DIM>> blocks)
        : blocks_(std::move(blocks)){};
    void apply(LinearAlgebra::SmallVector<DIM> k, Mat mat) const;
    std::vector<KPhaseShiftBlock<DIM>> blocks_;

   private:
};

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
    /// Function to compute the offdiagonal blocks that need to be
    /// KPhaseShift-ed.
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

    /// Set the function that extracts an element matrix for the given element,
    /// where the rows correspond the the CG basis functions and the columns to
    /// the DG ones. If the GlobalIndex of the cg or dg part does not include
    /// all unknowns, then the matrix should only have rows or columns for the
    /// included unknowns of the corresponding GlobalIndex.
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

};  // namespace DGMax

#endif  // HPGEM_KPHASESHIFT_H

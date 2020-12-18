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

};  // namespace DGMax

#endif  // HPGEM_KPHASESHIFT_H

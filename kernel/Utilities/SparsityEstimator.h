
#ifndef HPGEM_KERNEL_SPARSITYESTIMATOR_H
#define HPGEM_KERNEL_SPARSITYESTIMATOR_H

#include <vector>
#include "Utilities/Table2D.h"

namespace hpgem {

// Forward definitions
namespace Base {
class Element;
class MeshManipulatorBase;
}  // namespace Base
namespace Utilities {
class GlobalIndexing;
}  // namespace Utilities

namespace Utilities {
/// Computation for the estimate (upper bound) on the number of non zero entries
/// in a matrix.
class SparsityEstimator {
   public:
    /// Construct a SparsityEstimator with the same GlobalIndexing for rows and
    /// columns \param indexing The indexing for the basis functions
    SparsityEstimator(const GlobalIndexing& indexing)
        : SparsityEstimator(indexing, indexing) {}

    /// Construct a sparsity estimator with possibly different indices for the
    /// rows and columns. \param rowIndexing \param columnIndexing
    SparsityEstimator(const GlobalIndexing& rowIndexing,
                      const GlobalIndexing& columnIndexing);

    /// Compute the sparsity estimate
    /// \param includeFaceCoupling [in] Include coupling through face matrices
    ///
    /// \return Pair of vectors, each entry corresponding to a single row DoF in
    /// processor local ordering. For each row DoF the entry contains the number
    /// of column DoFs with which it can form a non zero entry. The first vector
    /// is for column DoFs owned by the current process, while second vector
    /// for those owned by different processes
    std::pair<std::vector<int>, std::vector<int>> computeSparsityEstimate(
        bool includeFaceCoupling = true) const;

    std::pair<std::vector<int>, std::vector<int>> computeSparsityEstimate(
        const Table2D<bool>& faceCoupling) const;

   private:
    /// The indexing for the basis functions on the rows of the matrix
    const GlobalIndexing& rowIndexing_;
    /// The indexing for the basis functions on the columns of the matrix
    const GlobalIndexing& columnIndexing_;
};
}  // namespace Utilities

}  // namespace hpgem

#endif  // HPGEM_KERNEL_SPARSITYESTIMATOR_H

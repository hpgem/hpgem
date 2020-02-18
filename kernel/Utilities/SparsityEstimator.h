
#ifndef HPGEM_SPARSITYESTIMATOR_H
#define HPGEM_SPARSITYESTIMATOR_H

#include "vector"

// Forward definitions
namespace Base
{
    class Element;
    class MeshManipulatorBase;
}
namespace Utilities
{
    class GlobalIndexing;
}

namespace Utilities
{
    /// Computation for the estimate (upper bound) on the number of non zero entries in a matrix.
    class SparsityEstimator
    {
    public:
        /// Construct a SparsityEstimator
        /// \param mesh  The mesh
        /// \param indexing The indexing for the basis functions
        SparsityEstimator(const Base::MeshManipulatorBase& mesh, const GlobalIndexing& indexing);

        /// Compute the sparsity estimate
        /// \param nonZeroPerRowOwned [out] Per local DoF the number of non zero columns from
        ///     locally owned basis functions.
        /// \param nonZeroPerRowNonOwned [out] Per local DoF the number of non zero columns
        ///     from non locally owned basis functions.
        void computeSparsityEstimate(std::vector<int>& nonZeroPerRowOwned, std::vector<int>& nonZeroPerRowNonOwned) const;

    private:
        /// Workspace variables for the computation
        struct Workspace;

        /// Add the DoFs that have support on an Element to the workspace.
        void addElementDoFs(const Base::Element* element, Workspace& workspace) const;


        /// The mesh for the sparisity estimate
        const Base::MeshManipulatorBase& mesh_;
        /// The indexing for the basis functions
        const GlobalIndexing& indexing_;
    };
}


#endif //HPGEM_SPARSITYESTIMATOR_H

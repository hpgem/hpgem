
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


        /// Write the DoF count for an Element, Face, etc. to the sparsity estimate vector
        /// \tparam GEOM The Element, Face, etc. type
        /// \param geom The Element, Face, etc. which supports the basis functions for which to write the DoF count.
        /// \param workspace The workspace with DoF counts
        /// \param nonZeroPerRowOwned see computeSparsityEstimate(std::vector<int>&, std::vector<int>&)
        /// \param nonZeroPerRowNonOwned see computeSparsityEstimate(std::vector<int>&, std::vector<int>&)
        template<typename GEOM>
        void writeDoFCount(const GEOM* geom, const Workspace& workspace,
                std::vector<int>& nonZeroPerRowOwned, std::vector<int>& nonZeroPerRowNonOwned) const;

        /// The mesh for the sparsity estimate
        const Base::MeshManipulatorBase& mesh_;
        /// The indexing for the basis functions
        const GlobalIndexing& indexing_;
    };
}


#endif //HPGEM_SPARSITYESTIMATOR_H

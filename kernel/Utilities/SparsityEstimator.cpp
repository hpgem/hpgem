
#include "SparsityEstimator.h"

#include "Logger.h"

#include <map>
#include <set>
#include <type_traits>
#include "Base/Element.h"
#include "Base/MeshManipulatorBase.h"
#include "Utilities/GlobalIndexing.h"
#include "Utilities/Table2D.h"

namespace hpgem {

namespace Utilities {
SparsityEstimator::SparsityEstimator(const GlobalIndexing& rowindexing,
                                     const GlobalIndexing& columnindexing)
    : rowIndexing_(rowindexing), columnIndexing_(columnindexing) {}

namespace Detail {
/// An interval (start, size) of global DoF-Indices
/// We assume that each interval is for a single unknown and all the DoFs in the
/// interval are associated with the same Element, Face, Edge or Node.
using DoFInterval = std::pair<std::size_t, std::size_t>;

/// Set of DoFIntervals
using DoFIntervals = std::set<DoFInterval>;

/// Enum to specify rows or columns of a matrix
enum MatrixSide { ROWS, COLUMNS };

/// For an element or face matrix, the global DoF-indices to which the rows or
/// columns correspond, grouped by the unknown.
///
/// The global DoF-indices differ for rows and columns:
///  - For columns we store the global DoF-indices, separated by ownership
///  - For rows we store the processor local global DoF-indices (and ignoring
///    rows owned by other processors).
class LocalMatrixSideDoFs {
   public:
    LocalMatrixSideDoFs(const Utilities::GlobalIndexing& indexing,
                        MatrixSide side)
        : indexing_(indexing),
          dimension_(indexing.getMesh()->dimension()),
          side_(side),
          isOnlyElementBased_(indexing.getTotalNumberOfUnknowns(), true),
          ownedIntervalsByUnknown_(indexing.getTotalNumberOfUnknowns()),
          nonOwnedIntervalsByUnknown_(indexing.getTotalNumberOfUnknowns())
    // By default assume all unknowns are only element based.
    {}

    /// Clear all values associated with the current local matrix.
    void clear() {
        for (auto& intervals : ownedIntervalsByUnknown_) {
            intervals.clear();
        }
        for (auto& intervals : nonOwnedIntervalsByUnknown_) {
            intervals.clear();
        }
    }

    /**
     * Add the global DoFs for the rows/columns of a FaceMatrix.
     * @param face The face for of the Face matrix.
     */
    void addLocalFaceMatrixDoFs(const Base::Face* face) {
        // DoFs that have support on the face will be added twice, but are
        // deduplicated by DoFIntervals.
        addLocalElementMatrixDoFs(face->getPtrElementLeft());
        if (face->isInternal()) {
            addLocalElementMatrixDoFs(face->getPtrElementRight());
        }
    }

    /**
     * Add the global DoFs for the rows/columns of an ElementMatrix.
     * @param element The element of the Element matrix.
     */
    void addLocalElementMatrixDoFs(const Base::Element* element) {
        addLocal(element);
        for (const Base::Face* face : element->getFacesList()) {
            addLocal(face);
        }
        for (const Base::Edge* edge : element->getEdgesList()) {
            addLocal(edge);
        }
        if (dimension_ > 1) {
            for (const Base::Node* node : element->getNodesList()) {
                addLocal(node);
            }
        }
    }

    const DoFIntervals& getOwnedIntervals(std::size_t unknown) {
        return ownedIntervalsByUnknown_[unknown];
    }

    const DoFIntervals& getNonOwnedIntervals(std::size_t unknown) {
        return nonOwnedIntervalsByUnknown_[unknown];
    }

    bool isOnlyElementBased(std::size_t unknown) {
        return isOnlyElementBased_[unknown];
    }

   private:
    /// Add the DoFs local to a geometric part (Element, Face, etc.)
    ///
    /// \tparam GEOM The type of geometric part
    /// \param geom The geometric part
    template <typename GEOM>
    void addLocal(const GEOM* geom) {
        for (const std::size_t& unknown : indexing_.getIncludedUnknowns()) {
            std::size_t size = geom->getLocalNumberOfBasisFunctions(unknown);
            if (size == 0) {
                continue;
            }

            if (!std::is_base_of<Base::Element, GEOM>()) {
                // Keep track of unknowns which have DoFs on the Faces Edges or
                // Nodes.
                isOnlyElementBased_[unknown] = false;
            }

            // For rows we only care about the owned rows.
            if (side_ == ROWS && !geom->isOwnedByCurrentProcessor()) {
                continue;
            }

            std::size_t start;  // start of the interval
            if (side_ == ROWS) {
                // Safe because we already excluded the non owned ones.
                start = indexing_.getProcessorLocalIndex(geom, unknown);
            } else {
                start = indexing_.getGlobalIndex(geom, unknown);
            }
            if (geom->isOwnedByCurrentProcessor()) {
                ownedIntervalsByUnknown_[unknown].emplace(start, size);
            } else {
                nonOwnedIntervalsByUnknown_[unknown].emplace(start, size);
            }
        }
    }

    ///// GLOBAL INFOMRATION ////

    const Utilities::GlobalIndexing& indexing_;
    /// The mesh dimension
    std::size_t dimension_;
    /// Which side of the matrix is this associated to
    MatrixSide side_;
    /// Whether the DoFs of each unknown are only associated with Elements
    /// (true), or also with Faces, Edges and Nodes (false). Only filled for the
    /// unknowns included in indexing_;
    std::vector<bool> isOnlyElementBased_;

    //// LOCAL INFORMATION ////
    ///////////////////////////

    /// The DoFIntervals for the current element/face matrix, that are owned by
    /// the current processor. Only filled for unknowns include in indexing_.
    std::vector<DoFIntervals> ownedIntervalsByUnknown_;
    /// Same as ownedIntervalsByUnknown_ but for the intervals that are owned by
    /// different processors.
    std::vector<DoFIntervals> nonOwnedIntervalsByUnknown_;
};

}  // namespace Detail

void SparsityEstimator::computeSparsityEstimate(
    std::vector<int>& nonZeroPerRowOwned,
    std::vector<int>& nonZeroPerRowNonOwned, bool includeFaceCoupling) const {
    Table2D<bool> faceCoupling(rowIndexing_.getNumberOfIncludedUnknowns(),
                               columnIndexing_.getNumberOfIncludedUnknowns(),
                               includeFaceCoupling);
    computeSparsityEstimate(nonZeroPerRowOwned, nonZeroPerRowNonOwned,
                            faceCoupling);
}
void SparsityEstimator::computeSparsityEstimate(
    std::vector<int>& nonZeroPerRowOwned,
    std::vector<int>& nonZeroPerRowNonOwned,
    const Table2D<bool>& faceCoupling) const {

    /*
     * ** APPROACH **
     * **************
     *
     * For each row-DoF we need to count the number of column-DoFs that it
     * interacts with through element and face matrices. The approach we take
     * here is to mimic the assembly procedure.
     *
     * We loop over all elements and faces, and for a corresponding
     * element/face-matrix we compute the row/column-DoFs that would be coupled
     * in such a matrix. These couplings are then combined so that for each
     * row-DoF that this processor owns, we get a list of column-DoFs to which
     * it couples. The number of these column-DoFs then gives the number of
     * non-zero entries.
     *
     * This procedure is complicated by a few things:
     *  - A combination of row and column DoF that couple could be part of
     *    only one element/face matrix, but also part of multiple ones. So we
     *    need to deduplicate them.
     *  - Face couplings are more expensive to compute and result in far larger
     *    non zero counts. So we skip it if possible.
     *  - On a distributed mesh we need to separately count the column-DoFs that
     *    are owned by the current processor, and those that are owned by the
     *    other processors.
     *  - On a distributed mesh we can not support face coupling of row-DoFs
     *    which are associated with Faces, Edges or Nodes (see explanation
     *    below).
     *  - For efficiency and to reduce the effects of the previous point, we
     *    allow restricting which unknowns couple through face matrices
     *    (faceCoupling table).
     */

    logger.assert_always(faceCoupling.getNumberOfRows() ==
                             rowIndexing_.getNumberOfIncludedUnknowns(),
                         "Different rows in face coupling from the row index");
    logger.assert_always(
        faceCoupling.getNumberOfColumns() ==
            columnIndexing_.getNumberOfIncludedUnknowns(),
        "Different number of columns in faceCoupling from the column index");

    logger(VERBOSE, "Computing sparsity estimate for mesh %",
           rowIndexing_.getMesh());
    const std::size_t totalNumberOfDoF =
        rowIndexing_.getNumberOfLocalBasisFunctions();

    nonZeroPerRowOwned.assign(totalNumberOfDoF, 0);
    nonZeroPerRowNonOwned.assign(totalNumberOfDoF, 0);

    if (rowIndexing_.getMesh() == nullptr || totalNumberOfDoF == 0) {
        // Nothing to do, result has already been resized.
        return;
    }

    Detail::LocalMatrixSideDoFs rowDoFs(rowIndexing_, Detail::ROWS);
    Detail::LocalMatrixSideDoFs columnDoFs(columnIndexing_, Detail::COLUMNS);

    std::map<Detail::DoFInterval, Detail::DoFIntervals> owned;
    std::map<Detail::DoFInterval, Detail::DoFIntervals> nonOwned;

    for (const Base::Element* element :
         rowIndexing_.getMesh()->getElementsList(Base::IteratorType::GLOBAL)) {
        // Face, Edge and Node based DoFs will be used in the element matrices
        // on the ghost layer. So we need to iterate over all elements and not
        // just the locally owned ones.
        rowDoFs.clear();
        columnDoFs.clear();
        rowDoFs.addLocalElementMatrixDoFs(element);
        columnDoFs.addLocalElementMatrixDoFs(element);
        // Everything is coupled, so for each
        for (std::size_t rowUnknown : rowIndexing_.getIncludedUnknowns()) {
            const Detail::DoFIntervals& rowIntervals =
                rowDoFs.getOwnedIntervals(rowUnknown);
            for (std::size_t colUnknown :
                 columnIndexing_.getIncludedUnknowns()) {
                const Detail::DoFIntervals& ownedColIntervals =
                    columnDoFs.getOwnedIntervals(colUnknown);
                const Detail::DoFIntervals& nonOwnedColIntervals =
                    columnDoFs.getNonOwnedIntervals(colUnknown);

                for (Detail::DoFInterval rowInterval : rowIntervals) {
                    owned[rowInterval].insert(ownedColIntervals.begin(),
                                              ownedColIntervals.end());
                    nonOwned[rowInterval].insert(nonOwnedColIntervals.begin(),
                                                 nonOwnedColIntervals.end());
                }
            }
        }
    }

    // If there is no face coupling anywhere (e.g. dg mass matrix) then we can
    // skip looping over the adjacent elements.
    bool anyFaceCoupling = false;
    for (std::size_t i = 0; i < faceCoupling.getSize(); ++i) {
        if (faceCoupling[i]) {
            anyFaceCoupling = true;
            break;
        }
    }

    bool distributedMesh = false;

    if (anyFaceCoupling) {
        for (const Base::Face* face : rowIndexing_.getMesh()->getFacesList()) {
            // Check for a distributed mesh
            if (face->isSubdomainBoundary()) {
                distributedMesh = true;
            }
            // Actual face matrix
            rowDoFs.clear();
            columnDoFs.clear();
            rowDoFs.addLocalFaceMatrixDoFs(face);
            columnDoFs.addLocalFaceMatrixDoFs(face);
            // Add the DoFs
            for (std::size_t ri = 0;
                 ri < rowIndexing_.getNumberOfIncludedUnknowns(); ++ri) {
                std::size_t rowUnknown = rowIndexing_.getIncludedUnknowns()[ri];
                const Detail::DoFIntervals& rowIntervals =
                    rowDoFs.getOwnedIntervals(rowUnknown);

                for (std::size_t ci = 0;
                     ci < columnIndexing_.getNumberOfIncludedUnknowns(); ++ci) {
                    if (!faceCoupling(ri, ci)) {
                        continue;
                    }
                    std::size_t colUnknown =
                        columnIndexing_.getIncludedUnknowns()[ci];

                    const Detail::DoFIntervals& ownedColIntervals =
                        columnDoFs.getOwnedIntervals(colUnknown);
                    const Detail::DoFIntervals& nonOwnedColIntervals =
                        columnDoFs.getNonOwnedIntervals(colUnknown);

                    for (Detail::DoFInterval rowInterval : rowIntervals) {
                        for (Detail::DoFInterval colInterval :
                             ownedColIntervals) {
                            owned[rowInterval].insert(colInterval);
                        }
                        for (Detail::DoFInterval colInterval :
                             nonOwnedColIntervals) {
                            nonOwned[rowInterval].insert(colInterval);
                        }
                    }
                }
            }
        }
    }

    /*
     * The combination of
     *  - distributed meshes
     *  - face coupling
     *  - DoFs associated with Faces, Edges or Nodes
     *  - Local computation of the sparsity estimate.
     * will form a problem.
     *
     * As example take the simple 1D mesh consisting of 3 elements and 4 face
     * (marked |).
     *
     * Mesh:      |-- e0 --|-- e1 --|-- e2 --|
     *            f0      f1        f2       f3
     * Owner:     0   0    0    0   1   1    1
     *
     * In a distributed mesh we only have a single layer of boundary elements,
     * so processor 1 effectively sees:
     * Mesh:      |-- e1 --|-- e2 --|
     *           f1        f2       f3
     * Owner:     0    0   1   1    1
     *
     * Let us consider the extend of a DoF that is associated with face f2 and
     * how it would couple to DoFs from other elements and faces:
     *  - Through the e1-element matrix it would couple to e1,     f1 and f2
     *  - Through the e2-element matrix it would couple to e2,     f2 and f3
     *  - Through the f1-face    matrix it would couple to e0, e1, f0, f1 and f2
     *  - Through the f2-face    matrix it would couple to e1, e2, f1, f2 and f3
     *  - Through the f3-face    matrix it would couple to e2,     f2 and f3
     *
     * Note that in the f1-face matrix there is coupling to unknown (for
     * processor 1) parts of the mesh, and this is to both an element (e0) and
     * face (f0).
     *
     * This restricts the options when constructing the sparsity estimate.
     *  - For rows with DoFs associated with Faces, Edges and Nodes there can
     *    not be any face coupling. As we miss information about the existence
     *    of elements, faces, edges and nodes on the far side of the ghost
     *    layer.
     *  - For rows with DoFs associated with Elements there is no restriction,
     *    they can be coupled to column DoFs, independent of the association of
     *    those column DoFs. As a face matrix for such a row-DoF could couple it
     *    to the column-DoFs on an element in the ghost layer. But we do have
     *    all the information about the DoFs on the ghost layer.
     *
     * At this moment I see two options of lifting this restriction:
     *  - Including a second layer of ghost elements.
     *  - Adding MPI communication to supply the necessary information.
     * Both of these require significant extra work, and we do not expect the
     * use of face matrices with Face/Edge/Node DoFs.
     */
    if (distributedMesh) {
        for (std::size_t i = 0; i < rowIndexing_.getNumberOfIncludedUnknowns();
             ++i) {
            std::size_t rowUnknown = rowIndexing_.getIncludedUnknowns()[i];
            for (std::size_t j = 0;
                 j < columnIndexing_.getNumberOfIncludedUnknowns(); ++j) {
                std::size_t colUnknown =
                    columnIndexing_.getIncludedUnknowns()[j];
                if (faceCoupling(i, j) && distributedMesh &&
                    !rowDoFs.isOnlyElementBased(i)) {
                    logger(ERROR,
                           "Face coupling of a row-DoF from an Face, Edge or "
                           "Node.");
                }
            }
        }
    }
    // Write the result
    for (auto const& o : owned) {
        const Detail::DoFInterval& rowInterval = o.first;
        std::size_t nonZero = 0;
        for (auto const& columns : o.second) {
            // Add the size of the DoFInterval, as we don't care about its
            // offset.
            nonZero += columns.second;
        }
        for (std::size_t ri = rowInterval.first;
             ri < rowInterval.first + rowInterval.second; ++ri) {
            nonZeroPerRowOwned[ri] = nonZero;
        }
    }
    // Same for the non owned part
    for (auto const& no : nonOwned) {
        const Detail::DoFInterval& rowInterval = no.first;
        std::size_t nonZero = 0;
        for (auto const& columns : no.second) {
            nonZero += columns.second;
        }
        for (std::size_t ri = rowInterval.first;
             ri < rowInterval.first + rowInterval.second; ++ri) {
            nonZeroPerRowNonOwned[ri] = nonZero;
        }
    }
}

}  // namespace Utilities
}  // namespace hpgem

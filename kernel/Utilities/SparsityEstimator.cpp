
#include "SparsityEstimator.h"

#include "Logger.h"

#include "map"
#include "Base/Element.h"
#include "Base/MeshManipulatorBase.h"
#include "Utilities/GlobalIndexing.h"
#include "Utilities/Table2D.h"

namespace hpgem {

namespace Utilities {
SparsityEstimator::SparsityEstimator(const GlobalIndexing& rowindexing,
                                     const GlobalIndexing& columnindexing)
    : rowIndexing_(rowindexing), columnIndexing_(columnindexing) {}

struct DoFRanges {
    /// Storage (unknown, dof, size) that the global indices (dof, dof+size-1)
    /// belong to unknown. A map is used for the inner storage to automatically
    /// deduplicate when the DoF indices for a Face/Edge/Node are added multiple
    /// times.
    std::vector<std::map<std::size_t, std::size_t>> dofs_;
};

struct SparsityEstimator::Workspace {

   private:
    // When counting DoFs in the columns that may interact with the DoFs of the
    // rows we need to group them based on two different properties:
    //  - Whether the pair of DoFs is part of element matrix or only of a face
    //    matrix. This is the same question as asking whether the support of the
    //    DoFs for the row and column share an element.
    //  - Whether the DoF of the column is owned by the current processor or
    //    not. This impacts whether it is in the owned or non owned part of the
    //    sparsity estimate.

    /// DoFs that have a shared supporting element and are locally owned.
    DoFRanges elementOwned;
    /// DoFs that have a shared supporting element and are not locally owned.
    DoFRanges elementNonOwned;
    /// DoFs without a shared supporting element and are locally owned.
    DoFRanges faceOwned;
    /// DoFs without a shared supporting element and are not locally owned.
    DoFRanges faceNonOwned;

    std::vector<bool> conformingColumns_;
    std::vector<bool> conformingRows_;

    /// The elements on which the DoF of the row has support.
    std::set<const Base::Element*> supportElements;
    /// The faces on which the DoF of the row has support.
    std::set<const Base::Face*> supportFaces;
    /// The edges on which the DoF of the row has support.
    std::set<const Base::Edge*> supportEdges;
    /// The nodes on which the DoF of the row has support.
    std::set<const Base::Node*> supportNodes;

    std::size_t dimension_;

   public:
    void init(std::size_t maxUnknown, std::size_t dimension) {
        dimension_ = dimension;
        // Ensure that we can directly address via the unknown of the column
        // index. The entries for unknowns not included in the column indexing
        // will not be used.
        elementOwned.dofs_.resize(maxUnknown);
        elementNonOwned.dofs_.resize(maxUnknown);
        faceOwned.dofs_.resize(maxUnknown);
        faceNonOwned.dofs_.resize(maxUnknown);
        // For checking which unknowns are conforming
        conformingColumns_.resize(maxUnknown);
        conformingRows_.resize(maxUnknown);
    }

    void clear() {
        // Clear the counted dofs, but leave the vectors intact as those
        // correspond to the number of unknowns.
        for (std::size_t i = 0; i < elementOwned.dofs_.size(); ++i) {
            elementOwned.dofs_[i].clear();
            elementNonOwned.dofs_[i].clear();
            faceOwned.dofs_[i].clear();
            faceNonOwned.dofs_[i].clear();
        }
        // Clear the support
        supportElements.clear();
        supportFaces.clear();
        supportEdges.clear();
        supportNodes.clear();
        // Not clearing the conforming property
    }

    void addSupport(const Base::Element* element) {
        supportElements.emplace(element);
        supportFaces.insert(element->getFacesList().begin(),
                            element->getFacesList().end());
        supportEdges.insert(element->getEdgesList().begin(),
                            element->getEdgesList().end());
        if (dimension_ > 1) {
            supportNodes.insert(element->getNodesList().begin(),
                                element->getNodesList().end());
        }
    }

    DoFRanges& getDoFRanges(const Base::Element* element) {
        bool owned = element->isOwnedByCurrentProcessor();
        bool support = supportElements.find(element) != supportElements.end();
        return get(support, owned);
    }

    DoFRanges& getDoFRanges(const Base::Face* face) {
        bool owned = face->isOwnedByCurrentProcessor();
        bool support = supportFaces.find(face) != supportFaces.end();
        return get(support, owned);
    }

    DoFRanges& getDoFRanges(const Base::Edge* edge) {
        bool owned = edge->isOwnedByCurrentProcessor();
        bool support = supportEdges.find(edge) != supportEdges.end();
        return get(support, owned);
    }

    DoFRanges& getDoFRanges(const Base::Node* node) {
        bool owned = node->isOwnedByCurrentProcessor();
        bool support = supportNodes.find(node) != supportNodes.end();
        return get(support, owned);
    }

    DoFRanges& get(bool support, bool owned) {
        if (support) {
            return owned ? elementOwned : elementNonOwned;
        } else {
            return owned ? faceOwned : faceNonOwned;
        }
    }

    const DoFRanges& get(bool support, bool owned) const {
        if (support) {
            return owned ? elementOwned : elementNonOwned;
        } else {
            return owned ? faceOwned : faceNonOwned;
        }
    }

    std::vector<bool>& getConformingRows() { return conformingRows_; }

    void setConformingColumn(std::size_t unknown) {
        conformingColumns_[unknown] = true;
    }

    std::vector<bool>& getConformingColumns() { return conformingColumns_; }
};

template <typename GEOM>
void SparsityEstimator::writeDoFCount(const GEOM* geom,
                                      const Workspace& workspace,
                                      std::vector<int>& nonZeroPerRowOwned,
                                      std::vector<int>& nonZeroPerRowNonOwned,
                                      const Table2D<bool>& faceCoupling) const {
    logger.assert_debug(
        geom->isOwnedByCurrentProcessor(),
        "Can not assign sparsity estimate for nonowned geometrical parts");
    for (std::size_t i = 0; i < rowIndexing_.getNumberOfIncludedUnknowns();
         ++i) {
        std::size_t rowUnknown = rowIndexing_.getIncludedUnknowns()[i];
        // Compute the number of entries
        std::size_t ownedSize = 0;
        std::size_t notOwnedSize = 0;

        for (std::size_t j = 0;
             j < columnIndexing_.getNumberOfIncludedUnknowns(); ++j) {
            std::size_t colUnknown = columnIndexing_.getIncludedUnknowns()[j];
            // Those with overlapping elements
            for (auto iter : workspace.get(true, true).dofs_[colUnknown]) {
                ownedSize += iter.second;
            }
            for (auto iter : workspace.get(true, false).dofs_[colUnknown]) {
                notOwnedSize += iter.second;
            }
            // Add the face coupling
            if (faceCoupling(i, j)) {
                for (auto iter : workspace.get(false, true).dofs_[colUnknown]) {
                    ownedSize += iter.second;
                }
                for (auto iter :
                     workspace.get(false, false).dofs_[colUnknown]) {
                    notOwnedSize += iter.second;
                }
                // TODO: Check for conforming-conforming coupling
            }
        }

        // Local placement
        std::size_t offset =
            rowIndexing_.getProcessorLocalIndex(geom, rowUnknown);
        std::size_t nbasis = geom->getLocalNumberOfBasisFunctions(rowUnknown);
        logger.assert_debug(offset + nbasis <= nonZeroPerRowOwned.size(),
                            "Indexing error index % is larger than size %",
                            offset + nbasis, nonZeroPerRowOwned.size());
        logger.assert_debug(offset + nbasis <= nonZeroPerRowNonOwned.size(),
                            "Indexing error index % is larger than size %",
                            offset + nbasis, nonZeroPerRowNonOwned.size());
        for (std::size_t b = 0; b < nbasis; ++b) {
            nonZeroPerRowOwned[b + offset] = ownedSize;
            nonZeroPerRowNonOwned[b + offset] = notOwnedSize;
        }
    }
}

bool hasLocalDoFs(const Base::Element* element,
                  const std::vector<std::size_t>& unknowns) {
    for (const std::size_t& unknown : unknowns) {
        if (element->getLocalNumberOfBasisFunctions(unknown) > 0) return true;
    }
    return false;
}

/// Test whether a geometrical object has DoF associated with it.
///
/// \tparam GEOM The type of geometry
/// \param geom The geometrical object
/// \param unknowns The unknows to check
/// \param conforming A vector to mark the conforming dofs
/// \return  Whether any of these unknowns has a DoF for the object.
template <typename GEOM>
bool hasLocalDoFs2(const GEOM* geom, const std::vector<std::size_t>& unknowns,
                   std::vector<bool>& conforming) {
    bool result = false;
    for (const std::size_t& unknown : unknowns) {
        if (geom->getLocalNumberOfBasisFunctions(unknown) > 0) {
            result = true;
            conforming[unknown] = true;
        }
    }
    return result;
}

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
    if (rowIndexing_.getMesh() == nullptr) {
        nonZeroPerRowOwned.clear();
        nonZeroPerRowOwned.resize(totalNumberOfDoF, 0);
        nonZeroPerRowNonOwned.clear();
        nonZeroPerRowNonOwned.resize(totalNumberOfDoF, 0);
        return;
    }

    Workspace workspace;
    workspace.init(columnIndexing_.getTotalNumberOfUnknowns(),
                   rowIndexing_.getMesh()->dimension());

    // Resize and initialize with error data
    nonZeroPerRowOwned.assign(totalNumberOfDoF, -1);
    nonZeroPerRowNonOwned.assign(totalNumberOfDoF, -1);

    // If there is no face coupling anywhere (e.g. dg mass matrix) then we can
    // skip looping over the adjacent elements.
    bool anyFaceCoupling = false;
    for (std::size_t i = 0; i < faceCoupling.getSize(); ++i) {
        if (faceCoupling[i]) {
            anyFaceCoupling = true;
            break;
        }
    }

    bool distributed = false;

    // To compute the number of non zero's in each row we assume
    // (pessimistically) that all face matrices are filled. For each DoF we then
    // compute the set of other DoFs that contribute to at least one face
    // matrix, and counting these gives the number of non zeros in the row.
    //
    // Counting DoFs is done by going over all elements where a basis function
    // has support. All basis functions with support on these elements can
    // result in non zero element matrix entries. For non zero face matrix
    // entries we also need to consider all basis functions that have support on
    // the element adjacent to the elements where the basis function has
    // support.
    //
    // For efficiency purposes we do not traverse each basis function
    // separately, instead we compute it for all basis functions of a single
    // element/face/etc.

    for (const Base::Element* element :
         rowIndexing_.getMesh()->getElementsList()) {
        if (!hasLocalDoFs(element, rowIndexing_.getIncludedUnknowns())) {
            // If there are no DoFs associated with this element, then then
            // there are no rows in the matrix for this element and there is
            // nothing to count.
            continue;
        }
        workspace.clear();
        workspace.addSupport(element);
        // Add local basis functions
        addElementDoFs(element, workspace);

        // Face matrix coupling with the element on the other side of the face
        if (anyFaceCoupling) {
            for (const Base::Face* face : element->getFacesList()) {
                if (!face->isInternal()) continue;
                const Base::Element* other = face->getPtrOtherElement(element);
                addElementDoFs(other, workspace);
            }
        }

        // With all the global indices counted, allocate them.
        writeDoFCount(element, workspace, nonZeroPerRowOwned,
                      nonZeroPerRowNonOwned, faceCoupling);
    }
    logger(VERBOSE, "Sparsity pattern for Element DoFs computed");
    // Faces
    for (const Base::Face* face : rowIndexing_.getMesh()->getFacesList()) {
        // Safe way to see if the mesh is distributed.
        if (face->getFaceType() == Geometry::FaceType::SUBDOMAIN_BOUNDARY ||
            face->getFaceType() == Geometry::FaceType::PERIODIC_SUBDOMAIN_BC) {
            distributed = true;
        }
        if (!face->isOwnedByCurrentProcessor() ||
            !hasLocalDoFs2(face, rowIndexing_.getIncludedUnknowns(),
                           workspace.getConformingRows())) {
            // If the face is not owned by the current processor, then it
            // does not have any rows in the matrix on this processor. The
            // number of non zero entries in the corresponding row does not
            // need to be estimated on this processor.
            continue;
        }
        workspace.clear();
        std::vector<const Base::Element*> faceElements = {
            face->getPtrElementLeft()};
        if (face->isInternal())
            faceElements.push_back(face->getPtrElementRight());

        for (const Base::Element* element : faceElements) {
            workspace.addSupport(element);
        }
        for (const Base::Element* element : faceElements) {
            // DoFs from supporting element
            addElementDoFs(element, workspace);
            // Face matrix coupling with the element on the other side of the
            // face
            if (!anyFaceCoupling) continue;
            for (const Base::Face* otherFace : element->getFacesList()) {
                if (otherFace == face)
                    continue;  // Both left & right element of the original face
                               // are already visited
                if (!otherFace->isInternal()) continue;
                const Base::Element* otherElement =
                    otherFace->getPtrOtherElement(element);
                // The basis function from the face may couple through any other
                // face of the elements adjacent to the face
                addElementDoFs(otherElement, workspace);
            }
        }
        // Store the information
        writeDoFCount(face, workspace, nonZeroPerRowOwned,
                      nonZeroPerRowNonOwned, faceCoupling);
    }
    logger(VERBOSE, "Sparsity pattern for Face DoFs computed");
    // Edges
    for (const Base::Edge* edge : rowIndexing_.getMesh()->getEdgesList()) {
        if (!edge->isOwnedByCurrentProcessor() ||
            !hasLocalDoFs2(edge, rowIndexing_.getIncludedUnknowns(),
                           workspace.getConformingRows())) {
            continue;
        }
        workspace.clear();
        for (const Base::Element* element : edge->getElements()) {
            workspace.addSupport(element);
        }
        for (const Base::Element* element : edge->getElements()) {
            // DoFs from supporting element
            addElementDoFs(element, workspace);
            // Face matrix coupling with the element on the other side of the
            // face
            if (!anyFaceCoupling) continue;
            for (const Base::Face* face : element->getFacesList()) {
                if (!face->isInternal()) continue;
                const Base::Element* otherElement =
                    face->getPtrOtherElement(element);
                addElementDoFs(otherElement, workspace);
            }
        }
        // Store the information
        writeDoFCount(edge, workspace, nonZeroPerRowOwned,
                      nonZeroPerRowNonOwned, faceCoupling);
    }
    logger(VERBOSE, "Sparsity pattern for Edge DoFs computed");
    // Nodes
    if (rowIndexing_.getMesh()->dimension() > 1) {
        for (const Base::Node* node : rowIndexing_.getMesh()->getNodesList()) {
            if (!node->isOwnedByCurrentProcessor() ||
                !hasLocalDoFs2(node, rowIndexing_.getIncludedUnknowns(),
                               workspace.getConformingRows())) {
                continue;
            }
            workspace.clear();
            for (const Base::Element* element : node->getElements()) {
                workspace.addSupport(element);
            }
            for (const Base::Element* element : node->getElements()) {
                // DoFs from the supporting element
                addElementDoFs(element, workspace);
                // Face matrix coupling with the element on the other side of
                // the face
                if (!anyFaceCoupling) continue;
                for (const Base::Face* face : element->getFacesList()) {
                    if (!face->isInternal()) continue;
                    const Base::Element* otherElement =
                        face->getPtrOtherElement(element);
                    addElementDoFs(otherElement, workspace);
                }
            }
            // Store the information
            writeDoFCount(node, workspace, nonZeroPerRowOwned,
                          nonZeroPerRowNonOwned, faceCoupling);
        }
        logger(VERBOSE, "Sparsity pattern for Node DoFs computed");
    }

    // The distributed mesh presents a problem when computing the sparsity
    // estimate for conforming DoFs that couple to another conforming DoF
    // through a face matrix. This problem is caused by that we only have
    // a single layer of ghost elements (thus ghost elements share at least a
    // node with an owned element).
    //
    // Consider now a conforming DoF on the subdomain boundary. Without loss of
    // generality we may assume that we own it (otherwise repeat the argument
    // for the owning processor). This DoF has support on some elements of the
    // layer of ghost elements. For this layer we have all the information and
    // we can thus compute the sparsity estimate for the corresponding element
    // matrices. For face matrices it will couple through all the faces of the
    // elements the DoF has support on. This includes all the faces of the ghost
    // elements, which will include faces that couple to elements that are
    // beyond the first layer of ghost elements. As we do not know of these
    // elements, we can not correctly compute the sparsity estimate of such a
    // face matrix.
    //
    // The only feasible way I can currently see around this is to add
    // another layer of ghost cells. Before this is added we just warn about
    // this incompatibility.
    for (std::size_t i = 0; i < rowIndexing_.getNumberOfIncludedUnknowns();
         ++i) {
        std::size_t rowUnknown = rowIndexing_.getIncludedUnknowns()[i];
        for (std::size_t j = 0;
             j < columnIndexing_.getNumberOfIncludedUnknowns(); ++j) {
            std::size_t colUnknown = columnIndexing_.getIncludedUnknowns()[j];
            if (faceCoupling(i, j) && distributed &&
                workspace.getConformingRows()[rowUnknown] &&
                workspace.getConformingColumns()[colUnknown]) {
                logger(ERROR,
                       "Face coupling between conforming DoFs on distributed "
                       "mesh is not supported");
            }
        }
    }
}

void SparsityEstimator::addElementDoFs(const Base::Element* element,
                                       Workspace& workspace) const {
    {
        // For the element
        DoFRanges& toChange = workspace.getDoFRanges(element);
        for (std::size_t unknown : columnIndexing_.getIncludedUnknowns()) {
            const std::size_t offset =
                columnIndexing_.getGlobalIndex(element, unknown);
            const std::size_t localBasisFunctions =
                element->getLocalNumberOfBasisFunctions(unknown);
            if (localBasisFunctions > 0) {
                toChange.dofs_[unknown][offset] = localBasisFunctions;
            }
        }
    }
    // For the faces
    for (const Base::Face* face : element->getFacesList()) {
        DoFRanges& toChange = workspace.getDoFRanges(face);
        for (std::size_t unknown : columnIndexing_.getIncludedUnknowns()) {
            const std::size_t offset =
                columnIndexing_.getGlobalIndex(face, unknown);
            const std::size_t localBasisFunctions =
                face->getLocalNumberOfBasisFunctions(unknown);
            if (localBasisFunctions > 0) {
                toChange.dofs_[unknown][offset] = localBasisFunctions;
                workspace.setConformingColumn(unknown);
            }
        }
    }
    // For the edges
    for (const Base::Edge* edge : element->getEdgesList()) {
        DoFRanges& toChange = workspace.getDoFRanges(edge);
        for (std::size_t unknown : columnIndexing_.getIncludedUnknowns()) {
            const std::size_t offset =
                columnIndexing_.getGlobalIndex(edge, unknown);
            const std::size_t localBasisFunctions =
                edge->getLocalNumberOfBasisFunctions(unknown);
            if (localBasisFunctions > 0) {
                toChange.dofs_[unknown][offset] = localBasisFunctions;
                workspace.setConformingColumn(unknown);
            }
        }
    }
    // For the nodes
    if (rowIndexing_.getMesh()->dimension() > 1) {
        // For the edges
        for (const Base::Node* node : element->getNodesList()) {
            DoFRanges& toChange = workspace.getDoFRanges(node);
            for (std::size_t unknown : columnIndexing_.getIncludedUnknowns()) {
                const std::size_t offset =
                    columnIndexing_.getGlobalIndex(node, unknown);
                const std::size_t localBasisFunctions =
                    node->getLocalNumberOfBasisFunctions(unknown);
                if (localBasisFunctions > 0) {
                    toChange.dofs_[unknown][offset] = localBasisFunctions;
                    workspace.setConformingColumn(unknown);
                }
            }
        }
    }
}
}  // namespace Utilities
}  // namespace hpgem

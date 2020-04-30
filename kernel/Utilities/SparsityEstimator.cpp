
#include "SparsityEstimator.h"

#include "Logger.h"

#include "map"
#include "Base/Element.h"
#include "Base/MeshManipulatorBase.h"
#include "Utilities/GlobalIndexing.h"

namespace Utilities
{
    SparsityEstimator::SparsityEstimator(
            const GlobalIndexing &rowindexing, const GlobalIndexing &columnindexing)
        : rowIndexing_ (rowindexing)
        , columnIndexing_ (columnindexing)
    {}

    struct SparsityEstimator::Workspace
    {
        /// DoFs owned by the current processor. Stored as (dof, size) where dof is the
        /// first dof of an element, face, etc.
        std::map<std::size_t, std::size_t> owned;
        // Same as #owned, but for the DoFs not owned by the current processor
        std::map<std::size_t, std::size_t> notOwned;

        /// Get a reference to either the owned or notOwned DoF set.
        std::map<std::size_t, std::size_t>& getDoFs(bool owning)
        {
            return owning ? owned : notOwned;
        }

        void clear()
        {
            owned.clear();
            notOwned.clear();
        }
    };

    template<typename GEOM>
    void SparsityEstimator::writeDoFCount(
            const GEOM* geom, const Workspace &workspace,
            std::vector<int> &nonZeroPerRowOwned, std::vector<int> &nonZeroPerRowNonOwned) const
    {
        logger.assert_debug(geom->isOwnedByCurrentProcessor(),
                "Can not assign sparsity estimate for nonowned geometrical parts");
        std::size_t ownedSize = 0;
        std::size_t notOwnedSize = 0;
        for (auto iter : workspace.owned)
            ownedSize += iter.second;
        for (auto iter : workspace.notOwned)
            notOwnedSize += iter.second;

        for (std::size_t unknown : rowIndexing_.getIncludedUnknowns())
        {
            std::size_t offset = rowIndexing_.getProcessorLocalIndex(geom, unknown);
            std::size_t nbasis = geom->getLocalNumberOfBasisFunctions(unknown);

            logger.assert_debug(offset + nbasis <= nonZeroPerRowOwned.size(),
                    "Indexing error index % is larger than size %",
                    offset + nbasis, nonZeroPerRowOwned.size());
            logger.assert_debug(offset + nbasis <= nonZeroPerRowNonOwned.size(),
                    "Indexing error index % is larger than size %",
                    offset + nbasis, nonZeroPerRowNonOwned.size());

            for (std::size_t i = 0; i < nbasis; ++i)
            {
                nonZeroPerRowOwned[i + offset] = ownedSize;
                nonZeroPerRowNonOwned[i + offset] = notOwnedSize;
            }
        }
    }

    /// Test whether a geometrical object has DoF associated with it.
    ///
    /// \tparam GEOM The type of geometry
    /// \param geom The geometrical object
    /// \param unknowns The unknows to check
    /// \return  Whether any of these unknowns has a DoF for the object.
    template<typename GEOM>
    bool hasLocalDoFs(const GEOM *geom, const std::vector<std::size_t>& unknowns)
    {
        for (const std::size_t& unknown : unknowns)
        {
            if (geom->getLocalNumberOfBasisFunctions(unknown > 0))
                return true;
        }
        return false;
    }

    void SparsityEstimator::computeSparsityEstimate(
            std::vector<int> &nonZeroPerRowOwned,
            std::vector<int> &nonZeroPerRowNonOwned,
            bool includeFaceCoupling) const
    {
        logger(VERBOSE, "Computing sparsity estimate for mesh %", rowIndexing_.getMesh());
        const std::size_t totalNumberOfDoF = rowIndexing_.getNumberOfLocalBasisFunctions();
        if (rowIndexing_.getMesh() == nullptr)
        {
            nonZeroPerRowOwned.clear();
            nonZeroPerRowOwned.resize(totalNumberOfDoF, 0);
            nonZeroPerRowNonOwned.clear();
            nonZeroPerRowNonOwned.resize(totalNumberOfDoF, 0);
            return;
        }

        Workspace workspace;

        // Resize and initialize with error data
        nonZeroPerRowOwned.assign(totalNumberOfDoF, -1);
        nonZeroPerRowNonOwned.assign(totalNumberOfDoF, -1);

        // To compute the number of non zero's in each row we assume (pessimistically)
        // that all face matrices are filled. For each DoF we then compute the set of
        // other DoFs that contribute to at least one face matrix, and counting these
        // gives the number of non zeros in the row.
        //
        // Counting DoFs is done by going over all elements where a basis function has
        // support. All basis functions with support on these elements can result in non
        // zero element matrix entries. For non zero face matrix entries we also need to
        // consider all basis functions that have support on the element adjacent to the
        // elements where the basis function has support.
        //
        // For efficiency purposes we do not traverse each basis function separately,
        // instead we compute it for all basis functions of a single element/face/etc.

        for (const Base::Element* element : rowIndexing_.getMesh()->getElementsList())
        {
            if (!hasLocalDoFs(element, rowIndexing_.getIncludedUnknowns()))
            {
                // If there are no DoFs associated with this element, then then
                // there are no rows in the matrix for this element and there is
                // nothing to count.
                continue;
            }
            workspace.clear();
            // Add local basis functions
            addElementDoFs(element, workspace);

            // Face matrix coupling with the element on the other side of the face
            if (includeFaceCoupling)
            {
                for (const Base::Face *face : element->getFacesList())
                {
                    if (!face->isInternal())
                        continue;
                    const Base::Element *other = face->getPtrOtherElement(element);
                    addElementDoFs(other, workspace);
                }
            }

            // With all the global indices counted, allocate them.
            writeDoFCount(element, workspace, nonZeroPerRowOwned, nonZeroPerRowNonOwned);
        }
        logger(VERBOSE, "Sparsity pattern for Element DoFs computed");
        // Faces
        for (const Base::Face* face : rowIndexing_.getMesh()->getFacesList())
        {
            if (!face->isOwnedByCurrentProcessor() || !hasLocalDoFs(face, rowIndexing_.getIncludedUnknowns()))
            {
                // If the face is not owned by the current processor, then it
                // does not have any rows in the matrix on this processor. The
                // number of non zero entries in the corresponding row does not
                // need to be estimated on this processor.
                continue;
            }
            workspace.clear();
            std::vector<const Base::Element*> faceElements = {face->getPtrElementLeft()};
            if (face->isInternal())
                faceElements.push_back(face->getPtrElementRight());
            for (const Base::Element* element : faceElements)
            {
                // DoFs from supporting element
                addElementDoFs(element, workspace);
                // Face matrix coupling with the element on the other side of the face
                if (!includeFaceCoupling)
                    continue;
                // Face coupling with conforming elements can couple DoFs that
                // only have shared support on the face between two elements.
                // For a DoF that is defined on a subdomain boundary face, it
                // has support on a ghost cell and all the faces of the ghost
                // cell. As there is only a single layer of ghost cells, we do
                // in general not have the element on the other side of the
                // faces of this ghost cell. Therefore, it is not possible to
                // know how many basis functions have support on this face.
                //
                // The only way I can currently see around this is to add
                // another layer of ghost cells.
                logger.assert_always(element->isOwnedByCurrentProcessor(),
                        "Face coupling with conforming basis functions on distributed meshes is not supported");
                for (const Base::Face* otherFace : element->getFacesList())
                {
                    if (otherFace == face)
                        continue; // Both left & right element of the original face are already visited
                    if (!otherFace->isInternal())
                        continue;
                    const Base::Element* otherElement = otherFace->getPtrOtherElement(element);
                    // The basis function from the face may couple through any other face of the elements adjacent to
                    // the face
                    addElementDoFs(otherElement, workspace);
                }
            }
            // Store the information
            writeDoFCount(face, workspace, nonZeroPerRowOwned, nonZeroPerRowNonOwned);
        }
        logger(VERBOSE, "Sparsity pattern for Face DoFs computed");
        // Edges
        for (const Base::Edge* edge : rowIndexing_.getMesh()->getEdgesList())
        {
            if (!edge->isOwnedByCurrentProcessor() || !hasLocalDoFs(edge, rowIndexing_.getIncludedUnknowns()))
            {
                continue;
            }
            workspace.clear();
            for (const Base::Element* element : edge->getElements())
            {
                // DoFs from supporting element
                addElementDoFs(element, workspace);
                // Face matrix coupling with the element on the other side of the face
                if (!includeFaceCoupling)
                    continue;
                // See remark for face based DoFs
                logger.assert_always(element->isOwnedByCurrentProcessor(),
                        "Face coupling with conforming basis functions on distributed meshes is not supported");
                for (const Base::Face* face : element->getFacesList())
                {
                    if (!face->isInternal())
                        continue;
                    const Base::Element* otherElement = face->getPtrOtherElement(element);
                    addElementDoFs(otherElement, workspace);
                }
            }
            // Store the information
            writeDoFCount(edge, workspace, nonZeroPerRowOwned, nonZeroPerRowNonOwned);
        }
        logger(VERBOSE, "Sparsity pattern for Edge DoFs computed");
        // Nodes
        if (rowIndexing_.getMesh()->dimension() > 1)
        {
            for (const Base::Node* node : rowIndexing_.getMesh()->getNodesList())
            {
                if (!node->isOwnedByCurrentProcessor() || !hasLocalDoFs(node, rowIndexing_.getIncludedUnknowns()))
                {
                    continue;
                }
                workspace.clear();
                for (const Base::Element* element : node->getElements())
                {
                    // DoFs from the supporting element
                    addElementDoFs(element, workspace);
                    // Face matrix coupling with the element on the other side of the face
                    if (!includeFaceCoupling)
                        continue;
                    // See remark for face based DoFs
                    logger.assert_always(element->isOwnedByCurrentProcessor(),
                            "Face coupling with conforming basis functions on distributed meshes is not supported");
                    for (const Base::Face* face : element->getFacesList())
                    {
                        if (!face->isInternal())
                            continue;
                        const Base::Element* otherElement = face->getPtrOtherElement(element);
                        addElementDoFs(otherElement, workspace);
                    }
                }
                // Store the information
                writeDoFCount(node, workspace, nonZeroPerRowOwned, nonZeroPerRowNonOwned);
            }
            logger(VERBOSE, "Sparsity pattern for Node DoFs computed");
        }
    }

    void SparsityEstimator::addElementDoFs(
            const Base::Element *element,
            Workspace &workspace) const
    {
        {
            // For the element
            std::map<std::size_t, std::size_t> &toChange = workspace.getDoFs(element->isOwnedByCurrentProcessor());
            for (std::size_t unknown : columnIndexing_.getIncludedUnknowns())
            {
                const std::size_t offset = columnIndexing_.getGlobalIndex(element, unknown);
                const std::size_t localBasisFunctions = element->getLocalNumberOfBasisFunctions(unknown);
                if (localBasisFunctions > 0)
                {
                    toChange[offset] = localBasisFunctions;
                }
            }
        }
        // For the faces
        for (const Base::Face* face : element->getFacesList())
        {
            std::map<std::size_t, std::size_t>& toChange = workspace.getDoFs(face->isOwnedByCurrentProcessor());
            for (std::size_t unknown : columnIndexing_.getIncludedUnknowns())
            {
                const std::size_t offset = columnIndexing_.getGlobalIndex(face, unknown);
                const std::size_t localBasisFunctions = face->getLocalNumberOfBasisFunctions(unknown);
                if (localBasisFunctions > 0)
                {
                    toChange[offset] = localBasisFunctions;
                }
            }
        }
        // For the edges
        for (const Base::Edge* edge : element->getEdgesList())
        {
            std::map<std::size_t, std::size_t>& toChange = workspace.getDoFs(edge->isOwnedByCurrentProcessor());
            for (std::size_t unknown : columnIndexing_.getIncludedUnknowns())
            {
                const std::size_t offset = columnIndexing_.getGlobalIndex(edge, unknown);
                const std::size_t localBasisFunctions = edge->getLocalNumberOfBasisFunctions(unknown);
                if (localBasisFunctions > 0)
                {
                    toChange[offset] = localBasisFunctions;
                }
            }
        }
        // For the nodes
        if (rowIndexing_.getMesh()->dimension() > 1)
        {
            // For the edges
            for (const Base::Node* node : element->getNodesList())
            {
                std::map<std::size_t, std::size_t>& toChange = workspace.getDoFs(node->isOwnedByCurrentProcessor());
                for (std::size_t unknown : columnIndexing_.getIncludedUnknowns())
                {
                    const std::size_t offset = columnIndexing_.getGlobalIndex(node, unknown);
                    const std::size_t localBasisFunctions = node->getLocalNumberOfBasisFunctions(unknown);
                    if (localBasisFunctions > 0)
                    {
                        toChange[offset] = localBasisFunctions;
                    }
                }
            }
        }
    }
}
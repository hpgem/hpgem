
#include "SparsityEstimator.h"

#include "Logger.h"

#include "Base/Element.h"
#include "Base/MeshManipulatorBase.h"
#include "Utilities/GlobalIndexing.h"

namespace Utilities
{
    SparsityEstimator::SparsityEstimator(const Base::MeshManipulatorBase &mesh, const GlobalIndexing &indexing)
        : mesh_ (mesh)
        , indexing_ (indexing)
    {}

    struct SparsityEstimator::Workspace
    {
        /// DoFs owned by the current processor
        std::set<std::size_t> owned;
        // DoFs not owned by the current processor
        std::set<std::size_t> notOwned;

        // Number of unknowns
        std::size_t nUnknowns;

        /// Get a reference to either the owned or notOwned DoF set.
        std::set<std::size_t>& getDoFs(bool owning)
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
        for (std::size_t unknown = 0; unknown < workspace.nUnknowns; ++unknown)
        {
            std::size_t offset = indexing_.getGlobalIndex(geom, unknown);
            std::size_t nbasis = geom->getLocalNumberOfBasisFunctions(unknown);
            for (std::size_t i = 0; i < nbasis; ++i)
            {
                nonZeroPerRowOwned[i + offset] = workspace.owned.size();
                nonZeroPerRowNonOwned[i + offset] = workspace.notOwned.size();
            }
        }
    }

    void SparsityEstimator::computeSparsityEstimate(
            std::vector<int> &nonZeroPerRowOwned,
            std::vector<int> &nonZeroPerRowNonOwned) const
    {
        const std::size_t totalNumberOfDoF = indexing_.getNumberOfLocalBasisFunctions();
        const std::size_t nUknowns = indexing_.getNumberOfUnknowns();
        Workspace workspace;
        workspace.nUnknowns = nUknowns;

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

        for (const Base::Element* element : mesh_.getElementsList())
        {
            workspace.clear();
            // Add local basis functions
            addElementDoFs(element, workspace);

            // Face matrix coupling with the element on the other side of the face
            for(const Base::Face* face : element->getFacesList())
            {
                if (!face->isInternal())
                    continue;
                const Base::Element* other = face->getPtrOtherElement(element);
                addElementDoFs(other, workspace);
            }

            // With all the global indices counted, allocate them.
            for (std::size_t unknown = 0; unknown < nUknowns; ++unknown)
            {
                writeDoFCount(element, workspace, nonZeroPerRowOwned, nonZeroPerRowNonOwned);
            }
        }
        // Faces
        for (const Base::Face* face : mesh_.getFacesList())
        {
            logger.assert_debug(face->isOwnedByCurrentProcessor(), "Face not owned by current processor");
            workspace.clear();
            std::vector<const Base::Element*> faceElements = {face->getPtrElementLeft()};
            if (face->isInternal())
                faceElements.push_back(face->getPtrElementRight());
            for (const Base::Element* element : faceElements)
            {
                // DoFs from supporting element
                addElementDoFs(element, workspace);
                // Face matrix coupling with the element on the other side of the face
                for (const Base::Face* otherFace : element->getFacesList())
                {
                    if (otherFace == face)
                        continue; // Both left & right element of the original face are already visited
                    if (!face->isInternal())
                        continue;
                    const Base::Element* otherElement = face->getPtrOtherElement(element);
                    // The basis function from the face may couple through any other face of the elements adjacent to
                    // the face
                    addElementDoFs(otherElement, workspace);
                }
            }
            // Store the information
            writeDoFCount(face, workspace, nonZeroPerRowOwned, nonZeroPerRowNonOwned);
        }
        // Edges
        for (const Base::Edge* edge : mesh_.getEdgesList())
        {
            logger.assert_debug(edge->isOwnedByCurrentProcessor(), "Non owned edge");
            workspace.clear();
            for (const Base::Element* element : edge->getElements())
            {
                // DoFs from supporting element
                addElementDoFs(element, workspace);
                // Face matrix coupling with the element on the other side of the face
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
        // Nodes
        if (mesh_.dimension() > 1)
        {
            for (const Base::Node* node : mesh_.getNodesList())
            {
                logger.assert_debug(node->isOwnedByCurrentProcessor(), "Non owned node");
                workspace.clear();
                for (const Base::Element* element : node->getElements())
                {
                    // DoFs from the supporting element
                    addElementDoFs(element, workspace);
                    // Face matrix coupling with the element on the other side of the face
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
        }
    }

    void SparsityEstimator::addElementDoFs(
            const Base::Element *element,
            Workspace &workspace) const
    {
        logger.assert_debug(element->isOwnedByCurrentProcessor(), "Element is not owned by processor");
        const std::size_t unknowns = element->getNumberOfUnknowns();
        {
            // For the element
            std::set<std::size_t> &toChange = workspace.getDoFs(element->isOwnedByCurrentProcessor());
            for (std::size_t unknown = 0; unknown < unknowns; ++unknown)
            {
                const std::size_t offset = indexing_.getGlobalIndex(element, unknown);
                const std::size_t localBasisFunctions = element->getLocalNumberOfBasisFunctions(unknown);
                for (std::size_t i = 0; i < localBasisFunctions; ++i)
                {
                    toChange.insert(offset + i);
                }
            }
        }
        // For the faces
        for (const Base::Face* face : element->getFacesList())
        {
            std::set<std::size_t>& toChange = workspace.getDoFs(face->isOwnedByCurrentProcessor());
            for (std::size_t unknown = 0; unknown < unknowns; ++unknown)
            {
                const std::size_t offset = indexing_.getGlobalIndex(face, unknown);
                const std::size_t localBasisFunctions = face->getLocalNumberOfBasisFunctions(unknown);
                for (std::size_t i = 0; i < localBasisFunctions; ++i)
                {
                    toChange.insert(offset + i);
                }
            }
        }
        // For the edges
        for (const Base::Edge* edge : element->getEdgesList())
        {
            std::set<std::size_t>& toChange = workspace.getDoFs(edge->isOwnedByCurrentProcessor());
            for (std::size_t unknown = 0; unknown < unknowns; ++unknown)
            {
                const std::size_t offset = indexing_.getGlobalIndex(edge, unknown);
                const std::size_t localBasisFunctions = edge->getLocalNumberOfBasisFunctions(unknown);
                for (std::size_t i = 0; i < localBasisFunctions; ++i)
                {
                    toChange.insert(offset + i);
                }
            }
        }
        // For the nodes
        if (mesh_.dimension() > 1)
        {
            // For the edges
            for (const Base::Node* node : element->getNodesList())
            {
                std::set<std::size_t>& toChange = workspace.getDoFs(node->isOwnedByCurrentProcessor());
                for (std::size_t unknown = 0; unknown < unknowns; ++unknown)
                {
                    const std::size_t offset = indexing_.getGlobalIndex(node, unknown);
                    const std::size_t localBasisFunctions = node->getLocalNumberOfBasisFunctions(unknown);
                    for (std::size_t i = 0; i < localBasisFunctions; ++i)
                    {
                        toChange.insert(offset + i);
                    }
                }
            }
        }
    }
}
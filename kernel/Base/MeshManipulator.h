/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2014, University of Twente
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

#ifndef HPGEM_KERNEL_MESHMANIPULATOR_H
#define HPGEM_KERNEL_MESHMANIPULATOR_H

#include <vector>
#include <fstream>

#include "Geometry/FaceGeometry.h"
#include "PhysicalElement.h"
#include "Mesh.h"
#include "GlobalNamespaceBase.h"
#include "BasisFunctionSet.h"
#include "MeshManipulatorBase.h"

namespace Base {
template <std::size_t DIM>
class MeshManipulator;
}

template <std::size_t DIM>
std::ostream& operator<<(std::ostream&, const Base::MeshManipulator<DIM>&);
namespace Geometry {
template <std::size_t DIM>
class PointPhysical;
template <std::size_t DIM>
class PointReference;
}  // namespace Geometry

namespace Base {
class BasisFunctionSet;
class OrientedBasisFunctionSet;
class Face;
template <std::size_t DIM>
class MeshMoverBase;
template <class V>
class LevelTree;
class Element;
struct ConfigurationData;
class Edge;

struct HalfFaceDescription {
    std::vector<std::size_t> nodesList;
    std::size_t elementNumber;
    std::size_t localFaceIndex;
};

template <std::size_t DIM>
class MeshManipulator : public MeshManipulatorBase {
   public:
    using CollectionOfBasisFunctionSets =
        Element::CollectionOfBasisFunctionSets;

    /*using ConstElementIterator = TreeIteratorConst<Element*>;
    using ElementIterator = TreeIterator<Element*>;

    using ConstFaceIterator = TreeIteratorConst<Face*>;
    using FaceIterator = TreeIterator<Face*>;*/

    // fixme: periodicity information is requested for legacy reasons, but
    // ignored in the current mesh reader
    MeshManipulator(const ConfigurationData* configData,
                    std::size_t numberOfElementMatrixes = 0,
                    std::size_t numberOfElementVectors = 0,
                    std::size_t numberOfFaceMatrixes = 0,
                    std::size_t numberOfFaceVectors = 0);

    MeshManipulator(const MeshManipulator& other);

    ~MeshManipulator() override;

    /// creates some cheap, easy to construct basis function set (monomials) to
    /// use as a placeholder Note that they usually result in very badly
    /// conditioned matrices.
    void useMonomialBasisFunctions(std::size_t order);

    /// \brief automatically creates the DG basis function set most appropriate
    /// for the shape of the element and sets that set as the basis function set
    /// to use \details This function takes the default conforming basis
    /// functions and cuts them off at element boundaries. The resulting basis
    /// functions have much better conditioning properties than monomials. Due
    /// to the nature of DG, this creates 'interior' basis functions only. (All
    /// basis functions are associated with the element) This function should be
    /// called after a mesh has been created to ensure basis functions exist for
    /// all types of elements needed. It is allowed to change the basis function
    /// during a computation, but you have to re-initialise the expansion
    /// coefficients, since they are expected to change such that the current
    /// solution stays the same. (Stored old solutions are not affected, but
    /// they may require a change back to the original basis functions before
    /// being usable again)
    void useDefaultDGBasisFunctions(std::size_t order);
    void useDefaultDGBasisFunctions(std::size_t order, std::size_t unknown);

    /// \brief automatically creates Nedelec DG basis functions for tetrahedra.
    /// \details This function should be called after a mesh has been created to
    /// ensure basis functions exist for all types of elements needed.
    void useNedelecDGBasisFunctions(std::size_t order);

    /// \brief automatically creates Ainsworth-Coyle DG basis functions for
    /// tetrahedra. \details This function should be called after a mesh has
    /// been created to ensure basis functions exist for all types of elements
    /// needed.
    void useAinsworthCoyleDGBasisFunctions(std::size_t order);

    /// \brief automatically creates conforming basis functions, even for mixed
    /// meshes \details For p=1, this creates a nodal basis function set
    /// associated with the mesh nodes. This function should be called after a
    /// mesh has been created to ensure basis functions exist for all types of
    /// elements needed. For higher polynomial orders it will add hierarchic
    /// basis functions associated with edges, faces and the interior of
    /// elements. For the basis functions to be truly conforming you will need
    /// some sort of global assembly structure that respects these associations
    /// (GlobalPETScMatrix does this). If you are building such an assembly
    /// structure yourself you may assume that the basisFunctions are indexed
    /// such that interior basis functions of the element are first. They are
    /// followed by the face basis functions(ordered by local face number
    /// primary and the different basis functions secondary) They are followed
    /// by the edge basis functions (using the same ordering) and finally the
    /// nodal basis functions (also using the same ordering) It is allowed to
    /// change the basis function during a computation, but you have to
    /// re-initialise the expansion coefficients, since they are expected to
    /// change such that the current solution stays the same. (Stored old
    /// solutions are not affected, but they may require a change back to the
    /// original basis functions before being usable again)
    void useDefaultConformingBasisFunctions(std::size_t order);
    void useDefaultConformingBasisFunctions(std::size_t order,
                                            std::size_t unknown);

    Element* addElement(const std::vector<std::size_t>& globalNodeIndexes,
                        std::size_t owner, bool owning);

    bool addFace(
        Element* leftElementPtr, std::size_t leftElementLocalFaceNo,
        Element* rightElementPtr, std::size_t rightElementLocalFaceNo,
        const Geometry::FaceType& faceType = Geometry::FaceType::WALL_BC);

    Edge* addEdge();

    void addNode();

    std::size_t getNumberOfElements(
        IteratorType part = IteratorType::LOCAL) const override {
        return theMesh_.getNumberOfElements(part);
    }

    std::size_t getNumberOfFaces(
        IteratorType part = IteratorType::LOCAL) const override {
        return theMesh_.getNumberOfFaces(part);
    }

    std::size_t getNumberOfEdges(
        IteratorType part = IteratorType::LOCAL) const override {
        return theMesh_.getNumberOfEdges(part);
    }

    std::size_t getNumberOfNodes(
        IteratorType part = IteratorType::LOCAL) const override {
        return theMesh_.getNumberOfNodes(part);
    }

    std::size_t getNumberOfNodeCoordinates() const {
        return theMesh_.getNumberOfNodeCoordinates();
    }

    /// *****************Iteration through the Elements*******************

    ConstElementIterator elementColBegin(
        IteratorType part = IteratorType::LOCAL) const override {
        return theMesh_.elementColBegin(part);
    }

    ConstElementIterator elementColEnd(
        IteratorType part = IteratorType::LOCAL) const override {
        return theMesh_.elementColEnd(part);
    }

    ElementIterator elementColBegin(
        IteratorType part = IteratorType::LOCAL) override {
        return theMesh_.elementColBegin(part);
    }

    ElementIterator elementColEnd(
        IteratorType part = IteratorType::LOCAL) override {
        return theMesh_.elementColEnd(part);
    }

    ConstFaceIterator faceColBegin(
        IteratorType part = IteratorType::LOCAL) const override {
        return theMesh_.faceColBegin(part);
    }

    ConstFaceIterator faceColEnd(
        IteratorType part = IteratorType::LOCAL) const override {
        return theMesh_.faceColEnd(part);
    }

    FaceIterator faceColBegin(
        IteratorType part = IteratorType::LOCAL) override {
        return theMesh_.faceColBegin(part);
    }

    FaceIterator faceColEnd(IteratorType part = IteratorType::LOCAL) override {
        return theMesh_.faceColEnd(part);
    }

    TreeIteratorConst<Edge*> edgeColBegin(
        IteratorType part = IteratorType::LOCAL) const override {
        return theMesh_.edgeColBegin(part);
    }

    TreeIteratorConst<Edge*> edgeColEnd(
        IteratorType part = IteratorType::LOCAL) const override {
        return theMesh_.edgeColEnd(part);
    }

    TreeIterator<Edge*> edgeColBegin(
        IteratorType part = IteratorType::LOCAL) override {
        return theMesh_.edgeColBegin(part);
    }

    TreeIterator<Edge*> edgeColEnd(
        IteratorType part = IteratorType::LOCAL) override {
        return theMesh_.edgeColEnd(part);
    }

    std::vector<Node*>::const_iterator nodeColBegin(
        IteratorType part = IteratorType::LOCAL) const override {
        return theMesh_.nodeColBegin(part);
    }

    std::vector<Node*>::const_iterator nodeColEnd(
        IteratorType part = IteratorType::LOCAL) const override {
        return theMesh_.nodeColEnd(part);
    }

    std::vector<Node*>::iterator nodeColBegin(
        IteratorType part = IteratorType::LOCAL) override {
        return theMesh_.nodeColBegin(part);
    }

    std::vector<Node*>::iterator nodeColEnd(
        IteratorType part = IteratorType::LOCAL) override {
        return theMesh_.nodeColEnd(part);
    }

    //  *****************Iteration through the Elements*******************

    /**
     * load a mesh that was generated and partitioned by the preprocessor
     */
    void readMesh(const std::string& filename);

#ifdef HPGEM_USE_QHULL
    /**
     * \brief create an unstructured triangular mesh
     * \details An iterative mesh generator based on "A Simple Mesh Generator in
     * Matlab" (Persson & Strang, 2004) The initial mesh uses the same structure
     * as createTriangularMesh. This algorithm will still function if the
     * bounding box of the domain is only known approximately If there is a very
     * large difference between the smallest desired edge length and the largest
     * desired edge length (ratio > ~4), it will usually help to overestimate
     * the domain size The domain has to be described implicitly by a distance
     * function, such that domainDescription(p) < 0 means p is inside the domain
     * and the gradient of domaindescription is normal to the boundary Local
     * element refinement is possible by providing desired relative edge lengths
     * of the output mesh. The scaling of this function has no effect on the
     * resulting mesh. If the edge scaling function returns NaN for some part of
     * the domain, it is expanded exponentially from the known parts This
     * routine cannot deal very well with concave(>pi) or very sharp
     * corners(~<pi/4), by placing fixed points on these location, performance
     * can be greatly improved Mesh quality is not guaraneed if growFactor is
     * much larger or smaller than 1
     * @param BottomLeft The bottom left corner of the bounding box of the
     * domain
     * @param TopRight The top right corner of the bounding box of the domain
     * @param TotalNoNodes The desired amount of nodes in the mesh
     * @param domainDescription A function that maps PointPhysicals to doubles,
     * such that negative numbers signify points inside the mesh
     * @param fixedPoints coordinates of point that MUST be in the mesh, no
     * matter what
     * @param relativeEdgeLength Allow r-refinement
     * @param growFactor specify how much larger than its neighbours an element
     * may be in areas where relativeEdgeLengths returns NaN
     */
    void createUnstructuredMesh(
        Geometry::PointPhysical<DIM> BottomLeft,
        Geometry::PointPhysical<DIM> TopRight, std::size_t TotalNoNodes,
        std::function<double(Geometry::PointPhysical<DIM>)> domainDescription,
        std::vector<Geometry::PointPhysical<DIM>> fixedPoints = {},
        std::function<double(Geometry::PointPhysical<DIM>)> relativeEdgeLength =
            [](Geometry::PointPhysical<DIM>) { return 1.; },
        double growFactor = 1.1);

    /**
     * \brief improve the mesh quality of an existing mesh
     * \details An iterative mesh generator based on "A Simple Mesh Generator in
     * Matlab" (Persson & Strang, 2004) The domain has to be described
     * implicitly by a distance function, such that domainDescription(p) < 0
     * means p is inside the domain and the gradient of domaindescription is
     * normal to the boundary Local element refinement is possible by providing
     * desired relative edge lengths of the output mesh. The scaling of this
     * function has no effect on the resulting mesh. If the edge scaling
     * function returns NaN for some part of the domain, it is expanded
     * exponentially from the known parts This routine cannot deal very well
     * with concave(>pi) or very sharp corners(~<pi/4), by placing fixed points
     * on these location, performance can be greatly improved If no implicit
     * description of the domain is available, fixing all boundary nodes usually
     * prevents the other nodes from escaping. In this case domainDescription
     * can return -1. for all p. \todo the current implementation throws away
     * all data, this behaviour should be replaced by an interpolation scheme
     * \todo this algoritm boils down to alternatingly doing delaunay
     * triangulations and moving nodes according to an optimally damped mass
     * spring system. This is currently done using a relatively crude
     * implementation. Once there is a proper coupling between mercury and
     * hpGEM, some thought should be given to the improvement of this algorithm
     * @param domainDescription A function that maps PointPhysicals to doubles,
     * such that negative numbers signify points inside the mesh
     * @param fixedPointIdxs pointIndexes of point that MUST remain in the same
     * location, no matter what
     * @param relativeEdgeLength Allow r-refinement
     * @param growFactor specify how much larger than its neighbours an element
     * may be in areas where relativeEdgeLengths returns NaN
     * @param isOnPeriodicBoundary A function that returns true if the point
     * passed is on a periodic boundary. This should only return true on the
     * 'master' half of the periodic boundary. Note that support for periodic
     * boundaries is not yet thoroughly tested and might be fragile. (Make sure
     * to allow a small tolerance/snapping distance.)
     * @param mapPeriodicNode A function that returns a point on the other
     * periodic boundary. You may assume isOnPeriodicBoundary returns true for
     * the argument of this function. To prevent infinite duplication you should
     * specify only half the periodic boundary as being periodic. The other half
     * will implicitly also become periodic because it receives periodic nodes
     * from this function
     * @param isOnOtherPeriodicBoundary A function that return true if the point
     * is on a periodic boundary, but isOnPeriodicBoundary returns false
     * @param safeNonPeriodicNode A function that maps the receiving end of the
     * periodic boundary to the interior of the domain so that non-periodic
     * nodes that stray near the periodic boundary can be safed
     * @param dontConnect specify a group of nodes that should not be connected
     * by elements, for example because they are part of a concave boundary.
     * Note that this can have unexpected effect if the nodes you specify are
     * not fixed
     */
    void updateMesh(
        std::function<double(Geometry::PointPhysical<DIM>)> domainDescription,
        std::vector<std::size_t> fixedPointIdxs = {},
        std::function<double(Geometry::PointPhysical<DIM>)> relativeEdgeLength =
            [](Geometry::PointPhysical<DIM>) { return 1.; },
        double growFactor = 1.1,
        std::function<bool(Geometry::PointPhysical<DIM>)> isOnPeriodicBoundary =
            [](Geometry::PointPhysical<DIM>) { return false; },
        std::function<
            Geometry::PointPhysical<DIM>(Geometry::PointPhysical<DIM>)>
            mapPeriodicNode = nullptr,
        std::function<bool(Geometry::PointPhysical<DIM>)>
            isOnOtherPeriodicBoundary =
                [](Geometry::PointPhysical<DIM>) { return false; },
        std::function<
            Geometry::PointPhysical<DIM>(Geometry::PointPhysical<DIM>)>
            safeNonPeriodicNode = nullptr,
        std::vector<std::size_t> dontConnect = {});
#endif

    //! Set MeshMoverBase object pointer, for moving meshes if needed
    void setMeshMover(const MeshMoverBase<DIM>* const meshMover);

    void move();

    // ********THESE SHOULD BE REPLACED by ITERABLE EDITIONS LATER**********

    //! Get const list of elements
    const LevelTree<Element*>& getElementsList(
        IteratorType part = IteratorType::LOCAL) const override {
        return theMesh_.getElementsList(part);
    }

    //! Get non-const list of elements
    LevelTree<Element*>& getElementsList(
        IteratorType part = IteratorType::LOCAL) override {
        return theMesh_.getElementsList(part);
    }

    //! Get const list of faces
    const LevelTree<Face*>& getFacesList(
        IteratorType part = IteratorType::LOCAL) const override {
        return theMesh_.getFacesList(part);
    }

    //! Get non-const list of faces
    LevelTree<Face*>& getFacesList(
        IteratorType part = IteratorType::LOCAL) override {
        return theMesh_.getFacesList(part);
    }

    const LevelTree<Edge*>& getEdgesList(
        IteratorType part = IteratorType::LOCAL) const override {
        return theMesh_.getEdgesList(part);
    }

    LevelTree<Edge*>& getEdgesList(
        IteratorType part = IteratorType::LOCAL) override {
        return theMesh_.getEdgesList(part);
    }

    const std::vector<Node*>& getNodesList(
        IteratorType part = IteratorType::LOCAL) const override {
        return theMesh_.getNodesList(part);
    }

    std::vector<Node*>& getNodesList(
        IteratorType part = IteratorType::LOCAL) override {
        return theMesh_.getNodesList(part);
    }

    const std::map<int, std::vector<Element*>>& getPullElements() override {
        return theMesh_.getPullElements();
    }

    const std::map<int, std::vector<Element*>>& getPushElements() override {
        return theMesh_.getPushElements();
    }
    // ************************************************************************

    //! Changes the default set of basisFunctions for this mesh and all of its
    //! elements. Ignores any conforming basisFunctionset that may be linked to
    //! faces/edges/...
    /// Using this to set the hpGEM provided conforming or DG basis functions is
    /// deprecated: The routines useDefaultDGBasisFunctionSet and
    /// useDefaultConformingBasisFunctionSet can do this more flexibly and also
    /// support mixed meshes
    void setDefaultBasisFunctionSet(BasisFunctionSet* bFSet);

    //! Adds vertex based degrees of freedom to the set of basisfunctions for
    //! this mesh and all of its vertices. This routine will assume that the
    //! first basisfunctionset connects to the first vertex (local numbering)
    //! and so on
    /// Using this to set the hpGEM provided conforming basis functions is
    /// deprecated: The routine useDefaultConformingBasis is more flexible and
    /// can also deal with mixed meshes
    void addVertexBasisFunctionSet(
        const std::vector<const BasisFunctionSet*>& bFsets);

    //! Adds face based degrees of freedom to the set of basisfunctions for this
    //! mesh and all of its faces. This routine will assume that all needed
    //! orientations are available in the collection of basisfunctionsets
    /// Using this to set the hpGEM provided conforming basis functions is
    /// deprecated: The routine useDefaultConformingBasis is more flexible and
    /// can also deal with mixed meshes
    void addFaceBasisFunctionSet(
        const std::vector<const OrientedBasisFunctionSet*>& bFsets);

    //! Adds edge based degrees of freedom to the set of basisfunctions for this
    //! mesh and all of its edges. This routine will assume that all needed
    //! orientations are available in the collection of basisfunctionsets
    /// Using this to set the hpGEM provided conforming basis functions is
    /// deprecated: The routine useDefaultConformingBasis is more flexible and
    /// can also deal with mixed meshes
    void addEdgeBasisFunctionSet(
        const std::vector<const OrientedBasisFunctionSet*>& bFsets);

    const std::vector<Geometry::PointPhysical<DIM>>& getNodeCoordinates()
        const {
        return theMesh_.getNodeCoordinates();
    }

    std::vector<Geometry::PointPhysical<DIM>>& getNodeCoordinates() {
        return theMesh_.getNodeCoordinates();
    }
    /**
     * Retrieves the Mesh as stored in this MeshManipulator
     */
    Mesh<DIM>& getMesh();
    const Mesh<DIM>& getMesh() const;

    /// given a PointPhysical inside the domain, finds an Element and a
    /// ReferencePoint such that Element::transform(PointReference) ==
    /// PointPhysical and Element::isInteriorPoint(PointReference)
    std::tuple<const Base::Element*, Geometry::PointReference<DIM>>
        physicalToReference(Geometry::PointPhysical<DIM>) const;

    /// add a PointPhysical where the solution is to be sampled more than once
    /// over the course of the simulation
    void addMeasurePoint(Geometry::PointPhysical<DIM>);

    /// returns all locations that were registered as interesting
    const std::vector<
        std::tuple<const Base::Element*, Geometry::PointReference<DIM>>>&
        getMeasurePoints() const;

    /// returns all locations that were registered as interesting
    std::vector<
        std::tuple<const Base::Element*, Geometry::PointReference<DIM>>>&
        getMeasurePoints();

    ///\brief placeholder refinement driver
    ///\details checks all elements on the current active level against the
    /// provided binary predicate. If it returns true it will refine the element
    /// according to the provided mapping. Note that the predicate should at
    /// least return false for elements that are already refined and for
    /// elements that have the wrong shape. Note also that you can do multiple
    /// passes on the same level if you want different refinement mappings for
    /// different elements. if you just want to refine every element on the
    /// level, you can omit the final argument. If the relative sizes of the
    /// refined elements don't suit your needs, use the mesh mover to relocate
    /// the nodes Note that doing so might render the coarse mesh useless, due
    /// to it using the old (incorrect) mapping
    // since a correct default for shouldRefine requires information from
    // refinementMapping it is set in the implementation
    void refine(const Geometry::RefinementMapping* refinementMapping,
                std::function<bool(const Element*)> shouldRefine = nullptr);

    //---------------------------------------------------------------------
   private:
    //! Construct the faces based on connectivity information about elements and
    //! nodes
    void faceFactory();

    //! Construct the faces based on connectivity information about elements and
    //! nodes
    void edgeFactory();

    // iterable should provide a begin() and an end() that return a
    // TreeEntry<Element*>
    // but concepts don't exist yet (will probably have to become an
    // Iterable<TreeEntry<Element*>> once they do)
    template <typename Iterable>
    std::tuple<const Base::Element*, Geometry::PointReference<DIM>>
        physicalToReference_detail(Geometry::PointPhysical<DIM>,
                                   Iterable elementContainer) const;

    Mesh<DIM> theMesh_;

    /// Pointer to MeshMoverBase, in order to move points in the mesh, when
    /// needed by user.
    const MeshMoverBase<DIM>* meshMover_;

    //! Collection of additional basis function set, if p-refinement is applied
    CollectionOfBasisFunctionSets collBasisFSet_;

    // when the mesh is updated, persistently store original node coordinates to
    // see if retriangulation is in order
    std::vector<Geometry::PointPhysical<DIM>> oldNodeLocations_;

    std::vector<std::tuple<const Base::Element*, Geometry::PointReference<DIM>>>
        measurePoints_;

    /// \brief Parse a double (possibly white space prefixed) from an input
    /// stream.
    ///
    /// This functions supports both the standard 3.2e+1 format and the
    /// hex-float format (e.g. 0x1.fp+4). Unfortunately the c++ specification is
    /// ambiguous on whether the latter should be parsable using istream and as
    /// a result the gnu's c++ and clang differ in their interpretation and
    /// implementation.
    /// \param stream The stream to read from.
    /// \return The resulting double
    double readDouble(std::istream& stream) const;
};

}  // namespace Base

#include "MeshManipulator_Impl.h"

#endif  // HPGEM_KERNEL_MESHMANIPULATOR_H

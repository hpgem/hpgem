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

//----------------------------------------------------------------
#ifndef HPGEM_KERNEL_ELEMENT_H
#define HPGEM_KERNEL_ELEMENT_H
//----------------------------------------------------------------
#include "Base/ElementData.h"
#include "Geometry/ElementGeometry.h"
#include "LinearAlgebra/MiddleSizeVector.h"
#include "FE/BasisFunctionSet.h"
#include "TreeEntry.h"

//#include "PhysGradientOfBasisFunction.h"//included after class definition due
// to cross dependencies (needed for templated function definition) #include
//"PhysicalElement.h"//included after class definition due to cross dependencies
//(needed for templated function definition)

#include <vector>
#include <iostream>
#include <memory>
#include <Integration/QuadratureRules/GaussQuadratureRule.h>
#include <Base/ElementBasisFunctions.h>
namespace hpgem {

namespace FE {
class BasisFunctionSet;
class BaseBasisFunction;
}  // namespace FE

namespace Base {
class Node;
class Edge;
class Face;

class ElementCacheData;
template <std::size_t DIM>
class PhysicalElement;

// class is final as a reminder that the virtual default destructor should be
// added once something inherits from this class
class Element final : public Geometry::ElementGeometry, public ElementData {
   public:
    using SolutionVector = LinearAlgebra::MiddleSizeVector;
    using CollectionOfBasisFunctionSets =
        std::vector<std::shared_ptr<const FE::BasisFunctionSet>>;

    template <std::size_t DIM>
    Element(const std::vector<std::size_t>& globalNodeIndexes,
            const CollectionOfBasisFunctionSets* basisFunctionSet,
            std::vector<Geometry::PointPhysical<DIM>>& allNodes,
            std::size_t numberOfUnkowns, std::size_t numberOfTimeLevels,
            std::size_t id, Zone* zone, std::size_t owner = 0,
            bool owning = true, std::size_t numberOfElementMatrices = 0,
            std::size_t numberOfElementVectors = 0);

    Element(const Element& other) = delete;
    Element& operator=(const Element& other) = delete;

    Element* copyWithoutFacesEdgesNodes();

    std::size_t getID() const;

    std::size_t getID();

    void setQuadratureRulesWithOrder(std::size_t quadrROrder);

    void setGaussQuadratureRule(
        QuadratureRules::GaussQuadratureRule* const quadR);

    void setDefaultBasisFunctionSet(std::size_t position) {
        // Clear all unknowns so that no old data is left from previous basis
        // function configurations.
        for (std::size_t unknown = 0; unknown < getNumberOfUnknowns();
             ++unknown) {
            basisFunctions_.clearBasisFunctionPosition(unknown);
        }
        for (std::size_t unknown = 0; unknown < getNumberOfUnknowns();
             ++unknown) {
            setDefaultBasisFunctionSet(position, unknown);
        }
    }

    void setDefaultBasisFunctionSet(std::size_t position, std::size_t unknown);

    void setVertexBasisFunctionSet(std::size_t position,
                                   std::size_t localIndex) {
        for (std::size_t unknown = 0; unknown < getNumberOfUnknowns();
             ++unknown) {
            setVertexBasisFunctionSet(position, localIndex, unknown);
        }
    }

    void setVertexBasisFunctionSet(std::size_t position, std::size_t localIndex,
                                   std::size_t unknown);

    void setEdgeBasisFunctionSet(std::size_t position, std::size_t localIndex) {
        for (std::size_t unknown = 0; unknown < getNumberOfUnknowns();
             ++unknown) {
            setEdgeBasisFunctionSet(position, localIndex, unknown);
        }
    }
    void setEdgeBasisFunctionSet(std::size_t position, std::size_t localIndex,
                                 std::size_t unknown);

    void setFaceBasisFunctionSet(std::size_t position, std::size_t localIndex) {
        for (std::size_t unknown = 0; unknown < getNumberOfUnknowns();
             ++unknown) {
            setFaceBasisFunctionSet(position, localIndex, unknown);
        }
    }
    void setFaceBasisFunctionSet(std::size_t position, std::size_t localIndex,
                                 std::size_t unknown);

    /// \brief Get a pointer to the quadrature rule used to do integration on
    /// this element.
    QuadratureRules::GaussQuadratureRule* getGaussQuadratureRule() const;

    // std::vector<Base::ElementCacheData>& getVecCacheData();

    /// \brief Get the value of the basis function (corresponding to index i) at
    /// the physical point corresponding to reference point p.
    template <std::size_t DIM>
    double basisFunction(std::size_t i,
                         const Geometry::PointReference<DIM>& p) const;

    template <std::size_t DIM>
    double basisFunction(std::size_t i, const Geometry::PointReference<DIM>& p,
                         std::size_t unknown) const;

    double basisFunction(std::size_t i,
                         QuadratureRules::GaussQuadratureRule* quadratureRule,
                         std::size_t quadraturePointIndex) const;
    double basisFunction(std::size_t i,
                         QuadratureRules::GaussQuadratureRule* quadratureRule,
                         std::size_t quadraturePointIndex,
                         std::size_t unknown) const;

    double basisFunction(
        std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule,
        std::size_t quadraturePointIndex,
        const Geometry::MappingReferenceToReference<1>* map) const;
    double basisFunction(std::size_t i,
                         QuadratureRules::GaussQuadratureRule* quadratureRule,
                         std::size_t quadraturePointIndex,
                         const Geometry::MappingReferenceToReference<1>* map,
                         std::size_t unknown) const;

    /// \brief returns the value of the i-th basisfunction at point p in ret.
    template <std::size_t DIM>
    void basisFunction(std::size_t i, const Geometry::PointReference<DIM>& p,
                       LinearAlgebra::SmallVector<DIM>& ret) const;

    template <std::size_t DIM>
    void basisFunction(std::size_t i, const Geometry::PointReference<DIM>& p,
                       LinearAlgebra::SmallVector<DIM>& ret,
                       std::size_t unknown) const;

    template <std::size_t DIM>
    void basisFunction(std::size_t i,
                       QuadratureRules::GaussQuadratureRule* quadratureRule,
                       std::size_t quadraturePointIndex,
                       LinearAlgebra::SmallVector<DIM>& ret) const;
    template <std::size_t DIM>
    void basisFunction(std::size_t i,
                       QuadratureRules::GaussQuadratureRule* quadratureRule,
                       std::size_t quadraturePointIndex,
                       LinearAlgebra::SmallVector<DIM>& ret,
                       std::size_t unknown) const;

    template <std::size_t DIM>
    void basisFunction(std::size_t i,
                       QuadratureRules::GaussQuadratureRule* quadratureRule,
                       std::size_t quadraturePointIndex,
                       const Geometry::MappingReferenceToReference<1>* map,
                       LinearAlgebra::SmallVector<DIM>& ret) const;

    template <std::size_t DIM>
    void basisFunction(std::size_t i,
                       QuadratureRules::GaussQuadratureRule* quadratureRule,
                       std::size_t quadraturePointIndex,
                       const Geometry::MappingReferenceToReference<1>* map,
                       LinearAlgebra::SmallVector<DIM>& ret,
                       std::size_t unknown) const;
    /// \brief Get the value of the derivative of the physical basisfunction
    /// (corresponding to index i) in the direction jDir at the physical point
    /// corresponding to reference point p.
    template <std::size_t DIM>
    double basisFunctionDeriv(std::size_t i, std::size_t jDir,
                              const Geometry::PointReference<DIM>& p) const;

    template <std::size_t DIM>
    double basisFunctionDeriv(std::size_t i, std::size_t jDir,
                              const Geometry::PointReference<DIM>& p,
                              std::size_t unknown) const;

    /// \brief Get the gradient of the physical basis function (corresponding to
    /// index i) at the physical point corresponding to reference point p.
    /// \details WARNING: contrary to some previous versions of hpGEM, this will
    /// NOT do any coordinate transformations!!
    template <std::size_t DIM>
    LinearAlgebra::SmallVector<DIM> basisFunctionDeriv(
        std::size_t i, const Geometry::PointReference<DIM>& p) const;
    template <std::size_t DIM>
    LinearAlgebra::SmallVector<DIM> basisFunctionDeriv(
        std::size_t i, const Geometry::PointReference<DIM>& p,
        std::size_t unknown) const;

    template <std::size_t DIM>
    LinearAlgebra::SmallVector<DIM> basisFunctionDeriv(
        std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule,
        std::size_t quadraturePointIndex) const;
    template <std::size_t DIM>
    LinearAlgebra::SmallVector<DIM> basisFunctionDeriv(
        std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule,
        std::size_t quadraturePointIndex, std::size_t unknown) const;

    template <std::size_t DIM>
    LinearAlgebra::SmallVector<DIM> basisFunctionDeriv(
        std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule,
        std::size_t quadraturePointIndex,
        const Geometry::MappingReferenceToReference<1>* map) const;
    template <std::size_t DIM>
    LinearAlgebra::SmallVector<DIM> basisFunctionDeriv(
        std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule,
        std::size_t quadraturePointIndex,
        const Geometry::MappingReferenceToReference<1>* map,
        std::size_t unknown) const;

    /// \brief Returns the curl of the i-th basisfunction at point p in ret.
    template <std::size_t DIM>
    LinearAlgebra::SmallVector<DIM> basisFunctionCurl(
        std::size_t i, const Geometry::PointReference<DIM>& p) const;
    /// \brief Returns the curl of the i-th basisfunction at point p in ret.
    template <std::size_t DIM>
    LinearAlgebra::SmallVector<DIM> basisFunctionCurl(
        std::size_t i, const Geometry::PointReference<DIM>& p,
        std::size_t unknown) const;

    template <std::size_t DIM>
    LinearAlgebra::SmallVector<DIM> basisFunctionCurl(
        std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule,
        std::size_t quadraturePointIndex) const;
    template <std::size_t DIM>
    LinearAlgebra::SmallVector<DIM> basisFunctionCurl(
        std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule,
        std::size_t quadraturePointIndex, std::size_t unknown) const;

    template <std::size_t DIM>
    LinearAlgebra::SmallVector<DIM> basisFunctionCurl(
        std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule,
        std::size_t quadraturePointIndex,
        const Geometry::MappingReferenceToReference<1>* map) const;
    template <std::size_t DIM>
    LinearAlgebra::SmallVector<DIM> basisFunctionCurl(
        std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule,
        std::size_t quadraturePointIndex,
        const Geometry::MappingReferenceToReference<1>* map,
        std::size_t unknown) const;

    /// \brief Returns the divergence of the i-th basisfunction at point p in
    /// ret.
    template <std::size_t DIM>
    double basisFunctionDiv(std::size_t i,
                            const Geometry::PointReference<DIM>& p) const;
    template <std::size_t DIM>
    double basisFunctionDiv(std::size_t i,
                            const Geometry::PointReference<DIM>& p,
                            std::size_t unknown) const;

    double basisFunctionDiv(
        std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule,
        std::size_t quadraturePointIndex) const;
    double basisFunctionDiv(
        std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule,
        std::size_t quadraturePointIndex, std::size_t unknown) const;

    double basisFunctionDiv(
        std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule,
        std::size_t quadraturePointIndex,
        const Geometry::MappingReferenceToReference<1>* map) const;
    double basisFunctionDiv(
        std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule,
        std::size_t quadraturePointIndex,
        const Geometry::MappingReferenceToReference<1>* map,
        std::size_t unknown) const;

    /// \brief Get the solution at the given timeLevel at the physical point
    /// corresponding to reference point p. \details This routine assumes the
    /// result of evaluating a basis function has to be transformed using the
    /// identity transformation
    template <std::size_t DIM>
    SolutionVector getSolution(std::size_t timeLevel,
                               const Geometry::PointReference<DIM>& p) const;

    /// \brief Get the gradient of the solution at the given timeLevel at the
    /// physical point corresponding to reference point p. \details returns a
    /// vector of gradients. This routine assumes the result of evaluating a
    /// gradient of a basis function has to be transformed using the identity
    /// transformation
    template <std::size_t DIM>
    std::vector<LinearAlgebra::SmallVector<DIM>> getSolutionGradient(
        std::size_t timeLevel, const Geometry::PointReference<DIM>& p) const;

    /// \brief Get the solution at the given timeLevel at the physical point
    /// corresponding to reference point p. \details uses the physical element
    /// for evaluation and transformation of the basis functions
    template <std::size_t DIM>
    SolutionVector getSolution(std::size_t timeLevel,
                               PhysicalElement<DIM>& element) const;

    /// \brief Get the gradient of the solution at the given timeLevel at the
    /// physical point corresponding to reference point p. \details returns a
    /// vector of gradients. Uses the physical element for evaluation and
    /// transformation of gradients of the basis functions
    template <std::size_t DIM>
    std::vector<LinearAlgebra::SmallVector<DIM>> getSolutionGradient(
        std::size_t timeLevel, PhysicalElement<DIM>& element) const;

    ///\todo not implemented
    void initialiseSolution(std::size_t timeLevel, std::size_t solutionId,
                            const SolutionVector& solution);

    void setFace(std::size_t localFaceNumber, Face* face);

    void setEdge(std::size_t localEdgeNumber, Edge* edge);

    void setNode(std::size_t localNodeNumber, Node* node);

    ///\deprecated Does not follow naming conventions, use
    /// getLocalNumberOfBasisFunctions instead
    std::size_t getLocalNrOfBasisFunctions() const {
        return getLocalNumberOfBasisFunctions();
    }

    std::size_t getLocalNrOfBasisFunctions(std::size_t unknown) const {
        return getLocalNumberOfBasisFunctions(unknown);
    }

    /// return the number of basis functions that are associated with this
    /// element only. This always includes functions with compact support on the
    /// interior of the element and DG basis function, but never include
    /// conforming basis functions that are nonzero on a face, edge or node
    std::size_t getLocalNumberOfBasisFunctions() const {
        return basisFunctions_.getNumberOfLocalBasisFunctions();
    }

    std::size_t getLocalNumberOfBasisFunctions(std::size_t unknown) const {
        return basisFunctions_.getNumberOfLocalBasisFunctions(unknown);
    }

    std::size_t getTotalLocalNumberOfBasisFunctions() const {
        return basisFunctions_.getTotalLocalNumberOfBasisFunctions();
    }

    Face* getFace(std::size_t localFaceNumber) const {
        logger.assert_debug(localFaceNumber < getNumberOfFaces(),
                            "Asked for face %, but there are only % faces",
                            localFaceNumber, getNumberOfFaces());
        return facesList_[localFaceNumber];
    }

    const std::vector<Face*> getFacesList() const { return facesList_; }

    Edge* getEdge(std::size_t localEdgeNumber) const {
        logger.assert_debug(localEdgeNumber < getNumberOfEdges(),
                            "Asked for edge %, but there are only % edges",
                            localEdgeNumber, getNumberOfEdges());
        return edgesList_[localEdgeNumber];
    }

    const std::vector<Edge*> getEdgesList() const { return edgesList_; }

    const Node* getNode(std::size_t localNodeNumber) const {
        logger.assert_debug(localNodeNumber < getNumberOfNodes(),
                            "Asked for node %, but there are only % nodes",
                            localNodeNumber, getNumberOfNodes());
        return nodesList_[localNodeNumber];
    }

    Node* getNode(std::size_t localNodeNumber) {
        logger.assert_debug(localNodeNumber < getNumberOfNodes(),
                            "Asked for node %, but there are only % nodes",
                            localNodeNumber, getNumberOfNodes());
        return nodesList_[localNodeNumber];
    }

    const std::vector<Node*> getNodesList() const { return nodesList_; }

    /// Compute the face id of a face that is adjacent to the element.
    std::size_t getLocalId(const Base::Face* face) const {
        for (std::size_t i = 0; i < facesList_.size(); ++i) {
            if (facesList_[i] == face) return i;
        }
        logger.assert_debug(false, "Not an face adjacent to the element");
        return std::numeric_limits<std::size_t>::max();
    }

    /// Compute the edge id of a edge that is adjacent to the element.
    std::size_t getLocalId(const Base::Edge* edge) const {
        for (std::size_t i = 0; i < edgesList_.size(); ++i) {
            if (edgesList_[i] == edge) return i;
        }
        logger.assert_debug(false, "Not an edge adjacent to the element");
        return std::numeric_limits<std::size_t>::max();
    }

    /// Compute the node id of a node that is adjacent to the element.
    std::size_t getLocalId(const Base::Node* node) const {
        for (std::size_t i = 0; i < nodesList_.size(); ++i) {
            if (nodesList_[i] == node) return i;
        }
        logger.assert_debug(false, "Not an node adjacent to the element");
        return std::numeric_limits<std::size_t>::max();
    }

    std::size_t getBasisFunctionOffset(const Base::Face* face,
                                       std::size_t unknown) const {
        return getFaceBasisFunctionOffset(getLocalId(face), unknown);
    }

    std::size_t getBasisFunctionOffset(const Base::Edge* edge,
                                       std::size_t unknown) const {
        return getEdgeBasisFunctionOffset(getLocalId(edge), unknown);
    }

    std::size_t getBasisFunctionOffset(const Base::Node* node,
                                       std::size_t unknown) const {
        return getNodeBasisFunctionOffset(getLocalId(node), unknown);
    }

    std::size_t getFaceBasisFunctionOffset(std::size_t localFaceId,
                                           std::size_t unknown) const {
        return basisFunctions_.getBasisFunctionOffset(unknown, 1 + localFaceId);
    }

    std::size_t getEdgeBasisFunctionOffset(std::size_t localEdgeId,
                                           std::size_t unknown) const {
        return basisFunctions_.getBasisFunctionOffset(
            unknown, 1 + getNumberOfFaces() + localEdgeId);
    }

    std::size_t getNodeBasisFunctionOffset(std::size_t localNodeId,
                                           std::size_t unknown) const {
        return basisFunctions_.getBasisFunctionOffset(
            unknown, 1 + getNumberOfFaces() + getNumberOfEdges() + localNodeId);
    }

    ///\deprecated Does not follow naming conventions, use getNumberOfFaces
    /// instead
    std::size_t getNrOfFaces() const { return getNumberOfFaces(); }

    ///\deprecated Does not follow naming conventions, use getNumberOfEdges
    /// instead
    std::size_t getNrOfEdges() const { return getNumberOfEdges(); }

    ///\deprecated Does not follow naming conventions, use getNumberOfNodes
    /// instead
    std::size_t getNrOfNodes() const { return getNumberOfNodes(); }

    std::size_t getNumberOfFaces() const { return facesList_.size(); }

    std::size_t getNumberOfEdges() const { return edgesList_.size(); }

    ///\todo This function overwrites the non-virtual function in
    /// ElementGeometry. Either remove this function, rename one of the two
    /// functions or make the one in ElementGeometry virtual
    std::size_t getNumberOfNodes() const { return nodesList_.size(); }

#ifndef NDEBUG
    const FE::BaseBasisFunction* getBasisFunction(std::size_t i) const;
#endif

    void setPositionInTree(const TreeEntry<Element*>* position) {
        logger.assert_debug(
            position->getData() == this,
            "Trying to set the position of another element as this element");
        positionInTheTree_ = position;
    }

    const TreeEntry<Element*>* getPositionInTree() const {
        return positionInTheTree_;
    }

    /// \brief Adjust the ownership flag of this Element.
    ///
    /// Note that this does not change the topology as adminstrated in Mesh
    /// and Submesh. The information there should be updated separately.
    /// \param owner The process owning this element.
    /// \param owned Whether this element is owned or not by the current
    /// processor.
    void setOwnedByCurrentProcessor(std::size_t owner, bool owned);
    bool isOwnedByCurrentProcessor() const;
    std::size_t getOwner() const;

    /// Output operator.
    friend std::ostream& operator<<(std::ostream& os, const Element& element);

   private:
    /// Constructor that copies the data and geometry of the given ElementData
    /// and ElementGeometry.
    Element(std::size_t owner, bool owned, const ElementData& otherData,
            const Geometry::ElementGeometry& otherGeometry);

    /// Quadrature rule used to do the integration on this element.
    QuadratureRules::GaussQuadratureRule* quadratureRule_;

    const TreeEntry<Element*>* positionInTheTree_;

    /// Identifier (index) of the element.
    std::size_t id_;

    /// Constant that describes a relation between the polynomial order of the
    /// basis function set and the accuracy of the quadrature rule.
    std::size_t orderCoeff_;

    std::vector<Face*> facesList_;
    std::vector<Edge*> edgesList_;
    std::vector<Node*> nodesList_;

    ElementBasisFunctions basisFunctions_;

    /// Whether the current processor owns this element.
    bool owned_;
    // MPI rank of the owning processor
    std::size_t owner_;
};
}  // namespace Base
}  // namespace hpgem
#include "PhysGradientOfBasisFunction.h"
#include "PhysicalElement.h"

namespace hpgem {

namespace Base {
/// \details The user does not need to worry about the construction of elements.
/// This is done by mesh-generators. For example the interface HpgemAPIBase can
/// be used to create meshes.
template <std::size_t DIM>
Element::Element(const std::vector<std::size_t>& globalNodeIndexes,
                 const CollectionOfBasisFunctionSets* basisFunctionSet,
                 std::vector<Geometry::PointPhysical<DIM>>& allNodes,
                 std::size_t numberOfUnknowns, std::size_t numberOfTimeLevels,
                 std::size_t id, Zone* zone, std::size_t owner, bool owned,
                 std::size_t numberOfElementMatrices,
                 std::size_t numberOfElementVectors)
    : ElementGeometry(globalNodeIndexes, allNodes),
      ElementData(numberOfTimeLevels, numberOfUnknowns, zone,
                  numberOfElementMatrices, numberOfElementVectors),
      quadratureRule_(nullptr),
      id_(id),
      basisFunctions_(basisFunctionSet, numberOfUnknowns),
      owned_(owned),
      owner_(owner) {
    logger.assert_debug(basisFunctionSet != nullptr,
                        "Invalid basis function set passed");
    logger(VERBOSE, "numberOfElementMatrices: %", numberOfElementMatrices);
    logger(VERBOSE, "numberOfElementVectors: %", numberOfElementVectors);

    orderCoeff_ = 2;  // for safety
    basisFunctions_.validatePositions();

    for (std::size_t i = 0; i < numberOfUnknowns; ++i) {
        setNumberOfBasisFunctions(basisFunctions_.getNumberOfBasisFunctions(i),
                                  i);
    }

    setQuadratureRulesWithOrder(orderCoeff_ *
                                basisFunctions_.getMaximumOrder());
    facesList_.assign(getReferenceGeometry()->getNumberOfCodim1Entities(),
                      nullptr);
    if (getReferenceGeometry()->getNumberOfCodim3Entities() > 0) {
        edgesList_.assign(getReferenceGeometry()->getNumberOfCodim2Entities(),
                          nullptr);
    }
    nodesList_.assign(getReferenceGeometry()->getNumberOfNodes(), nullptr);
}

template <std::size_t DIM>
double Element::basisFunctionDeriv(
    std::size_t i, std::size_t jDir,
    const Geometry::PointReference<DIM>& p) const {
    logger.assert_debug(
        i < getNumberOfBasisFunctions(),
        "Asked for basis function %, but there are only % basis functions", i,
        getNumberOfBasisFunctions());
    logger.assert_debug(
        (jDir < p.size()),
        "Error in BasisFunctionSet.EvalDeriv: invalid derivative direction!");

    const FE::BasisFunctionSet* subSet;
    std::size_t subIndex;
    std::tie(subSet, subIndex) = basisFunctions_.getBasisFunctionSetAndIndex(i);
    return subSet->evalDeriv(subIndex, jDir, p);
}

template <std::size_t DIM>
double Element::basisFunctionDeriv(std::size_t i, std::size_t jDir,
                                   const Geometry::PointReference<DIM>& p,
                                   std::size_t unknown) const {
    logger.assert_debug(
        i < getNumberOfBasisFunctions(unknown),
        "Asked for basis function %, but there are only % basis functions", i,
        getNumberOfBasisFunctions(unknown));
    logger.assert_debug(
        (jDir < p.size()),
        "Error in BasisFunctionSet.EvalDeriv: invalid derivative direction!");

    const FE::BasisFunctionSet* subSet;
    std::size_t subIndex;
    std::tie(subSet, subIndex) =
        basisFunctions_.getBasisFunctionSetAndIndex(i, unknown);
    return subSet->evalDeriv(subIndex, jDir, p);
}

template <std::size_t DIM>
double Element::basisFunction(std::size_t i,
                              const Geometry::PointReference<DIM>& p) const {
    logger.assert_debug(
        i < getNumberOfBasisFunctions(),
        "Asked for basis function %, but there are only % basis functions", i,
        getNumberOfBasisFunctions());
    const FE::BasisFunctionSet* subSet;
    std::size_t subIndex;
    std::tie(subSet, subIndex) = basisFunctions_.getBasisFunctionSetAndIndex(i);
    return subSet->eval(subIndex, p);
}

template <std::size_t DIM>
double Element::basisFunction(std::size_t i,
                              const Geometry::PointReference<DIM>& p,
                              std::size_t unknown) const {
    logger.assert_debug(
        i < getNumberOfBasisFunctions(unknown),
        "Asked for basis function %, but there are only % basis functions", i,
        getNumberOfBasisFunctions(unknown));
    const FE::BasisFunctionSet* subSet;
    std::size_t subIndex;
    std::tie(subSet, subIndex) =
        basisFunctions_.getBasisFunctionSetAndIndex(i, unknown);
    return subSet->eval(subIndex, p);
}

template <std::size_t DIM>
void Element::basisFunction(std::size_t i,
                            const Geometry::PointReference<DIM>& p,
                            LinearAlgebra::SmallVector<DIM>& ret) const {
    logger.assert_debug(
        i < getNumberOfBasisFunctions(),
        "Asked for basis function %, but there are only % basis functions", i,
        getNumberOfBasisFunctions());
    const FE::BasisFunctionSet* subSet;
    std::size_t subIndex;
    std::tie(subSet, subIndex) = basisFunctions_.getBasisFunctionSetAndIndex(i);
    subSet->eval(subIndex, p, ret);
}

template <std::size_t DIM>
void Element::basisFunction(std::size_t i,
                            const Geometry::PointReference<DIM>& p,
                            LinearAlgebra::SmallVector<DIM>& ret,
                            std::size_t unknown) const {
    logger.assert_debug(
        i < getNumberOfBasisFunctions(unknown),
        "Asked for basis function %, but there are only % basis functions", i,
        getNumberOfBasisFunctions(unknown));
    const FE::BasisFunctionSet* subSet;
    std::size_t subIndex;
    std::tie(subSet, subIndex) =
        basisFunctions_.getBasisFunctionSetAndIndex(i, unknown);
    subSet->eval(subIndex, p, ret);
}

template <std::size_t DIM>
LinearAlgebra::SmallVector<DIM> Element::basisFunctionCurl(
    std::size_t i, const Geometry::PointReference<DIM>& p) const {
    logger.assert_debug(
        i < getNumberOfBasisFunctions(),
        "Asked for basis function %, but there are only % basis functions", i,
        getNumberOfBasisFunctions());
    const FE::BasisFunctionSet* subSet;
    std::size_t subIndex;
    std::tie(subSet, subIndex) = basisFunctions_.getBasisFunctionSetAndIndex(i);
    return subSet->evalCurl(subIndex, p);
}

template <std::size_t DIM>
LinearAlgebra::SmallVector<DIM> Element::basisFunctionCurl(
    std::size_t i, const Geometry::PointReference<DIM>& p,
    std::size_t unknown) const {
    logger.assert_debug(
        i < getNumberOfBasisFunctions(unknown),
        "Asked for basis function %, but there are only % basis functions", i,
        getNumberOfBasisFunctions(unknown));
    const FE::BasisFunctionSet* subSet;
    std::size_t subIndex;
    std::tie(subSet, subIndex) =
        basisFunctions_.getBasisFunctionSetAndIndex(i, unknown);
    return subSet->evalCurl(subIndex, p);
}

template <std::size_t DIM>
double Element::basisFunctionDiv(std::size_t i,
                                 const Geometry::PointReference<DIM>& p) const {
    logger.assert_debug(
        i < getNumberOfBasisFunctions(),
        "Asked for basis function %, but there are only % basis functions", i,
        getNumberOfBasisFunctions());
    const FE::BasisFunctionSet* subSet;
    std::size_t subIndex;
    std::tie(subSet, subIndex) = basisFunctions_.getBasisFunctionSetAndIndex(i);
    return subSet->evalDiv(subIndex, p);
}

template <std::size_t DIM>
double Element::basisFunctionDiv(std::size_t i,
                                 const Geometry::PointReference<DIM>& p,
                                 std::size_t unknown) const {
    logger.assert_debug(
        i < getNumberOfBasisFunctions(unknown),
        "Asked for basis function %, but there are only % basis functions", i,
        getNumberOfBasisFunctions(unknown));
    const FE::BasisFunctionSet* subSet;
    std::size_t subIndex;
    std::tie(subSet, subIndex) =
        basisFunctions_.getBasisFunctionSetAndIndex(i, unknown);
    return subSet->evalDiv(subIndex, p);
}

template <std::size_t DIM>
LinearAlgebra::SmallVector<DIM> Element::basisFunctionDeriv(
    std::size_t i, const Geometry::PointReference<DIM>& p) const {
    logger.assert_debug(
        i < getNumberOfBasisFunctions(),
        "Asked for basis function %, but there are only % basis functions", i,
        getNumberOfBasisFunctions());
    const FE::BasisFunctionSet* subSet;
    std::size_t subIndex;
    std::tie(subSet, subIndex) = basisFunctions_.getBasisFunctionSetAndIndex(i);
    return subSet->evalDeriv(subIndex, p);
}

template <std::size_t DIM>
LinearAlgebra::SmallVector<DIM> Element::basisFunctionDeriv(
    std::size_t i, const Geometry::PointReference<DIM>& p,
    std::size_t unknown) const {
    logger.assert_debug(
        i < getNumberOfBasisFunctions(unknown),
        "Asked for basis function %, but there are only % basis functions", i,
        getNumberOfBasisFunctions(unknown));
    const FE::BasisFunctionSet* subSet;
    std::size_t subIndex;
    std::tie(subSet, subIndex) =
        basisFunctions_.getBasisFunctionSetAndIndex(i, unknown);
    return subSet->evalDeriv(subIndex, p);
}

template <std::size_t DIM>
Element::SolutionVector Element::getSolution(
    std::size_t timeIntegrationVectorId,
    const Geometry::PointReference<DIM>& p) const {
    std::size_t numberOfUnknowns = ElementData::getNumberOfUnknowns();
    std::vector<std::size_t> numberOfBasisFunctions =
        std::vector<std::size_t>(numberOfUnknowns, 0);
    SolutionVector solution(numberOfUnknowns);

    LinearAlgebra::MiddleSizeVector data =
        ElementData::getTimeIntegrationVector(timeIntegrationVectorId);

    std::size_t iVB = 0;
    for (std::size_t iV = 0; iV < numberOfUnknowns; ++iV) {
        numberOfBasisFunctions[iV] = ElementData::getNumberOfBasisFunctions(iV);
        for (std::size_t iB = 0; iB < numberOfBasisFunctions[iV]; ++iB) {
            iVB = convertToSingleIndex(iB, iV);
            solution[iV] += data(iVB) * basisFunction(iB, p);
        }
    }
    return solution;
}

template <std::size_t DIM>
std::vector<LinearAlgebra::SmallVector<DIM>> Element::getSolutionGradient(
    std::size_t timeIntegrationVectorId,
    const Geometry::PointReference<DIM>& p) const {
    std::size_t numberOfUnknowns = ElementData::getNumberOfUnknowns();
    std::vector<std::size_t> numberOfBasisFunctions =
        std::vector<std::size_t>(numberOfUnknowns, 0);
    std::vector<LinearAlgebra::SmallVector<DIM>> solution(numberOfUnknowns);
    auto jacobean = getReferenceToPhysicalMap()->calcJacobian(p);
    jacobean = jacobean.transpose();

    LinearAlgebra::MiddleSizeVector data =
        ElementData::getTimeIntegrationVector(timeIntegrationVectorId);

    std::size_t iVB = 0;
    for (std::size_t iV = 0; iV < numberOfUnknowns; ++iV) {
        numberOfBasisFunctions[iV] = ElementData::getNumberOfBasisFunctions(iV);
        for (std::size_t iB = 0; iB < numberOfBasisFunctions[iV]; ++iB) {
            iVB = convertToSingleIndex(iB, iV);
            auto derivative = basisFunctionDeriv(iB, p);
            jacobean.solve(derivative);
            solution[iV] += data(iVB) * derivative;
        }
    }
    return solution;
}

template <std::size_t DIM>
Element::SolutionVector Element::getSolution(
    std::size_t timeIntegrationVectorId, PhysicalElement<DIM>& element) const {
    logger.assert_debug(element.getElement() == this,
                        "Cannot find the solution in a different element!");
    std::size_t numberOfUnknowns = ElementData::getNumberOfUnknowns();
    std::vector<std::size_t> numberOfBasisFunctions =
        std::vector<std::size_t>(numberOfUnknowns, 0);
    SolutionVector solution(numberOfUnknowns);

    const LinearAlgebra::MiddleSizeVector& data =
        ElementData::getTimeIntegrationVector(timeIntegrationVectorId);

    std::size_t iVB = 0;
    for (std::size_t iV = 0; iV < numberOfUnknowns; ++iV) {
        numberOfBasisFunctions[iV] = ElementData::getNumberOfBasisFunctions(iV);
        for (std::size_t iB = 0; iB < numberOfBasisFunctions[iV]; ++iB) {
            iVB = convertToSingleIndex(iB, iV);
            solution[iV] += data(iVB) * element.basisFunction(iB);
        }
    }
    return solution;
}

template <std::size_t DIM>
std::vector<LinearAlgebra::SmallVector<DIM>> Element::getSolutionGradient(
    std::size_t timeIntegrationVectorId, PhysicalElement<DIM>& element) const {
    logger.assert_debug(
        element.getElement() == this,
        "Cannot find the gradient of the solution in a different element!");
    std::size_t numberOfUnknowns = ElementData::getNumberOfUnknowns();
    std::vector<std::size_t> numberOfBasisFunctions =
        std::vector<std::size_t>(numberOfUnknowns, 0);
    std::vector<LinearAlgebra::SmallVector<DIM>> solution(numberOfUnknowns);

    LinearAlgebra::MiddleSizeVector data =
        ElementData::getTimeIntegrationVector(timeIntegrationVectorId);

    std::size_t iVB = 0;
    for (std::size_t iV = 0; iV < numberOfUnknowns; ++iV) {
        numberOfBasisFunctions[iV] = ElementData::getNumberOfBasisFunctions(iV);
        for (std::size_t iB = 0; iB < numberOfBasisFunctions[iV]; ++iB) {
            iVB = convertToSingleIndex(iB, iV);
            solution[iV] += data(iVB) * element.basisFunctionDeriv(iB);
        }
    }
    return solution;
}

template <std::size_t DIM>
void Element::basisFunction(
    std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule,
    std::size_t quadraturePointIndex,
    LinearAlgebra::SmallVector<DIM>& ret) const {
    logger.assert_debug(
        i < getNumberOfBasisFunctions(),
        "Asked for basis function %, but there are only % basis functions", i,
        getNumberOfBasisFunctions());
    const FE::BasisFunctionSet* subSet;
    std::size_t subIndex;
    std::tie(subSet, subIndex) = basisFunctions_.getBasisFunctionSetAndIndex(i);
    quadratureRule->eval(subSet, subIndex, quadraturePointIndex, ret);
}

template <std::size_t DIM>
void Element::basisFunction(
    std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule,
    std::size_t quadraturePointIndex, LinearAlgebra::SmallVector<DIM>& ret,
    std::size_t unknown) const {
    logger.assert_debug(
        i < getNumberOfBasisFunctions(unknown),
        "Asked for basis function %, but there are only % basis functions", i,
        getNumberOfBasisFunctions(unknown));
    const FE::BasisFunctionSet* subSet;
    std::size_t subIndex;
    std::tie(subSet, subIndex) =
        basisFunctions_.getBasisFunctionSetAndIndex(i, unknown);
    quadratureRule->eval(subSet, subIndex, quadraturePointIndex, ret);
}

template <std::size_t DIM>
LinearAlgebra::SmallVector<DIM> Element::basisFunctionDeriv(
    std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule,
    std::size_t quadraturePointIndex) const {
    logger.assert_debug(
        i < getNumberOfBasisFunctions(),
        "Asked for basis function %, but there are only % basis functions", i,
        getNumberOfBasisFunctions());
    const FE::BasisFunctionSet* subSet;
    std::size_t subIndex;
    std::tie(subSet, subIndex) = basisFunctions_.getBasisFunctionSetAndIndex(i);
    return quadratureRule->evalGrad(subSet, subIndex, quadraturePointIndex);
}

template <std::size_t DIM>
LinearAlgebra::SmallVector<DIM> Element::basisFunctionDeriv(
    std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule,
    std::size_t quadraturePointIndex, std::size_t unknown) const {
    logger.assert_debug(
        i < getNumberOfBasisFunctions(unknown),
        "Asked for basis function %, but there are only % basis functions", i,
        getNumberOfBasisFunctions(unknown));
    const FE::BasisFunctionSet* subSet;
    std::size_t subIndex;
    std::tie(subSet, subIndex) =
        basisFunctions_.getBasisFunctionSetAndIndex(i, unknown);
    return quadratureRule->evalGrad(subSet, subIndex, quadraturePointIndex);
}

template <std::size_t DIM>
LinearAlgebra::SmallVector<DIM> Element::basisFunctionCurl(
    std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule,
    std::size_t quadraturePointIndex) const {
    logger.assert_debug(
        i < getNumberOfBasisFunctions(),
        "Asked for basis function %, but there are only % basis functions", i,
        getNumberOfBasisFunctions());
    const FE::BasisFunctionSet* subSet;
    std::size_t subIndex;
    std::tie(subSet, subIndex) = basisFunctions_.getBasisFunctionSetAndIndex(i);
    return quadratureRule->evalCurl<DIM>(subSet, subIndex,
                                         quadraturePointIndex);
}

template <std::size_t DIM>
LinearAlgebra::SmallVector<DIM> Element::basisFunctionCurl(
    std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule,
    std::size_t quadraturePointIndex, std::size_t unknown) const {
    logger.assert_debug(
        i < getNumberOfBasisFunctions(unknown),
        "Asked for basis function %, but there are only % basis functions", i,
        getNumberOfBasisFunctions(unknown));
    const FE::BasisFunctionSet* subSet;
    std::size_t subIndex;
    std::tie(subSet, subIndex) =
        basisFunctions_.getBasisFunctionSetAndIndex(i, unknown);
    return quadratureRule->evalCurl<DIM>(subSet, subIndex,
                                         quadraturePointIndex);
}

template <std::size_t DIM>
void Element::basisFunction(
    std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule,
    std::size_t quadraturePointIndex,
    const Geometry::MappingReferenceToReference<1>* map,
    LinearAlgebra::SmallVector<DIM>& ret) const {
    logger.assert_debug(
        i < getNumberOfBasisFunctions(),
        "Asked for basis function %, but there are only % basis functions", i,
        getNumberOfBasisFunctions());
    const FE::BasisFunctionSet* subSet;
    std::size_t subIndex;
    std::tie(subSet, subIndex) = basisFunctions_.getBasisFunctionSetAndIndex(i);
    quadratureRule->eval(subSet, subIndex, quadraturePointIndex, map, ret);
}

template <std::size_t DIM>
void Element::basisFunction(
    std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule,
    std::size_t quadraturePointIndex,
    const Geometry::MappingReferenceToReference<1>* map,
    LinearAlgebra::SmallVector<DIM>& ret, std::size_t unknown) const {
    logger.assert_debug(
        i < getNumberOfBasisFunctions(unknown),
        "Asked for basis function %, but there are only % basis functions", i,
        getNumberOfBasisFunctions(unknown));
    const FE::BasisFunctionSet* subSet;
    std::size_t subIndex;
    std::tie(subSet, subIndex) =
        basisFunctions_.getBasisFunctionSetAndIndex(i, unknown);
    quadratureRule->eval(subSet, subIndex, quadraturePointIndex, map, ret);
}

template <std::size_t DIM>
LinearAlgebra::SmallVector<DIM> Element::basisFunctionDeriv(
    std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule,
    std::size_t quadraturePointIndex,
    const Geometry::MappingReferenceToReference<1>* map) const {
    logger.assert_debug(
        i < getNumberOfBasisFunctions(),
        "Asked for basis function %, but there are only % basis functions", i,
        getNumberOfBasisFunctions());
    const FE::BasisFunctionSet* subSet;
    std::size_t subIndex;
    std::tie(subSet, subIndex) = basisFunctions_.getBasisFunctionSetAndIndex(i);
    return quadratureRule->evalGrad(subSet, subIndex, quadraturePointIndex,
                                    map);
}

template <std::size_t DIM>
LinearAlgebra::SmallVector<DIM> Element::basisFunctionDeriv(
    std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule,
    std::size_t quadraturePointIndex,
    const Geometry::MappingReferenceToReference<1>* map,
    std::size_t unknown) const {
    logger.assert_debug(
        i < getNumberOfBasisFunctions(unknown),
        "Asked for basis function %, but there are only % basis functions", i,
        getNumberOfBasisFunctions(unknown));
    const FE::BasisFunctionSet* subSet;
    std::size_t subIndex;
    std::tie(subSet, subIndex) =
        basisFunctions_.getBasisFunctionSetAndIndex(i, unknown);
    return quadratureRule->evalGrad(subSet, subIndex, quadraturePointIndex,
                                    map);
}

template <std::size_t DIM>
LinearAlgebra::SmallVector<DIM> Element::basisFunctionCurl(
    std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule,
    std::size_t quadraturePointIndex,
    const Geometry::MappingReferenceToReference<1>* map) const {
    logger.assert_debug(
        i < getNumberOfBasisFunctions(),
        "Asked for basis function %, but there are only % basis functions", i,
        getNumberOfBasisFunctions());
    const FE::BasisFunctionSet* subSet;
    std::size_t subIndex;
    std::tie(subSet, subIndex) = basisFunctions_.getBasisFunctionSetAndIndex(i);
    return quadratureRule->evalCurl<DIM>(subSet, subIndex, quadraturePointIndex,
                                         map);
}

template <std::size_t DIM>
LinearAlgebra::SmallVector<DIM> Element::basisFunctionCurl(
    std::size_t i, QuadratureRules::GaussQuadratureRule* quadratureRule,
    std::size_t quadraturePointIndex,
    const Geometry::MappingReferenceToReference<1>* map,
    std::size_t unknown) const {
    logger.assert_debug(
        i < getNumberOfBasisFunctions(unknown),
        "Asked for basis function %, but there are only % basis functions", i,
        getNumberOfBasisFunctions(unknown));
    const FE::BasisFunctionSet* subSet;
    std::size_t subIndex;
    std::tie(subSet, subIndex) =
        basisFunctions_.getBasisFunctionSetAndIndex(i, unknown);
    return quadratureRule->evalCurl<DIM>(subSet, subIndex, quadraturePointIndex,
                                         map);
}
}  // namespace Base

}  // namespace hpgem

#endif  // HPGEM_KERNEL_ELEMENT_H

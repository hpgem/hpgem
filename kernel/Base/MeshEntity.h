/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2021, University of Twente
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
#ifndef HPGEM_MESHENTITY_H
#define HPGEM_MESHENTITY_H

#include <functional>
#include <memory>
#include <type_traits>

namespace hpgem {
namespace Base {

// Predeclaration to prevent circular dependency issues
class Edge;
class Element;
class Face;
class Node;

/**
 * Visitor for the MeshEntity type, could either be const or non const with
 * respect to the visited entities.
 * @tparam constify Whether the MeshEntities are const or not
 */
template <bool constify = false>
class MeshEntityVisitorT {
   public:
    virtual ~MeshEntityVisitorT() = default;
    /// Helper to add const
    template <typename T>
    using Constified =
        typename std::conditional<constify, typename std::add_const<T>::type,
                                  T>::type;

    virtual void visit(Constified<Base::Element>& element) = 0;
    virtual void visit(Constified<Base::Edge>& edge) = 0;
    virtual void visit(Constified<Base::Face>& face) = 0;
    virtual void visit(Constified<Base::Node>& node) = 0;
};

using MeshEntityVisitor = MeshEntityVisitorT<false>;
using ConstMeshEntityVisitor = MeshEntityVisitorT<true>;

/**
 * Entity of the Mesh, that is an Element, Face, Edge or Node.
 */
class MeshEntity {
   protected:
    // Deletion of MeshEntities should be governed by the actual instances.
    ~MeshEntity() = default;

   public:
    /// Enumeration of the possible mesh entities. Defines a bijection between
    /// the EntityTypes and the numbers 0 to N-1 (with N the number of
    /// EntityTypes). The ordering of the numbers in the bijection is
    /// unspecified.
    enum class EntityType : std::size_t {
        ELEMENT = 0,
        FACE = 1,
        EDGE = 2,
        NODE = 3
    };

    constexpr static std::size_t NUM_MESH_ENTITY_TYPES = 4;
    /**
     * Convert an entity type to a unique number from 0 to
     * NUM_MESH_ENTITY_TYPES-1.
     * @param type The type
     */
    constexpr static std::size_t entityTypeId(EntityType type) {
        return static_cast<std::size_t>(type);
    }

    virtual void accept(MeshEntityVisitor& vistor) = 0;
    virtual void accept(ConstMeshEntityVisitor& vistor) const = 0;

    /**
     * @return The identifier for this MeshEntity
     */
    virtual std::size_t getID() const = 0;

    virtual EntityType getType() const = 0;

    /**
     * In a distributed (MPI) problem each mesh entity has a single owning
     * process.
     * @return Whether the current process owns the MeshEntity.
     */
    virtual bool isOwnedByCurrentProcessor() const = 0;
    /**
     * Each MeshEntity is associated with a single element. In a distributed
     * setting this information is not always available due to the finite size
     * of the ghost layer. In such a setting three guarantees are given:
     * 1. If the owning element can not be determined the function will error.
     * 2. An owning element is always available if this MeshEntity is adjacent
     *    to an element that is owned by the current process.
     *
     * @return The owning element
     */
    virtual const Base::Element* getOwningElement() const = 0;

    /**
     * The number of basis functions per unknown that are associated with this
     * MeshEntity.
     *
     * Legacy behaviour: With multiple unknowns this will crash if the value is
     * different for the differnt unknowns.
     * @return
     */
    virtual std::size_t getLocalNumberOfBasisFunctions() const = 0;
    /**
     * The number of basis functions for a specific unknown that are associated
     * with this MeshEntity.
     * @param unknown The unknown to ask the number of basis functions of.
     * @return The number of basis functions
     */
    virtual std::size_t getLocalNumberOfBasisFunctions(
        std::size_t unknown) const = 0;
    /**
     * @return The number of basis functions for all unknowns that are
     * associated with this MeshEntity.
     */
    virtual std::size_t getTotalLocalNumberOfBasisFunctions() const = 0;
};

}  // namespace Base
}  // namespace hpgem

#endif  // HPGEM_MESHENTITY_H

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
#ifndef HPGEM_IDTYPES_H
#define HPGEM_IDTYPES_H

#include <memory>
#include <iostream>

namespace Preprocessor {

namespace {

/**
 * Wrapper around std::size_t to create type safe identifiers that
 * can't be mixed.
 * @tparam T A label to be used to distinguish between different ids.
 */
template <typename T>
struct Id {
    Id() = default;
    constexpr explicit Id(std::size_t id) : id(id){};

    std::size_t id;

    // Comparison operators
    inline bool operator==(const Id& other) const { return id == other.id; }
    inline bool operator!=(const Id& other) const { return id != other.id; }
    inline bool operator<(const Id& other) const { return id < other.id; }
    inline bool operator<=(const Id& other) const { return id <= other.id; }
    inline bool operator>(const Id& other) const { return id > other.id; }
    inline bool operator>=(const Id& other) const { return id >= other.id; }

    // Increment/decrement operators, primarily used in for-loops
    inline Id& operator++() {
        ++id;
        return *this;
    }

    Id operator++(int) {
        Id temp = *this;
        ++id;
        return temp;
    }

    inline Id& operator--() {
        --id;
        return *this;
    }

    Id operator--(int) {
        Id temp = *this;
        --id;
        return temp;
    }
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const Id<T>& id) {
    os << id.id;
    return os;
}

/**
 * Wrapper around std::size_t to create type safe identifiers that
 * can't be mixed, but with an implicit conversion to/from std::size_t. This
 * allows much easier looping and array indexing at the cost of safety.
 *
 * @tparam T A tag to distinguish different identifiers.
 */
template <typename T>
struct ImplicitId {
    ImplicitId() = default;
    // NOLINTNEXTLINE(google-explicit-constructor)
    constexpr ImplicitId(std::size_t id) : id(id){};

    std::size_t id;

    // Comparison operators
    inline bool operator==(const ImplicitId& other) const {
        return id == other.id;
    }
    inline bool operator!=(const ImplicitId& other) const {
        return id != other.id;
    }
    inline bool operator<(const ImplicitId& other) const {
        return id < other.id;
    }
    inline bool operator<=(const ImplicitId& other) const {
        return id <= other.id;
    }
    inline bool operator>(const ImplicitId& other) const {
        return id > other.id;
    }
    inline bool operator>=(const ImplicitId& other) const {
        return id >= other.id;
    }

    // Increment/decrement operators, primarily used in for-loops
    inline ImplicitId& operator++() {
        ++id;
        return *this;
    }

    ImplicitId operator++(int) {
        ImplicitId temp = *this;
        ++id;
        return temp;
    }

    inline ImplicitId& operator--() {
        --id;
        return *this;
    }

    ImplicitId operator--(int) {
        ImplicitId temp = *this;
        --id;
        return temp;
    }

    // NOLINTNEXTLINE(google-explicit-constructor)
    operator std::size_t() { return id; }
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const ImplicitId<T>& id) {
    os << id.id;
    return os;
}

struct CoordLabel {};
struct EntityGIdLabel {};
struct EntityLIdLabel {};
struct RegionIdLabel;

}  // namespace

/**
 * Identitifer for the physical coordinates in a mesh.
 */
using CoordId = Id<CoordLabel>;

/**
 * Identifier for a Mesh Entity. The combination of this id and the dimension of
 * the entity will uniquely identify the entity in a mesh.
 */
using EntityGId = Id<EntityGIdLabel>;
/**
 * Local identifier for a MeshEntity. For example, this ID type is used to point
 * to a node/edge/face of an element.
 */
using EntityLId = ImplicitId<EntityLIdLabel>;

/**
 * Identifier for the regions in a mesh. The combination of this id and the mesh
 * will uniquely identify the region.
 */
using RegionId = Id<RegionIdLabel>;

/**
 * Identifier for the lack of a region
 */
constexpr static RegionId NO_REGION_ID = RegionId (std::numeric_limits<std::size_t>::max());

}  // namespace Preprocessor

#endif  // HPGEM_REGIONMETA_H

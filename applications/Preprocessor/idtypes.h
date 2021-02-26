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

/// \brief Id for a physical coordinate
///
/// A numeric ID for physical coordinates in a mesh, wrapped in a class to
/// prevent accidental switching with other types of numeric indices.
struct CoordId {
    CoordId() = default;
    explicit CoordId(std::size_t id) : id(id) {}

    /// The actual id
    std::size_t id;

    inline bool operator==(const CoordId& other) const {
        return id == other.id;
    }
    inline bool operator!=(const CoordId& other) const {
        return id != other.id;
    }
    inline bool operator<(const CoordId& other) const { return id < other.id; }
    inline bool operator<=(const CoordId& other) const {
        return id <= other.id;
    }
    inline bool operator>(const CoordId& other) const { return id > other.id; }
    inline bool operator>=(const CoordId& other) const {
        return id >= other.id;
    }
    inline CoordId& operator++() {
        ++id;
        return *this;
    }
    CoordId operator++(int) {
        CoordId temp = *this;
        ++id;
        return temp;
    }

};

std::ostream& operator<<(std::ostream& os, const CoordId& coord) {
    os << coord.id;
    return os;
}


/// \brief Mesh-global ID for an entity
///
/// A numeric ID for a entity in a mesh. This id, combined with the dimension of
/// the entity will uniquely define an entity in the mesh.
struct EntityGId {
    EntityGId() = default;
    explicit EntityGId(std::size_t id) : id(id) {}

    /// The actual id
    std::size_t id;

    inline bool operator==(const EntityGId& other) const {
        return id == other.id;
    }
    inline bool operator!=(const EntityGId& other) const {
        return id != other.id;
    }
    inline bool operator<(const EntityGId& other) const {
        return id < other.id;
    }
    inline bool operator<=(const EntityGId& other) const {
        return id <= other.id;
    }
    inline bool operator>(const EntityGId& other) const {
        return id > other.id;
    }
    inline bool operator>=(const EntityGId& other) const {
        return id >= other.id;
    }
    inline EntityGId& operator++() {
        ++id;
        return *this;
    }
    EntityGId operator++(int) {
        EntityGId temp = *this;
        ++id;
        return temp;
    }
};

std::ostream& operator<<(std::ostream& os, const EntityGId& entityGId) {
    os << entityGId.id;
    return os;
}

/// \brief Local ID for an entity
///
/// A numeric ID for an entity, but in a local fashion. For example, this ID may
/// an index in the nodes/edges/faces of an element.
struct EntityLId {
    EntityLId() = default;
    EntityLId(std::size_t id) : id(id) {}

    /// The actual id
    std::size_t id;

    // Allow implicit conversion
    operator std::size_t() { return id; }

    inline bool operator==(const EntityLId& other) const {
        return id == other.id;
    }
    inline bool operator!=(const EntityLId& other) const {
        return id != other.id;
    }
    inline bool operator<(const EntityLId& other) const {
        return id < other.id;
    }
    inline bool operator<=(const EntityLId& other) const {
        return id <= other.id;
    }
    inline bool operator>(const EntityLId& other) const {
        return id > other.id;
    }
    inline bool operator>=(const EntityLId& other) const {
        return id >= other.id;
    }
    inline EntityLId& operator++() {
        ++id;
        return *this;
    }
    EntityLId operator++(int) {
        EntityLId temp = *this;
        ++id;
        return temp;
    }
};

std::ostream& operator<<(std::ostream& os, const EntityLId& entityLId) {
    os << entityLId.id;
    return os;
}

}  // namespace Preprocessor

#endif  // HPGEM_IDTYPES_H

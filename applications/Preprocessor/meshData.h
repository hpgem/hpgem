/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2017, University of Twente
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

#ifndef HPGEM_APP_MESHDATA_H
#define HPGEM_APP_MESHDATA_H

#include <cstddef>
#include "mesh.h"

using namespace hpgem;


namespace Preprocessor {

// todo: conversion operators

template <typename dataType, std::size_t meshDimension,
          std::size_t associatedDimension>
class MeshData {
   public:
    MeshData(Mesh<meshDimension>* mesh) : mesh(mesh) {
        data_.resize(mesh->template getNumberOfEntities<associatedDimension>());
    }
    MeshData(const MeshData&) = default;

    /**
     * MeshData gives the additional guarantee that moved-from MeshDatas are
     * still associated with their mesh after moving
     */
    MeshData(MeshData&& other) noexcept
        : mesh(other.mesh), data_(std::move(other.data_)) {}

    MeshData& operator=(const MeshData&) = default;
    MeshData& operator=(MeshData&& other) noexcept {
        mesh = other.mesh;
        data_ = std::move(other.data_);
        return *this;
    }
    ~MeshData() = default;

    dataType& operator[](std::size_t i) {
        logger.assert_debug(
            i < mesh->template getNumberOfEntities<associatedDimension>(),
            "There is not enough data in the mesh");
        data_.resize(mesh->template getNumberOfEntities<associatedDimension>());
        return data_[i];
    }

    dataType operator[](std::size_t i) const {
        logger.assert_debug(
            i < mesh->template getNumberOfEntities<associatedDimension>(),
            "There is not enough data in the mesh");
        if (i < data_.size()) return data_[i];
        return dataType{};
    }

    dataType& operator[](
        const MeshEntity<associatedDimension, meshDimension>& entity) {
        logger.assert_debug(mesh == entity.getMesh(),
                            "The entity does not belong to this mesh");
        return (*this)[entity.getGlobalIndex()];
    }

    dataType operator[](
        const MeshEntity<associatedDimension, meshDimension>& entity) const {
        logger.assert_debug(mesh == entity.getMesh(),
                            "The entity does not belong to this mesh");
        return (*this)[entity.getGlobalIndex()];
    }

    dataType* data() { return data_.data(); }

    std::size_t size() { return data_.size(); }

    // todo: think about appropriate iterator types for the results of begin/end
   private:
    Mesh<meshDimension>* mesh;
    std::vector<dataType> data_;
};

template <typename dataType, std::size_t meshDimension,
          std::size_t associatedDimension>
class PartialData {
   public:
    PartialData(Mesh<meshDimension>* mesh) : mesh(mesh) {}
    PartialData(const PartialData&) = default;

    /**
     * PartialData gives the additional guarantee that moved-from PartialDatas
     * are still associated with their mesh after moving
     */
    PartialData(PartialData&& other) noexcept
        : mesh(other.mesh), data_(std::move(other.data_)) {}

    PartialData& operator=(const PartialData&) = default;
    PartialData& operator=(PartialData&& other) noexcept {
        mesh = other.mesh;
        data_ = std::move(other.data_);
        return *this;
    }
    ~PartialData() = default;

    dataType& operator[](std::size_t i) {
        logger.assert_debug(
            i < mesh->template getNumberOfEntities<associatedDimension>(),
            "There is not enough data in the mesh");
        return data_[i];
    }

    dataType operator[](std::size_t i) const {
        logger.assert_debug(
            i < mesh->template getNumberOfEntities<associatedDimension>(),
            "There is not enough data in the mesh");
        try {
            return data_.at(i);
        } catch (std::out_of_range&) {
            return dataType{};
        }
    }

    dataType& operator[](
        const MeshEntity<associatedDimension, meshDimension>& entity) {
        logger.assert_debug(mesh == entity.getMesh(),
                            "The entity does not belong to this mesh");
        return (*this)[entity.getGlobalIndex()];
    }

    dataType operator[](
        const MeshEntity<associatedDimension, meshDimension>& entity) const {
        logger.assert_debug(mesh == entity.getMesh(),
                            "The entity does not belong to this mesh");
        return (*this)[entity.getGlobalIndex()];
    }

    // todo: think about appropriate iterator types for the results of begin/end
   private:
    Mesh<meshDimension>* mesh;
    std::map<std::size_t, dataType> data_{};
};
}  // namespace Preprocessor

#endif  // HPGEM_APP_MESHDATA_H

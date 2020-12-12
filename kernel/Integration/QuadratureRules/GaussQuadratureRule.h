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

#ifndef HPGEM_KERNEL_GAUSSQUADRATURERULE_H
#define HPGEM_KERNEL_GAUSSQUADRATURERULE_H

#include <string>
#include <cstring>
#include "Geometry/Mappings/ConcatenatedMapping.h"
namespace hpgem {
namespace Geometry {
// forward declaration
class ReferenceGeometry;
class PointReferenceBase;
template <int CODIM>
class MappingReferenceToReference;
}  // namespace Geometry

namespace Base {
class BasisFunctionSet;
}

/// helper class intended to store face to element maps in a way that is usable
/// by std::map (includes special treatment for concatenated mapping)
class faceMapContainer {
   public:
    faceMapContainer(const Geometry::MappingReferenceToReference<1>* map) {
        if (typeid(*map) == typeid(Geometry::ConcatenatedMapping)) {
            isConcatenated_ = true;
            auto concatenatedMap =
                static_cast<const Geometry::ConcatenatedMapping*>(map);
            data_.indirect_ = concatenatedMap->getSubMaps();
        } else {
            isConcatenated_ = false;
            data_.direct_ = map;
        }
    }

    bool operator<(const faceMapContainer& RHS) const {
        if (!isConcatenated_) {
            if (RHS.isConcatenated_) {
                return true;
            }
            return data_.direct_ < RHS.data_.direct_;

        } else {
            if (RHS.isConcatenated_) {
                return data_.indirect_ < RHS.data_.indirect_;
            }
            return false;
        }
    }

   private:
    // isConcatenated decides which field in the union is active
    // if it is ever changed, also make sure to write to the union to activate
    // the appropriate field
    bool isConcatenated_;
    union pointerContainer {
        const Geometry::MappingReferenceToReference<1>* direct_;
        std::pair<const Geometry::MappingReferenceToReference<0>*,
                  const Geometry::MappingReferenceToReference<1>*>
            indirect_;

        pointerContainer() { std::memset(this, 0, sizeof(pointerContainer)); }
    } data_;
};

namespace QuadratureRules {
class GaussQuadratureRule {
   public:
    GaussQuadratureRule() = default;
    GaussQuadratureRule(const GaussQuadratureRule& other) = delete;

    virtual ~GaussQuadratureRule() = default;

    //! Return the name of the quadrature.
    virtual std::string getName() const = 0;

    //! Return the dimension.
    virtual std::size_t dimension() const = 0;

    //! Return the order of the quadrature.
    virtual std::size_t order() const = 0;

    //! Return the number of points used in the quadrature.
    virtual std::size_t getNumberOfPoints() const = 0;

    //! Return the number of points used in the quadrature.
    //! \deprecated use getNumberOfPoints instead.
    std::size_t nrOfPoints() const { return getNumberOfPoints(); }

    //! Return the weight attached to the function value of the requested point
    //! number.
    virtual double weight(std::size_t) const = 0;

    //! Return the coordinates of the point with the given index.
    virtual const Geometry::PointReferenceBase& getPoint(std::size_t) const = 0;

    //! Each rule also knows which ReferenceGeometry it is meant for.
    virtual Geometry::ReferenceGeometry* forReferenceGeometry() const = 0;

    /// tell the quadrature rule that a pointer to a basis function set is no
    /// longer suitable for quick lookup
    void unregisterBasisFunctionSet(const Base::BasisFunctionSet* set) {
        basisFunctionValues_.erase(set);
        basisFunctionGrads_.erase(set);
        basisFunctionCurls_.erase(set);
        basisFunctionDivs_.erase(set);
        faceBasisFunctionValues_.erase(set);
        faceBasisFunctionGrads_.erase(set);
        faceBasisFunctionCurls_.erase(set);
        faceBasisFunctionDivs_.erase(set);
    }

    /// pre-evaluate a set of basisfunctions to speed up computation
    double eval(const Base::BasisFunctionSet* set,
                std::size_t basisFunctionIndex,
                std::size_t quadraturePointIndex);

    /// pre-evaluate a set of basisfunctions to speed up computation. First maps
    /// the quadrature points to an element using the provided mapping
    double eval(const Base::BasisFunctionSet* set,
                std::size_t basisFunctionIndex,
                std::size_t quadraturePointIndex,
                const Geometry::MappingReferenceToReference<1>* map);

    /// pre-evaluate a set of basisfunctions to speed up computation
    template <std::size_t DIM>
    void eval(const Base::BasisFunctionSet* set, std::size_t basisFunctionIndex,
              std::size_t quadraturePointIndex,
              LinearAlgebra::SmallVector<DIM>& result);

    /// pre-evaluate a set of basisfunctions to speed up computation. First maps
    /// the quadrature points to an element using the provided mapping
    template <std::size_t DIM>
    void eval(const Base::BasisFunctionSet* set, std::size_t basisFunctionIndex,
              std::size_t quadraturePointIndex,
              const Geometry::MappingReferenceToReference<1>* map,
              LinearAlgebra::SmallVector<DIM>& result);

    /// pre-evaluate the derivative of a set of basisfunctions to speed up
    /// computation the result of this function should be assigned back into a
    /// smallvector of appropriate size
    const LinearAlgebra::MiddleSizeVector& evalGrad(
        const Base::BasisFunctionSet* set, std::size_t basisFunctionIndex,
        std::size_t quadraturePointIndex);

    /// pre-evaluate the derivative of a set of basisfunctions to speed up
    /// computation. First maps the quadrature points to an element using the
    /// provided mapping the result of this function should be assigned back
    /// into a smallvector of appropriate size
    const LinearAlgebra::MiddleSizeVector& evalGrad(
        const Base::BasisFunctionSet* set, std::size_t basisFunctionIndex,
        std::size_t quadraturePointIndex,
        const Geometry::MappingReferenceToReference<1>* map);

    /// pre-evaluate the curl of a set of basisfunctions to speed up computation
    template <std::size_t DIM>
    LinearAlgebra::SmallVector<DIM> evalCurl(const Base::BasisFunctionSet* set,
                                             std::size_t basisFunctionIndex,
                                             std::size_t quadraturePointIndex) {
        logger(ERROR, "Curl only valid for 2D and 3D not in %", dimension());
        LinearAlgebra::SmallVector<DIM> dummy;
        return dummy;
    }

    /// pre-evaluate the curl of a set of basisfunctions to speed up
    /// computation. First maps the quadrature points to an element using the
    /// provided mapping
    template <std::size_t DIM>
    LinearAlgebra::SmallVector<DIM> evalCurl(
        const Base::BasisFunctionSet* set, std::size_t basisFunctionIndex,
        std::size_t quadraturePointIndex,
        const Geometry::MappingReferenceToReference<1>* map) {
        logger(ERROR, "Curl only valid for 2D and 3D not in %", dimension());
        LinearAlgebra::SmallVector<DIM> dummy;
        return dummy;
    }

    /// pre-evaluate the divergence of a set of basisfunctions to speed up
    /// computation
    double evalDiv(const Base::BasisFunctionSet* set,
                   std::size_t basisFunctionIndex,
                   std::size_t quadraturePointIndex);

    /// pre-evaluate the divergence of a set of basisfunctions to speed up
    /// computation. First maps the quadrature points to an element using the
    /// provided mapping
    double evalDiv(const Base::BasisFunctionSet* set,
                   std::size_t basisFunctionIndex,
                   std::size_t quadraturePointIndex,
                   const Geometry::MappingReferenceToReference<1>* map);

   private:
    // this uses maps because if you can swap two vectors of vectors efficiently
    // you lose the benefit of consecutive memory access anyway and there are
    // relatively few entries created, but a lot of lookups so the allocation
    // overhead is relatively minor
    std::map<const Base::BasisFunctionSet*, std::vector<std::vector<double>>>
        basisFunctionValues_;
    std::map<const Base::BasisFunctionSet*,
             std::vector<std::vector<LinearAlgebra::MiddleSizeVector>>>
        basisFunctionVectorValues_;
    std::map<const Base::BasisFunctionSet*,
             std::vector<std::vector<LinearAlgebra::MiddleSizeVector>>>
        basisFunctionGrads_;
    std::map<const Base::BasisFunctionSet*,
             std::vector<std::vector<LinearAlgebra::SmallVector<3>>>>
        basisFunctionCurls_;
    std::map<const Base::BasisFunctionSet*,
             std::vector<std::vector<LinearAlgebra::SmallVector<2>>>>
        basisFunctionCurls2D_;
    std::map<const Base::BasisFunctionSet*, std::vector<std::vector<double>>>
        basisFunctionDivs_;
    std::map<const Base::BasisFunctionSet*,
             std::map<faceMapContainer, std::vector<std::vector<double>>>>
        faceBasisFunctionValues_;
    std::map<
        const Base::BasisFunctionSet*,
        std::map<faceMapContainer,
                 std::vector<std::vector<LinearAlgebra::MiddleSizeVector>>>>
        faceBasisFunctionVectorValues_;
    std::map<
        const Base::BasisFunctionSet*,
        std::map<faceMapContainer,
                 std::vector<std::vector<LinearAlgebra::MiddleSizeVector>>>>
        faceBasisFunctionGrads_;
    std::map<const Base::BasisFunctionSet*,
             std::map<faceMapContainer,
                      std::vector<std::vector<LinearAlgebra::SmallVector<3>>>>>
        faceBasisFunctionCurls_;
    std::map<const Base::BasisFunctionSet*,
             std::map<faceMapContainer,
                      std::vector<std::vector<LinearAlgebra::SmallVector<2>>>>>
        faceBasisFunctionCurls2D_;
    std::map<const Base::BasisFunctionSet*,
             std::map<faceMapContainer, std::vector<std::vector<double>>>>
        faceBasisFunctionDivs_;
};
}  // namespace QuadratureRules
}  // namespace hpgem
#include "FE/BasisFunctionSet.h"

namespace hpgem {

template <std::size_t DIM>
inline void QuadratureRules::GaussQuadratureRule::eval(
    const Base::BasisFunctionSet* set, std::size_t basisFunctionIndex,
    std::size_t quadraturePointIndex, LinearAlgebra::SmallVector<DIM>& result) {
    logger.assert_debug(set != nullptr, "Invalid basis function set passed");
    logger.assert_debug(dimension() == DIM,
                        "Quadrature rule not for this dimension");
    logger.assert_debug(quadraturePointIndex < getNumberOfPoints(),
                        "Asked for point %, but this rule only has % points",
                        quadraturePointIndex, getNumberOfPoints());
    logger.assert_debug(basisFunctionIndex < set->size(),
                        "Asked for basis function %, but the provided basis "
                        "function set only has % points",
                        basisFunctionIndex, set->size());
    try {
        result = basisFunctionVectorValues_.at(
            set)[quadraturePointIndex][basisFunctionIndex];
    } catch (std::out_of_range&) {
        // we store smallVectors as middleSizeVectors so we dont have to
        // template the quadrature rule, but this means we have to silence the
        // efficiency warning efficiency is not a big issue here since we only
        // do a heap allocation once per basis function per quadrature point for
        // the entire computation
        auto oldWarn = loggerOutput->onWarn;
        loggerOutput->onWarn = [](std::string, std::string) {};
        set->registerDestructorListener(
            [set, this] { unregisterBasisFunctionSet(set); });
        basisFunctionVectorValues_[set].resize(getNumberOfPoints());
        for (std::size_t i = 0; i < getNumberOfPoints(); ++i) {
            basisFunctionVectorValues_[set][i].resize(set->size());
            const Geometry::PointReference<DIM>& point = getPoint(i);
            for (std::size_t j = 0; j < set->size(); ++j) {
                // borrow result for transferring the information; set to the
                // appropriate value later
                set->eval(j, point, result);
                basisFunctionVectorValues_[set][i][j] = result;
            }
        }
        loggerOutput->onWarn = oldWarn;
        result = basisFunctionVectorValues_[set][quadraturePointIndex]
                                           [basisFunctionIndex];
    }
}

template <std::size_t DIM>
inline void QuadratureRules::GaussQuadratureRule::eval(
    const Base::BasisFunctionSet* set, std::size_t basisFunctionIndex,
    std::size_t quadraturePointIndex,
    const Geometry::MappingReferenceToReference<1>* map,
    LinearAlgebra::SmallVector<DIM>& result) {
    logger.assert_debug(set != nullptr, "Invalid basis function set passed");
    logger.assert_debug(dimension() == DIM - 1,
                        "Quadrature rule not for this dimension");
    logger.assert_debug(map != nullptr,
                        "Invalid coordinate transformation passed");
    logger.assert_debug(quadraturePointIndex < getNumberOfPoints(),
                        "Asked for point %, but this rule only has % points",
                        quadraturePointIndex, getNumberOfPoints());
    logger.assert_debug(basisFunctionIndex < set->size(),
                        "Asked for basis function %, but the provided basis "
                        "function set only has % points",
                        basisFunctionIndex, set->size());
    auto containedMap = faceMapContainer(map);
    try {
        result = faceBasisFunctionVectorValues_.at(set).at(
            containedMap)[quadraturePointIndex][basisFunctionIndex];
    } catch (std::out_of_range&) {
        // we store smallVectors as middleSizeVectors so we dont have to
        // template the quadrature rule, but this means we have to silence the
        // efficiency warning efficiency is not a big issue here since we only
        // do a heap allocation once per basis function per quadrature point for
        // the entire computation
        auto oldWarn = loggerOutput->onWarn;
        loggerOutput->onWarn = [](std::string, std::string) {};
        set->registerDestructorListener(
            [set, this] { unregisterBasisFunctionSet(set); });
        faceBasisFunctionVectorValues_[set][containedMap].resize(
            getNumberOfPoints());
        for (std::size_t i = 0; i < getNumberOfPoints(); ++i) {
            faceBasisFunctionVectorValues_[set][containedMap][i].resize(
                set->size());
            const Geometry::PointReference<DIM - 1>& facePoint = getPoint(i);
            const Geometry::PointReference<DIM>& point =
                map->transform(facePoint);
            for (std::size_t j = 0; j < set->size(); ++j) {
                // borrow result for transferring the information; set to the
                // appropriate value later
                set->eval(j, point, result);
                faceBasisFunctionVectorValues_[set][containedMap][i][j] =
                    result;
            }
        }
        loggerOutput->onWarn = oldWarn;
        result = faceBasisFunctionVectorValues_[set][containedMap]
                                               [quadraturePointIndex]
                                               [basisFunctionIndex];
    }
}

template <>
inline LinearAlgebra::SmallVector<3>
    QuadratureRules::GaussQuadratureRule::evalCurl(
        const Base::BasisFunctionSet* set, std::size_t basisFunctionIndex,
        std::size_t quadraturePointIndex) {
    logger.assert_debug(set != nullptr, "Invalid basis function set passed");
    logger.assert_debug(quadraturePointIndex < getNumberOfPoints(),
                        "Asked for point %, but this rule only has % points",
                        quadraturePointIndex, getNumberOfPoints());
    logger.assert_debug(basisFunctionIndex < set->size(),
                        "Asked for basis function %, but the provided basis "
                        "function set only has % points",
                        basisFunctionIndex, set->size());
    try {
        return basisFunctionCurls_.at(
            set)[quadraturePointIndex][basisFunctionIndex];
    } catch (std::out_of_range&) {
        set->registerDestructorListener(
            [set, this] { unregisterBasisFunctionSet(set); });
        basisFunctionCurls_[set].resize(getNumberOfPoints());
        for (std::size_t i = 0; i < getNumberOfPoints(); ++i) {
            basisFunctionCurls_[set][i].resize(set->size());
            switch (dimension()) {
                case 3: {
                    const Geometry::PointReference<3>& point3D = getPoint(i);
                    for (std::size_t j = 0; j < set->size(); ++j) {
                        basisFunctionCurls_[set][i][j] =
                            set->evalCurl(j, point3D);
                    }
                } break;
                default:
                    logger(ERROR, "curl is only defined in R^3");
            }
        }
        return basisFunctionCurls_[set][quadraturePointIndex]
                                  [basisFunctionIndex];
    }
}

template <>
inline LinearAlgebra::SmallVector<2>
    QuadratureRules::GaussQuadratureRule::evalCurl(
        const Base::BasisFunctionSet* set, std::size_t basisFunctionIndex,
        std::size_t quadraturePointIndex) {
    logger.assert_debug(set != nullptr, "Invalid basis function set passed");
    logger.assert_debug(quadraturePointIndex < getNumberOfPoints(),
                        "Asked for point %, but this rule only has % points",
                        quadraturePointIndex, getNumberOfPoints());
    logger.assert_debug(basisFunctionIndex < set->size(),
                        "Asked for basis function %, but the provided basis "
                        "function set only has % points",
                        basisFunctionIndex, set->size());
    try {
        return basisFunctionCurls2D_.at(
            set)[quadraturePointIndex][basisFunctionIndex];
    } catch (std::out_of_range&) {
        set->registerDestructorListener(
            [set, this] { unregisterBasisFunctionSet(set); });
        basisFunctionCurls2D_[set].resize(getNumberOfPoints());
        for (std::size_t i = 0; i < getNumberOfPoints(); ++i) {
            basisFunctionCurls2D_[set][i].resize(set->size());
            switch (dimension()) {
                case 2: {
                    const Geometry::PointReference<2>& point2D = getPoint(i);
                    for (std::size_t j = 0; j < set->size(); ++j) {
                        basisFunctionCurls2D_[set][i][j] =
                            set->evalCurl(j, point2D);
                    }
                } break;
                default:
                    logger(ERROR, "curl is only defined in R^2");
            }
        }
        return basisFunctionCurls2D_[set][quadraturePointIndex]
                                    [basisFunctionIndex];
    }
}

template <>
inline LinearAlgebra::SmallVector<3>
    QuadratureRules::GaussQuadratureRule::evalCurl(
        const Base::BasisFunctionSet* set, std::size_t basisFunctionIndex,
        std::size_t quadraturePointIndex,
        const Geometry::MappingReferenceToReference<1>* map) {
    logger.assert_debug(set != nullptr, "Invalid basis function set passed");
    logger.assert_debug(map != nullptr,
                        "Invalid coordinate transformation passed");
    logger.assert_debug(quadraturePointIndex < getNumberOfPoints(),
                        "Asked for point %, but this rule only has % points",
                        quadraturePointIndex, getNumberOfPoints());
    logger.assert_debug(basisFunctionIndex < set->size(),
                        "Asked for basis function %, but the provided basis "
                        "function set only has % points",
                        basisFunctionIndex, set->size());
    auto containedMap = faceMapContainer(map);
    try {
        return faceBasisFunctionCurls_.at(set).at(
            containedMap)[quadraturePointIndex][basisFunctionIndex];
    } catch (std::out_of_range&) {
        set->registerDestructorListener(
            [set, this] { unregisterBasisFunctionSet(set); });
        faceBasisFunctionCurls_[set][containedMap].resize(getNumberOfPoints());
        for (std::size_t i = 0; i < getNumberOfPoints(); ++i) {
            faceBasisFunctionCurls_[set][containedMap][i].resize(set->size());
            switch (dimension()) {
                case 2: {
                    const Geometry::PointReference<2>& facePoint2D =
                        getPoint(i);
                    const Geometry::PointReference<3>& point3D =
                        map->transform(facePoint2D);
                    for (std::size_t j = 0; j < set->size(); ++j) {
                        faceBasisFunctionCurls_[set][containedMap][i][j] =
                            set->evalCurl(j, point3D);
                    }
                } break;
                default:
                    logger(ERROR, "curl is only defined in R^3", dimension());
            }
        }
        return faceBasisFunctionCurls_[set][containedMap][quadraturePointIndex]
                                      [basisFunctionIndex];
    }
}

template <>
inline LinearAlgebra::SmallVector<2>
    QuadratureRules::GaussQuadratureRule::evalCurl(
        const Base::BasisFunctionSet* set, std::size_t basisFunctionIndex,
        std::size_t quadraturePointIndex,
        const Geometry::MappingReferenceToReference<1>* map) {
    logger.assert_debug(set != nullptr, "Invalid basis function set passed");
    logger.assert_debug(map != nullptr,
                        "Invalid coordinate transformation passed");
    logger.assert_debug(quadraturePointIndex < getNumberOfPoints(),
                        "Asked for point %, but this rule only has % points",
                        quadraturePointIndex, getNumberOfPoints());
    logger.assert_debug(basisFunctionIndex < set->size(),
                        "Asked for basis function %, but the provided basis "
                        "function set only has % points",
                        basisFunctionIndex, set->size());
    auto containedMap = faceMapContainer(map);
    try {
        return faceBasisFunctionCurls2D_.at(set).at(
            containedMap)[quadraturePointIndex][basisFunctionIndex];
    } catch (std::out_of_range&) {
        set->registerDestructorListener(
            [set, this] { unregisterBasisFunctionSet(set); });
        faceBasisFunctionCurls2D_[set][containedMap].resize(
            getNumberOfPoints());
        for (std::size_t i = 0; i < getNumberOfPoints(); ++i) {
            faceBasisFunctionCurls2D_[set][containedMap][i].resize(set->size());
            switch (dimension()) {
                case 1: {
                    const Geometry::PointReference<1>& facePoint1D =
                        getPoint(i);
                    const Geometry::PointReference<2>& point2D =
                        map->transform(facePoint1D);
                    for (std::size_t j = 0; j < set->size(); ++j) {
                        faceBasisFunctionCurls2D_[set][containedMap][i][j] =
                            set->evalCurl(j, point2D);
                    }
                } break;
                default:
                    logger(ERROR, "curl is only defined in R^2", dimension());
            }
        }
        return faceBasisFunctionCurls2D_[set][containedMap]
                                        [quadraturePointIndex]
                                        [basisFunctionIndex];
    }
}

}  // namespace hpgem

#endif  // HPGEM_KERNEL_GAUSSQUADRATURERULE_H

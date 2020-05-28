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

#include "GaussQuadratureRule.h"
#include "Geometry/PointReference.h"

double QuadratureRules::GaussQuadratureRule::eval(
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
        return basisFunctionValues_.at(
            set)[quadraturePointIndex][basisFunctionIndex];
    } catch (std::out_of_range&) {
        set->registerQuadratureRule(this);
        basisFunctionValues_[set].resize(getNumberOfPoints());
        for (std::size_t i = 0; i < getNumberOfPoints(); ++i) {
            basisFunctionValues_[set][i].resize(set->size());
            switch (dimension()) {
                case 1: {
                    const Geometry::PointReference<1>& point1D = getPoint(i);
                    for (std::size_t j = 0; j < set->size(); ++j) {
                        basisFunctionValues_[set][i][j] = set->eval(j, point1D);
                    }
                } break;
                case 2: {
                    const Geometry::PointReference<2>& point2D = getPoint(i);
                    for (std::size_t j = 0; j < set->size(); ++j) {
                        basisFunctionValues_[set][i][j] = set->eval(j, point2D);
                    }
                } break;
                case 3: {
                    const Geometry::PointReference<3>& point3D = getPoint(i);
                    for (std::size_t j = 0; j < set->size(); ++j) {
                        basisFunctionValues_[set][i][j] = set->eval(j, point3D);
                    }
                } break;
                case 4: {
                    const Geometry::PointReference<4>& point4D = getPoint(i);
                    for (std::size_t j = 0; j < set->size(); ++j) {
                        basisFunctionValues_[set][i][j] = set->eval(j, point4D);
                    }
                } break;
                default:
                    logger(ERROR,
                           "the dimension of the quadrature rule is unsuitable "
                           "for direct evaluation of the basis function, "
                           "please also provide a face to element map");
            }
        }
        return basisFunctionValues_[set][quadraturePointIndex]
                                   [basisFunctionIndex];
    }
}

double QuadratureRules::GaussQuadratureRule::eval(
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
        return faceBasisFunctionValues_.at(set).at(
            containedMap)[quadraturePointIndex][basisFunctionIndex];
    } catch (std::out_of_range&) {
        set->registerQuadratureRule(this);
        faceBasisFunctionValues_[set][containedMap].resize(getNumberOfPoints());
        for (std::size_t i = 0; i < getNumberOfPoints(); ++i) {
            faceBasisFunctionValues_[set][containedMap][i].resize(set->size());
            switch (dimension()) {
                case 0: {
                    const Geometry::PointReference<0>& facePoint0D =
                        getPoint(i);
                    const Geometry::PointReference<1>& point1D =
                        map->transform(facePoint0D);
                    for (std::size_t j = 0; j < set->size(); ++j) {
                        faceBasisFunctionValues_[set][containedMap][i][j] =
                            set->eval(j, point1D);
                    }
                } break;
                case 1: {
                    const Geometry::PointReference<1>& facePoint1D =
                        getPoint(i);
                    const Geometry::PointReference<2>& point2D =
                        map->transform(facePoint1D);
                    for (std::size_t j = 0; j < set->size(); ++j) {
                        faceBasisFunctionValues_[set][containedMap][i][j] =
                            set->eval(j, point2D);
                    }
                } break;
                case 2: {
                    const Geometry::PointReference<2>& facePoint2D =
                        getPoint(i);
                    const Geometry::PointReference<3>& point3D =
                        map->transform(facePoint2D);
                    for (std::size_t j = 0; j < set->size(); ++j) {
                        faceBasisFunctionValues_[set][containedMap][i][j] =
                            set->eval(j, point3D);
                    }
                } break;
                case 3: {
                    const Geometry::PointReference<3>& facePoint3D =
                        getPoint(i);
                    const Geometry::PointReference<4>& point4D =
                        map->transform(facePoint3D);
                    for (std::size_t j = 0; j < set->size(); ++j) {
                        faceBasisFunctionValues_[set][containedMap][i][j] =
                            set->eval(j, point4D);
                    }
                } break;
                default:
                    logger(ERROR, "hpGEM does not support faces of dimension %",
                           dimension());
            }
        }
        return faceBasisFunctionValues_[set][containedMap][quadraturePointIndex]
                                       [basisFunctionIndex];
    }
}

const LinearAlgebra::MiddleSizeVector&
    QuadratureRules::GaussQuadratureRule::evalGrad(
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
        return basisFunctionGrads_.at(
            set)[quadraturePointIndex][basisFunctionIndex];
    } catch (std::out_of_range&) {
        // we store smallVectors as middleSizeVectors so we dont have to
        // template the quadrature rule, but this means we have to silence the
        // efficiency warning efficiency is not a big issue here since we only do
        // a heap allocation once per basis function per quadrature point for the
        // entire computation
        auto oldWarn = loggerOutput->onWarn;
        loggerOutput->onWarn = [](std::string, std::string) {};
        set->registerQuadratureRule(this);
        basisFunctionGrads_[set].resize(getNumberOfPoints());
        for (std::size_t i = 0; i < getNumberOfPoints(); ++i) {
            basisFunctionGrads_[set][i].resize(set->size());
            switch (dimension()) {
                case 1: {
                    const Geometry::PointReference<1>& point1D = getPoint(i);
                    for (std::size_t j = 0; j < set->size(); ++j) {
                        basisFunctionGrads_[set][i][j] =
                            set->evalDeriv(j, point1D);
                    }
                } break;
                case 2: {
                    const Geometry::PointReference<2>& point2D = getPoint(i);
                    for (std::size_t j = 0; j < set->size(); ++j) {
                        basisFunctionGrads_[set][i][j] =
                            set->evalDeriv(j, point2D);
                    }
                } break;
                case 3: {
                    const Geometry::PointReference<3>& point3D = getPoint(i);
                    for (std::size_t j = 0; j < set->size(); ++j) {
                        basisFunctionGrads_[set][i][j] =
                            set->evalDeriv(j, point3D);
                    }
                } break;
                case 4: {
                    const Geometry::PointReference<4>& point4D = getPoint(i);
                    for (std::size_t j = 0; j < set->size(); ++j) {
                        basisFunctionGrads_[set][i][j] =
                            set->evalDeriv(j, point4D);
                    }
                } break;
                default:
                    logger(ERROR,
                           "the dimension of the quadrature rule is unsuitable "
                           "for direct evaluation of the basis function, "
                           "please also provide a face to element map");
            }
        }
        loggerOutput->onWarn = oldWarn;
        return basisFunctionGrads_[set][quadraturePointIndex]
                                  [basisFunctionIndex];
    }
}

const LinearAlgebra::MiddleSizeVector&
    QuadratureRules::GaussQuadratureRule::evalGrad(
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
        return faceBasisFunctionGrads_.at(set).at(
            containedMap)[quadraturePointIndex][basisFunctionIndex];
    } catch (std::out_of_range&) {
        // we store smallVectors as middleSizeVectors so we dont have to
        // template the quadrature rule, but this means we have to silence the
        // efficiency warning efficiency is not a big issue here since we only do
        // a heap allocation once per basis function per quadrature point for the
        // entire computation
        auto oldWarn = loggerOutput->onWarn;
        loggerOutput->onWarn = [](std::string, std::string) {};
        set->registerQuadratureRule(this);
        faceBasisFunctionGrads_[set][containedMap].resize(getNumberOfPoints());
        for (std::size_t i = 0; i < getNumberOfPoints(); ++i) {
            faceBasisFunctionGrads_[set][containedMap][i].resize(set->size());
            switch (dimension()) {
                case 0: {
                    const Geometry::PointReference<0>& facePoint0D =
                        getPoint(i);
                    const Geometry::PointReference<1>& point1D =
                        map->transform(facePoint0D);
                    for (std::size_t j = 0; j < set->size(); ++j) {
                        faceBasisFunctionGrads_[set][containedMap][i][j] =
                            set->evalDeriv(j, point1D);
                    }
                } break;
                case 1: {
                    const Geometry::PointReference<1>& facePoint1D =
                        getPoint(i);
                    const Geometry::PointReference<2>& point2D =
                        map->transform(facePoint1D);
                    for (std::size_t j = 0; j < set->size(); ++j) {
                        faceBasisFunctionGrads_[set][containedMap][i][j] =
                            set->evalDeriv(j, point2D);
                    }
                } break;
                case 2: {
                    const Geometry::PointReference<2>& facePoint2D =
                        getPoint(i);
                    const Geometry::PointReference<3>& point3D =
                        map->transform(facePoint2D);
                    for (std::size_t j = 0; j < set->size(); ++j) {
                        faceBasisFunctionGrads_[set][containedMap][i][j] =
                            set->evalDeriv(j, point3D);
                    }
                } break;
                case 3: {
                    const Geometry::PointReference<3>& facePoint3D =
                        getPoint(i);
                    const Geometry::PointReference<4>& point4D =
                        map->transform(facePoint3D);
                    for (std::size_t j = 0; j < set->size(); ++j) {
                        faceBasisFunctionGrads_[set][containedMap][i][j] =
                            set->evalDeriv(j, point4D);
                    }
                } break;
                default:
                    logger(ERROR, "hpGEM does not support faces of dimension %",
                           dimension());
            }
        }
        loggerOutput->onWarn = oldWarn;
        return faceBasisFunctionGrads_[set][containedMap][quadraturePointIndex]
                                      [basisFunctionIndex];
    }
}

LinearAlgebra::SmallVector<3> QuadratureRules::GaussQuadratureRule::evalCurl(
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
        set->registerQuadratureRule(this);
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

LinearAlgebra::SmallVector<2> QuadratureRules::GaussQuadratureRule::evalCurl2D(
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
        set->registerQuadratureRule(this);
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

LinearAlgebra::SmallVector<3> QuadratureRules::GaussQuadratureRule::evalCurl(
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
        set->registerQuadratureRule(this);
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

LinearAlgebra::SmallVector<2> QuadratureRules::GaussQuadratureRule::evalCurl2D(
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
        set->registerQuadratureRule(this);
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

double QuadratureRules::GaussQuadratureRule::evalDiv(
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
        return basisFunctionDivs_.at(
            set)[quadraturePointIndex][basisFunctionIndex];
    } catch (std::out_of_range&) {
        set->registerQuadratureRule(this);
        basisFunctionDivs_[set].resize(getNumberOfPoints());
        for (std::size_t i = 0; i < getNumberOfPoints(); ++i) {
            basisFunctionDivs_[set][i].resize(set->size());
            switch (dimension()) {
                case 2: {
                    const Geometry::PointReference<2>& point2D = getPoint(i);
                    for (std::size_t j = 0; j < set->size(); ++j) {
                        basisFunctionDivs_[set][i][j] =
                            set->evalDiv(j, point2D);
                    }
                } break;
                case 3: {
                    const Geometry::PointReference<3>& point3D = getPoint(i);
                    for (std::size_t j = 0; j < set->size(); ++j) {
                        basisFunctionDivs_[set][i][j] =
                            set->evalDiv(j, point3D);
                    }
                } break;
                case 4: {
                    const Geometry::PointReference<4>& point4D = getPoint(i);
                    for (std::size_t j = 0; j < set->size(); ++j) {
                        basisFunctionDivs_[set][i][j] =
                            set->evalDiv(j, point4D);
                    }
                } break;
                default:
                    logger(ERROR,
                           "the dimension of the quadrature rule is unsuitable "
                           "for direct evaluation of a divergence, please also "
                           "provide a face to element map");
            }
        }
        return basisFunctionDivs_[set][quadraturePointIndex]
                                 [basisFunctionIndex];
    }
}

double QuadratureRules::GaussQuadratureRule::evalDiv(
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
        return faceBasisFunctionDivs_.at(set).at(
            containedMap)[quadraturePointIndex][basisFunctionIndex];
    } catch (std::out_of_range&) {
        set->registerQuadratureRule(this);
        faceBasisFunctionDivs_[set][containedMap].resize(getNumberOfPoints());
        for (std::size_t i = 0; i < getNumberOfPoints(); ++i) {
            faceBasisFunctionDivs_[set][containedMap][i].resize(set->size());
            switch (dimension()) {
                case 1: {
                    const Geometry::PointReference<1>& facePoint1D =
                        getPoint(i);
                    const Geometry::PointReference<2>& point2D =
                        map->transform(facePoint1D);
                    for (std::size_t j = 0; j < set->size(); ++j) {
                        faceBasisFunctionDivs_[set][containedMap][i][j] =
                            set->evalDiv(j, point2D);
                    }
                } break;
                case 2: {
                    const Geometry::PointReference<2>& facePoint2D =
                        getPoint(i);
                    const Geometry::PointReference<3>& point3D =
                        map->transform(facePoint2D);
                    for (std::size_t j = 0; j < set->size(); ++j) {
                        faceBasisFunctionDivs_[set][containedMap][i][j] =
                            set->evalDiv(j, point3D);
                    }
                } break;
                case 3: {
                    const Geometry::PointReference<3>& facePoint3D =
                        getPoint(i);
                    const Geometry::PointReference<4>& point4D =
                        map->transform(facePoint3D);
                    for (std::size_t j = 0; j < set->size(); ++j) {
                        faceBasisFunctionDivs_[set][containedMap][i][j] =
                            set->evalDiv(j, point4D);
                    }
                } break;
                case 0:
                    logger(ERROR, "dimension 1 does not have a div-operator");
                    break;
                default:
                    logger(ERROR, "hpGEM does not support faces of dimension %",
                           dimension());
            }
        }
        return faceBasisFunctionDivs_[set][containedMap][quadraturePointIndex]
                                     [basisFunctionIndex];
    }
}

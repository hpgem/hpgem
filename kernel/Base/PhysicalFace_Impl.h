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

#include "CoordinateTransformation.h"

namespace hpgem {

namespace Base {

template <std::size_t DIM>
inline double PhysicalFace<DIM>::basisFunction(std::size_t i) {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    if (i < nLeftBasisFunctions[0]) {
        return left.basisFunction(i);
    }
    logger.assert_debug(isInternal_, "basis function index out of bounds");
    double value = right.basisFunction(i - nLeftBasisFunctions[0]);
    if (requiresTransformation) {
        value = transform_[0]->transform(value, rightVectorTransform);
        if (std::abs(getUnitNormalVector()[2]) > 0.5) {
            // Assume negative mirror character
            value *= -1.00;;
        }
    }
    return value;
}

template <std::size_t DIM>
inline double PhysicalFace<DIM>::basisFunction(std::size_t i,
                                               std::size_t unknown) {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    if (i < nLeftBasisFunctions[unknown]) {
        return left.basisFunction(i, unknown);
    }
    logger.assert_debug(isInternal_, "basis function index out of bounds");
    double value =
        right.basisFunction(i - nLeftBasisFunctions[unknown], unknown);
    if (requiresTransformation) {
        value = transform_[unknown]->transform(value, rightVectorTransform);
        if (std::abs(getUnitNormalVector()[2]) > 0.5) {
            // Assume negative mirror character
            value *= -1.00;;
        }
    }
    return value;
}

template <std::size_t DIM>
inline double PhysicalFace<DIM>::basisFunction(Side side, std::size_t i) {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    if (side == Side::LEFT) {
        return left.basisFunction(i);
    }
    logger.assert_debug(isInternal_,
                        "cannot find the right element for a boundary face");
    double value = right.basisFunction(i);
    if (requiresTransformation) {
        value = transform_[0]->transform(value, rightVectorTransform);
        if (std::abs(getUnitNormalVector()[2]) > 0.5) {
            // Assume negative mirror character
            value *= -1.00;;
        }
    }
    return value;
}

template <std::size_t DIM>
inline double PhysicalFace<DIM>::basisFunction(Side side, std::size_t i,
                                               std::size_t unknown) {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    if (side == Side::LEFT) {
        return left.basisFunction(i, unknown);
    }
    logger.assert_debug(isInternal_,
                        "cannot find the right element for a boundary face");
    double value = right.basisFunction(i, unknown);
    if (requiresTransformation) {
        value = transform_[unknown]->transform(value, rightVectorTransform);
        if (std::abs(getUnitNormalVector()[2]) > 0.5) {
            // Assume negative mirror character
            value *= -1.00;;
        }
    }
    return value;
}

template <std::size_t DIM>
inline const LinearAlgebra::SmallVector<DIM>&
    PhysicalFace<DIM>::basisFunctionDeriv(std::size_t i, std::size_t unknown) {
    return basisFunctionDeriv_[unknown][i];
}

template <std::size_t DIM>
inline const LinearAlgebra::SmallVector<DIM>&
    PhysicalFace<DIM>::basisFunctionDeriv(Side side, std::size_t i,
                                          std::size_t unknown) {
    if (side == Side::RIGHT) {
        i += nLeftBasisFunctions[unknown];
    }
    return basisFunctionDeriv(i, unknown);
}

template <std::size_t DIM>
inline const LinearAlgebra::SmallVector<DIM>&
    PhysicalFace<DIM>::basisFunctionNormal(std::size_t i) {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    if (hasBasisFunctionNormal[0]) {
        return basisFunctionNormal_[0][i];
    }
    hasBasisFunctionNormal[0] = true;
    for (std::size_t j = 0; j < face_->getNumberOfBasisFunctions(); ++j) {
        basisFunctionNormal_[0][j] = getNormalVector() * basisFunction(j);
        if (j >= nLeftBasisFunctions[0]) {
            basisFunctionNormal_[0][j] *= -1.;
        }
    }
    return basisFunctionNormal_[0][i];
}

template <std::size_t DIM>
inline const LinearAlgebra::SmallVector<DIM>&
    PhysicalFace<DIM>::basisFunctionNormal(std::size_t i, std::size_t unknown) {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    logger.assert_debug(
        unknown < face_->getPtrElementLeft()->getNumberOfUnknowns(),
        "Unknown % does not exist", unknown);
    if (hasBasisFunctionNormal[unknown]) {
        return basisFunctionNormal_[unknown][i];
    }
    hasBasisFunctionNormal[unknown] = true;

    for (std::size_t j = 0; j < face_->getNumberOfBasisFunctions(unknown);
         ++j) {
        basisFunctionNormal_[unknown][j] = getNormalVector() * basisFunction(j);
        if (j >= nLeftBasisFunctions[unknown]) {
            basisFunctionNormal_[unknown][j] *= -1.;
        }
    }
    return basisFunctionNormal_[unknown][i];
}

template <std::size_t DIM>
inline const LinearAlgebra::SmallVector<DIM>&
    PhysicalFace<DIM>::basisFunctionNormal(Side side, std::size_t i) {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    if (side == Side::LEFT) {
        return basisFunctionNormal(i);
    }
    logger.assert_debug(isInternal_,
                        "cannot find the right element for a boundary face");
    return basisFunctionNormal(i + nLeftBasisFunctions[0]);
}

template <std::size_t DIM>
inline const LinearAlgebra::SmallVector<DIM>&
    PhysicalFace<DIM>::basisFunctionNormal(Side side, std::size_t i,
                                           std::size_t unknown) {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    if (side == Side::LEFT) {
        return basisFunctionNormal(i, unknown);
    }
    logger.assert_debug(isInternal_,
                        "cannot find the right element for a boundary face");
    return basisFunctionNormal(i + nLeftBasisFunctions[unknown], unknown);
}

template <std::size_t DIM>
inline const LinearAlgebra::SmallVector<DIM>&
    PhysicalFace<DIM>::basisFunctionUnitNormal(std::size_t i) {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    if (hasBasisFunctionUnitNormal[0]) {
        return basisFunctionUnitNormal_[0][i];
    }
    hasBasisFunctionUnitNormal[0] = true;
    for (std::size_t j = 0; j < face_->getNumberOfBasisFunctions(); ++j) {
        basisFunctionUnitNormal_[0][j] =
            getUnitNormalVector() * basisFunction(j);
        if (j >= nLeftBasisFunctions[0]) {
            basisFunctionUnitNormal_[0][j] *= -1.;
        }
    }
    return basisFunctionUnitNormal_[0][i];
}

template <std::size_t DIM>
inline const LinearAlgebra::SmallVector<DIM>&
    PhysicalFace<DIM>::basisFunctionUnitNormal(std::size_t i,
                                               std::size_t unknown) {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");

    logger.assert_debug(
        unknown < face_->getPtrElementLeft()->getNumberOfUnknowns(),
        "Unknown % does not exist", unknown);
    if (hasBasisFunctionUnitNormal[unknown]) {
        return basisFunctionUnitNormal_[unknown][i];
    }
    hasBasisFunctionUnitNormal[unknown] = true;

    for (std::size_t j = 0; j < face_->getNumberOfBasisFunctions(unknown);
         ++j) {
        basisFunctionUnitNormal_[unknown][j] =
            getUnitNormalVector() * basisFunction(j);
        if (j >= nLeftBasisFunctions[unknown]) {
            basisFunctionUnitNormal_[unknown][j] *= -1.;
        }
    }
    return basisFunctionUnitNormal_[unknown][i];
}

template <std::size_t DIM>
inline const LinearAlgebra::SmallVector<DIM>&
    PhysicalFace<DIM>::basisFunctionUnitNormal(Side side, std::size_t i) {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    if (side == Side::LEFT) {
        return basisFunctionUnitNormal(i);
    }
    logger.assert_debug(isInternal_,
                        "cannot find the right element for a boundary face");
    return basisFunctionUnitNormal(i + nLeftBasisFunctions[0]);
}

template <std::size_t DIM>
inline const LinearAlgebra::SmallVector<DIM>&
    PhysicalFace<DIM>::basisFunctionUnitNormal(Side side, std::size_t i,
                                               std::size_t unknown) {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    if (side == Side::LEFT) {
        return basisFunctionUnitNormal(i, unknown);
    }
    logger.assert_debug(isInternal_,
                        "cannot find the right element for a boundary face");
    return basisFunctionUnitNormal(i + nLeftBasisFunctions[unknown], unknown);
}

template <std::size_t DIM>
inline void PhysicalFace<DIM>::basisFunction(
    std::size_t i, LinearAlgebra::SmallVector<DIM>& result) {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    if (i < nLeftBasisFunctions[0]) {
        left.basisFunction(i, result);
        return;
    }
    logger.assert_debug(isInternal_, "basis function index out of bounds");
    right.basisFunction(i - nLeftBasisFunctions[0], result);
    if (requiresTransformation) {
        result = transform_[0]->transform(result, rightVectorTransform);
        if (std::abs(getUnitNormalVector()[2]) > 0.5) {
            // Assume negative mirror character
            result *= -1.00;;
        }
    }
    return;
}
template <std::size_t DIM>
inline void PhysicalFace<DIM>::basisFunction(
    std::size_t i, LinearAlgebra::SmallVector<DIM>& result,
    std::size_t unknown) {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    if (i < nLeftBasisFunctions[unknown]) {
        left.basisFunction(i, result, unknown);
        return;
    }
    logger.assert_debug(isInternal_, "basis function index out of bounds");
    right.basisFunction(i - nLeftBasisFunctions[unknown], result, unknown);
    if (requiresTransformation) {
        result = transform_[unknown]->transform(result, rightVectorTransform);
        if (std::abs(getUnitNormalVector()[2]) > 0.5) {
            // Assume negative mirror character
            result *= -1.00;;
        }
    }
    return;
}

template <std::size_t DIM>
inline void PhysicalFace<DIM>::basisFunction(
    Side side, std::size_t i, LinearAlgebra::SmallVector<DIM>& result) {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    if (side == Side::LEFT) {
        left.basisFunction(i, result);
    } else {
        logger.assert_debug(
            isInternal_, "cannot find the right element for a boundary face");
        right.basisFunction(i, result);
        if (requiresTransformation) {
            result = transform_[0]->transform(result, rightVectorTransform);
            if (std::abs(getUnitNormalVector()[2]) > 0.5) {
                // Assume negative mirror character
                result *= -1.00;;
            }
        }
    }
}

template <std::size_t DIM>
inline void PhysicalFace<DIM>::basisFunction(
    Side side, std::size_t i, LinearAlgebra::SmallVector<DIM>& result,
    std::size_t unknown) {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    if (side == Side::LEFT) {
        left.basisFunction(i, result, unknown);
    } else {
        logger.assert_debug(
            isInternal_, "cannot find the right element for a boundary face");
        right.basisFunction(i, result, unknown);
        if (requiresTransformation) {
            result =
                transform_[unknown]->transform(result, rightVectorTransform);
            if (std::abs(getUnitNormalVector()[2]) > 0.5) {
                // Assume negative mirror character
                result *= -1.00;;
            }
        }
    }
}

template <std::size_t DIM>
inline const LinearAlgebra::SmallVector<DIM>&
    PhysicalFace<DIM>::basisFunctionCurl(std::size_t i, std::size_t unknown) {
    return basisFunctionCurl_[unknown][i];
}

template <std::size_t DIM>
inline const LinearAlgebra::SmallVector<DIM>&
    PhysicalFace<DIM>::basisFunctionCurl(Side side, std::size_t i,
                                         std::size_t unknown) {
    if (side == Side::RIGHT) {
        i += nLeftBasisFunctions[unknown];
    }
    return basisFunctionCurl(i, unknown);
}

template <std::size_t DIM>
inline const double& PhysicalFace<DIM>::basisFunctionDiv(std::size_t i,
                                                         std::size_t unknown) {
    return basisFunctionDiv_[unknown][i];
}

template <std::size_t DIM>
inline const double& PhysicalFace<DIM>::basisFunctionDiv(Side side,
                                                         std::size_t i,
                                                         std::size_t unknown) {
    if (side == Side::RIGHT) {
        i += nLeftBasisFunctions[unknown];
    }
    return basisFunctionDiv(i, unknown);
}

template <std::size_t DIM>
inline void PhysicalFace<DIM>::basisFunctionNormalCross(
    std::size_t i, LinearAlgebra::SmallVector<DIM>& result)  // Needed for DGMax
{
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    if (hasVectorBasisFunctionNormal[0]) {
        result = vectorBasisFunctionNormal_[0][i];
        return;
    }
    hasVectorBasisFunctionNormal[0] = true;
    for (std::size_t j = 0; j < face_->getNumberOfBasisFunctions(); ++j) {
        basisFunction(j, result);
        vectorBasisFunctionNormal_[0][j] =
            LinearAlgebra::SmallMatrix<DIM, DIM - 1>{
                {getNormalVector(), result}}
                .computeWedgeStuffVector();
        if (j >= nLeftBasisFunctions[0]) {
            vectorBasisFunctionNormal_[0][j] *= -1.;
        }
    }
    result = vectorBasisFunctionNormal_[0][i];
    return;
}

template <std::size_t DIM>
inline void PhysicalFace<DIM>::basisFunctionNormalCross(
    std::size_t i, LinearAlgebra::SmallVector<DIM>& result,
    std::size_t unknown)  // Needed for DGMax
{
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");

    logger.assert_debug(
        unknown < face_->getPtrElementLeft()->getNumberOfUnknowns(),
        "Unknown % does not exist", unknown);
    if (hasVectorBasisFunctionNormal[unknown]) {
        result = vectorBasisFunctionNormal_[unknown][i];
        return;
    }
    hasVectorBasisFunctionNormal[unknown] = true;

    for (std::size_t j = 0; j < face_->getNumberOfBasisFunctions(unknown);
         ++j) {
        basisFunction(j, result, unknown);
        vectorBasisFunctionNormal_[unknown][j] =
            LinearAlgebra::SmallMatrix<DIM, DIM - 1>{
                {getNormalVector(), result}}
                .computeWedgeStuffVector();
        if (j >= nLeftBasisFunctions[unknown]) {
            vectorBasisFunctionNormal_[unknown][j] *= -1.;
        }
    }
    result = vectorBasisFunctionNormal_[unknown][i];
    return;
}

template <std::size_t DIM>
inline void PhysicalFace<DIM>::basisFunctionNormalCross(
    Side side, std::size_t i, LinearAlgebra::SmallVector<DIM>& result) {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    if (side == Side::LEFT) {
        return basisFunctionNormal(i, result);
    }
    logger.assert_debug(isInternal_,
                        "cannot find the right element for a boundary face");
    return basisFunctionNormal(i + nLeftBasisFunctions[0], result);
}

template <std::size_t DIM>
inline void PhysicalFace<DIM>::basisFunctionNormalCross(
    Side side, std::size_t i, LinearAlgebra::SmallVector<DIM>& result,
    std::size_t unknown) {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    if (side == Side::LEFT) {
        return basisFunctionNormal(i, result, unknown);
    }
    logger.assert_debug(isInternal_,
                        "cannot find the right element for a boundary face");
    return basisFunctionNormal(i + nLeftBasisFunctions[unknown], result,
                               unknown);
}

template <std::size_t DIM>
inline void PhysicalFace<DIM>::basisFunctionUnitNormalCross(
    std::size_t i, LinearAlgebra::SmallVector<DIM>& result) {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    if (hasVectorBasisFunctionUnitNormal[0]) {
        result = vectorBasisFunctionUnitNormal_[0][i];
        return;
    }
    hasVectorBasisFunctionUnitNormal[0] = true;
    for (std::size_t j = 0; j < face_->getNumberOfBasisFunctions(); ++j) {
        basisFunction(j, result);
        switch (DIM) {
            case 2: {
                LinearAlgebra::SmallVector<DIM> tangentialUnitVector;
                tangentialUnitVector =
                    LinearAlgebra::SmallMatrix<DIM, DIM - 1>{
                        {getUnitNormalVector()}}
                        .computeWedgeStuffVector();
                vectorBasisFunctionUnitNormal_[0][j][0] =
                    tangentialUnitVector * result;
                vectorBasisFunctionUnitNormal_[0][j][1] = 0.0;
            } break;
            default:
                vectorBasisFunctionUnitNormal_[0][j] =
                    LinearAlgebra::SmallMatrix<DIM, DIM - 1>{
                        {getUnitNormalVector(), result}}
                        .computeWedgeStuffVector();
        }
        if (j >= nLeftBasisFunctions[0]) {
            vectorBasisFunctionUnitNormal_[0][j] *= -1.;  // current
        }
    }
    result = vectorBasisFunctionUnitNormal_[0][i];
    return;
}

template <std::size_t DIM>
inline void PhysicalFace<DIM>::basisFunctionUnitNormalCross(
    std::size_t i, LinearAlgebra::SmallVector<DIM>& result,
    std::size_t unknown) {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    logger.assert_debug(
        unknown < face_->getPtrElementLeft()->getNumberOfUnknowns(),
        "Unknown % does not exist", unknown);
    if (hasVectorBasisFunctionUnitNormal[unknown]) {
        result = vectorBasisFunctionUnitNormal_[unknown][i];
        return;
    }
    hasVectorBasisFunctionUnitNormal[unknown] = true;
    for (std::size_t j = 0; j < face_->getNumberOfBasisFunctions(unknown);
         ++j) {
        basisFunction(j, result, unknown);
        switch (DIM) {
            case 2: {
                LinearAlgebra::SmallVector<DIM> tangentialUnitVector;
                tangentialUnitVector =
                    LinearAlgebra::SmallMatrix<DIM, DIM - 1>{
                        {getUnitNormalVector()}}
                        .computeWedgeStuffVector();
                vectorBasisFunctionUnitNormal_[unknown][j][0] =
                    tangentialUnitVector * result;
                vectorBasisFunctionUnitNormal_[unknown][j][1] = 0.0;
                break;
            }
            default:
                vectorBasisFunctionUnitNormal_[unknown][j] =
                    LinearAlgebra::SmallMatrix<DIM, DIM - 1>{
                        {getUnitNormalVector(), result}}
                        .computeWedgeStuffVector();
        }
        if (j >= nLeftBasisFunctions[unknown]) {
            vectorBasisFunctionUnitNormal_[unknown][j] *= -1.;  // current
        }
    }
    result = vectorBasisFunctionUnitNormal_[unknown][i];
    return;
}

template <std::size_t DIM>
inline void PhysicalFace<DIM>::basisFunctionUnitNormalCross(
    Side side, std::size_t i, LinearAlgebra::SmallVector<DIM>& result) {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    if (side == Side::LEFT) {
        return basisFunctionUnitNormal(i, result);
    }
    logger.assert_debug(isInternal_,
                        "cannot find the right element for a boundary face");
    return basisFunctionUnitNormal(i + nLeftBasisFunctions[0], result);
}
template <std::size_t DIM>
inline void PhysicalFace<DIM>::basisFunctionUnitNormalCross(
    Side side, std::size_t i, LinearAlgebra::SmallVector<DIM>& result,
    std::size_t unknown) {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    if (side == Side::LEFT) {
        return basisFunctionUnitNormal(i, result, unknown);
    }
    logger.assert_debug(isInternal_,
                        "cannot find the right element for a boundary face");
    return basisFunctionUnitNormal(i + nLeftBasisFunctions[unknown], result,
                                   unknown);
}

template <std::size_t DIM>
inline const LinearAlgebra::MiddleSizeVector& PhysicalFace<DIM>::getSolution(
    Side side) {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    if (side == Side::LEFT) {
        return left.getSolution();
    }
    logger.assert_debug(isInternal_,
                        "cannot find the right element for a boundary face");
    return right.getSolution();
}

template <std::size_t DIM>
inline void PhysicalFace<DIM>::getSolution(
    Side side, std::vector<LinearAlgebra::SmallVector<DIM>>& result) {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    if (side == Side::LEFT) {
        return left.getSolution(result);
    }
    logger.assert_debug(isInternal_,
                        "cannot find the right element for a boundary face");
    return right.getSolution(result);
}

template <std::size_t DIM>
inline const std::vector<LinearAlgebra::SmallVector<DIM>>&
    PhysicalFace<DIM>::getSolutionDeriv(Side side) {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    if (side == Side::LEFT) {
        return left.getSolutionDeriv();
    }
    logger.assert_debug(isInternal_,
                        "cannot find the right element for a boundary face");
    return right.getSolutionDeriv();
}

template <std::size_t DIM>
inline const std::vector<LinearAlgebra::SmallVector<DIM>>&
    PhysicalFace<DIM>::getSolutionCurl(Side side) {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    if (side == Side::LEFT) {
        return left.getSolutionCurl();
    }
    logger.assert_debug(isInternal_,
                        "cannot find the right element for a boundary face");
    return right.getSolutionCurl();
}

template <std::size_t DIM>
inline std::vector<LinearAlgebra::SmallVector<DIM>>
    PhysicalFace<DIM>::getSolutionNormal(Side side) {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    if (side == Side::LEFT) {
        return getNormalVector() * left.getSolution();
    }
    logger.assert_debug(isInternal_,
                        "cannot find the right element for a boundary face");
    return -getNormalVector() * right.getSolution();
}

template <std::size_t DIM>
inline std::vector<LinearAlgebra::SmallVector<DIM>>
    PhysicalFace<DIM>::getSolutionUnitNormal(Side side) {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    if (side == Side::LEFT) {
        return getUnitNormalVector() * left.getSolution();
    }
    logger.assert_debug(isInternal_,
                        "cannot find the right element for a boundary face");
    return -getUnitNormalVector() * right.getSolution();
}

template <std::size_t DIM>
inline void PhysicalFace<DIM>::getSolutionNormal(
    Side side, LinearAlgebra::SmallVector<DIM>& result) {
    logger(ERROR, "Not supported by Element yet");
}

template <std::size_t DIM>
inline void PhysicalFace<DIM>::getSolutionUnitNormal(
    Side side, LinearAlgebra::SmallVector<DIM>& result) {
    logger(ERROR, "Not supported by Element yet");
}

template <std::size_t DIM>
inline const Geometry::PointReference<DIM - 1>&
    PhysicalFace<DIM>::getPointReference() {
    logger.assert_debug(hasPointReference,
                        "Need a location to evaluate the data");
    return pointReference_;
}

template <std::size_t DIM>
inline const Geometry::PointPhysical<DIM>&
    PhysicalFace<DIM>::getPointPhysical() {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    return left.getPointPhysical();
}

template <std::size_t DIM>
inline const LinearAlgebra::SmallVector<DIM>&
    PhysicalFace<DIM>::getNormalVector() {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    if (hasNormal) {
        return normal;
    }
    hasNormal = true;
    normal = face_->getNormalVector(pointReference_);
    return normal;
}

template <std::size_t DIM>
inline const LinearAlgebra::SmallVector<DIM>&
    PhysicalFace<DIM>::getUnitNormalVector() {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    if (hasUnitNormal) {
        return unitNormal;
    }
    hasUnitNormal = true;
    unitNormal = getNormalVector() / getRelativeSurfaceArea();
    return unitNormal;
}

template <std::size_t DIM>
inline double PhysicalFace<DIM>::getRelativeSurfaceArea() {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    if (hasNormalNorm) {
        return normalNorm;
    }
    hasNormalNorm = true;
    normalNorm = getNormalVector().l2Norm();
    return normalNorm;
}

template <std::size_t DIM>
inline FaceMatrix& PhysicalFace<DIM>::getResultMatrix() {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    logger.assert_debug(
        hasFaceMatrix,
        "Matrix has already been requested for this face/point combination");
    hasFaceMatrix = false;
    return resultMatrix;
}

template <std::size_t DIM>
inline LinearAlgebra::MiddleSizeMatrix& PhysicalFace<DIM>::getResultMatrix(
    Side iSide, Side jSide) {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    if (iSide == Side::LEFT) {
        if (jSide == Side::LEFT) {
            return left.getResultMatrix();
        }
        logger.assert_debug(
            isInternal_, "cannot find the right element for a boundary face");
        logger.assert_debug(hasLeftRightMatrix,
                            "Matrix has already been requested for this "
                            "face/point combination");
        hasLeftRightMatrix = false;
        return leftRightMatrix;

    } else {
        logger.assert_debug(
            isInternal_, "cannot find the right element for a boundary face");
        if (jSide == Side::LEFT) {
            logger.assert_debug(hasRightLeftMatrix,
                                "Matrix has already been requested for this "
                                "face/point combination");
            hasRightLeftMatrix = false;
            return rightLeftMatrix;
        }
        return right.getResultMatrix();
    }
}

template <std::size_t DIM>
inline LinearAlgebra::MiddleSizeVector& PhysicalFace<DIM>::getResultVector() {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    logger.assert_debug(
        hasFaceVector,
        "Vector has already been requested for this face/point combination");
    hasFaceVector = false;
    return resultVector;
}

template <std::size_t DIM>
inline LinearAlgebra::MiddleSizeVector& PhysicalFace<DIM>::getResultVector(
    Side side) {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");
    if (side == Side::LEFT) {
        return left.getResultVector();
    }
    logger.assert_debug(isInternal_,
                        "cannot find the right element for a boundary face");
    return right.getResultVector();
}

template <std::size_t DIM>
inline bool PhysicalFace<DIM>::isInternal() {
    return isInternal_;
}

template <std::size_t DIM>
inline const Face* PhysicalFace<DIM>::getFace() {
    logger.assert_debug(hasFace, "Need a location to evaluate the data");
    return face_;
}

template <std::size_t DIM>
inline const CoordinateTransformation<DIM>* PhysicalFace<DIM>::getTransform() {
    for (std::size_t i = 1; i < transform_.size(); ++i)
        logger.assert_debug(
            transform_[0] == transform_[i],
            "Different unknowns have different coordinate transformation, call "
            "getTransformation for a given unknown");
    return transform_[0].get();
}

template <std::size_t DIM>
inline const CoordinateTransformation<DIM>* PhysicalFace<DIM>::getTransform(
    std::size_t unknown) {
    return transform_[unknown].get();
}

template <std::size_t DIM>
inline void PhysicalFace<DIM>::setPointReference(
    const Geometry::PointReference<DIM - 1>& point) {
    pointReference_ = point;
    hasPointReference = true;
    if (hasFace) {
        left.setPointReference(face_->mapRefFaceToRefElemL(point));
        if (isInternal_) {
            right.setPointReference(face_->mapRefFaceToRefElemR(point));
        }
    }
    // even if they are already computed, the information is now out of date
    std::size_t unknowns = face_->getPtrElementLeft()->getNumberOfUnknowns();
    hasBasisFunctionNormal.assign(unknowns, false);
    hasBasisFunctionUnitNormal.assign(unknowns, false);
    hasVectorBasisFunctionNormal.assign(unknowns, false);
    hasVectorBasisFunctionUnitNormal.assign(unknowns, false);

    resetLazyCaches(unknowns);

    hasSolutionNormal = false;
    hasSolutionUnitNormal = false;
    hasVectorSolutionNormal = false;
    hasVectorSolutionUnitNormal = false;
    hasNormal = false;
    hasUnitNormal = false;
    hasNormalNorm = false;
    if (!hasFaceMatrix) {
        resultMatrix *= 0;
    }
    if (!hasFaceVector) {
        resultVector *= 0;
    }
    if (!hasLeftRightMatrix) {
        leftRightMatrix *= 0;
    }
    if (!hasRightLeftMatrix) {
        rightLeftMatrix *= 0;
    }
    hasLeftRightMatrix = true;
    hasRightLeftMatrix = true;
    hasFaceMatrix = true;
    hasFaceVector = true;
}

template <std::size_t DIM>
inline void PhysicalFace<DIM>::setFace(const Face* face) {
    logger.assert_debug(isInternal_ == face->isInternal(),
                        "This face is not supported by this physical face");
    std::size_t numberOfUnknowns =
        face->getPtrElementLeft()->getNumberOfUnknowns();
    if (!hasFace) {
        std::size_t leftCoefficients =
            face->getPtrElementLeft()->getTotalNumberOfBasisFunctions();
        nLeftBasisFunctions.resize(numberOfUnknowns);
        basisFunctionNormal_.resize(numberOfUnknowns);
        vectorBasisFunctionNormal_.resize(numberOfUnknowns);
        basisFunctionUnitNormal_.resize(numberOfUnknowns);
        vectorBasisFunctionUnitNormal_.resize(numberOfUnknowns);
        std::size_t rightCoefficients = 0;
        for (std::size_t i = 0; i < numberOfUnknowns; ++i) {
            std::size_t numBasis =
                face->getPtrElementLeft()->getNumberOfBasisFunctions(i);
            basisFunctionNormal_[i].resize(numBasis);
            vectorBasisFunctionNormal_[i].resize(numBasis);
            basisFunctionUnitNormal_[i].resize(numBasis);
            vectorBasisFunctionUnitNormal_[i].resize(numBasis);
            nLeftBasisFunctions[i] = numBasis;
        }
        if (isInternal_) {
            rightCoefficients =
                face->getPtrElementRight()->getTotalNumberOfBasisFunctions();
            // We assume that in the left and right element the number of
            // unknowns is equal.
            for (std::size_t i = 0;
                 i < face->getPtrElementRight()->getNumberOfUnknowns(); ++i) {
                std::size_t numBasis =
                    face->getPtrElementLeft()->getNumberOfBasisFunctions(i) +
                    face->getPtrElementRight()->getNumberOfBasisFunctions(i);
                basisFunctionNormal_[i].resize(numBasis);
                vectorBasisFunctionNormal_[i].resize(numBasis);
                basisFunctionUnitNormal_[i].resize(numBasis);
                vectorBasisFunctionUnitNormal_[i].resize(numBasis);
            }
        }
        resultMatrix.resize(leftCoefficients, rightCoefficients);
        leftRightMatrix.resize(leftCoefficients, rightCoefficients);
        rightLeftMatrix.resize(rightCoefficients, leftCoefficients);
        resultVector.resize(leftCoefficients + rightCoefficients);
    }
    face_ = face;
    hasFace = true;
    left.setElement(face->getPtrElementLeft());
    if (isInternal_) {
        right.setElement(face->getPtrElementRight());
    }
    if (hasPointReference) {
        left.setPointReference(face_->mapRefFaceToRefElemL(pointReference_));
        if (isInternal_) {
            right.setPointReference(
                face_->mapRefFaceToRefElemR(pointReference_));
        }
    }
    // even if they are already computed, the information is now out of date
    hasBasisFunctionNormal.assign(numberOfUnknowns, false);
    hasBasisFunctionUnitNormal.assign(numberOfUnknowns, false);
    hasVectorBasisFunctionNormal.assign(numberOfUnknowns, false);
    hasVectorBasisFunctionUnitNormal.assign(numberOfUnknowns, false);

    resetLazyCaches(numberOfUnknowns);

    hasSolutionNormal = false;
    hasSolutionUnitNormal = false;
    hasVectorSolutionNormal = false;
    hasVectorSolutionUnitNormal = false;
    hasNormal = false;
    hasUnitNormal = false;
    hasNormalNorm = false;
    if (!hasFaceMatrix) {
        resultMatrix *= 0;
    }
    if (!hasFaceVector) {
        resultVector *= 0;
    }
    if (!hasLeftRightMatrix) {
        leftRightMatrix *= 0;
    }
    if (!hasRightLeftMatrix) {
        rightLeftMatrix *= 0;
    }
    hasLeftRightMatrix = true;
    hasRightLeftMatrix = true;
    hasFaceMatrix = true;
    hasFaceVector = true;
    updateLeftRightTransform();
}

template <std::size_t DIM>
inline void PhysicalFace<DIM>::setTransform(
    std::shared_ptr<Base::CoordinateTransformation<DIM>> transform,
    std::size_t unknown) {
    if (transform_.size() <= unknown) {
        // We should not need to resize when we know the exact number of
        // unknowns.
        logger.assert_debug(!hasFace, "Resizing with a face");
        transform_.resize(unknown + 1);
    }
    transform_[unknown] = transform;
    left.setTransformation(transform, unknown);
    if (isInternal_) {
        right.setTransformation(transform, unknown);
    }

    hasBasisFunctionNormal.assign(hasBasisFunctionNormal.size(), false);
    hasBasisFunctionUnitNormal.assign(hasBasisFunctionUnitNormal.size(), false);
    hasVectorBasisFunctionNormal.assign(hasVectorBasisFunctionNormal.size(),
                                        false);
    hasVectorBasisFunctionUnitNormal.assign(
        hasVectorBasisFunctionUnitNormal.size(), false);

    resetLazyCaches(basisFunctionDeriv_.size());

    hasSolutionNormal = false;
    hasSolutionUnitNormal = false;
    hasVectorSolutionNormal = false;
    hasVectorSolutionUnitNormal = false;

    if (!hasFaceMatrix) {
        resultMatrix *= 0;
    }
    if (!hasFaceVector) {
        resultVector *= 0;
    }
    if (!hasLeftRightMatrix) {
        leftRightMatrix *= 0;
    }
    if (!hasRightLeftMatrix) {
        rightLeftMatrix *= 0;
    }
    hasLeftRightMatrix = true;
    hasRightLeftMatrix = true;
    hasFaceMatrix = true;
    hasFaceVector = true;
}

template <std::size_t DIM>
inline void PhysicalFace<DIM>::setQuadratureRule(
    QuadratureRules::GaussQuadratureRule* rule) {
    mapToLeftElement = face_->refFaceToRefElemMapL();
    left.setQuadratureRule(rule, mapToLeftElement.get());
    if (isInternal()) {
        mapToRightElement = face_->refFaceToRefElemMapR();
        right.setQuadratureRule(rule, mapToRightElement.get());
    }
    quadratureRule_ = rule;
    setQuadraturePointIndex(0);
}

template <std::size_t DIM>
inline void PhysicalFace<DIM>::setQuadraturePointIndex(std::size_t index) {
    setPointReference(quadratureRule_->getPoint(index));
    // setPointReference tells the element that it is using a pointReference
    // so we have to set the quadrature rule back. If this turns out to be slow
    // resort to code duplication to prevent double work
    left.setQuadratureRule(quadratureRule_, mapToLeftElement.get());
    left.setQuadraturePointIndex(index);
    if (isInternal()) {
        right.setQuadratureRule(quadratureRule_, mapToRightElement.get());
        right.setQuadraturePointIndex(index);
    }
}

template <std::size_t DIM>
void PhysicalFace<DIM>::resetLazyCaches(std::size_t currentUnknowns) {
    basisFunctionDeriv_.reset(currentUnknowns);
    basisFunctionCurl_.reset(currentUnknowns);
    basisFunctionDiv_.reset(currentUnknowns);
}

template <std::size_t DIM>
void PhysicalFace<DIM>::computeBasisFunctionDeriv(
    std::vector<LinearAlgebra::SmallVector<DIM>>& values, std::size_t unknown) {
    logger.assert_debug(hasPointReference && hasFace,
                        "Need a location to evaluate the data");

    std::size_t numBasisFunctions = face_->getNumberOfBasisFunctions(unknown);
    std::size_t nLeft = nLeftBasisFunctions[unknown];
    values.resize(numBasisFunctions);
    for (std::size_t i = 0; i < numBasisFunctions; ++i) {
        if (i < nLeft) {
            values[i] = left.basisFunctionDeriv(i, unknown);
        } else {
            values[i] = right.basisFunctionDeriv(i - nLeft, unknown);
            if (requiresTransformation) {
                values[i] = transform_[unknown]->transformDeriv(
                    values[i], rightVectorTransform);
                if (std::abs(getUnitNormalVector()[2]) > 0.5) {
                    // Assume negative mirror character
                    values[i] *= -1.00;;
                }
            }
        }
    }
}

template <std::size_t DIM>
void PhysicalFace<DIM>::computeBasisFunctionCurl(
    std::vector<LinearAlgebra::SmallVector<DIM>>& values, std::size_t unknown) {
    std::size_t numBasisFunctions = face_->getNumberOfBasisFunctions(unknown);
    std::size_t nLeft = nLeftBasisFunctions[unknown];
    values.resize(numBasisFunctions);
    for (std::size_t i = 0; i < numBasisFunctions; ++i) {
        if (i < nLeft) {
            values[i] = left.basisFunctionCurl(i, unknown);
        } else {
            values[i] = right.basisFunctionCurl(i - nLeft, unknown);
            if (requiresTransformation) {
                values[i] = transform_[unknown]->transformCurl(
                    values[i], rightVectorTransform);
                if (std::abs(getUnitNormalVector()[2]) > 0.5) {
                    // Assume negative mirror character
                    values[i] *= -1.00;;
                }
            }
        }
    }
}

template <std::size_t DIM>
void PhysicalFace<DIM>::computeBasisFunctionDiv(std::vector<double>& values,
                                                std::size_t unknown) {
    std::size_t numBasisFunctions = face_->getNumberOfBasisFunctions(unknown);
    std::size_t nLeft = nLeftBasisFunctions[unknown];
    values.resize(numBasisFunctions);

    for (std::size_t i = 0; i < numBasisFunctions; ++i) {
        if (i < nLeft) {
            values[i] = left.basisFunctionDiv(i, unknown);
        } else {
            values[i] = right.basisFunctionDiv(i - nLeft, unknown);
            if (requiresTransformation) {
                values[i] = transform_[unknown]->transformDiv(
                    values[i], rightVectorTransform);
                if (std::abs(getUnitNormalVector()[2]) > 0.5) {
                    // Assume negative mirror character
                    values[i] *= -1.00;;
                }
            }
        }
    }
}

template <std::size_t DIM>
void PhysicalFace<DIM>::updateLeftRightTransform() {
    if (!face_->isInternal() || DIM == 1) {
        requiresTransformation = false;
        return;
    }
    // NOTE: face->getFaceType() is not very reliable, otherwise a very good
    // optimization would be to check if it is a periodic (subdomain) boundary
    // and only compute the transform if so.

    // We need to compute how the affine transformation from the right face to
    // the left one. We proceed as follows:
    //
    // 1. Compute the DIM-1 physical vectors in real space that correspond to
    //    the unit vectors e_i on the reference face.
    // 2. Add a physical normal vector for one the outer for the other the
    //    inner so that they should transform onto each other. Giving DIM
    //    independent physical vectors for both faces.
    // 3. Compute the rotation matrix that transforms the right set of vectors
    //    onto that of the left Set.
    //
    // Optionally one could compute the offset of the affine transformation, but
    // that is not needed.

    std::array<LinearAlgebra::SmallVector<DIM>, DIM> leftPoints, rightPoints;

    bool identicalPoints = true;
    for (std::size_t i = 0; i < DIM; ++i) {
        // We map the zero vector and unit vectors e_i. This abuses the
        // knowledge that these vectors are inside the reference shape for each
        // reference shape.
        Geometry::PointReference<DIM - 1> facePoint;
        if (i > 0) {
            facePoint[i - 1] = 1.0;
        }

        // Compute the physical coordinates for both the left and right face
        leftPoints[i] =
            left.getElement()
                ->referenceToPhysical(
                    face_->refFaceToRefElemMapL()->transform(facePoint))
                .getCoordinates();
        rightPoints[i] =
            right.getElement()
                ->referenceToPhysical(
                    face_->refFaceToRefElemMapR()->transform(facePoint))
                .getCoordinates();

        // Small optimization, most faces are non-periodic. If all left and
        // right points match we can directly conclude no transform is needed.
        identicalPoints &=
            (leftPoints[i] - rightPoints[i]).l2NormSquared() < 1e-24;
    }
    if (identicalPoints) {
        requiresTransformation = false;
        return;
    }

    // Compute the direction vectors relative to the left and right face.
    LinearAlgebra::SmallMatrix<DIM, DIM> leftMatrix =
        computeDirectionVectors(leftPoints);

    LinearAlgebra::SmallMatrix<DIM, DIM> rightMatrix =
        computeDirectionVectors(rightPoints);

    const auto& normal = getUnitNormalVector();
    if (std::abs(normal[2]) > 0.5) {
        for (std::size_t i = 0; i < DIM; ++i) {
            // Needed as the computeDirectionVectors assumes a transformation
            // with det(R) = 1, while our current case has -1.
            rightMatrix(DIM - 1, i) *= -1;
        }
    }

    // Computes the transformation matrix to map the directions on the right
    // face to those on the left face.
    leftMatrix.solve(rightMatrix);

    // Check if transformation is needed by comparing the matrix to the identity
    // matrix.
    requiresTransformation = false;
    for (std::size_t i = 0; i < DIM; ++i) {
        for (std::size_t j = 0; j < DIM; ++j) {
            double val = i == j ? 1.0 : 0.0;
            requiresTransformation |= std::abs(rightMatrix(i, j) - val) > 1e-8;
        }
    }
    if (requiresTransformation) {
        if (std::abs(rightMatrix(1, 1) - 1) > 5e-2) {
            logger(INFO, "WRONG!");
        }
        // Check that it is a proper rotation matrix, no scaling and no improper
        // rotations.
        logger.assert_debug(
            std::abs(rightMatrix.determinant() - 1.0) < 1e-8 ||
                std::abs(rightMatrix.determinant() + 1.0) < 1e-8,
            "Not a good transformation matrix.");
        rightVectorTransform.setJacobian(rightMatrix);
    }
}

template <std::size_t DIM>
LinearAlgebra::SmallMatrix<DIM, DIM> PhysicalFace<DIM>::computeDirectionVectors(
    const std::array<LinearAlgebra::SmallVector<DIM>, DIM>& points) {
    LinearAlgebra::SmallMatrix<DIM, DIM - 1> temp;
    LinearAlgebra::SmallMatrix<DIM, DIM> result;
    // Use the first vector as origin and compute direction vectors to the other
    // DIM-1 vectors.
    for (std::size_t i = 0; i < DIM - 1; ++i) {
        for (std::size_t j = 0; j < DIM; ++j) {
            double value = points[i + 1][j] - points[0][j];
            temp(j, i) = value;
            result(i, j) = value;
        }
    }
    // ALl vectors are in the plane of the face, compute a normal vector to the
    // face as DIM-th independent direction.
    LinearAlgebra::SmallVector<DIM> normal = temp.computeWedgeStuffVector();
    for (std::size_t j = 0; j < DIM; ++j) {
        result(DIM - 1, j) = normal[j];
    }
    return result;
}

}  // namespace Base
}  // namespace hpgem
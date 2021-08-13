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

#include "H1ConformingTransformation.h"

namespace hpgem {

namespace Base {
template <std::size_t DIM>
PhysicalElement<DIM>::PhysicalElement()
    : hasPointReference(false),
      hasElement(false),
      hasQuadratureRule(false)  // other data will get initialized when we have
                                // more info
{
    std::shared_ptr<Base::CoordinateTransformation<DIM> > transform(
        new H1ConformingTransformation<DIM>{});
    transform_.push_back(transform);
    hasElementMatrix = false;
    hasElementVector = false;
}

template <std::size_t DIM>
inline double PhysicalElement<DIM>::basisFunction(std::size_t i) {
    logger.assert_debug(hasPointReference && hasElement,
                        "Need a location to evaluate the data");
    if (hasFunctionValue[0]) {
        return basisFunctionValue[0][i];
    }
    hasFunctionValue[0] = true;
    for (std::size_t j = 0; j < theElement_->getNumberOfBasisFunctions(); ++j) {
        if (hasQuadratureRule) {
            if (doesMapQuadraturePointFromFace) {
                basisFunctionValue[0][j] = transform_[0]->transform(
                    theElement_->basisFunction(j, quadratureRule_,
                                               quadraturePointIndex_,
                                               faceToElementMap_),
                    *this);
            } else {
                basisFunctionValue[0][j] = transform_[0]->transform(
                    theElement_->basisFunction(j, quadratureRule_,
                                               quadraturePointIndex_),
                    *this);
            }
        } else {
            basisFunctionValue[0][j] = transform_[0]->transform(
                theElement_->basisFunction(j, pointReference_), *this);
        }
    }
    return basisFunctionValue[0][i];
}

template <std::size_t DIM>
inline double PhysicalElement<DIM>::basisFunction(std::size_t i,
                                                  std::size_t unknown) {
    logger.assert_debug(hasPointReference && hasElement,
                        "Need a location to evaluate the data");
    logger.assert_debug(unknown < theElement_->getNumberOfUnknowns(),
                        "Unknown % does not exist", unknown);
    if (hasFunctionValue[unknown]) {
        return basisFunctionValue[unknown][i];
    }
    hasFunctionValue[unknown] = true;
    for (std::size_t j = 0; j < theElement_->getNumberOfBasisFunctions(unknown);
         ++j) {
        if (hasQuadratureRule) {
            if (doesMapQuadraturePointFromFace) {
                basisFunctionValue[unknown][j] = transform_[unknown]->transform(
                    theElement_->basisFunction(j, quadratureRule_,
                                               quadraturePointIndex_,
                                               faceToElementMap_, unknown),
                    *this);
            } else {
                basisFunctionValue[unknown][j] = transform_[unknown]->transform(
                    theElement_->basisFunction(j, quadratureRule_,
                                               quadraturePointIndex_, unknown),
                    *this);
            }
        } else {
            basisFunctionValue[unknown][j] = transform_[unknown]->transform(
                theElement_->basisFunction(j, pointReference_, unknown), *this);
        }
    }
    return basisFunctionValue[unknown][i];
}

template <std::size_t DIM>
inline const LinearAlgebra::SmallVector<DIM>&
    PhysicalElement<DIM>::basisFunctionDeriv(std::size_t i) {
    logger.assert_debug(hasPointReference && hasElement,
                        "Need a location to evaluate the data");
    if (hasFunctionDeriv[0]) {
        return basisFunctionDeriv_[0][i];
    }
    hasFunctionDeriv[0] = true;
    for (std::size_t j = 0; j < theElement_->getNumberOfBasisFunctions(); ++j) {
        if (hasQuadratureRule) {
            if (doesMapQuadraturePointFromFace) {
                basisFunctionDeriv_[0][j] = transform_[0]->transformDeriv(
                    theElement_->basisFunctionDeriv<DIM>(j, quadratureRule_,
                                                         quadraturePointIndex_,
                                                         faceToElementMap_),
                    *this);
            } else {
                basisFunctionDeriv_[0][j] = transform_[0]->transformDeriv(
                    theElement_->basisFunctionDeriv<DIM>(j, quadratureRule_,
                                                         quadraturePointIndex_),
                    *this);
            }
        } else {
            basisFunctionDeriv_[0][j] = transform_[0]->transformDeriv(
                theElement_->basisFunctionDeriv(j, pointReference_), *this);
        }
    }
    return basisFunctionDeriv_[0][i];
}

template <std::size_t DIM>
inline const LinearAlgebra::SmallVector<DIM>&
    PhysicalElement<DIM>::basisFunctionDeriv(std::size_t i,
                                             std::size_t unknown) {
    logger.assert_debug(hasPointReference && hasElement,
                        "Need a location to evaluate the data");
    logger.assert_debug(unknown < theElement_->getNumberOfUnknowns(),
                        "Unknown % does not exist", unknown);
    if (hasFunctionDeriv[unknown]) {
        return basisFunctionDeriv_[unknown][i];
    }
    hasFunctionDeriv[unknown] = true;
    for (std::size_t j = 0; j < theElement_->getNumberOfBasisFunctions(unknown);
         ++j) {
        if (hasQuadratureRule) {
            if (doesMapQuadraturePointFromFace) {
                basisFunctionDeriv_[unknown][j] =
                    transform_[unknown]->transformDeriv(
                        theElement_->basisFunctionDeriv<DIM>(
                            j, quadratureRule_, quadraturePointIndex_,
                            faceToElementMap_, unknown),
                        *this);
            } else {
                basisFunctionDeriv_[unknown][j] =
                    transform_[unknown]->transformDeriv(
                        theElement_->basisFunctionDeriv<DIM>(
                            j, quadratureRule_, quadraturePointIndex_, unknown),
                        *this);
            }
        } else {
            basisFunctionDeriv_[unknown][j] =
                transform_[unknown]->transformDeriv(
                    theElement_->basisFunctionDeriv(j, pointReference_,
                                                    unknown),
                    *this);
        }
    }
    return basisFunctionDeriv_[unknown][i];
}

template <std::size_t DIM>
inline void PhysicalElement<DIM>::basisFunction(
    std::size_t i, LinearAlgebra::SmallVector<DIM>& result) {
    logger.assert_debug(hasPointReference && hasElement,
                        "Need a location to evaluate the data");
    if (hasVectorFunctionValue[0]) {
        result = vectorBasisFunctionValue[0][i];
    } else {
        hasVectorFunctionValue[0] = true;
        for (std::size_t j = 0; j < theElement_->getNumberOfBasisFunctions();
             ++j) {
            if (hasQuadratureRule) {
                if (doesMapQuadraturePointFromFace) {
                    theElement_->basisFunction(
                        j, quadratureRule_, quadraturePointIndex_,
                        faceToElementMap_, vectorBasisFunctionValue[0][j]);
                    vectorBasisFunctionValue[0][j] = transform_[0]->transform(
                        vectorBasisFunctionValue[0][j], *this);
                } else {
                    theElement_->basisFunction(j, quadratureRule_,
                                               quadraturePointIndex_,
                                               vectorBasisFunctionValue[0][j]);
                    vectorBasisFunctionValue[0][j] = transform_[0]->transform(
                        vectorBasisFunctionValue[0][j], *this);
                }
            } else {
                theElement_->basisFunction(j, pointReference_,
                                           vectorBasisFunctionValue[0][j]);
                vectorBasisFunctionValue[0][j] = transform_[0]->transform(
                    vectorBasisFunctionValue[0][j], *this);
            }
        }
        result = vectorBasisFunctionValue[0][i];
    }
}

template <std::size_t DIM>
inline void PhysicalElement<DIM>::basisFunction(
    std::size_t i, LinearAlgebra::SmallVector<DIM>& result,
    std::size_t unknown) {
    logger.assert_debug(hasPointReference && hasElement,
                        "Need a location to evaluate the data");
    logger.assert_debug(unknown < theElement_->getNumberOfUnknowns(),
                        "Unknown % does not exist", unknown);
    if (hasVectorFunctionValue[unknown]) {
        result = vectorBasisFunctionValue[unknown][i];
    } else {
        hasVectorFunctionValue[unknown] = true;
        for (std::size_t j = 0;
             j < theElement_->getNumberOfBasisFunctions(unknown); ++j) {
            if (hasQuadratureRule) {
                if (doesMapQuadraturePointFromFace) {
                    theElement_->basisFunction(
                        j, quadratureRule_, quadraturePointIndex_,
                        faceToElementMap_, vectorBasisFunctionValue[unknown][j],
                        unknown);
                    vectorBasisFunctionValue[unknown][j] =
                        transform_[unknown]->transform(
                            vectorBasisFunctionValue[unknown][j], *this);
                } else {
                    theElement_->basisFunction(
                        j, quadratureRule_, quadraturePointIndex_,
                        vectorBasisFunctionValue[unknown][j], unknown);
                    vectorBasisFunctionValue[unknown][j] =
                        transform_[unknown]->transform(
                            vectorBasisFunctionValue[unknown][j], *this);
                }
            } else {
                theElement_->basisFunction(j, pointReference_,
                                           vectorBasisFunctionValue[unknown][j],
                                           unknown);
                vectorBasisFunctionValue[unknown][j] =
                    transform_[unknown]->transform(
                        vectorBasisFunctionValue[unknown][j], *this);
            }
        }
        result = vectorBasisFunctionValue[unknown][i];
    }
}

template <std::size_t DIM>
inline const LinearAlgebra::SmallVector<DIM>&
    PhysicalElement<DIM>::basisFunctionCurl(std::size_t i) {
    logger.assert_debug(hasPointReference && hasElement,
                        "Need a location to evaluate the data");
    if (hasFunctionCurl[0]) {
        return basisFunctionCurl_[0][i];
    }
    hasFunctionCurl[0] = true;
    for (std::size_t j = 0; j < theElement_->getNumberOfBasisFunctions(); ++j) {
        if (hasQuadratureRule) {
            if (doesMapQuadraturePointFromFace) {
                basisFunctionCurl_[0][j] = transform_[0]->transformCurl(
                    theElement_->basisFunctionCurl<DIM>(j, quadratureRule_,
                                                        quadraturePointIndex_,
                                                        faceToElementMap_),
                    *this);
            } else {
                basisFunctionCurl_[0][j] = transform_[0]->transformCurl(
                    theElement_->basisFunctionCurl<DIM>(j, quadratureRule_,
                                                        quadraturePointIndex_),
                    *this);
            }
        } else {
            basisFunctionCurl_[0][j] = transform_[0]->transformCurl(
                theElement_->basisFunctionCurl(j, pointReference_), *this);
        }
    }
    return basisFunctionCurl_[0][i];
}

template <std::size_t DIM>
inline const LinearAlgebra::SmallVector<DIM>&
    PhysicalElement<DIM>::basisFunctionCurl(std::size_t i,
                                            std::size_t unknown) {
    logger.assert_debug(hasPointReference && hasElement,
                        "Need a location to evaluate the data");
    logger.assert_debug(unknown < theElement_->getNumberOfUnknowns(),
                        "Unknown % does not exist", unknown);
    if (hasFunctionCurl[unknown]) {
        return basisFunctionCurl_[unknown][i];
    }
    hasFunctionCurl[unknown] = true;
    for (std::size_t j = 0; j < theElement_->getNumberOfBasisFunctions(unknown);
         ++j) {
        if (hasQuadratureRule) {
            if (doesMapQuadraturePointFromFace) {
                basisFunctionCurl_[unknown][j] =
                    transform_[unknown]->transformCurl(
                        theElement_->basisFunctionCurl<DIM>(
                            j, quadratureRule_, quadraturePointIndex_,
                            faceToElementMap_, unknown),
                        *this);
            } else {
                basisFunctionCurl_[unknown][j] =
                    transform_[unknown]->transformCurl(
                        theElement_->basisFunctionCurl<DIM>(
                            j, quadratureRule_, quadraturePointIndex_, unknown),
                        *this);
            }
        } else {
            basisFunctionCurl_[unknown][j] = transform_[unknown]->transformCurl(
                theElement_->basisFunctionCurl(j, pointReference_, unknown),
                *this);
        }
    }
    return basisFunctionCurl_[unknown][i];
}

template <std::size_t DIM>
inline const double& PhysicalElement<DIM>::basisFunctionDiv(std::size_t i) {
    logger.assert_debug(hasPointReference && hasElement,
                        "Need a location to evaluate the data");
    if (hasFunctionDiv[0]) {
        return basisFunctionDiv_[0][i];
    }
    hasFunctionDiv[0] = true;
    for (std::size_t j = 0; j < theElement_->getNumberOfBasisFunctions(); ++j) {
        if (hasQuadratureRule) {
            if (doesMapQuadraturePointFromFace) {
                basisFunctionDiv_[0][j] = transform_[0]->transformDiv(
                    theElement_->basisFunctionDiv(j, quadratureRule_,
                                                  quadraturePointIndex_,
                                                  faceToElementMap_),
                    *this);
            } else {
                basisFunctionDiv_[0][j] = transform_[0]->transformDiv(
                    theElement_->basisFunctionDiv(j, quadratureRule_,
                                                  quadraturePointIndex_),
                    *this);
            }
        } else {
            basisFunctionDiv_[0][j] = transform_[0]->transformDiv(
                theElement_->basisFunctionDiv(j, pointReference_), *this);
        }
    }
    return basisFunctionDiv_[0][i];
}

template <std::size_t DIM>
inline const double& PhysicalElement<DIM>::basisFunctionDiv(
    std::size_t i, std::size_t unknown) {
    logger.assert_debug(hasPointReference && hasElement,
                        "Need a location to evaluate the data");
    if (hasFunctionDiv[unknown]) {
        return basisFunctionDiv_[unknown][i];
    }
    hasFunctionDiv[unknown] = true;
    for (std::size_t j = 0; j < theElement_->getNumberOfBasisFunctions(unknown);
         ++j) {
        if (hasQuadratureRule) {
            if (doesMapQuadraturePointFromFace) {
                basisFunctionDiv_[unknown][j] =
                    transform_[unknown]->transformDiv(
                        theElement_->basisFunctionDiv(
                            j, quadratureRule_, quadraturePointIndex_,
                            faceToElementMap_, unknown),
                        *this);
            } else {
                basisFunctionDiv_[unknown][j] =
                    transform_[unknown]->transformDiv(
                        theElement_->basisFunctionDiv(
                            j, quadratureRule_, quadraturePointIndex_, unknown),
                        *this);
            }
        } else {
            basisFunctionDiv_[unknown][j] = transform_[unknown]->transformDiv(
                theElement_->basisFunctionDiv(j, pointReference_, unknown),
                *this);
        }
    }
    return basisFunctionDiv_[unknown][i];
}

template <std::size_t DIM>
inline const LinearAlgebra::MiddleSizeVector&
    PhysicalElement<DIM>::getSolution() {
    logger.assert_debug(hasPointReference && hasElement,
                        "Need a location to evaluate the data");
    if (hasSolution) {
        return solution;
    }
    hasSolution = true;

    std::size_t numberOfUnknowns = getNumberOfUnknowns();
    std::vector<std::size_t> numberOfBasisFunctions =
            std::vector<std::size_t>(numberOfUnknowns, 0);

    solution.resize(numberOfUnknowns);
    solution.set(0.0);

    const LinearAlgebra::MiddleSizeVector& data =
        theElement_->getTimeIntegrationVector(0);

    std::size_t iVb = 0;
    for (std::size_t iV = 0; iV < numberOfUnknowns; ++iV) {
        numberOfBasisFunctions[iV] =
            theElement_->getNumberOfBasisFunctions(iV);
        for (std::size_t iB = 0; iB < numberOfBasisFunctions[iV]; ++iB) {
            iVb = theElement_->convertToSingleIndex(iB, iV);
            solution[iV] += data(iVb) * basisFunction(iB);
        }
    }
    return solution;
}

template <std::size_t DIM>
inline const std::vector<LinearAlgebra::SmallVector<DIM> >&
    PhysicalElement<DIM>::getSolutionDeriv() {
    logger.assert_debug(hasPointReference && hasElement,
                        "Need a location to evaluate the data");
    if (hasSolutionDeriv) {
        return solutionDeriv;
    }
    hasSolutionDeriv = true;

    std::size_t numberOfUnknowns = getNumberOfUnknowns();
    std::vector<std::size_t> numberOfBasisFunctions =
        std::vector<std::size_t>(numberOfUnknowns, 0);
    solutionDeriv.resize(numberOfUnknowns);

    const LinearAlgebra::MiddleSizeVector& data =
        theElement_->getTimeIntegrationVector(0);

    std::size_t iVB = 0;
    for (std::size_t iV = 0; iV < numberOfUnknowns; ++iV) {
        solutionDeriv[iV].set(0.0);
        numberOfBasisFunctions[iV] = theElement_->getNumberOfBasisFunctions(iV);
        for (std::size_t iB = 0; iB < numberOfBasisFunctions[iV]; ++iB) {
            iVB = convertToSingleIndex(iB, iV);
            solutionDeriv[iV] += data(iVB) * basisFunctionDeriv(iB);
        }
    }
    return solutionDeriv;
}

template <std::size_t DIM>
inline void PhysicalElement<DIM>::getSolution(
    std::vector<LinearAlgebra::SmallVector<DIM> >& result) {
    logger(ERROR, "not supported by element yet");
}

template <std::size_t DIM>
inline const std::vector<LinearAlgebra::SmallVector<DIM> >&
    PhysicalElement<DIM>::getSolutionCurl() {
    logger(ERROR, "not supported by element yet");
    return std::vector<LinearAlgebra::SmallVector<DIM> >();
}

template <std::size_t DIM>
inline const LinearAlgebra::MiddleSizeVector&
    PhysicalElement<DIM>::getSolutionDiv() {
    logger(ERROR, "not supported by element yet");
    // just inlining a default-constructed vector complains about stack return
    // despite being dead code
    static LinearAlgebra::MiddleSizeVector dummy;
    return dummy;
}

template <std::size_t DIM>
inline const Geometry::PointReference<DIM>&
    PhysicalElement<DIM>::getPointReference() {
    logger.assert_debug(hasPointReference,
                        "Need a location to evaluate the data");
    return pointReference_;
}

template <std::size_t DIM>
inline const Geometry::PointPhysical<DIM>&
    PhysicalElement<DIM>::getPointPhysical() {
    logger.assert_debug(hasPointReference && hasElement,
                        "Need a location to evaluate the data");
    if (hasPointPhysical) {
        return pointPhysical;
    }
    hasPointPhysical = true;
    pointPhysical = theElement_->referenceToPhysical(pointReference_);
    return pointPhysical;
}

template <std::size_t DIM>
inline const Geometry::Jacobian<DIM, DIM>& PhysicalElement<DIM>::getJacobian()
    const {
    logger.assert_debug(hasPointReference && hasElement,
                        "Need a location to evaluate the data");
    if (hasJacobian) {
        return jacobian;
    }
    hasJacobian = true;
    jacobian = theElement_->calcJacobian(pointReference_);
    return jacobian;
}

template <std::size_t DIM>
inline const Geometry::Jacobian<DIM, DIM>&
    PhysicalElement<DIM>::getInverseTransposeJacobian() {
    logger.assert_debug(hasPointReference && hasElement,
                        "Need a location to evaluate the data");
    if (hasInverseTransposeJacobian) {
        return inverseTransposeJacobian;
    }
    hasInverseTransposeJacobian = true;
    inverseTransposeJacobian = getTransposeJacobian().inverse();
    return inverseTransposeJacobian;
}

template <std::size_t DIM>
inline const Geometry::Jacobian<DIM, DIM>&
    PhysicalElement<DIM>::getTransposeJacobian() const {
    logger.assert_debug(hasPointReference && hasElement,
                        "Need a location to evaluate the data");
    if (hasTransposeJacobian) {
        return transposeJacobian;
    }
    hasTransposeJacobian = true;
    transposeJacobian = getJacobian().transpose();
    return transposeJacobian;
}

template <std::size_t DIM>
inline double PhysicalElement<DIM>::getJacobianAbsDet() {
    logger.assert_debug(hasPointReference && hasElement,
                        "Need a location to evaluate the data");
    if (hasJacobianAbsDet) {
        return jacobianAbsDet;
    }
    hasJacobianAbsDet = true;
    jacobianAbsDet = std::abs(getJacobianDet());
    return jacobianAbsDet;
}

template <std::size_t DIM>
inline double PhysicalElement<DIM>::getJacobianDet() const {
    logger.assert_debug(hasPointReference && hasElement,
                        "Need a location to evaluate the data");
    if (hasJacobianDet) {
        return jacobianDet;
    }
    hasJacobianDet = true;
    jacobianDet = getJacobian().determinant();
    return jacobianDet;
}

template <std::size_t DIM>
inline LinearAlgebra::MiddleSizeMatrix&
    PhysicalElement<DIM>::getResultMatrix() {
    logger.assert_debug(hasPointReference && hasElement,
                        "Need a location to evaluate the data");
    logger.assert_debug(hasElementMatrix,
                        "Can only provide the matrix once per coordinate");
    hasElementMatrix = false;
    return resultMatrix;
}

template <std::size_t DIM>
inline LinearAlgebra::MiddleSizeVector&
    PhysicalElement<DIM>::getResultVector() {
    logger.assert_debug(hasPointReference && hasElement,
                        "Need a location to evaluate the data");
    logger.assert_debug(hasElementVector,
                        "Can only provide the vector once per coordinate");
    hasElementVector = false;
    return resultVector;
}

template <std::size_t DIM>
inline const Base::Element* PhysicalElement<DIM>::getElement() {
    logger.assert_debug(hasElement, "Need a location to evaluate the data");
    return theElement_;
}

template <std::size_t DIM>
inline const Base::CoordinateTransformation<DIM>*
    PhysicalElement<DIM>::getTransformation() {
    for (std::size_t i = 1; i < transform_.size(); ++i)
        logger.assert_debug(
            transform_[0] == transform_[i],
            "Different unknowns have different coordinate transformation, call "
            "getTransformation for a given unknown");
    return transform_[0].get();
}

template <std::size_t DIM>
inline const Base::CoordinateTransformation<DIM>*
    PhysicalElement<DIM>::getTransformation(std::size_t unknown) {
    logger.assert_debug(unknown < theElement_->getNumberOfUnknowns(),
                        "Unknown % does not exist", unknown);
    return transform_[unknown].get();
}

template <std::size_t DIM>
inline void PhysicalElement<DIM>::setPointReference(
    const Geometry::PointReference<DIM>& point) {
    pointReference_ = point;
    hasQuadratureRule = false;
    hasPointReference = true;
    // even if they are already computed, the information is now out of date
    std::size_t unknowns = theElement_->getNumberOfUnknowns();
    hasFunctionValue.assign(unknowns, false);
    hasVectorFunctionValue.assign(unknowns, false);
    hasFunctionDeriv.assign(unknowns, false);
    hasFunctionCurl.assign(unknowns, false);
    hasFunctionDiv.assign(unknowns, false);
    hasSolution = false;
    hasVectorSolution = false;
    hasSolutionDeriv = false;
    hasSolutionCurl = false;
    hasSolutionDiv = false;
    hasPointPhysical = false;
    hasJacobian = false;
    hasTransposeJacobian = false;
    hasInverseTransposeJacobian = false;
    hasJacobianDet = false;
    hasJacobianAbsDet = false;
    if (!hasElementMatrix) {
        resultMatrix *= 0;
    }
    if (!hasElementVector) {
        resultVector *= 0;
    }
    hasElementMatrix = true;
    hasElementVector = true;
}

template <std::size_t DIM>
inline void PhysicalElement<DIM>::setElement(const Element* element) {
    if (!hasElement || element->getTotalNumberOfBasisFunctions() !=
                           theElement_->getTotalNumberOfBasisFunctions()) {
        basisFunctionValue.resize(element->getNumberOfUnknowns());
        vectorBasisFunctionValue.resize(element->getNumberOfUnknowns());
        basisFunctionDeriv_.resize(element->getNumberOfUnknowns());
        basisFunctionCurl_.resize(element->getNumberOfUnknowns());
        basisFunctionDiv_.resize(element->getNumberOfUnknowns());
        std::size_t numberOfEntries = 0;
        for (std::size_t i = 0; i < element->getNumberOfUnknowns(); ++i) {
            numberOfEntries += element->getNumberOfBasisFunctions(i);
            basisFunctionValue[i].resize(element->getNumberOfBasisFunctions(i));
            vectorBasisFunctionValue[i].resize(
                element->getNumberOfBasisFunctions(i));
            basisFunctionDeriv_[i].resize(
                element->getNumberOfBasisFunctions(i));
            basisFunctionCurl_[i].resize(element->getNumberOfBasisFunctions(i));
            basisFunctionDiv_[i].resize(element->getNumberOfBasisFunctions(i));
        }
        resultMatrix.resize(numberOfEntries, numberOfEntries);
        resultVector.resize(numberOfEntries);
    }
    theElement_ = element;
    hasElement = true;
    // even if they are already computed, the information is now out of date
    std::size_t unknowns = theElement_->getNumberOfUnknowns();

    hasFunctionValue.assign(unknowns, false);
    hasVectorFunctionValue.assign(unknowns, false);
    hasFunctionDeriv.assign(unknowns, false);
    hasFunctionCurl.assign(unknowns, false);
    hasFunctionDiv.assign(unknowns, false);
    hasSolution = false;
    hasVectorSolution = false;
    hasSolutionDeriv = false;
    hasSolutionCurl = false;
    hasSolutionDiv = false;
    hasPointPhysical = false;
    hasJacobian = false;
    hasTransposeJacobian = false;
    hasInverseTransposeJacobian = false;
    hasJacobianDet = false;
    hasJacobianAbsDet = false;
    if (!hasElementMatrix) {
        resultMatrix *= 0;
    }
    if (!hasElementVector) {
        resultVector *= 0;
    }
    hasElementMatrix = true;
    hasElementVector = true;
}

template <std::size_t DIM>
inline void PhysicalElement<DIM>::setTransformation(
    std::shared_ptr<Base::CoordinateTransformation<DIM> >& transform,
    std::size_t unknown) {
    if (transform_.size() <= unknown) {
        // We should not need to resize when we know the exact number of
        // unknowns.
        logger.assert_debug(!hasElement,
                            "Resizing the transform with an element");
        // Grow the transform vector automatically as we do not know how big it
        // should be.
        transform_.resize(unknown + 1);
    }
    transform_[unknown] = transform;
    // even if they are already computed, the information is now out of date
    // Note that if there is no element, the setElement will resize them
    // to the correct size. When there is an element this will just clear
    // the current values.
    hasFunctionValue.assign(hasFunctionValue.size(), false);
    hasVectorFunctionValue.assign(hasVectorFunctionValue.size(), false);
    hasFunctionDeriv.assign(hasFunctionDeriv.size(), false);
    hasFunctionCurl.assign(hasFunctionCurl.size(), false);
    hasFunctionDiv.assign(hasFunctionDiv.size(), false);
    hasSolution = false;
    hasVectorSolution = false;
    hasSolutionDeriv = false;
    hasSolutionCurl = false;
    hasSolutionDiv = false;
    if (!hasElementMatrix) {
        resultMatrix *= 0;
    }
    if (!hasElementVector) {
        resultVector *= 0;
    }
    hasElementMatrix = true;
    hasElementVector = true;
}

template <std::size_t DIM>
inline void PhysicalElement<DIM>::setQuadratureRule(
    QuadratureRules::GaussQuadratureRule* rule) {
    quadratureRule_ = rule;
    hasQuadratureRule = true;
    doesMapQuadraturePointFromFace = false;
    setQuadraturePointIndex(0);
}

template <std::size_t DIM>
inline void PhysicalElement<DIM>::setQuadratureRule(
    QuadratureRules::GaussQuadratureRule* rule,
    const Geometry::MappingReferenceToReference<1>* map) {
    quadratureRule_ = rule;
    faceToElementMap_ = map;
    hasQuadratureRule = true;
    doesMapQuadraturePointFromFace = true;
    setQuadraturePointIndex(0);
}

template <std::size_t DIM>
inline void PhysicalElement<DIM>::setQuadraturePointIndex(std::size_t index) {
    logger.assert_debug(hasQuadratureRule,
                        "It only makes sense to index into a quadrature rule "
                        "if you set one first");
    quadraturePointIndex_ = index;
    if (doesMapQuadraturePointFromFace) {
        Geometry::PointReference<DIM - 1> pointOnFace =
            quadratureRule_->getPoint(index);
        setPointReference(faceToElementMap_->transform(pointOnFace));
    } else {
        setPointReference(quadratureRule_->getPoint(index));
    }
    hasQuadratureRule = true;
}
}  // namespace Base
}  // namespace hpgem
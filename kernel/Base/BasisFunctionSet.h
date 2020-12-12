/*
This file forms part of hpGEM. This package has been developed over a number of
years by various people at the University of Twente and a full list of
contributors can be found at http://hpgem.org/about-the-code/team

This code is distributed using BSD 3-Clause License. A copy of which can found
below.


Copyright (c) 2014, University of Twente
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef HPGEM_KERNEL_BASISFUNCTIONSET_H
#define HPGEM_KERNEL_BASISFUNCTIONSET_H

#include <vector>
#include "Logger.h"
#include "BaseBasisFunction.h"
namespace hpgem {
namespace LinearAlgebra {
template <std::size_t DIM>
class SmallVector;
}

namespace Geometry {
template <std::size_t DIM>
class PointReference;

template <int CODIM>
class MappingReferenceToReference;
}  // namespace Geometry

namespace QuadratureRules {
class GaussQuadratureRule;
}

namespace Base {
class BaseBasisFunction;

class BasisFunctionSet {
   public:
    using BaseBasisFunctions = std::vector<BaseBasisFunction *>;  // check again

    explicit BasisFunctionSet(std::size_t order);

    // BasisFunctionSets should not be copied, therefore the copy constructor is
    // deleted.
    BasisFunctionSet(const BasisFunctionSet &other) = delete;

    virtual ~BasisFunctionSet();

    std::size_t size() const;

    std::size_t getOrder() const;

    void addBasisFunction(BaseBasisFunction *bf);

    template <std::size_t DIM>
    double eval(std::size_t i, const Geometry::PointReference<DIM> &p) const;

    ///\evaluate basis function i at a point needed for a quadrature rule, where
    /// the quadrature rule is meant for integrating over an element
    double eval(std::size_t i,
                QuadratureRules::GaussQuadratureRule *elementQuadratureRule,
                std::size_t quadraturePointIndex) const;

    ///\evaluate basis function i at a point needed for a quadrature rule, where
    /// the quadrature rule is meant for integrating over a face
    double eval(
        std::size_t i, QuadratureRules::GaussQuadratureRule *faceQuadratureRule,
        std::size_t quadraturePointIndex,
        const Geometry::MappingReferenceToReference<1> *faceToElementMap) const;

    ///\brief returns the value of the i-th basisfunction at point p in ret
    template <std::size_t DIM>
    void eval(std::size_t i, const Geometry::PointReference<DIM> &p,
              LinearAlgebra::SmallVector<DIM> &ret) const;

    ///\evaluate basis function i at a point needed for a quadrature rule, where
    /// the quadrature rule is meant for integrating over an element
    template <std::size_t DIM>
    void eval(std::size_t i,
              QuadratureRules::GaussQuadratureRule *elementQuadratureRule,
              std::size_t quadraturePointIndex,
              LinearAlgebra::SmallVector<DIM> &result) const;

    ///\evaluate basis function i at a point needed for a quadrature rule, where
    /// the quadrature rule is meant for integrating over a face
    template <std::size_t DIM>
    void eval(std::size_t i,
              QuadratureRules::GaussQuadratureRule *faceQuadratureRule,
              std::size_t quadraturePointIndex,
              const Geometry::MappingReferenceToReference<1> *faceToElementMap,
              LinearAlgebra::SmallVector<DIM> &result) const;

    template <std::size_t DIM>
    double evalDeriv(std::size_t i, std::size_t jDir,
                     const Geometry::PointReference<DIM> &p) const;

    template <std::size_t DIM>
    LinearAlgebra::SmallVector<DIM> evalDeriv(
        std::size_t i, const Geometry::PointReference<DIM> &p) const;

    ///\evaluate the gradient of basis function i at a point needed for a
    /// quadrature rule, where the quadrature rule is meant for integrating over
    /// an element
    template <std::size_t DIM>
    LinearAlgebra::SmallVector<DIM> evalDeriv(
        std::size_t i,
        QuadratureRules::GaussQuadratureRule *elementQuadratureRule,
        std::size_t quadraturePointIndex) const;

    ///\evaluate the gradient of basis function i at a point needed for a
    /// quadrature rule, where the quadrature rule is meant for integrating over
    /// a face
    template <std::size_t DIM>
    LinearAlgebra::SmallVector<DIM> evalDeriv(
        std::size_t i, QuadratureRules::GaussQuadratureRule *faceQuadratureRule,
        std::size_t quadraturePointIndex,
        const Geometry::MappingReferenceToReference<1> *faceToElementMap) const;

    ///\brief returns the curl of the i-th basisfunction at point p in ret
    template <std::size_t DIM>
    LinearAlgebra::SmallVector<DIM> evalCurl(
        std::size_t i, const Geometry::PointReference<DIM> &p) const;

    ///\evaluate the curl of basis function i at a point needed for a quadrature
    /// rule, where the quadrature rule is meant for integrating over an element
    template <std::size_t DIM>
    LinearAlgebra::SmallVector<DIM> evalCurl(
        std::size_t i,
        QuadratureRules::GaussQuadratureRule *elementQuadratureRule,
        std::size_t quadraturePointIndex) const;

    ///\evaluate the curl of basis function i at a point needed for a quadrature
    /// rule, where the quadrature rule is meant for integrating over a face
    template <std::size_t DIM>
    LinearAlgebra::SmallVector<DIM> evalCurl(
        std::size_t i, QuadratureRules::GaussQuadratureRule *faceQuadratureRule,
        std::size_t quadraturePointIndex,
        const Geometry::MappingReferenceToReference<1> *faceToElementMap) const;

    ///\brief returns the divergence of the i-th basisfunction at point p in ret
    template <std::size_t DIM>
    double evalDiv(std::size_t i, const Geometry::PointReference<DIM> &p) const;

    ///\evaluate the divergence of basis function i at a point needed for a
    /// quadrature rule, where the quadrature rule is meant for integrating over
    /// an element
    double evalDiv(std::size_t i,
                   QuadratureRules::GaussQuadratureRule *elementQuadratureRule,
                   std::size_t quadraturePointIndex) const;

    ///\evaluate the divergence of basis function i at a point needed for a
    /// quadrature rule, where the quadrature rule is meant for integrating over
    /// a face
    double evalDiv(
        std::size_t i, QuadratureRules::GaussQuadratureRule *faceQuadratureRule,
        std::size_t quadraturePointIndex,
        const Geometry::MappingReferenceToReference<1> *faceToElementMap) const;

    const BaseBasisFunction *operator[](std::size_t i) const {
        logger.assert_debug(
            i < size(),
            "Asked for basis function %, but there are only % basis functions",
            i, size());
        return vecOfBasisFcn_[i];
    }

    /// iterators (for range-based for loop)
    BaseBasisFunctions::const_iterator begin() const {
        return vecOfBasisFcn_.begin();
    }

    BaseBasisFunctions::iterator begin() { return vecOfBasisFcn_.begin(); }

    BaseBasisFunctions::const_iterator end() const {
        return vecOfBasisFcn_.end();
    }

    BaseBasisFunctions::iterator end() { return vecOfBasisFcn_.end(); }

    /// Register a function that will be called when this BasisFunctionSet is
    /// deconstructed.
    void registerDestructorListener(std::function<void()> callback) const {
        destructorListeners_.push_back(callback);
    }

   private:
    std::size_t order_;
    BaseBasisFunctions vecOfBasisFcn_;

    mutable std::vector<std::function<void()>> destructorListeners_;
};
}  // namespace Base
}  // namespace hpgem
#include "Integration/QuadratureRules/GaussQuadratureRule.h"

// reopening namespace because gcc has a bug where it doesn't accept template
// specializations in an enclosing namespace
namespace hpgem {
namespace Base {
template <std::size_t DIM>
double BasisFunctionSet::eval(std::size_t i,
                              const Geometry::PointReference<DIM> &p) const {
    logger.assert_debug(
        i < size(),
        "Asked for basis function %, but there are only % basis functions", i,
        size());
    return vecOfBasisFcn_[i]->eval(p);
}

// DIM = 1 doesn't have e.g. evalDeriv2 so we have to manually implement the 4
// special cases
template <std::size_t DIM>
double BasisFunctionSet::evalDeriv(
    std::size_t i, std::size_t jDir,
    const Geometry::PointReference<DIM> &p) const {
    logger(ERROR, "Only supports points of dimension up to 4");
    return 0.;
}

template <>
inline double BasisFunctionSet::evalDeriv(
    std::size_t i, std::size_t jDir,
    const Geometry::PointReference<1> &p) const {
    logger.assert_debug(
        i < size(),
        "Asked for basis function %, but there are only % basis functions", i,
        size());
    logger.assert_debug(
        (jDir < 1),
        "Error in BasisFunctionSet.EvalDeriv: invalid derivative direction!");

    switch (jDir) {
        case 0:
            return vecOfBasisFcn_[i]->evalDeriv0(p);
        default:
            return 0.;
    }
}

template <>
inline double BasisFunctionSet::evalDeriv(
    std::size_t i, std::size_t jDir,
    const Geometry::PointReference<2> &p) const {
    logger.assert_debug(
        i < size(),
        "Asked for basis function %, but there are only % basis functions", i,
        size());
    logger.assert_debug(
        (jDir < 2),
        "Error in BasisFunctionSet.EvalDeriv: invalid derivative direction!");

    switch (jDir) {
        case 0:
            return vecOfBasisFcn_[i]->evalDeriv0(p);
        case 1:
            return vecOfBasisFcn_[i]->evalDeriv1(p);
        default:
            return 0.;
    }
}

template <>
inline double BasisFunctionSet::evalDeriv(
    std::size_t i, std::size_t jDir,
    const Geometry::PointReference<3> &p) const {
    logger.assert_debug(
        i < size(),
        "Asked for basis function %, but there are only % basis functions", i,
        size());
    logger.assert_debug(
        (jDir < 3),
        "Error in BasisFunctionSet.EvalDeriv: invalid derivative direction!");

    switch (jDir) {
        case 0:
            return vecOfBasisFcn_[i]->evalDeriv0(p);
        case 1:
            return vecOfBasisFcn_[i]->evalDeriv1(p);
        case 2:
            return vecOfBasisFcn_[i]->evalDeriv2(p);
        default:
            return 0.;
    }
}

template <>
inline double BasisFunctionSet::evalDeriv(
    std::size_t i, std::size_t jDir,
    const Geometry::PointReference<4> &p) const {
    logger.assert_debug(
        i < size(),
        "Asked for basis function %, but there are only % basis functions", i,
        size());
    logger.assert_debug(
        (jDir < 4),
        "Error in BasisFunctionSet.EvalDeriv: invalid derivative direction!");

    switch (jDir) {
        case 0:
            return vecOfBasisFcn_[i]->evalDeriv0(p);
        case 1:
            return vecOfBasisFcn_[i]->evalDeriv1(p);
        case 2:
            return vecOfBasisFcn_[i]->evalDeriv2(p);
        case 3:
            return vecOfBasisFcn_[i]->evalDeriv3(p);
        default:
            return 0.;
    }
}

template <std::size_t DIM>
void BasisFunctionSet::eval(std::size_t i,
                            const Geometry::PointReference<DIM> &p,
                            LinearAlgebra::SmallVector<DIM> &ret) const {
    logger.assert_debug(
        i < size(),
        "Asked for basis function %, but there are only % basis functions", i,
        size());
    vecOfBasisFcn_[i]->eval(p, ret);
}

template <std::size_t DIM>
LinearAlgebra::SmallVector<DIM> BasisFunctionSet::evalDeriv(
    std::size_t i, const Geometry::PointReference<DIM> &p) const {
    logger.assert_debug(
        i < size(),
        "Asked for basis function %, but there are only % basis functions", i,
        size());
    return vecOfBasisFcn_[i]->evalDeriv(p);
}

template <std::size_t DIM>
LinearAlgebra::SmallVector<DIM> BasisFunctionSet::evalCurl(
    std::size_t i, const Geometry::PointReference<DIM> &p) const {
    logger.assert_debug(
        i < size(),
        "Asked for basis function %, but there are only % basis functions", i,
        size());
    return vecOfBasisFcn_[i]->evalCurl(p);
}

template <std::size_t DIM>
double BasisFunctionSet::evalDiv(std::size_t i,
                                 const Geometry::PointReference<DIM> &p) const {
    logger.assert_debug(
        i < size(),
        "Asked for basis function %, but there are only % basis functions", i,
        size());
    return vecOfBasisFcn_[i]->evalDiv(p);
}

template <std::size_t DIM>
inline void BasisFunctionSet::eval(
    std::size_t i, QuadratureRules::GaussQuadratureRule *elementQuadratureRule,
    std::size_t quadraturePointIndex,
    LinearAlgebra::SmallVector<DIM> &result) const {
    elementQuadratureRule->eval(this, i, quadraturePointIndex, result);
}

template <std::size_t DIM>
inline void BasisFunctionSet::eval(
    std::size_t i, QuadratureRules::GaussQuadratureRule *faceQuadratureRule,
    std::size_t quadraturePointIndex,
    const Geometry::MappingReferenceToReference<1> *faceToElementMap,
    LinearAlgebra::SmallVector<DIM> &result) const {
    faceQuadratureRule->eval(this, i, quadraturePointIndex, faceToElementMap,
                             result);
}

template <std::size_t DIM>
inline LinearAlgebra::SmallVector<DIM> BasisFunctionSet::evalDeriv(
    std::size_t i, QuadratureRules::GaussQuadratureRule *elementQuadratureRule,
    std::size_t quadraturePointIndex) const {
    return elementQuadratureRule->evalGrad(this, i, quadraturePointIndex);
}

template <std::size_t DIM>
inline LinearAlgebra::SmallVector<DIM> BasisFunctionSet::evalDeriv(
    std::size_t i, QuadratureRules::GaussQuadratureRule *faceQuadratureRule,
    std::size_t quadraturePointIndex,
    const Geometry::MappingReferenceToReference<1> *faceToElementMap) const {
    return faceQuadratureRule->evalGrad(this, i, quadraturePointIndex,
                                        faceToElementMap);
}

template <std::size_t DIM>
inline LinearAlgebra::SmallVector<DIM> BasisFunctionSet::evalCurl(
    std::size_t i, QuadratureRules::GaussQuadratureRule *elementQuadratureRule,
    std::size_t quadraturePointIndex) const {
    logger(ERROR, "only implemented for DIM = 3");
    return LinearAlgebra::SmallVector<DIM>();
}

template <>
inline LinearAlgebra::SmallVector<2> BasisFunctionSet::evalCurl(
    std::size_t i, QuadratureRules::GaussQuadratureRule *elementQuadratureRule,
    std::size_t quadraturePointIndex) const {
    return elementQuadratureRule->evalCurl2D(this, i, quadraturePointIndex);
}

template <>
inline LinearAlgebra::SmallVector<3> BasisFunctionSet::evalCurl(
    std::size_t i, QuadratureRules::GaussQuadratureRule *elementQuadratureRule,
    std::size_t quadraturePointIndex) const {
    return elementQuadratureRule->evalCurl(this, i, quadraturePointIndex);
}

template <std::size_t DIM>
inline LinearAlgebra::SmallVector<DIM> BasisFunctionSet::evalCurl(
    std::size_t i, QuadratureRules::GaussQuadratureRule *faceQuadratureRule,
    std::size_t quadraturePointIndex,
    const Geometry::MappingReferenceToReference<1> *faceToElementMap) const {
    logger(ERROR, "only implemented for DIM = 3");
    return LinearAlgebra::SmallVector<DIM>();
}

template <>
inline LinearAlgebra::SmallVector<3> BasisFunctionSet::evalCurl(
    std::size_t i, QuadratureRules::GaussQuadratureRule *faceQuadratureRule,
    std::size_t quadraturePointIndex,
    const Geometry::MappingReferenceToReference<1> *faceToElementMap) const {
    return faceQuadratureRule->evalCurl(this, i, quadraturePointIndex,
                                        faceToElementMap);
}

template <>
inline LinearAlgebra::SmallVector<2> BasisFunctionSet::evalCurl(
    std::size_t i, QuadratureRules::GaussQuadratureRule *faceQuadratureRule,
    std::size_t quadraturePointIndex,
    const Geometry::MappingReferenceToReference<1> *faceToElementMap) const {
    return faceQuadratureRule->evalCurl2D(this, i, quadraturePointIndex,
                                          faceToElementMap);
}

}  // namespace Base
}  // namespace hpgem

#endif  // HPGEM_KERNEL_BASISFUNCTIONSET_H

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

#ifndef ELEMENTINTEGRAL_IMPL_HPP_
#define ELEMENTINTEGRAL_IMPL_HPP_

#include "Logger.h"
#include "Base/Element.h"
#include "ElementIntegrandBase.h"
#include "QuadratureRules/GaussQuadratureRule.h"
#include "Geometry/ReferenceGeometry.h"
#include "ElementIntegral.h"
#include "LinearAlgebra/Axpy.h"
#include "Base/DoNotScaleIntegrands.h"

namespace Integration {
/*!
\param[in]  el        the Element to be integrated on,
\param[in]  rule      the GaussQuadratureRule to use,
\param[in]  integrand a function/functor with operator(Element, PointReference,
ResultType&), \return    a reference to the variable with result storage.

This function integrates the function in integrand over the element el with
the given Gauss Quadrature rule.
Note that one has the possibility to leave the rule argument
away, in which case the default for the ReferenceGeometry of
the passed element will be used.
\deprecated Please use integrate(Element*, std::function<...>, const
QuadratureRulesT) wherever possible.
 */
template <std::size_t DIM>
template <typename ReturnTrait1>
ReturnTrait1 ElementIntegral<DIM>::integrate(
    const Base::Element* el, ElementIntegrandBase<ReturnTrait1, DIM>* integrand,
    QuadratureRules::GaussQuadratureRule* qdrRule) {
    logger.assert_debug(el != nullptr, "Invalid element detected");
    logger.assert_debug(integrand != nullptr, "Invalid integrand detected");
    // quadrature rule is allowed to be equal to nullptr!
    std::function<ReturnTrait1(Base::PhysicalElement<DIM>&)> integrandFun =
        [=](Base::PhysicalElement<DIM>& el) -> ReturnTrait1 {
        ReturnTrait1 result;
        integrand->elementIntegrand(el, result);
        return result;
    };
    return integrate(el, integrandFun, qdrRule);
}

/*!
\param[in]  el        the Element to be integrated on
\param[in]  qdrRule      the GaussQuadratureRule to use
\param[in]  integrandFun a function with parameters Element, PointReference,
ResultType& \return    a reference to the variable with result storage

This function integrates the function in integrand over the element el with
the given Gauss Quadrature rule.
Note that one has the possibility to leave the rule argument
away, in which case the default for the ReferenceGeometry of
the passed element will be used.
*/
template <std::size_t DIM>
template <typename FunctionType>
std::result_of_t<FunctionType(Base::PhysicalElement<DIM>&)>
    ElementIntegral<DIM>::integrate(
        const Base::Element* el, FunctionType integrandFun,
        QuadratureRules::GaussQuadratureRule* qdrRule) {
    using ReturnType =
        std::result_of_t<FunctionType(Base::PhysicalElement<DIM>&)>;

    logger.assert_debug(el != nullptr, "Invalid element detected");
    element_.setElement(el);
    // quadrature rule is allowed to be equal to nullptr!
    QuadratureRules::GaussQuadratureRule* qdrRuleLoc =
        (qdrRule == nullptr ? el->getGaussQuadratureRule() : qdrRule);

    // check whether the GaussQuadratureRule is actually for the element's
    // ReferenceGeometry
    logger.assert_debug(
        (qdrRuleLoc->forReferenceGeometry() == el->getReferenceGeometry()),
        "ElementIntegral: wrong geometry.");

    // value returned by the integrand
    ReturnType value, result;

    // number of Gauss quadrature points
    std::size_t numberOfPoints = qdrRuleLoc->getNumberOfPoints();
    logger.assert_debug(
        numberOfPoints > 0,
        "Did not get any points from qdrRuleLoc->getNumberOfPoints");

    element_.setQuadratureRule(qdrRuleLoc);
    // element_.setPointReference(qdrRuleLoc->getPoint(0));

    // first Gauss point
    // first we calculate the jacobian, then compute the function value on one
    // of the reference points and finally we multiply this value with a weight
    // and the jacobian and save it in result.

    result = integrandFun(element_);
    // We use the same quadrature rule for all unknowns.
    result *=
        (qdrRuleLoc->weight(0) *
         element_.getTransformation(0)->getIntegrandScaleFactor(element_));

    // next Gauss points, again calculate the jacobian, value at gauss point and
    // add this value multiplied with jacobian and weight to result.
    for (std::size_t i = 1; i < numberOfPoints; ++i) {
        element_.setQuadraturePointIndex(i);
        // element_.setPointReference(qdrRuleLoc->getPoint(i));
        value = integrandFun(element_);

        // axpy: Y = alpha * X + Y
        LinearAlgebra::axpy(
            qdrRuleLoc->weight(i) *
                element_.getTransformation(0)->getIntegrandScaleFactor(
                    element_),
            value, result);
    }
    return result;
}

/// \param[in] ptrQdrRule A pointer to a quadrature rule used for the
/// integration. \param[in] integrandFunction A function that is integrated on
/// the reference element. It takes as input argument a reference point and
/// returns an object of the class <IntegrandType>.
/*!
 \details This function computes the integral of a function \f$ f_{ref}\f$ on a
 reference element \f$ E_{ref} \f$, so it returns the following value \f[
 \int_{E_{ref}} f_{ref}(\xi) \,d\xi, \f] where \f$ f_{ref}:E_{ref}\rightarrow
 R\f$, with \f$ R \f$ a linear function space (e.g. space of vectors/matrices).
 In many cases the weak formulation is based on integrals on a physical element
 \f$ E_{phys} \f$ of the form given below \f[ \int_{E_{phys}} f_{phys}(x)
 \,dx.\f] Let \f$ \phi:E_{ref}\rightarrow E_{phys} \f$ be the mapping from the
 reference element to the physical element, let \f$ J\f$ be the Jacobian of \f$
 \phi\f$ and \f$ |J| \f$ the reference-to-physical element scale (usually the
 absolute value of the determinant of J). Then we can write \f[ \int_{E_{phys}}
 f_{phys}(x) \,dx = \int_{E_{ref}} f_{phys}(\phi(\xi)) |J| \,d\xi, \f] so \f$
 f_{ref}(\xi) = f_{phys}(\phi(\xi)) |J| \f$. In some cases it is more
 advantageous to compute \f$ f_{ref}(\xi) \f$ instead of \f$ f_{phys}(x)\f$.

 NOTE: do not mix up gradients of pyhsical and reference basis functions with
 integrals on physical and reference elements. If \f$ f_{phys}(x) \f$ contains a
 (physical) gradient of a physical basis function then so does \f$ f_{ref}(\xi)
 = f_{phys}(\phi(\xi)) |J| \f$. The difference is the input argument (reference
 point \f$ \xi \f$ instead of physical point \f$ x \f$ ) and the scaling \f$ |J|
 \f$. (Ofcourse it is possible to rewrite the gradient of a physical basis
 function in terms of the gradient of the corresponding reference basis
 function). This routine assumes the user takes responsibility of setting the
 correct element in the physical element
 */
template <std::size_t DIM>
template <typename IntegrandType>
IntegrandType ElementIntegral<DIM>::referenceElementIntegral(
    const QuadratureRules::GaussQuadratureRule* ptrQdrRule,
    std::function<IntegrandType()> integrandFunction) {
    std::size_t numberOfPoints = ptrQdrRule->getNumberOfPoints();
    std::size_t iPoint = 0;  // Index for the quadrature points.

    const Geometry::PointReference<DIM>& pRef0 = ptrQdrRule->getPoint(iPoint);
    element_.setQuadratureRule(ptrQdrRule);
    // element_.setPointReference(pRef0);
    IntegrandType integral(ptrQdrRule->weight(iPoint) * integrandFunction());
    for (iPoint = 1; iPoint < numberOfPoints; iPoint++) {
        element_.setQuadraturePointIndex(iPoint);
        // element_.setPointReference(ptrQdrRule->getPoint(iPoint));
        LinearAlgebra::axpy(ptrQdrRule->weight(iPoint), integrandFunction(),
                            integral);
    }
    return integral;
}

//! \brief Construct an ElementIntegral with cache on.
template <std::size_t DIM>
ElementIntegral<DIM>::ElementIntegral() {}

//! \brief Class destructor
template <std::size_t DIM>
ElementIntegral<DIM>::~ElementIntegral() {}

template <std::size_t DIM>
void ElementIntegral<DIM>::setTransformation(
    std::shared_ptr<Base::CoordinateTransformation<DIM> > transform,
    std::size_t unknown) {
    element_.setTransformation(transform, unknown);
}

template <std::size_t DIM>
Base::CoordinateTransformation<DIM>& ElementIntegral<DIM>::getTransformation(
    std::size_t unknown) {
    return element_.getTransformation(unknown);
}

template <std::size_t DIM>
Base::PhysicalElement<DIM>& ElementIntegral<DIM>::getPhysicalElement() {
    return element_;
}
}  // namespace Integration

#endif /* ELEMENTINTEGRAL_IMPL_HPP_ */

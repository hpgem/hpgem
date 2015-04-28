/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, University of Twente
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef ELEMENTINTEGRAL_IMPL_HPP_
#define ELEMENTINTEGRAL_IMPL_HPP_

#include "Base/ShortTermStorageElementH1.h"
#include "Logger.h"
#include "ElementIntegrandBase.h"
#include "QuadratureRules/GaussQuadratureRule.h"
#include "Geometry/ReferenceGeometry.h"
#include "ElementIntegral.h"
#include "LinearAlgebra/Axpy.h"

namespace Integration
{
    /*!
    \param[in]  el        the Element to be integrated on,
    \param[in]  rule      the GaussQuadratureRule to use,
    \param[in]  integrand a function/functor with operator(Element, PointReference, ResultType&),
    \return    a reference to the variable with result storage.

    This function integrates the function in integrand over the element el with
    the given Gauss Quadrature rule.
    Note that one has the possibility to leave the rule argument
    away, in which case the default for the ReferenceGeometry of
    the passed element will be used.
    \deprecated Please use integrate(Element*, std::function<...>, const QuadratureRulesT) wherever possible.
     */
    template<typename ReturnTrait1>
    ReturnTrait1 ElementIntegral::integrate(Base::Element* el, ElementIntegrandBase<ReturnTrait1>* integrand, const QuadratureRulesT * const qdrRule)
    {
        logger.assert(el!=nullptr, "Invalid element detected");
        logger.assert(integrand!=nullptr, "Invalid integrand detected");
        //quadrature rule is allowed to be equal to nullptr!
        std::function<ReturnTrait1(const Base::Element*, const PointReferenceT&)> integrandFun = [=](const Base::Element* el, const PointReferenceT& p)-> ReturnTrait1
        {   
            ReturnTrait1 result;
            integrand -> elementIntegrand(el, p, result);
            return result;
        };
        return integrate(el, integrandFun, qdrRule);
    }
    
    /*!
    \param[in]  el        the Element to be integrated on
    \param[in]  qdrRule      the GaussQuadratureRule to use
    \param[in]  integrandFun a function with parameters Element, PointReference, ResultType&
    \return    a reference to the variable with result storage

    This function integrates the function in integrand over the element el with
    the given Gauss Quadrature rule.
    Note that one has the possibility to leave the rule argument
    away, in which case the default for the ReferenceGeometry of
    the passed element will be used.
    */
    template<typename ReturnType>
    ReturnType ElementIntegral::integrate(Base::Element* el, std::function<ReturnType(const Base::Element*, const PointReferenceT&)> integrandFun, const QuadratureRulesT * const qdrRule)
    {
        logger.assert(el!=nullptr, "Invalid element detected");
        //quadrature rule is allowed to be equal to nullptr!
        if (localElement_ == nullptr)
        {
            localElement_ = new Base::ShortTermStorageElementH1(el->getGaussQuadratureRule()->dimension());
            if (useCache_)
            {
                localElement_->cacheOn();
            }
            else
            {
                localElement_->cacheOff();
            }
        }
        *localElement_ = *el;
        const QuadratureRulesT * const qdrRuleLoc = (qdrRule == nullptr ? localElement_->getGaussQuadratureRule() : qdrRule);
        
        // check whether the GaussQuadratureRule is actually for the element's ReferenceGeometry
        logger.assert((qdrRuleLoc->forReferenceGeometry() == localElement_->getReferenceGeometry()), "ElementIntegral: wrong geometry.");
        
        // value returned by the integrand
        ReturnType value, result;
        
        // number of Gauss quadrature points
        std::size_t nrOfPoints = qdrRuleLoc->nrOfPoints();
        logger.assert(nrOfPoints > 0, "Did not get any points from qdrRuleLoc->nrOfPoints");
        
        // Initialize Gauss quadrature point
        const Geometry::PointReference& p0 = qdrRuleLoc->getPoint(0);
        
        // first Gauss point
        // first we calculate the jacobian, then compute the function value on one of
        // the reference points and finally we multiply this value with a weight and
        // the jacobian and save it in result.
        
        Geometry::Jacobian jac = localElement_->calcJacobian(p0);
        result = integrandFun(localElement_, p0);
        result *= (qdrRuleLoc->weight(0) * std::abs(jac.determinant()));
        
        // next Gauss points, again calculate the jacobian, value at gauss point and
        // add this value multiplied with jacobian and weight to result.
        for (std::size_t i = 1; i < nrOfPoints; ++i)
        {
            
            const Geometry::PointReference& p = qdrRuleLoc->getPoint(i);
            jac = localElement_->calcJacobian(p);
            value = integrandFun(localElement_, p);
            
            //axpy: Y = alpha * X + Y
            LinearAlgebra::axpy(qdrRuleLoc->weight(i) * std::abs(jac.determinant()), value, result);
        }
        return result;
    }
    
    /// \param[in] ptrQdrRule A pointer to a quadrature rule used for the integration.
    /// \param[in] integrandFunction A function that is integrated on the reference element. It takes as input argument a reference point and returns an object of the class <IntegrandType>.
    /*!
     \details This function computes the integral of a function \f$ f_{ref}\f$ on a reference element \f$ E_{ref} \f$, so it returns the following value
     \f[ \int_{E_{ref}} f_{ref}(\xi) \,d\xi, \f]
     where \f$ f_{ref}:E_{ref}\rightarrow R\f$, with \f$ R \f$ a linear function space (e.g. space of vectors/matrices). In many cases the weak formulation is based on integrals on a physical element \f$ E_{phys} \f$ of the form given below
     \f[ \int_{E_{phys}} f_{phys}(x) \,dx.\f]
     Let \f$ \phi:E_{ref}\rightarrow E_{phys} \f$ be the mapping from the reference element to the physical element, let \f$ J\f$ be the Jacobian of \f$ \phi\f$ and \f$ |J| \f$ the reference-to-physical element scale (usually the absolute value of the determinant of J). Then we can write
     \f[ \int_{E_{phys}} f_{phys}(x) \,dx = \int_{E_{ref}} f_{phys}(\phi(\xi)) |J| \,d\xi, \f]
     so \f$ f_{ref}(\xi) = f_{phys}(\phi(\xi)) |J| \f$. In some cases it is more advantageous to compute \f$ f_{ref}(\xi) \f$ instead of \f$ f_{phys}(x)\f$.
     
     NOTE: do not mix up gradients of pyhsical and reference basis functions with integrals on physical and reference elements. If \f$ f_{phys}(x) \f$ contains a (physical) gradient of a physical basis function then so does \f$ f_{ref}(\xi) = f_{phys}(\phi(\xi)) |J| \f$. The difference is the input argument (reference point \f$ \xi \f$ instead of physical point \f$ x \f$ ) and the scaling \f$ |J| \f$. (Ofcourse it is possible to rewrite the gradient of a physical basis function in terms of the gradient of the corresponding reference basis function).
     */
    template<typename IntegrandType>
    IntegrandType ElementIntegral::referenceElementIntegral(const QuadratureRules::GaussQuadratureRule *ptrQdrRule, std::function<IntegrandType(const Geometry::PointReference &)> integrandFunction)
    {
        std::size_t numOfPoints = ptrQdrRule->nrOfPoints();
        std::size_t iPoint = 0; // Index for the quadrature points.
        
        const Geometry::PointReference& pRef0 = ptrQdrRule->getPoint(iPoint);
        IntegrandType integral(ptrQdrRule->weight(iPoint) * integrandFunction(pRef0));
        for (iPoint = 1; iPoint < numOfPoints; iPoint++)
        {
            const Geometry::PointReference& pRef = ptrQdrRule->getPoint(iPoint);
            LinearAlgebra::axpy(ptrQdrRule->weight(iPoint), integrandFunction(pRef), integral);
        }
        return integral;
    }

}

#endif /* ELEMENTINTEGRAL_IMPL_HPP_ */

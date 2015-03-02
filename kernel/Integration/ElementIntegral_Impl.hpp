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

#include <cassert>

#include "Base/ShortTermStorageElementH1.hpp"

#include "ElementIntegrandBase.hpp"
#include "QuadratureRules/GaussQuadratureRule.hpp"
#include "Geometry/ReferenceGeometry.hpp"
#include "Base/TestErrorDebug.hpp"
#include "ElementIntegral.hpp"

namespace Integration
{

  template <typename ReturnTrait1>
  ReturnTrait1 ElementIntegral::integrate(ElementT* el,
                                  ElementIntegrandBase<ReturnTrait1>* integrand,
                                  const QuadratureRulesT * const qdrRule)
  {
      //(@dducks) this lambda definition is not allowed inside the integrate(), why not?
      std::function<ReturnTrait1(const ElementT*, const  PointReferenceT&) > integrandFun = 
      [=](const ElementT* el, const PointReferenceT& p)-> ReturnTrait1{
          ReturnTrait1 result;
          integrand -> elementIntegrand(el, p, result);
          return result;
      };
      return integrate(el, integrandFun, qdrRule);
  }

  template<typename ReturnType>
  ReturnType ElementIntegral::integrate(ElementT* el,
                 std::function<ReturnType(const ElementT*, const  PointReferenceT&) > integrandFun,
                 const QuadratureRulesT * const qdrRule)
  {
    if (localElement_ == nullptr)
    {
      localElement_=new Base::ShortTermStorageElementH1(el->getGaussQuadratureRule()->dimension());
      if (useCache_)
      {
        localElement_->cacheOn();
      } else
      {
        localElement_->cacheOff();
      }
    }
    *localElement_=*el;
    const QuadratureRulesT * const qdrRuleLoc = (qdrRule == nullptr ? localElement_->getGaussQuadratureRule() : qdrRule);

    // check whether the GaussQuadratureRule is actually for the element's ReferenceGeometry
    assert((qdrRuleLoc->forReferenceGeometry() == localElement_->getReferenceGeometry()));

    // value returned by the integrand
    ReturnType value,result;

    // number of Gauss quadrature points
        std::size_t nrOfPoints = qdrRuleLoc->nrOfPoints();
    assert(nrOfPoints > 0);

    // Initialize Gauss quadrature point
    Geometry::PointReference p = qdrRuleLoc->getPoint(0);

    // first Gauss point
    // first we calculate the jacobian, then compute the function value on one of
    // the reference points and finally we multiply this value with a weight and
    // the jacobian and save it in result.
    Geometry::Jacobian jac(qdrRuleLoc->dimension(), qdrRuleLoc->dimension());
    localElement_->calcJacobian(p, jac);
    result = integrandFun(localElement_, p);
    result *= (qdrRuleLoc->weight(0) * std::abs(jac.determinant()));

    // next Gauss points, again calculate the jacobian, value at gauss point and
    // add this value multiplied with jacobian and weight to result.
        for (std::size_t i = 1; i < nrOfPoints; ++i)
    {
      p = qdrRuleLoc->getPoint(i);
      localElement_->calcJacobian(p, jac);
      value = integrandFun(localElement_, p);

      //axpy: Y = alpha * X + Y
      result.axpy(qdrRuleLoc->weight(i) * std::abs(jac.determinant()), value);
    }
    return result;
  }
}



#endif /* ELEMENTINTEGRAL_IMPL_HPP_ */

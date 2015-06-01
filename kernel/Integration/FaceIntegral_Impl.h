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
#ifndef FACEINTEGRAL_IMPL_HPP_
#define FACEINTEGRAL_IMPL_HPP_
#include "Base/ShortTermStorageFaceH1.h"

#include "QuadratureRules/GaussQuadratureRule.h"
#include "Logger.h"
#include "Base/L2Norm.h"
#include "FaceIntegrandBase.h"
#include "LinearAlgebra/Axpy.h"

namespace Integration
{
    
    template<typename ReturnTrait1>
    ReturnTrait1 FaceIntegral::integrate(Base::Face* fa, FaceIntegrandBase<ReturnTrait1>* integrand, const QuadratureRules::GaussQuadratureRule* qdrRule)
    {
        logger.assert(fa!=nullptr, "Invalid face detected");
        logger.assert(integrand!=nullptr, "Invalid integrand detected");
        //quadrature rule is allowed to be equal to nullptr!
        std::function<ReturnTrait1(const Base::Face*, const LinearAlgebra::NumericalVector&, const Geometry::PointReference&)> integrandFunc = [=](const Base::Face* face, const LinearAlgebra::NumericalVector& n, const Geometry::PointReference& p)
        {   
            ReturnTrait1 result;
            integrand->faceIntegrand(face,n,p,result);
            return result;
        };
        return integrate(fa, integrandFunc, qdrRule);
    }
    
    template<typename ReturnTrait1>
    ReturnTrait1 FaceIntegral::integrate(Base::Face* fa, std::function<ReturnTrait1(const Base::Face*, const LinearAlgebra::NumericalVector&, const Geometry::PointReference&)> integrandFunc, const QuadratureRulesT* const qdrRule)
    {
        logger.assert(fa!=nullptr, "Invalid face detected");
        //quadrature rule is allowed to be equal to nullptr!
        if (localFace_ == nullptr)
        {
            localFace_ = new Base::ShortTermStorageFaceH1(fa->getGaussQuadratureRule()->dimension() + 1);
        }
        *localFace_ = *fa;
        const QuadratureRulesT * const qdrRuleLoc = (qdrRule == nullptr ? localFace_->getGaussQuadratureRule() : qdrRule);
        
        // check whether the GaussIntegrationRule is actually for the
        // Element's ReferenceGeometry
        logger.assert((qdrRuleLoc->forReferenceGeometry() == localFace_->getReferenceGeometry()), "FaceIntegral: " + qdrRuleLoc->getName() + " rule is not for THIS ReferenceGeometry!");
        
        // value returned by the integrand
        ReturnTrait1 value, result;
        
        // number of Gauss quadrature points
        std::size_t nrOfPoints = qdrRuleLoc->nrOfPoints();
        
        // Gauss quadrature point
        const Geometry::PointReference& p0 = qdrRuleLoc->getPoint(0);
        
        LinearAlgebra::NumericalVector normal = localFace_->getNormalVector(p0);
        
        // first Gauss point;
        result = integrandFunc(localFace_, normal, p0);
        result *= (qdrRuleLoc->weight(0) * Base::L2Norm(normal));
        
        // next Gauss points
        for (std::size_t i = 1; i < nrOfPoints; ++i)
        {
            const Geometry::PointReference& p = qdrRuleLoc->getPoint(i);
            normal = localFace_->getNormalVector(p);
            value = integrandFunc(localFace_, normal, p);
            
            //Y = alpha * X + Y
            LinearAlgebra::axpy(qdrRuleLoc->weight(i) * Base::L2Norm(normal), value, result);
            
        }
        return result;
    } // function
    
    /// \param[in] ptrQdrRule A pointer to a quadrature rule used for the integration.
    /// \param[in] integrandFunction A function that is integrated on the reference face. It takes as input argument a reference point and returns an object of the class <IntegrandType>.
    /*!
     \details This function computes the integral of a function \f$ f_{ref}\f$ on a reference face \f$ F_{ref} \f$, so it returns the following value
     \f[ \int_{F_{ref}} f_{ref}(\xi) \,d\xi, \f]
     where \f$ f_{ref}:F_{ref}\rightarrow R\f$, with \f$ R \f$ a linear function space (e.g. space of vectors/matrices). In many cases the weak formulation is based on integrals on a physical face \f$ F_{phys} \f$ of the form given below
     \f[ \int_{F_{phys}} f_{phys}(x) \,dx.\f]
     Let \f$ \phi:F_{ref}\rightarrow F_{phys} \f$ be the mapping from the reference face to the physical face, let \f$ J\f$ be the Jacobian of \f$ \phi\f$ and \f$ |J| \f$ the reference-to-physical face scale. Then we can write
     \f[ \int_{F_{phys}} f_{phys}(x) \,dx = \int_{F_{ref}} f_{phys}(\phi(\xi)) |J| \,d\xi, \f]
     so \f$ f_{ref}(\xi) = f_{phys}(\phi(\xi)) |J| \f$. In some cases it is more advantageous to compute \f$ f_{ref}(\xi) \f$ instead of \f$ f_{phys}(x)\f$.
     
     NOTE: do not mix up gradients of pyhsical and reference basis functions with integrals on physical and reference faces. If \f$ f_{phys}(x) \f$ contains a (physical) gradient of a physical basis function then so does \f$ f_{ref}(\xi) = f_{phys}(\phi(\xi)) |J| \f$. The difference is the input argument (reference point \f$ \xi \f$ instead of physical point \f$ x \f$ ) and the scaling \f$ |J| \f$. (Ofcourse it is possible to rewrite the gradient of a physical basis function in terms of the gradient of the corresponding reference basis function).
     */
    template<typename IntegrandType>
    IntegrandType FaceIntegral::referenceFaceIntegral(const QuadratureRules::GaussQuadratureRule *ptrQdrRule, std::function<IntegrandType(const Geometry::PointReference &)> integrandFunction) const
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

#endif /* FACEINTEGRAL_IMPL_HPP_ */

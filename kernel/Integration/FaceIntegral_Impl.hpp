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
#include "Base/ShortTermStorageFaceH1.hpp"

#include "QuadratureRules/GaussQuadratureRule.hpp"
#include "Base/TestErrorDebug.hpp"
#include "Base/L2Norm.hpp"
#include "FaceIntegrandBase.hpp"

namespace Integration{
    
        template <class ReturnTrait1>
        void FaceIntegral::integrate(Base::Face* fa, FaceIntegrandBase<ReturnTrait1>* integrand, ReturnTrait1& result, const QuadratureRules::GaussQuadratureRule* qdrRule)
        {
            std::function<void(const Base::Face*, const LinearAlgebra::NumericalVector&, const Geometry::PointReference&, ReturnTrait1&)> integrandFunc = 
            [=](const Base::Face* face, const LinearAlgebra::NumericalVector& n, const Geometry::PointReference& p, ReturnTrait1& ret){integrand->faceIntegrand(face,n,p,ret);};
            integrate(fa,integrandFunc,result,qdrRule);
        }

template <class ReturnTrait1>
    void
    FaceIntegral::integrate(FaceT*                                           fa,
                                 std::function<void(const Base::Face*, const LinearAlgebra::NumericalVector&, const Geometry::PointReference&, ReturnTrait1&)> integrandFunc,
                                  ReturnTrait1&   result,
                                 const QuadratureRulesT* const                    qdrRule )
    {
		if(localFace_==NULL){
			localFace_=new Base::ShortTermStorageFaceH1(fa->getGaussQuadratureRule()->dimension()+1);
		}
		*localFace_=*fa;
        const QuadratureRulesT* const qdrRuleLoc = (qdrRule==NULL? localFace_->getGaussQuadratureRule(): qdrRule);

            // check whether the GaussIntegrationRule is actually for the
            // Element's ReferenceGeometry
        TestErrorDebug((qdrRuleLoc->forReferenceGeometry() == localFace_->getReferenceGeometry()),
                       "FaceIntegral: " + qdrRuleLoc->getName() + " rule is not for THIS ReferenceGeometry!");

            // value returned by the integrand
         ReturnTrait1 value(result);

            // number of Gauss quadrature points
        unsigned int nrOfPoints = qdrRuleLoc->nrOfPoints();

            // Gauss quadrature point
        Geometry::PointReference p(qdrRuleLoc->dimension());

        //if (!useCache_)//caching of transformation data is delegated to ShortTermStorageBase
        //{
            LinearAlgebra::NumericalVector Normal(qdrRuleLoc->dimension()+1);

                // first Gauss point
            qdrRuleLoc->getPoint(0, p);
            localFace_->getNormalVector(p, Normal);
            integrandFunc(localFace_, Normal, p, result);
            result *= (qdrRuleLoc->weight(0) * Base::L2Norm(Normal));

                // next Gauss points
            for (unsigned int i = 1; i < nrOfPoints; ++i)
            {
                qdrRuleLoc->getPoint(i, p);
                localFace_->getNormalVector(p, Normal);
                integrandFunc(localFace_, Normal, p, value);

                 //Y = alpha * X + Y
                result.axpy(qdrRuleLoc->weight(i) * Base::L2Norm(Normal),value);

            }
        //}
        /*else  // useCache_
        {///\TODO long term cache dont work well with short term cache
                // get vector of cache data
            VecCacheT vecCache = fa->getVecCacheData();

                // Calculate the cache
            if ((vecCache.size()!=nrOfPoints) || recomputeCache_)
            {
                std::cout << qdrRuleLoc->getName() << " ";
                std::cout << "FaceIntegral: filling up the cache (" << nrOfPoints << "points)!\n";

                vecCache.resize(nrOfPoints,qdrRuleLoc->dimension()+1);
                for (unsigned int i=0; i<nrOfPoints; ++i)
                {
                    qdrRuleLoc->getPoint(i, p);
                    vecCache[i](*fa,p);//this non-intuitive bit of notation computes and stores the outward-pointing normal vector and its norm
                }
            }

                // first Gauss point
            qdrRuleLoc->getPoint(0, p);
            integrand->faceIntegrand(fa, vecCache[0].Normal, p, result);
            result *= (qdrRuleLoc->weight(0) * vecCache[0].L2Normal);

                // next Gauss point(s)
            for (unsigned int i = 1; i < nrOfPoints; ++i)
            {
                qdrRuleLoc->getPoint(i, p);
                integrand->faceIntegrand(fa, vecCache[i].Normal, p, value);

                    //Y = alpha * X + Y
                result.axpy(qdrRuleLoc->weight(i) * vecCache[i].L2Normal,value);

            } // for integration points
        } // if cached data (else)*/
    } // function

}



#endif /* FACEINTEGRAL_IMPL_HPP_ */

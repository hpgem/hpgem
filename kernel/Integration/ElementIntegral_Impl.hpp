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
#include "Base/ShortTermStorageElementH1.hpp"

namespace Integration{

template <class ReturnTrait1>
    void
    ElementIntegral::integrate(ElementT* el,
                                    ElementIntegrandBase<ReturnTrait1>* integrand,
                                    ReturnTrait1& result,
                                    const QuadratureRulesT* const qdrRule)
    {
		if(localElement_==NULL){
			localElement_=new Base::ShortTermStorageElementH1(el->getGaussQuadratureRule()->dimension());
			if(useCache_){
				localElement_->cacheOn();
			}else{
				localElement_->cacheOff();
			}
		}
		*localElement_=*el;
        const QuadratureRulesT* const qdrRuleLoc = (qdrRule==NULL? localElement_->getGaussQuadratureRule(): qdrRule);

            // check whether the GaussQuadratureRule is actually for the element's ReferenceGeometry
        TestErrorDebug((qdrRuleLoc->forReferenceGeometry() == localElement_->getReferenceGeometry()),
                       "ElementIntegral: " + qdrRuleLoc->getName() + " rule is not for " + localElement_->getReferenceGeometry()->getName() + "!");

            // value returned by the integrand
        ReturnTrait1 value(result);

            // number of Gauss quadrature points
        unsigned int nrOfPoints = qdrRuleLoc->nrOfPoints();
        TestErrorDebug(nrOfPoints>0,"ElementIntegral: Invalid quadrature rule!");

            // Gauss quadrature point
            Geometry::PointReference p(qdrRuleLoc->dimension());

        //if (!useCache_)//caching of transformation data is delegated to ShortTermStorageBase
        //{
                // first Gauss point
            Geometry::Jacobian jac(qdrRuleLoc->dimension(),qdrRuleLoc->dimension());
            qdrRuleLoc->getPoint(0, p);
            localElement_->calcJacobian(p, jac);
            integrand->elementIntegrand(localElement_, p, result);
            result *= (qdrRuleLoc->weight(0) * std::abs(jac.determinant()));

//            cout <<"Result = "<<result<<endl;

                // next Gauss point(s)
            for (unsigned int i = 1; i < nrOfPoints; ++i)
            {
                qdrRuleLoc->getPoint(i, p);

                localElement_->calcJacobian(p, jac);

//                cout << "p="<<p<<endl;

                integrand->elementIntegrand(localElement_, p, value);

                    //Base::Axpy(qdrRuleLoc->weight(i) * std::abs(jac.determinant()), value, result);

                    //Y = alpha * X + Y

                    //cout<<"std::abs(jac.determinant()="<<std::abs(jac.determinant())<<endl;

//                cout <<"qdrRuleLoc->weight(i)="<<qdrRuleLoc->weight(i)<<endl;

//                cout <<"value="<<value<<endl;

                result.axpy(qdrRuleLoc->weight(i) * std::abs(jac.determinant()), value);

//                cout <<"Result2 = "<<result<<endl;
//                cout <<"*******************************************"<<endl;
            }
        //}
        /*else
        {

                // get vector of cache data
            VecCacheT& vecCache = el->getVecCacheData();

                // Calculate the cache
            if ((vecCache.size()!=nrOfPoints) || recomputeCache_)
            {
                std::cout << qdrRuleLoc->getName() << " ";
                std::cout << "ElementIntegral: filling up the cache (" << nrOfPoints << "points)!\n";

                vecCache.resize(nrOfPoints);
                for (unsigned int i=0; i<nrOfPoints; ++i)
                {
                    qdrRuleLoc->getPoint(i, p);
                    vecCache[i](el,p);//this non-intuitive bit of notation computes and stores the jacobian and its determinant
                }
            }

                // first Gauss point
            qdrRuleLoc->getPoint(0, p);
            integrand->elementIntegrand(el, p, result);
            result *= (qdrRuleLoc->weight(0) * vecCache[0].absDetJac_);

                // next Gauss point(s)
            for (unsigned int i = 1; i < nrOfPoints; ++i)
            {
                qdrRuleLoc->getPoint(i, p);
                integrand->elementIntegrand(el, p, value);
                    //Y = alpha * X + Y
                result.axpy(qdrRuleLoc->weight(i) * vecCache[i].absDetJac_,value);
                    //Base::Axpy(qdrRuleLoc->weight(i) * vecCache[i].absDetJac_, value, result);
            }

        }*/
    }



}



#endif /* ELEMENTINTEGRAL_IMPL_HPP_ */

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
#include "ElementIntegral.hpp"
#include "Geometry/ReferenceGeometry.hpp"
#include "Base/ElementCacheData.hpp"

namespace Integration
{
    
        //! \brief Construct an ElementIntegral with cache on.
    ElementIntegral::ElementIntegral(bool useCache):
        useCache_(useCache)
    {
        localElement_=nullptr;
    }
    
        //! \brief Class destructor
    ElementIntegral::~ElementIntegral()
    {
        if(localElement_!=nullptr)
            delete localElement_;
    }
        //! \brief Start caching (geometry) information now.
    void
    ElementIntegral::cacheOn()
    {
        useCache_ = true;
        if(localElement_!=nullptr){
            localElement_->cacheOn();
        }
    }
    
        //! \brief Stop using cache.
    void
    ElementIntegral::cacheOff()
    {
        useCache_ = false;
        if(localElement_!=nullptr){
        	localElement_->cacheOff();
        }
    }
    
        //! \brief Set recompute the cache ON.
    void
    ElementIntegral::recomputeCacheOn()
    {
        if(localElement_!=nullptr){
        	localElement_->recomputeCacheOn();
        }
    }
    
        //! \brief Set recompute the cache OFF.
    void
    ElementIntegral::recomputeCacheOff()
    {
        if(localElement_!=nullptr){
        	localElement_->recomputeCacheOff();
        }
    }
    
    /*!
     \param[in]  el        the \c Element to be integrated on,
     \param[in]  rule      the GaussQuadratureRule to use,
     \param[in]  integrand a function/functor with operator(Element, PointReference<DIM>,ResultType&),
     \param[out] result    a reference to the variable with result storage.
     
     Note that one now has the possibility to leave the \c rule argument
     away, in which case the default for the \c ReferenceGeometry of
     the passed element will be used. */
    
    
//    template <unsigned int DIM>
//    template <template<unsigned int> class OBJ, typename IntegrandT>
//    void
//    ElementIntegral<DIM>::integrate(ElementT* el,
//                                    IntegrandT& integrand,
//                                    typename ReturnTrait1<IntegrandT>::ReturnType& result,
//                                    OBJ<DIM>*                objPtr,
//                                    const QuadratureRulesT* const qdrRule
//                                    )                                    
//    {
    //        const QuadratureRulesT* const qdrRuleLoc = (qdrRule==nullptr? el->getGaussQuadratureRule(): qdrRule);
//            
//            // check whether the GaussQuadratureRule is actually for the element's ReferenceGeometry
//        logger.assert((qdrRuleLoc->forReferenceGeometry() == el->getReferenceGeometry()),
//                       "ElementIntegral: " + qdrRuleLoc->getName() + " rule is not for " + el->getReferenceGeometry()->getName() + "!");
//        
//            // value returned by the integrand
//        typename ReturnTrait1<IntegrandT>::ReturnType value(result);
//
//            // number of Gauss quadrature points
//        unsigned int nrOfPoints = qdrRuleLoc->nrOfPoints();
//        logger.assert(nrOfPoints>0,"ElementIntegral: Invalid quadrature rule!");
//        
//            // Gauss quadrature point
//            Geometry::PointReference<DIM> p;
//        
//        if (!useCache_)
//        {
//                // first Gauss point
//            Geometry::Jacobian<DIM,DIM> jac;
//            qdrRuleLoc->getPoint(0, p);
//            el->calcJacobian(p, jac);
//            (objPtr->*integrand)(el, p, result);
//            result *= (qdrRuleLoc->weight(0) * std::abs(jac.determinant()));
//            
//            cout <<"Result = "<<result<<endl;
//            
//                // next Gauss point(s)
//            for (unsigned int i = 1; i < nrOfPoints; ++i)
//            {
//                qdrRuleLoc->getPoint(i, p);
//                
//                el->calcJacobian(p, jac);
//                
//                cout << "p="<<p<<endl;
//                
//                (objPtr->*integrand)(el, p, value);
//                
//                    //Base::Axpy(qdrRuleLoc->weight(i) * std::abs(jac.determinant()), value, result);
//                
//                    //Y = alpha * X + Y
//                
//                cout<<"std::abs(jac.determinant()="<<std::abs(jac.determinant())<<endl;
//                                                            
//                cout <<"qdrRuleLoc->weight(i)="<<qdrRuleLoc->weight(i)<<endl;
//                
//                cout <<"value="<<value<<endl;
//                
//                result.axpy(qdrRuleLoc->weight(i) * std::abs(jac.determinant()), value);
//
//                cout <<"Result2 = "<<result<<endl;
//                cout <<"*******************************************"<<endl;
//            }
//        }
//        else
//        {
//            
//                // get vector of cache data
//            VecCacheT& vecCache = el->getVecCacheData();
//            
//                // Calculate the cache
//            if ((vecCache.size()!=nrOfPoints) || recomputeCache_)
//            {
//                std::cout << qdrRuleLoc->getName() << " ";
//                std::cout << "ElementIntegral: filling up the cache (" << nrOfPoints << "points)!\n";
//                
//                vecCache.resize(nrOfPoints);
//                for (unsigned int i=0; i<nrOfPoints; ++i)
//                {
//                    qdrRuleLoc->getPoint(i, p);
//                    vecCache[i](el,p);
//                }
//            }
//            
//                // first Gauss point
//            qdrRuleLoc->getPoint(0, p);
//            (objPtr->*integrand)(el, p, result);
//            result *= (qdrRuleLoc->weight(0) * vecCache[0].absDetJac_);
//            
//                // next Gauss point(s)
//            for (unsigned int i = 1; i < nrOfPoints; ++i)
//            {
//                qdrRuleLoc->getPoint(i, p);
//                (objPtr->*integrand)(el, p, value);
//                    //Y = alpha * X + Y
//                result.axpy(qdrRuleLoc->weight(i) * vecCache[i].absDetJac_,value);
//                    //Base::Axpy(qdrRuleLoc->weight(i) * vecCache[i].absDetJac_, value, result);
//            }
//            
//        }
//    }
}

void Integration::ElementIntegral::setStorageWrapper(Base::ShortTermStorageElementBase* transform) {
    //if(localElement_!=nullptr)
	delete localElement_;
    localElement_=transform;
    
    if(useCache_){
            localElement_->cacheOn();
        
    }else{
            localElement_->cacheOff();
        //std::cout<<"Working storage Wrapper"<<std::endl;
    }
}
    //! \brief AXPY operation, i.e. Y = alpha * X + Y, for various data type
    //        template <typename T>
    //        void Axpy(double alpha, const T& x, T& y)
    //        {
    //            y += alpha * x;
    //        }
    //
    //            //! \brief AXPY operation for Matrix
    //        inline void Axpy(double a, const LinearAlgebra::Matrix& x, LinearAlgebra::Matrix& y)
    //        {
    //            y.axpy(a,x);
    //        }
    //
    //            //! \brief AXPY operation for LinearAlgebra::NumericalVector
    //        inline void Axpy(double alpha,
    //                         const LinearAlgebra::NumericalVector& x,
    //                         LinearAlgebra::NumericalVector& y)
    //        {
    //            if (y.size()) // if not empty
    //                y += (alpha*x);
    //        }
    //
    //            //! \brief AXPY operation for std::vector<double>
    //        inline void Axpy(double alpha,
    //                         const std::vector<double>& x,
    //                         std::vector<double>& y)
    //        {
    //            const unsigned int neq = y.size();
    //            if (neq) // if not empty
    //            {
    //                for (unsigned int i=0;i<neq; ++i)
    //                {
    //                    y[i] += alpha * x[i];
    //                }
    //            }
    //        }

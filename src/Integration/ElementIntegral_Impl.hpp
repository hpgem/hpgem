//
//  ElementIntegral_Impl.hpp
//  
//
//  Created by Shavarsh Nurijanyan on 2/26/13.
//
//

#ifndef _ElementIntegral_Impl_hpp
#define _ElementIntegral_Impl_hpp

namespace Integration
{
    template<unsigned int DIM>
    class ElementIntegral;
    
        //! \brief Construct an ElementIntegral with cache on.
    template<unsigned int DIM>
    ElementIntegral<DIM>::ElementIntegral(bool useCache):
        useCache_(useCache),
        recomputeCache_(false)
    {
        
    }
    
        //! \brief Class destructor
    template<unsigned int DIM>
    ElementIntegral<DIM>::~ElementIntegral()
    {
        
    }
        //! \brief Start caching (geometry) information now.
    template<unsigned int DIM>
    void
    ElementIntegral<DIM>::cacheOn()
    {
        useCache_ = true;
        recomputeCache_ = false;
    }
    
        //! \brief Stop using cache.
    template<unsigned int DIM>
    void
    ElementIntegral<DIM>::cacheOff()
    {
        useCache_ = true;
        recomputeCache_ = false;
    }
    
        //! \brief Set recompute the cache ON.
    template<unsigned int DIM>
    void
    ElementIntegral<DIM>::recomputeCacheOn()
    {
        recomputeCache_ = true;
    }
    
        //! \brief Set recompute the cache OFF.
    template<unsigned int DIM>
    void
    ElementIntegral<DIM>::recomputeCacheOff()
    {
        recomputeCache_ = false;
    }
    
    /*!
     \param[in]  el        the \c Element to be integrated on,
     \param[in]  rule      the GaussQuadratureRule to use,
     \param[in]  integrand a function/functor with operator(Element, PointReference<DIM>,ResultType&),
     \param[out] result    a reference to the variable with result storage.
     
     Note that one now has the possibility to leave the \c rule argument
     away, in which case the default for the \c ReferenceGeometry of
     the passed element will be used. */
    template <unsigned int DIM>
    template <typename IntegrandT>
    void
    ElementIntegral<DIM>::integrate(ElementT* el,
                                    IntegrandT& integrand,
                                    typename ReturnTrait1<IntegrandT>::ReturnType& result,
                                    const QuadratureRulesT* const qdrRule)
    {
        const QuadratureRulesT* const qdrRuleLoc = (qdrRule==NULL? el->getGaussQuadratureRule(): qdrRule);
            
            // check whether the GaussQuadratureRule is actually for the element's ReferenceGeometry
        TestErrorDebug((qdrRuleLoc->forReferenceGeometry() == el->getReferenceGeometry()),
                       "ElementIntegral: " + qdrRuleLoc->getName() + " rule is not for " + el->getReferenceGeometry()->getName() + "!");
        
            // value returned by the integrand
        typename ReturnTrait1<IntegrandT>::ReturnType value(result);

            // number of Gauss quadrature points
        unsigned int nrOfPoints = qdrRuleLoc->nrOfPoints();
        TestErrorDebug(nrOfPoints>0,"ElementIntegral: Invalid quadrature rule!");
        
            // Gauss quadrature point
            Geometry::PointReference<DIM> p;
        
        if (!useCache_)
        {
                // first Gauss point
            Geometry::Jacobian<DIM,DIM> jac;
            qdrRuleLoc->getPoint(0, p);
            el->calcJacobian(p, jac);
            integrand(el, p, result);
            result *= (qdrRuleLoc->weight(0) * std::abs(jac.determinant()));
            
//            cout <<"Result = "<<result<<endl;
            
                // next Gauss point(s)
            for (unsigned int i = 1; i < nrOfPoints; ++i)
            {
                qdrRuleLoc->getPoint(i, p);
                
                el->calcJacobian(p, jac);
                
//                cout << "p="<<p<<endl;
                
                integrand(el, p, value);
                
                    //Base::Axpy(qdrRuleLoc->weight(i) * std::abs(jac.determinant()), value, result);
                
                    //Y = alpha * X + Y
                
                    //cout<<"std::abs(jac.determinant()="<<std::abs(jac.determinant())<<endl;
                                                            
//                cout <<"qdrRuleLoc->weight(i)="<<qdrRuleLoc->weight(i)<<endl;
                
//                cout <<"value="<<value<<endl;
                
                result.axpy(qdrRuleLoc->weight(i) * std::abs(jac.determinant()), value);

//                cout <<"Result2 = "<<result<<endl;
//                cout <<"*******************************************"<<endl;
            }
        }
        else
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
                    vecCache[i](el,p);
                }
            }
            
                // first Gauss point
            qdrRuleLoc->getPoint(0, p);
            integrand(el, p, result);
            result *= (qdrRuleLoc->weight(0) * vecCache[0].absDetJac_);
            
                // next Gauss point(s)
            for (unsigned int i = 1; i < nrOfPoints; ++i)
            {
                qdrRuleLoc->getPoint(i, p);
                integrand(el, p, value);
                    //Y = alpha * X + Y
                result.axpy(qdrRuleLoc->weight(i) * vecCache[i].absDetJac_,value);
                    //Base::Axpy(qdrRuleLoc->weight(i) * vecCache[i].absDetJac_, value, result);
            }
            
        }
    }

    
    template <unsigned int DIM>
    template <typename OBJ, typename IntegrandT>
    void
    ElementIntegral<DIM>::integrate(ElementT* el,
                                    IntegrandT& integrand,
                                    typename ReturnTrait1<IntegrandT>::ReturnType& result,
                                    OBJ*                objPtr,
                                    const QuadratureRulesT* const qdrRule
                                    )
    {
        const QuadratureRulesT* const qdrRuleLoc = (qdrRule==NULL? el->getGaussQuadratureRule(): qdrRule);
        
            // check whether the GaussQuadratureRule is actually for the element's ReferenceGeometry
        TestErrorDebug((qdrRuleLoc->forReferenceGeometry() == el->getReferenceGeometry()),
                       "ElementIntegral: " + qdrRuleLoc->getName() + " rule is not for " + el->getReferenceGeometry()->getName() + "!");
        
            // value returned by the integrand
        typename ReturnTrait1<IntegrandT>::ReturnType value(result);
        
            // number of Gauss quadrature points
        unsigned int nrOfPoints = qdrRuleLoc->nrOfPoints();
        TestErrorDebug(nrOfPoints>0,"ElementIntegral: Invalid quadrature rule!");
        
            // Gauss quadrature point
        Geometry::PointReference<DIM> p;
        
        if (!useCache_)
        {
                // first Gauss point
            Geometry::Jacobian<DIM,DIM> jac;
            qdrRuleLoc->getPoint(0, p);
            el->calcJacobian(p, jac);
            (objPtr->*integrand)(el, p, result);
            result *= (qdrRuleLoc->weight(0) * std::abs(jac.determinant()));
            
//            cout <<"Result = "<<result<<endl;
            
                // next Gauss point(s)
            for (unsigned int i = 1; i < nrOfPoints; ++i)
            {
                qdrRuleLoc->getPoint(i, p);
                
                el->calcJacobian(p, jac);
                
//                cout << "p="<<p<<endl;
                
                (objPtr->*integrand)(el, p, value);
                
                    //Base::Axpy(qdrRuleLoc->weight(i) * std::abs(jac.determinant()), value, result);
                
                    //Y = alpha * X + Y
                
//                cout<<"std::abs(jac.determinant()="<<std::abs(jac.determinant())<<endl;
                
//                cout <<"qdrRuleLoc->weight(i)="<<qdrRuleLoc->weight(i)<<endl;
                
//                cout <<"value="<<value<<endl;
                
                result.axpy(qdrRuleLoc->weight(i) * std::abs(jac.determinant()), value);
                
//                cout <<"Result2 = "<<result<<endl;
//                cout <<"*******************************************"<<endl;
            }
        }
        else
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
                    vecCache[i](el,p);
                }
            }
            
                // first Gauss point
            qdrRuleLoc->getPoint(0, p);
            (objPtr->*integrand)(el, p, result);
            result *= (qdrRuleLoc->weight(0) * vecCache[0].absDetJac_);
            
                // next Gauss point(s)
            for (unsigned int i = 1; i < nrOfPoints; ++i)
            {
                qdrRuleLoc->getPoint(i, p);
                (objPtr->*integrand)(el, p, value);
                    //Y = alpha * X + Y
                result.axpy(qdrRuleLoc->weight(i) * vecCache[i].absDetJac_,value);
                    //Base::Axpy(qdrRuleLoc->weight(i) * vecCache[i].absDetJac_, value, result);
            }
            
        }
    }
    
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
//        const QuadratureRulesT* const qdrRuleLoc = (qdrRule==NULL? el->getGaussQuadratureRule(): qdrRule);
//            
//            // check whether the GaussQuadratureRule is actually for the element's ReferenceGeometry
//        TestErrorDebug((qdrRuleLoc->forReferenceGeometry() == el->getReferenceGeometry()),
//                       "ElementIntegral: " + qdrRuleLoc->getName() + " rule is not for " + el->getReferenceGeometry()->getName() + "!");
//        
//            // value returned by the integrand
//        typename ReturnTrait1<IntegrandT>::ReturnType value(result);
//
//            // number of Gauss quadrature points
//        unsigned int nrOfPoints = qdrRuleLoc->nrOfPoints();
//        TestErrorDebug(nrOfPoints>0,"ElementIntegral: Invalid quadrature rule!");
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
#endif

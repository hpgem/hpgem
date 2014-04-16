/*
 * ElementIntegral_Impl.hpp
 *
 *  Created on: Feb 3, 2014
 *      Author: brinkf
 */

#ifndef ELEMENTINTEGRAL_IMPL_HPP_
#define ELEMENTINTEGRAL_IMPL_HPP_

namespace Integration{

template <class ReturnTrait1>
    void
    ElementIntegral::integrate(ElementT* el,
                                    ElementIntegrandBase<ReturnTrait1>* integrand,
                                    ReturnTrait1& result,
                                    const QuadratureRulesT* const qdrRule)
    {
        const QuadratureRulesT* const qdrRuleLoc = (qdrRule==NULL? el->getGaussQuadratureRule(): qdrRule);

            // check whether the GaussQuadratureRule is actually for the element's ReferenceGeometry
        TestErrorDebug((qdrRuleLoc->forReferenceGeometry() == el->getReferenceGeometry()),
                       "ElementIntegral: " + qdrRuleLoc->getName() + " rule is not for " + el->getReferenceGeometry()->getName() + "!");

            // value returned by the integrand
        ReturnTrait1 value(result);

            // number of Gauss quadrature points
        unsigned int nrOfPoints = qdrRuleLoc->nrOfPoints();
        TestErrorDebug(nrOfPoints>0,"ElementIntegral: Invalid quadrature rule!");

            // Gauss quadrature point
            Geometry::PointReference p(qdrRuleLoc->dimension());

        if (!useCache_)
        {
                // first Gauss point
            Geometry::Jacobian jac(qdrRuleLoc->dimension(),qdrRuleLoc->dimension());
            qdrRuleLoc->getPoint(0, p);
            el->calcJacobian(p, jac);
            integrand->elementIntegrand(el, p, result);
            result *= (qdrRuleLoc->weight(0) * std::abs(jac.determinant()));

//            cout <<"Result = "<<result<<endl;

                // next Gauss point(s)
            for (unsigned int i = 1; i < nrOfPoints; ++i)
            {
                qdrRuleLoc->getPoint(i, p);

                el->calcJacobian(p, jac);

//                cout << "p="<<p<<endl;

                integrand->elementIntegrand(el, p, value);

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

        }
    }



}



#endif /* ELEMENTINTEGRAL_IMPL_HPP_ */

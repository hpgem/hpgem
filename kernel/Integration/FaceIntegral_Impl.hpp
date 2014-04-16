/*
 * FaceIntegral_Impl.hpp
 *
 *  Created on: Feb 3, 2014
 *      Author: brinkf
 */

#ifndef FACEINTEGRAL_IMPL_HPP_
#define FACEINTEGRAL_IMPL_HPP_

namespace Integration{

template <class ReturnTrait1>
    void
    FaceIntegral::integrate(FaceT*                                           fa,
                                 FaceIntegrandBase<ReturnTrait1>*                                      integrand,
                                  ReturnTrait1&   result,
                                 const QuadratureRulesT* const                    qdrRule )
    {
        const QuadratureRulesT* const qdrRuleLoc = (qdrRule==NULL? fa->getGaussQuadratureRule(): qdrRule);

            // check whether the GaussIntegrationRule is actually for the
            // Element's ReferenceGeometry
        TestErrorDebug((qdrRuleLoc->forReferenceGeometry() == fa->getReferenceGeometry()),
                       "FaceIntegral: " + qdrRuleLoc->getName() + " rule is not for THIS ReferenceGeometry!");

            // value returned by the integrand
         ReturnTrait1 value(result);

            // number of Gauss quadrature points
        unsigned int nrOfPoints = qdrRuleLoc->nrOfPoints();

            // Gauss quadrature point
        Geometry::PointReference p(qdrRuleLoc->dimension());

        if (!useCache_)
        {
            LinearAlgebra::NumericalVector Normal(qdrRuleLoc->dimension()+1);

                // first Gauss point
            qdrRuleLoc->getPoint(0, p);
            fa->getNormalVector(p, Normal);
            integrand->faceIntegrand(fa, Normal, p, result);
            result *= (qdrRuleLoc->weight(0) * Base::L2Norm(Normal));

                // next Gauss points
            for (unsigned int i = 1; i < nrOfPoints; ++i)
            {
                qdrRuleLoc->getPoint(i, p);
                fa->getNormalVector(p, Normal);
                integrand->faceIntegrand(fa, Normal, p, value);

                 //Y = alpha * X + Y
                result.axpy(qdrRuleLoc->weight(i) * Base::L2Norm(Normal),value);

            }
        }
        else  // useCache_
        {
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
        } // if cached data (else)
    } // function

}



#endif /* FACEINTEGRAL_IMPL_HPP_ */

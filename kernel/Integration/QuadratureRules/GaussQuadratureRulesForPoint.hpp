/*
 * GaussQuadratureRulesForPoint.h
 *
 *  Created on: Feb 7, 2014
 *      Author: brinkf
 */

#ifndef GAUSSQUADRATURERULESFORPOINT_H_
#define GAUSSQUADRATURERULESFORPOINT_H_

#include "GaussQuadratureRule.hpp"
#include "Base/TestErrorDebug.hpp"
#include "Geometry/ReferencePoint.hpp"
#include <cstdint>

namespace QuadratureRules {

///'Quadrature rule' for a single point is always exact (and just does an evaluation)
class Cn0_inf_1: public QuadratureRules::GaussQuadratureRule {
public:
    static Cn0_inf_1& Instance()
        {
            static Cn0_inf_1 theInstance;
            return theInstance;
        }

    virtual std::string             getName() const{return name_;}

    virtual unsigned int            order() const{return INT32_MAX;}

    virtual unsigned int            dimension() const{return 0;}

    virtual unsigned int            nrOfPoints() const{return 1;}

    virtual double                  weight(unsigned int i) const{TestErrorDebug(i==0,"Cn0_inf_1: This quadrature rule only has one point!");return weight_[0];}

    virtual void                    getPoint(unsigned int i, Geometry::PointReference& p) const{TestErrorDebug(i==0,"Cn0_inf_1: This quadrature rule only has one point!");}

    virtual Geometry::ReferenceGeometry*     forReferenceGeometry() const{return refGeoPtr_;}

private:
    Cn0_inf_1():
    	refGeoPtr_(&Geometry::ReferencePoint::Instance()),
    	name_("Cn0_inf_1")
	{
    	weight_[0]=1;
    	refGeoPtr_->addGaussQuadratureRule(this);
	}
    Cn0_inf_1(const Cn0_inf_1&):refGeoPtr_(NULL){}
    virtual ~Cn0_inf_1(){}
private:
    const std::string           name_;
    double                      weight_[1];
    Geometry::ReferenceGeometry* const   refGeoPtr_;
};

} /* namespace QuadratureRules */

#endif /* GAUSSQUADRATURERULESFORPOINT_H_ */

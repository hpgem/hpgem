/*
 * ElementCacheData.cpp
 *
 *  Created on: Jan 31, 2014
 *      Author: brinkf
 */

#include "ElementCacheData.hpp"
#include "Element.hpp"

void Base::ElementCacheData::operator ()(const Element* el, const Geometry::PointReference& p){
    Geometry::Jacobian jac(p.size(),p.size());
    el->calcJacobian(p, jac);
    absDetJac_ = std::abs(jac.determinant());

}



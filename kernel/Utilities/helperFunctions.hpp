/*
 * helperFunctions.hpp
 *
 *  Created on: Mar 5, 2014
 *      Author: brinkf
 */

#ifndef HELPERFUNCTIONS_HPP_
#define HELPERFUNCTIONS_HPP_

namespace Geometry{
class PointReference;
}

namespace Utilities{

double LobattoPolynomial(int degree, double x);

double LobattoPolynomialDerivative(int degree, double x);

double LegendrePolynomial(int degree, double x);

double LegendrePolynomialDerivative(int degree, double x);

double baricentric_3D(int node, const Geometry::PointReference& p);

double baricentric_2D(int node, const Geometry::PointReference& p);

double baricentricDeriv(int node, int direction);
}


#endif /* HELPERFUNCTIONS_HPP_ */

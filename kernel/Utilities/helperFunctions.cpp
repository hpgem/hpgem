#include "helperFunctions.hpp"

#include <math.h>
#include "Geometry/PointReference.hpp"
namespace Utilities {

	double LobattoPolynomial(int degree, double x) {
		switch (degree) {
		case 0:
			return -2 * sqrt(3. / 2.);
		case 1:
			return -2 * sqrt(5. / 2.) * x;
		case 2:
			return -2 * sqrt(7. / 2.) * (1.25 * x * x - .25);
		case 3:
			return -2 * sqrt(9. / 2.) * (1.75 * x * x * x - .75 * x);
		default:
			throw "lobatto polynomials of this order have not been implemented";
		}
	}

	double LobattoPolynomialDerivative(int degree, double x) {
		switch (degree) {
		case 0:
			return 0;
		case 1:
			return -2 * sqrt(5. / 2.);
		case 2:
			return -2 * sqrt(7. / 2.) * (2.5 * x);
		case 3:
			return -2 * sqrt(9. / 2.) * (5.25 * x * x - .75);
		default:
			throw "derivatives of lobatto polynomials of this order have not been implemented";
		}
	}

	double LegendrePolynomial(int degree, double x) {
		switch (degree) {
		case 0:
			return 1;
		case 1:
			return x;
		default:
			return (2. * degree - 1.) / double(degree) * x * LegendrePolynomial(degree - 1, x)
					- (degree - 1.) / double(degree) * LegendrePolynomial(degree - 2, x);
		}
	}

	double LegendrePolynomialDerivative(int degree, double x) {
		switch (degree) {
		case 0:
			return 0;
		case 1:
			return 1;
		default:
			return (2. * degree - 1.) / double(degree) * (x * LegendrePolynomialDerivative(degree - 1, x) + LegendrePolynomial(degree - 1, x))
					- (degree - 1.) / double(degree) * LegendrePolynomialDerivative(degree - 2, x);
		}
	}

	double baricentric_3D(int node, const Geometry::PointReference& p) {
		if (node == 0) {
			return 1 - p[0] - p[1] - p[2];
		} else {
			return p[node - 1];
		}
	}

	double baricentric_2D(int node, const Geometry::PointReference& p) {
		if (node == 0) {
			return 1 - p[0] - p[1];
		} else {
			return p[node - 1];
		}
	}

	double baricentricDeriv(int node, int direction) {
		if (node == 0) {
			return -1;
		} else if (node == direction + 1) {
			return 1;
		} else {
			return 0;
		}
	}

}

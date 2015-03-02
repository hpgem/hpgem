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

#include "helperFunctions.hpp"

#include <cmath>
#include "Geometry/PointReference.hpp"
namespace Utilities {

	double LobattoPolynomial(int degree, double x) {
		switch (degree) {
		case 0:
			return std::sqrt(3. / 2.) * -2.;
		case 1:
			return std::sqrt(5. / 2.) * x * -2.;
		case 2:
			return std::sqrt(7. / 2.) * (5 * x * x - 1) / -2.;
		case 3:
			return std::sqrt(9. / 2.) * (7 * x * x - 3) * x / -2.;
		case 4://higher order cases currently implemented only to support automatic testing of quadrature rules
			return std::sqrt(11./ 2.) * ((21 * x * x - 14) * x * x + 1) / -4.;
		case 5:
			return std::sqrt(13./ 2.) * ((33 * x * x - 30 * x * x) + 5) * x / -4.;
		case 6:
			return std::sqrt(15./ 2.) * (((429 * x * x - 495) * x * x + 135) * x * x - 5) / -32.;
		case 7:
			return std::sqrt(17./ 2.) * (((715 * x * x - 1001) * x * x + 385) * x * x - 35) * x / -32.;
		case 8:
			return std::sqrt(19./ 2.) * ((((2431 * x * x - 4004) * x * x + 2002) * x * x - 308) * x * x + 7) / -64.;
		case 9:
			throw "cannot find this one";
		default:
			throw "lobatto polynomials of this order have not been implemented";
		}
	}

	double LobattoPolynomialDerivative(int degree, double x) {
		switch (degree) {
		case 0:
			return 0;
		case 1:
			return -2 * std::sqrt(5. / 2.);
		case 2:
			return -2 * std::sqrt(7. / 2.) * (2.5 * x);
		case 3:
			return -2 * std::sqrt(9. / 2.) * (5.25 * x * x - .75);
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

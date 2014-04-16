/*
 * 020BasisFunction_UnitTest.cpp
 *
 *  Created on: Apr 15, 2014
 *      Author: brinkf
 */

//naming convention: <Digit><ClassName>_UnitTest.cpp where <Digit> is a number that will make sure
//the unit tests are ordered such that the first failing unit test indicate the culprit class and
//other 'unit' tests may assume correct execution of all prior unit tests

//I can only check that basisFunctions return finite numbers for some of the points that you can pass to them
//actually validating the correct implementation of the promised basisfunction requires computing said basisfunction
//however, the 'unit' tests for the quadrature rules and most of the self tests rely quite heavily on the correct functionality
//of the basisfunctions so some assurance can be gained from the entire test suite.
//Testing derivatives is much easier since it relies on using the same numerical approximation for all basis-function over and over again
//(but of course this is not as accurate as the actual derivative should be) -FB

#include "Base/AssembleBasisFunctionSet.hpp"
#include "cassert"

int main() {

	// 1D

	Base::BasisFunctionSet all1DbasisFunctions(5);//WARNING: this breaks the ordering of the unit tests, but it is basically the only way to collect all basisfunctions in an indexable way
	Base::AssembleBasisFunctionSet_1D_Ord5_A0(all1DbasisFunctions);
	Geometry::PointReference point1D(1);
	LinearAlgebra::NumericalVector ret(1);
	for(int i=0;i<all1DbasisFunctions.size();++i){
		const Base::BaseBasisFunction* test = all1DbasisFunctions[i];
		for(point1D[0]=-1.5;point1D[0]<1.51;point1D[0]+=0.1){
			test->eval(point1D,ret);
			assert(("eval",(test->eval(point1D)-ret[0])<1e-12));

			point1D[0]+=-1.e-8;
			double x0=test->eval(point1D);
			point1D[0]+=2.e-8;
			double x1=test->eval(point1D);

			point1D[0]+=-1e-8;
			test->evalDeriv(point1D,ret);
			assert(("derivative",fabs(ret[0]-5.e7*(x1-x0))<1e-5));
			assert(("derivative",fabs(test->evalDeriv0(point1D)-5.e7*(x1-x0))<1e-5));
		}
	}

	// 2D

	Base::BasisFunctionSet all2DbasisFunctions(5);//WARNING: this breaks the ordering of the unit tests, but it is basically the only way to collect all basisfunctions in an indexable way
	Base::AssembleBasisFunctionSet_2D_Ord5_A1(all2DbasisFunctions);
	Geometry::PointReference point2D(2);
	for(int i=0;i<all2DbasisFunctions.size();++i){
		const Base::BaseBasisFunction* test = all2DbasisFunctions[i];
		for(point2D[0]=-1.5;point2D[0]<1.51;point2D[0]+=0.1){
			for(point2D[1]=-1.5;point2D[1]<1.51;point2D[1]+=0.1){
				test->eval(point2D,ret);
				assert(("eval",(test->eval(point2D)-ret[0])<1e-12));

				point2D[0]+=-1.e-8;
				double x0=test->eval(point2D);
				point2D[0]+=2.e-8;
				double x1=test->eval(point2D);

				point2D[0]+=-1e-8;
				ret.resize(2);
				test->evalDeriv(point2D,ret);
				assert(("derivative",fabs(ret[0]-5.e7*(x1-x0))<1e-5));
				assert(("derivative",fabs(test->evalDeriv0(point2D)-5.e7*(x1-x0))<1e-5));

				point2D[1]+=-1.e-8;
				x0=test->eval(point2D);
				point2D[1]+=2.e-8;
				x1=test->eval(point2D);

				point2D[1]+=-1e-8;
				assert(("derivative",fabs(ret[1]-5.e7*(x1-x0))<1e-5));
				assert(("derivative",fabs(test->evalDeriv1(point2D)-5.e7*(x1-x0))<1e-5));

				ret.resize(1);
			}
		}
	}

	//3D

	Base::BasisFunctionSet all3DbasisFunctions(5);//WARNING: this breaks the ordering of the unit tests, but it is basically the only way to collect all basisfunctions in an indexable way
	Base::AssembleBasisFunctionSet_3D_Ord5_A1(all3DbasisFunctions);
	Geometry::PointReference point3D(3);
	for(int i=0;i<all3DbasisFunctions.size();++i){
		const Base::BaseBasisFunction* test = all3DbasisFunctions[i];
		for(point3D[0]=-1.5;point3D[0]<1.51;point3D[0]+=0.15){
			for(point3D[1]=-1.5;point3D[1]<1.51;point3D[1]+=0.15){
				for(point3D[2]=-1.5;point3D[2]<1.51;point3D[2]+=0.15){
					test->eval(point3D,ret);
					assert(("eval",(test->eval(point3D)-ret[0])<1e-12));

					point3D[0]+=-1.e-8;
					double x0=test->eval(point3D);
					point3D[0]+=2.e-8;
					double x1=test->eval(point3D);

					point3D[0]+=-1e-8;
					ret.resize(3);
					test->evalDeriv(point3D,ret);
					assert(("derivative",fabs(ret[0]-5.e7*(x1-x0))<1e-5));
					assert(("derivative",fabs(test->evalDeriv0(point3D)-5.e7*(x1-x0))<1e-5));

					point3D[1]+=-1.e-8;
					x0=test->eval(point3D);
					point3D[1]+=2.e-8;
					x1=test->eval(point3D);

					point3D[1]+=-1e-8;
					assert(("derivative",fabs(ret[1]-5.e7*(x1-x0))<1e-5));
					assert(("derivative",fabs(test->evalDeriv1(point3D)-5.e7*(x1-x0))<1e-5));

					point3D[2]+=-1.e-8;
					x0=test->eval(point3D);
					point3D[2]+=2.e-8;
					x1=test->eval(point3D);

					point3D[2]+=-1e-8;
					assert(("derivative",fabs(ret[2]-5.e7*(x1-x0))<1e-5));
					assert(("derivative",fabs(test->evalDeriv2(point3D)-5.e7*(x1-x0))<1e-5));

					ret.resize(1);
				}
			}
		}
	}

	return 0;
}
